#include "BoundaryConditions.h"
#include "DoubleObstacle.h"
#include "DrivingForce.h"
#include "InterfaceProperties.h"
#include "InterfaceRegularization.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <ctime>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <tuple>
#include <vector>

using namespace openphase;

namespace {
volatile std::sig_atomic_t stop_requested = 0;
void request_stop(int) { stop_requested = 1; }

constexpr int NX = 100, NY = 2400;
constexpr double H = 0.01, YMIN = -12.0;
constexpr std::array<int, 3> LO{225, 1125, 0};
constexpr std::array<int, 3> HI{1125, 2325, 0};
constexpr std::array<const char*, 2> PAIR_NAMES{"short_L7", "long_L10"};

double number(const char* value, const char* name) {
    size_t used = 0; const double x = std::stod(value, &used);
    if (used != std::string(value).size() || !std::isfinite(x))
        throw std::runtime_error(std::string("invalid ") + name);
    return x;
}
int integer(const char* value, const char* name) {
    size_t used = 0; const long x = std::stol(value, &used);
    if (used != std::string(value).size() || x < 0 || x > 2147483647L)
        throw std::runtime_error(std::string("invalid ") + name);
    return static_cast<int>(x);
}

struct Point { double x, y; };
double triangle_negative_area(const std::array<Point,3>& p,
                              const std::array<double,3>& f,
                              double full_area) {
    const auto [mn, mx] = std::minmax_element(f.begin(), f.end());
    if (*mx <= 0.0 && *mn < 0.0) return full_area;
    if (!(*mn <= 0.0 && *mx >= 0.0) || (*mn == 0.0 && *mx == 0.0)) return 0.0;
    std::array<Point,6> poly{}; int n = 0;
    for (int i=0;i<3;++i) {
        const int j=(i+1)%3;
        if (f[i] <= 0.0) poly[n++]=p[i];
        if (f[i]*f[j] < 0.0) {
            const double t=-f[i]/(f[j]-f[i]);
            poly[n++]={p[i].x+t*(p[j].x-p[i].x), p[i].y+t*(p[j].y-p[i].y)};
        }
    }
    if (n < 3) return 0.0;
    double twice=0.0;
    for (int i=0;i<n;++i) {
        const auto& a=poly[i]; const auto& b=poly[(i+1)%n];
        twice += a.x*b.y-a.y*b.x;
    }
    return 0.5*std::abs(twice);
}

struct PairMetric {
    double area{}, central_f{}, diffuse_integral{}, edge_min{std::numeric_limits<double>::infinity()};
    std::vector<double> roots;
};
struct Sample {
    int step{}; double raw{}, physical{}, wall{};
    double min_phi{std::numeric_limits<double>::infinity()};
    double max_phi{-std::numeric_limits<double>::infinity()};
    double partition_linf{}, x_swap_linf{};
    int negative_f_outside{};
    std::array<double,3> phase_volume{}, phase_volume_ratio{};
    std::array<int,3> phase_exists{};
    std::array<uint64_t,3> phase_entries{};
    std::array<PairMetric,2> pair;
};

double fnode(const PhaseField& phase, const size_t ids[3], int i, int j) {
    const NodePF& n=phase.Fields(i,j,0);
    return n.get_value(ids[0])-std::max(n.get_value(ids[1]),n.get_value(ids[2]));
}

PairMetric measure_pair(const PhaseField& phase, const size_t ids[3], int pidx) {
    const int lo=LO[pidx], hi=HI[pidx], ny=hi-lo;
    PairMetric m;
    // Same reconstruction as the frozen reader: constant edge extension,
    // cell-centre samples, SW-NE diagonal, and clipped linear F=0 triangles.
    auto val=[&](int r,int c) {
        const int jj = (r==0) ? lo : (r==ny+1 ? hi-1 : lo+r-1);
        const int ii = (c==0) ? 0 : (c==NX+1 ? NX-1 : c-1);
        return fnode(phase,ids,ii,jj);
    };
    auto xpos=[](int c) { return c==0 ? 0.0 : (c==NX+1 ? 1.0 : (c-0.5)*H); };
    auto ypos=[&](int r) { return r==0 ? 0.0 : (r==ny+1 ? ny*H : (r-0.5)*H); };
    for (int r=0;r<=ny;++r) for (int c=0;c<=NX;++c) {
        const double x0=xpos(c), x1=xpos(c+1), y0=ypos(r), y1=ypos(r+1);
        const double a=(x1-x0)*(y1-y0)*0.5;
        const Point sw{x0,y0}, se{x1,y0}, ne{x1,y1}, nw{x0,y1};
        const double fsw=val(r,c), fse=val(r,c+1), fne=val(r+1,c+1), fnw=val(r+1,c);
        m.area += triangle_negative_area({sw,se,ne},{fsw,fse,fne},a);
        m.area += triangle_negative_area({sw,ne,nw},{fsw,fne,fnw},a);
    }
    const int center=lo+ny/2;
    m.central_f=0.25*(fnode(phase,ids,49,center-1)+fnode(phase,ids,50,center-1)
                     +fnode(phase,ids,49,center)+fnode(phase,ids,50,center));
    for(int i=0;i<NX;++i) {
        m.edge_min=std::min(m.edge_min,std::min(fnode(phase,ids,i,lo),fnode(phase,ids,i,hi-1)));
    }
    for(int j=lo;j<hi;++j) for(int i=0;i<NX;++i) {
        const NodePF& n=phase.Fields(i,j,0);
        m.diffuse_integral += (n.get_value(ids[1])+n.get_value(ids[2]))*H*H;
    }
    std::vector<double> line(static_cast<size_t>(ny));
    for(int j=lo;j<hi;++j) line[static_cast<size_t>(j-lo)]=0.5*(fnode(phase,ids,49,j)+fnode(phase,ids,50,j));
    for(int k=0;k<ny-1;++k) {
        const double a0=line[static_cast<size_t>(k)], b0=line[static_cast<size_t>(k+1)];
        const double y=YMIN+(lo+k+0.5)*H;
        if (a0==0.0) m.roots.push_back(y);
        if (a0*b0<0.0) m.roots.push_back(y-H*a0/(b0-a0));
    }
    if (!line.empty() && line.back()==0.0) m.roots.push_back(YMIN+(hi-0.5)*H);
    std::sort(m.roots.begin(),m.roots.end());
    m.roots.erase(std::unique(m.roots.begin(),m.roots.end()),m.roots.end());
    return m;
}

Sample inspect(const PhaseField& phase,const size_t ids[3],int step,double dt,double factor,double wall) {
    if(!std::isfinite(dt) || dt<=0.0 || !std::isfinite(factor) || factor<=0.0 || !std::isfinite(wall))
        throw std::runtime_error("nonfinite or nonpositive observation clock");
    Sample s; s.step=step; s.raw=step*dt; s.physical=s.raw/factor; s.wall=wall;
    for(int p=0;p<3;++p) {
        const auto& grain=phase.FieldsProperties[ids[p]];
        if(!std::isfinite(grain.Volume) || !std::isfinite(grain.VolumeRatio))
            throw std::runtime_error("nonfinite active-volume metadata");
        s.phase_volume[static_cast<size_t>(p)]=grain.Volume;
        s.phase_volume_ratio[static_cast<size_t>(p)]=grain.VolumeRatio;
        s.phase_exists[static_cast<size_t>(p)]=grain.Exist?1:0;
    }
    std::array<bool,NY> inside{};
    for(int j=LO[0];j<HI[1];++j) inside[static_cast<size_t>(j)]=true;
    for(int j=0;j<NY;++j) for(int i=0;i<NX;++i) {
        const NodePF& n=phase.Fields(i,j,0);
        std::array<double,3> v{}; double sum=0.0;
        for(int p=0;p<3;++p) { v[p]=n.get_value(ids[p]);
            if(!std::isfinite(v[p])) throw std::runtime_error("nonfinite phase-field value");
            sum+=v[p];
            s.phase_entries[static_cast<size_t>(p)] += n.present(ids[p]) ? 1 : 0;
            s.min_phi=std::min(s.min_phi,v[p]); s.max_phi=std::max(s.max_phi,v[p]); }
        s.partition_linf=std::max(s.partition_linf,std::abs(sum-1.0));
        s.x_swap_linf=std::max(s.x_swap_linf,std::abs(v[0]-phase.Fields(NX-1-i,j,0).get_value(ids[0])));
        s.x_swap_linf=std::max(s.x_swap_linf,std::abs(v[1]-phase.Fields(NX-1-i,j,0).get_value(ids[2])));
        s.x_swap_linf=std::max(s.x_swap_linf,std::abs(v[2]-phase.Fields(NX-1-i,j,0).get_value(ids[1])));
        const double f=v[0]-std::max(v[1],v[2]);
        if(!inside[static_cast<size_t>(j)] && f<0.0) ++s.negative_f_outside;
    }
    s.pair[0]=measure_pair(phase,ids,0); s.pair[1]=measure_pair(phase,ids,1);
    if(!std::isfinite(s.raw) || !std::isfinite(s.physical) || !std::isfinite(s.min_phi) ||
       !std::isfinite(s.max_phi) || !std::isfinite(s.partition_linf) || !std::isfinite(s.x_swap_linf))
        throw std::runtime_error("nonfinite scalar observation");
    for(const auto& pair:s.pair) {
        if(!std::isfinite(pair.area) || !std::isfinite(pair.central_f) ||
           !std::isfinite(pair.diffuse_integral) || !std::isfinite(pair.edge_min))
            throw std::runtime_error("nonfinite pair observation");
        for(double root:pair.roots) if(!std::isfinite(root))
            throw std::runtime_error("nonfinite crossing observation");
    }
    for(int p=0;p<3;++p) {
        if(s.phase_exists[static_cast<size_t>(p)] &&
           std::abs(s.phase_volume_ratio[static_cast<size_t>(p)]-1.0)>1.0e-12)
            throw std::runtime_error("active VolumeRatio differs from one");
        if(!s.phase_exists[static_cast<size_t>(p)] && s.phase_volume[static_cast<size_t>(p)]!=0.0)
            throw std::runtime_error("absent phase has nonzero volume");
    }
    return s;
}

std::vector<double> capture(const PhaseField& phase,const size_t ids[3]) {
    std::vector<double> out(static_cast<size_t>(3)*NY*NX);
    size_t k=0;
    for(int p=0;p<3;++p) for(int j=0;j<NY;++j) for(int i=0;i<NX;++i) {
        const double value=phase.Fields(i,j,0).get_value(ids[p]);
        if(!std::isfinite(value)) throw std::runtime_error("nonfinite phase field before archive write");
        out[k++]=value;
    }
    return out;
}
void set_initial_stable_volumes(PhaseField& phase,const size_t ids[3]) {
    for(int p=0;p<3;++p) { auto& g=phase.FieldsProperties[ids[p]]; if(g.Volume>0.0) {
        g.Stage=GrainStages::Stable; g.RefVolume=g.Volume; g.MAXVolume=g.Volume;
    }}
}
void load_2d(PhaseField& phase,const size_t ids[3],const std::string& path) {
    std::ifstream in(path,std::ios::binary); if(!in) throw std::runtime_error("cannot open parent field");
    std::vector<double> a(static_cast<size_t>(3)*NY*NX);
    in.read(reinterpret_cast<char*>(a.data()),static_cast<std::streamsize>(a.size()*sizeof(double)));
    if(!in || in.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("wrong input byte length");
    size_t k=0; for(int p=0;p<3;++p) for(int j=0;j<NY;++j) for(int i=0;i<NX;++i) {
        if(!std::isfinite(a[k])) throw std::runtime_error("nonfinite imported field");
        phase.Fields(i,j,0).set_value(ids[p],a[k++]);
    }
}

struct HaloSource { int i, j; double value; };
struct HaloStep {
    uint64_t source_cells{};
    uint64_t halo_cells{};
    uint64_t capillary_entries_total{};
    uint64_t capillary_entries_in_halo{};
    uint64_t scaled_entries{};
    std::array<uint64_t,3> pair_entries{};
    uint64_t nonzero_value2{};
    double max_source_product{};
    double max_H{};
    double max_multiplier{1.0};
    uint64_t limiter_candidate_nodes{};
    uint64_t limiter_changed_nodes{};
    uint64_t limiter_changed_entries{};
    double max_limiter_delta{};
};
struct HaloTotals {
    uint64_t steps{};
    uint64_t source_cells_sum{};
    uint64_t halo_cells_sum{};
    uint64_t capillary_entries_total{};
    uint64_t capillary_entries_in_halo{};
    uint64_t scaled_entries{};
    std::array<uint64_t,3> pair_entries{};
    uint64_t nonzero_value2{};
    uint64_t limiter_candidate_nodes{};
    uint64_t limiter_changed_nodes{};
    uint64_t limiter_changed_entries{};
    double max_limiter_delta{};
    double max_source_product{};
    double max_H{};
    double max_multiplier{1.0};
    uint64_t max_source_cells{};
    uint64_t max_halo_cells{};
};
struct DotSnapshot {
    int i, j;
    size_t a, b;
    double value1, value2;
};
struct HaloModel {
    double chi{};
    int radius{};
    std::vector<double> raw;
    std::vector<double> h;
    std::vector<double> work;
    std::vector<HaloSource> sources;
    HaloStep last;
    HaloTotals total;
    std::vector<DotSnapshot> limiter_snapshot;
};

int pair_slot(size_t a, size_t b, const size_t ids[3]) {
    if ((a==ids[0]&&b==ids[1])||(a==ids[1]&&b==ids[0])) return 0;
    if ((a==ids[0]&&b==ids[2])||(a==ids[2]&&b==ids[0])) return 1;
    if ((a==ids[1]&&b==ids[2])||(a==ids[2]&&b==ids[1])) return 2;
    return -1;
}

void build_square_max(const std::vector<double>& raw, int nx, int ny, int radius,
                      std::vector<double>& out, std::vector<double>& work,
                      std::vector<HaloSource>* saved_sources=nullptr) {
    if (nx<=0||ny<=0||radius<0||raw.size()!=static_cast<size_t>(nx)*ny)
        throw std::runtime_error("invalid square-max filter dimensions");
    out.assign(raw.size(),0.0);
    work.assign(raw.size(),0.0);
    std::vector<HaloSource> local_sources;
    auto& sources=saved_sources ? *saved_sources : local_sources;
    sources.clear();
    for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) {
        const double v=raw[static_cast<size_t>(j)*nx+i];
        if(!std::isfinite(v)||v<0.0) throw std::runtime_error("invalid max-filter source value");
        if(v>0.0) sources.push_back({i,j,v});
    }
    // Exact separable square max filter. This remains O(nx*ny) even when
    // a full ternary footprint occupies much of the domain.
    for(int j=0;j<ny;++j) {
        std::deque<int> q; int right=-1;
        for(int i=0;i<nx;++i) {
            const int want_right=std::min(nx-1,i+radius);
            while(right<want_right) {
                ++right;
                const double v=raw[static_cast<size_t>(j)*nx+right];
                while(!q.empty() && raw[static_cast<size_t>(j)*nx+q.back()]<=v) q.pop_back();
                if(v>0.0) q.push_back(right);
            }
            const int left=std::max(0,i-radius);
            while(!q.empty()&&q.front()<left) q.pop_front();
            work[static_cast<size_t>(j)*nx+i]=q.empty()?0.0:raw[static_cast<size_t>(j)*nx+q.front()];
        }
    }
    for(int i=0;i<nx;++i) {
        std::deque<int> q; int bottom=-1;
        for(int j=0;j<ny;++j) {
            const int want_bottom=std::min(ny-1,j+radius);
            while(bottom<want_bottom) {
                ++bottom;
                const double v=work[static_cast<size_t>(bottom)*nx+i];
                while(!q.empty() && work[static_cast<size_t>(q.back())*nx+i]<=v) q.pop_back();
                if(v>0.0) q.push_back(bottom);
            }
            const int top=std::max(0,j-radius);
            while(!q.empty()&&q.front()<top) q.pop_front();
            out[static_cast<size_t>(j)*nx+i]=q.empty()?0.0:work[static_cast<size_t>(q.front())*nx+i];
        }
    }
}

void self_test_halo_math() {
    const int nx=5,ny=4,r=2;
    std::vector<double> raw(static_cast<size_t>(nx)*ny,0.0), h;
    raw[0]=0.25; raw[static_cast<size_t>(2)*nx+2]=0.75;
    raw[static_cast<size_t>(3)*nx+4]=0.5;
    std::vector<HaloSource> sources;
    std::vector<double> work;
    build_square_max(raw,nx,ny,r,h,work,&sources);
    if(sources.size()!=3) throw std::runtime_error("square-max fixture source count failed");
    for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) {
        double want=0.0;
        for(int sj=0;sj<ny;++sj) for(int si=0;si<nx;++si)
            if(std::abs(i-si)<=r&&std::abs(j-sj)<=r)
                want=std::max(want,raw[static_cast<size_t>(sj)*nx+si]);
        if(h[static_cast<size_t>(j)*nx+i]!=want)
            throw std::runtime_error("square-max differs from brute-force fixture");
    }
    raw.assign(raw.size(),0.0);
    build_square_max(raw,nx,ny,r,h,work,&sources);
    if(!sources.empty()||std::any_of(h.begin(),h.end(),[](double v){return v!=0.0;}))
        throw std::runtime_error("globally binary H=0 fixture failed");
    NodeAB<double,double> pair;
    pair.add_asym1(0,1,2.0);
    pair.add_asym1(1,2,-3.0);
    const auto before=pair.get_asym_pair(0,1);
    if(before.second!=0.0) throw std::runtime_error("asymmetric capillary value2 fixture failed");
    const double multiplier=1.0+27.0*0.125;
    for(auto it=pair.begin();it!=pair.end();++it) {
        if(it->value2!=0.0) throw std::runtime_error("fixture has unexpected value2 channel");
        it->value1*=multiplier;
    }
    const auto after=pair.get_asym_pair(0,1);
    if(after.first!=before.first*multiplier||after.second!=0.0||
       pair.get_asym1(1,0)!=-after.first||pair.size()!=2)
        throw std::runtime_error("antisymmetric value1-only scaling fixture failed");
}

void apply_halo(PhaseField& phase, const size_t ids[3], HaloModel& model) {
    HaloStep s{};
    if(model.chi<=0.0) { model.last=s; ++model.total.steps; return; }
    std::fill(model.raw.begin(),model.raw.end(),0.0);
    for(int j=0;j<NY;++j) for(int i=0;i<NX;++i) {
        const NodePF& node=phase.Fields(i,j,0);
        const double a=node.get_value(ids[0]), b=node.get_value(ids[1]), c=node.get_value(ids[2]);
        if(a>0.0&&b>0.0&&c>0.0) {
            const double value=27.0*a*b*c;
            if(!std::isfinite(value)||value<0.0) throw std::runtime_error("invalid ternary product");
            model.raw[static_cast<size_t>(j)*NX+i]=value;
            ++s.source_cells;
            s.max_source_product=std::max(s.max_source_product,value);
        }
    }
    build_square_max(model.raw,NX,NY,model.radius,model.h,model.work,&model.sources);
    for(double h:model.h) if(h>0.0) { ++s.halo_cells; s.max_H=std::max(s.max_H,h); }
    s.max_multiplier=1.0+model.chi*s.max_H;
    for(int j=0;j<NY;++j) for(int i=0;i<NX;++i) {
        const double h=model.h[static_cast<size_t>(j)*NX+i];
        NodeAB<double,double>& cap=phase.FieldsDot(i,j,0);
        if(cap.size()==0) continue;
        const double multiplier=(h>0.0)?1.0+model.chi*h:1.0;
        for(auto it=cap.begin();it!=cap.end();++it) {
            ++s.capillary_entries_total;
            if(it->value2!=0.0) ++s.nonzero_value2;
            if(it->value2!=0.0)
                throw std::runtime_error("native capillary insertion has nonzero value2; refusing to scale only value1");
            if(h>0.0) {
                ++s.capillary_entries_in_halo;
                const int slot=pair_slot(it->indexA,it->indexB,ids);
                if(slot>=0) ++s.pair_entries[static_cast<size_t>(slot)];
                it->value1*=multiplier;
                ++s.scaled_entries;
            }
        }
    }
    model.last=s;
    ++model.total.steps;
    model.total.source_cells_sum+=s.source_cells;
    model.total.halo_cells_sum+=s.halo_cells;
    model.total.capillary_entries_total+=s.capillary_entries_total;
    model.total.capillary_entries_in_halo+=s.capillary_entries_in_halo;
    model.total.scaled_entries+=s.scaled_entries;
    for(size_t p=0;p<3;++p) model.total.pair_entries[p]+=s.pair_entries[p];
    model.total.nonzero_value2+=s.nonzero_value2;
    model.total.max_source_cells=std::max(model.total.max_source_cells,s.source_cells);
    model.total.max_halo_cells=std::max(model.total.max_halo_cells,s.halo_cells);
    model.total.max_source_product=std::max(model.total.max_source_product,s.max_source_product);
    model.total.max_H=std::max(model.total.max_H,s.max_H);
    model.total.max_multiplier=std::max(model.total.max_multiplier,s.max_multiplier);
}

void snapshot_dots_before_normalize(PhaseField& phase, const size_t ids[3], double dt,
                                    HaloModel& halo) {
    halo.limiter_snapshot.clear();
    for(int j=0;j<NY;++j) for(int i=0;i<NX;++i) {
        const NodePF& node=phase.Fields(i,j,0);
        const NodeAB<double,double>& dot=phase.FieldsDot(i,j,0);
        if(dot.size()==0) continue;
        if(dot.size()>3) throw std::runtime_error("more than three pair entries at a three-phase cell");
        bool candidate=false;
        for(auto it=dot.cbegin();it!=dot.cend();++it)
            halo.limiter_snapshot.push_back({i,j,it->indexA,it->indexB,it->value1,it->value2});
        for(auto a=node.cbegin();a!=node.cend();++a) {
            double d=0.0;
            for(auto b=node.cbegin();b!=node.cend();++b) if(a!=b)
                d+=dot.get_asym1(a->index,b->index)+dot.get_asym2(a->index,b->index);
            const double after=a->value+dt*d;
            if(after<0.0||after>1.0) candidate=true;
        }
        if(candidate) ++halo.last.limiter_candidate_nodes;
    }
    halo.total.limiter_candidate_nodes+=halo.last.limiter_candidate_nodes;
}

void compare_dots_after_normalize(PhaseField& phase, HaloModel& halo) {
    int cell_i=-1,cell_j=-1; bool cell_changed=false;
    auto finish_cell=[&](){ if(cell_changed) ++halo.last.limiter_changed_nodes; };
    for(const auto& old:halo.limiter_snapshot) {
        if(old.i!=cell_i||old.j!=cell_j) {
            finish_cell(); cell_i=old.i; cell_j=old.j; cell_changed=false;
        }
        const auto now=phase.FieldsDot(old.i,old.j,0).get_asym_pair(old.a,old.b);
        if(now.first!=old.value1||now.second!=old.value2) {
            ++halo.last.limiter_changed_entries;
            cell_changed=true;
        }
        halo.last.max_limiter_delta=std::max(halo.last.max_limiter_delta,
            std::max(std::abs(now.first-old.value1),std::abs(now.second-old.value2)));
    }
    finish_cell();
    halo.total.limiter_changed_entries+=halo.last.limiter_changed_entries;
    halo.total.limiter_changed_nodes+=halo.last.limiter_changed_nodes;
    halo.total.max_limiter_delta=std::max(halo.total.max_limiter_delta,halo.last.max_limiter_delta);
}

void step_once(PhaseField& phase,InterfaceProperties& ip,DoubleObstacle& obstacle,
               BoundaryConditions& boundary,DrivingForce& drive,InterfaceRegularization& reg,
               const size_t ids[3],bool regularized,double dt,HaloModel& halo,bool track_limiter) {
    ip.Set(phase,boundary);
    if(regularized) { obstacle.CalculatePhaseFieldIncrements(phase,ip,reg); reg.Average(phase,boundary); reg.MergePhaseFieldIncrements(phase,ip); }
    else obstacle.CalculatePhaseFieldIncrements(phase,ip);
    apply_halo(phase,ids,halo);
    drive.Clear();
    for(int j=0;j<NY;++j) for(int i=0;i<NX;++i) {
        const NodePF& n=phase.Fields(i,j,0); if(!n.wide_interface()) continue;
        const double u0=n.get_value(ids[0]),u1=n.get_value(ids[1]),u2=n.get_value(ids[2]);
        if(u0>0.0&&u1>0.0) drive.Force(i,j,0).set_all(ids[0],ids[1],30.0,0.0,0.0,std::sqrt(u0*u1));
        if(u0>0.0&&u2>0.0) drive.Force(i,j,0).set_all(ids[0],ids[2],30.0,0.0,0.0,std::sqrt(u0*u2));
    }
    drive.SetBoundaryConditions(boundary); drive.MergePhaseFieldIncrements(phase,ip);
    if(track_limiter) snapshot_dots_before_normalize(phase,ids,dt,halo);
    phase.NormalizeIncrements(boundary,dt);
    if(track_limiter) compare_dots_after_normalize(phase,halo);
    phase.MergeIncrements(boundary,dt);
}

struct Schedule { std::set<int> scalar; std::map<int,std::vector<std::string>> fields; };
Schedule read_schedule(const std::string& path) {
    std::ifstream in(path); if(!in) throw std::runtime_error("missing frozen schedule file");
    Schedule s; std::string line; std::getline(in,line);
    while(std::getline(in,line)) { if(line.empty()) continue; std::istringstream row(line);
        char kind; int step; std::string tag; row>>kind>>step>>tag;
        if(!row || step<0) throw std::runtime_error("invalid schedule row");
        if(kind=='S') s.scalar.insert(step); else if(kind=='F') s.fields[step].push_back(tag); else throw std::runtime_error("bad schedule kind");
    }
    if(s.scalar.empty()||*s.scalar.begin()!=0) throw std::runtime_error("scalar schedule must start at zero");
    return s;
}

void write_field_header(std::ofstream& out) {
    const char magic[8]={'O','P','H','2','D','F','1','\0'}; const uint32_t nx=NX,ny=NY,np=3,count=0; const double h=H;
    out.write(magic,8); out.write(reinterpret_cast<const char*>(&nx),4); out.write(reinterpret_cast<const char*>(&ny),4);
    out.write(reinterpret_cast<const char*>(&np),4); out.write(reinterpret_cast<const char*>(&count),4); out.write(reinterpret_cast<const char*>(&h),8);
}
void finish_header(std::ofstream& out,uint32_t count) { out.seekp(20,std::ios::beg); out.write(reinterpret_cast<const char*>(&count),4); out.seekp(0,std::ios::end); out.flush(); if(!out) throw std::runtime_error("field header update failed"); }
void write_observation(std::ofstream& out,const Sample& s) {
    if(!std::isfinite(s.raw) || !std::isfinite(s.physical) || !std::isfinite(s.wall))
        throw std::runtime_error("nonfinite scalar row before write");
    out<<std::setprecision(17)<<s.step<<'\t'<<s.raw<<'\t'<<s.physical<<'\t'<<s.wall<<'\t'<<s.min_phi<<'\t'<<s.max_phi<<'\t'<<s.partition_linf<<'\t'<<s.x_swap_linf<<'\t'<<s.negative_f_outside;
    for(int p=0;p<3;++p) out<<'\t'<<s.phase_volume[static_cast<size_t>(p)]
        <<'\t'<<s.phase_volume_ratio[static_cast<size_t>(p)]
        <<'\t'<<s.phase_exists[static_cast<size_t>(p)]
        <<'\t'<<s.phase_entries[static_cast<size_t>(p)];
    for(const auto& p:s.pair) {
        out<<'\t'<<p.area<<'\t'<<p.central_f<<'\t'<<p.diffuse_integral<<'\t'<<p.edge_min
           <<'\t'<<p.roots.size()<<'\t'<<'[';
        for(size_t i=0;i<p.roots.size();++i) { if(i) out<<','; out<<p.roots[i]; }
        out<<']';
    }
    out<<'\n'; out.flush(); if(!out) throw std::runtime_error("observation write failed");
}

class FieldWriter {
public:
    FieldWriter(const std::string& bin,const std::string& tagfile):bin_(bin,std::ios::binary|std::ios::trunc),tags_(tagfile,std::ios::trunc) {
        if(!bin_||!tags_) throw std::runtime_error("cannot create field output"); write_field_header(bin_);
        tags_<<"step\ttag\n";
    }
    void add(int step,double raw,const std::string& tag,const std::vector<double>& values) {
        if(values.size()!=static_cast<size_t>(3)*NY*NX) throw std::runtime_error("field snapshot size mismatch");
        if(!std::isfinite(raw) || !std::all_of(values.begin(),values.end(),[](double value){return std::isfinite(value);}))
            throw std::runtime_error("nonfinite field snapshot/time before archive write");
        auto& seen=labels_[step]; if(std::find(seen.begin(),seen.end(),tag)==seen.end()) { seen.push_back(tag); tags_<<step<<'\t'<<tag<<'\n'; tags_.flush(); }
        if(saved_.count(step)) return;
        const uint64_t st=static_cast<uint64_t>(step); bin_.write(reinterpret_cast<const char*>(&st),8); bin_.write(reinterpret_cast<const char*>(&raw),8);
        bin_.write(reinterpret_cast<const char*>(values.data()),static_cast<std::streamsize>(values.size()*sizeof(double)));
        bin_.flush(); if(!bin_) throw std::runtime_error("field snapshot write failed");
        saved_.insert(step); ++count_; finish_header(bin_,count_);
    }
    void close() { finish_header(bin_,count_); bin_.close(); tags_.close(); }
private:
    std::ofstream bin_,tags_; uint32_t count_{}; std::set<int> saved_;
    std::map<int,std::vector<std::string>> labels_;
};

void write_scalar_header(std::ofstream& o) {
    o<<"step\traw_time_W2\tphysical_time_W2\twall_from_start_s\tfield_min\tfield_max\tpartition_linf\tx_phase_swap_linf\tnegative_f_outside";
    for(int p=0;p<3;++p) o<<"\tphase"<<p<<"_volume\tphase"<<p<<"_volume_ratio\tphase"<<p<<"_exists\tphase"<<p<<"_entry_count";
    for(const char* name:PAIR_NAMES) o<<'\t'<<name<<"_area_W2\t"<<name<<"_central_F\t"<<name<<"_diffuse_integral_W2\t"<<name<<"_edge_F_min\t"<<name<<"_root_count\t"<<name<<"_roots_y_W_json";
    o<<'\n';
}

void write_halo_header(std::ofstream& o) {
    o << "step\tchi\tradius\tsource_cells\thalo_cells\tcapillary_entries_total\tcapillary_entries_in_halo"
         "\tscaled_entries\tpair01\tpair02\tpair12\tnonzero_value2"
         "\tmax_source_product\tmax_H\tmax_multiplier\tlimiter_candidate_nodes"
         "\tlimiter_changed_nodes\tlimiter_changed_entries\tmax_limiter_delta\n";
}
void write_halo_row(std::ofstream& o,int step,const HaloModel& h) {
    const auto& s=h.last;
    o<<std::setprecision(17)<<step<<'\t'<<h.chi<<'\t'<<h.radius<<'\t'
     <<s.source_cells<<'\t'<<s.halo_cells<<'\t'<<s.capillary_entries_total<<'\t'<<s.capillary_entries_in_halo<<'\t'
     <<s.scaled_entries<<'\t'<<s.pair_entries[0]<<'\t'<<s.pair_entries[1]<<'\t'
     <<s.pair_entries[2]<<'\t'<<s.nonzero_value2<<'\t'<<s.max_source_product<<'\t'
     <<s.max_H<<'\t'<<s.max_multiplier<<'\t'<<s.limiter_candidate_nodes<<'\t'
     <<s.limiter_changed_nodes<<'\t'<<s.limiter_changed_entries<<'\t'<<s.max_limiter_delta<<'\n';
    if(!o) throw std::runtime_error("halo metrics write failed");
}

void run_full(const std::string& route,const std::string& opi_path,const std::string& input_path,
              const std::string& summary_path,const std::string& fields_path,const std::string& tags_path,
              const std::string& schedule_path,double factor,double t99_short,double t99_long,double max_seconds,
              bool preflight_only,double chi,int radius,bool track_limiter,const std::string& halo_path) {
    const bool regularized=route=="regularized"; if(!regularized&&route!="native") throw std::runtime_error("bad route");
    const auto begin=std::chrono::steady_clock::now(); const std::clock_t cpu_begin=std::clock();
    auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();};
    Settings settings(opi_path); RunTimeControl rtc(settings,opi_path); PhaseField phase(settings,opi_path);
    DoubleObstacle obstacle(settings,opi_path); InterfaceProperties ip(settings,opi_path); BoundaryConditions boundary(settings,opi_path);
    DrivingForce drive(settings,opi_path); InterfaceRegularization reg(settings,opi_path);
    const size_t ids[3]={phase.AddGrainInfo(0),phase.AddGrainInfo(1),phase.AddGrainInfo(2)};
    if(phase.Grid.Nx!=NX||phase.Grid.Ny!=NY||phase.Grid.Active()!=2||std::abs(phase.Grid.dx-H)>1e-14||std::abs(phase.Grid.Eta-0.08)>1e-14)
        throw std::runtime_error("H24 OPI dimensions/eta differ");
    const double dt=rtc.dt; if(!std::isfinite(dt)||dt<=0.0) throw std::runtime_error("raw dt must be positive");
    if(!std::isfinite(chi)||chi<0.0||radius<0||radius>std::max(NX,NY))
        throw std::runtime_error("invalid halo parameters");
    if(!std::isfinite(factor) || factor<=0.0) throw std::runtime_error("physical clock factor must be positive");
    load_2d(phase,ids,input_path); phase.FinalizeInitialization(boundary); set_initial_stable_volumes(phase,ids); phase.Finalize(boundary);
    const double init_wall=elapsed();
    HaloModel halo; halo.chi=chi; halo.radius=radius;
    halo.raw.assign(static_cast<size_t>(NX)*NY,0.0);
    halo.h.assign(static_cast<size_t>(NX)*NY,0.0);
    halo.work.assign(static_cast<size_t>(NX)*NY,0.0);
    const Schedule sched=read_schedule(schedule_path); const int terminal=*sched.scalar.rbegin();
    if(preflight_only) {
        if(terminal!=0) throw std::runtime_error("preflight schedule must contain only step zero");
    } else if(terminal!=rtc.nSteps) throw std::runtime_error("OPI nSteps differs from nearest-step terminal");
    std::ofstream summary(summary_path,std::ios::trunc); if(!summary) throw std::runtime_error("cannot create scalar output"); write_scalar_header(summary);
    std::ofstream haloout(halo_path,std::ios::trunc); if(!haloout) throw std::runtime_error("cannot create halo output"); write_halo_header(haloout);
    FieldWriter fields(fields_path,tags_path);
    double stepper_s=0.0,observer_s=0.0,export_s=0.0; uint64_t observer_count=0;
    std::array<double,2> previous_area{}, previous_central{}; bool have_previous=false;
    std::array<std::array<bool,3>,2> event_found{};
    const std::array<double,2> targets{t99_short,t99_long};
    std::vector<double> previous_phi;
    int previous_step=0;
    auto save_fixed=[&](int step,double raw,const std::vector<double>& values) {
        auto it=sched.fields.find(step); if(it==sched.fields.end()) return;
        const auto e=std::chrono::steady_clock::now(); for(const auto& tag:it->second) fields.add(step,raw,tag,values);
        export_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-e).count();
    };
    auto process_scalar=[&](int step) {
        const auto ob0=std::chrono::steady_clock::now();
        Sample s=inspect(phase,ids,step,dt,factor,elapsed());
        if(s.min_phi < -1e-12 || s.max_phi > 1.0+1e-12 || s.partition_linf>2e-12 || s.x_swap_linf>2e-11)
            throw std::runtime_error("saved observer invariant failed");
        std::vector<double> current=capture(phase,ids); const double raw=step*dt;
        if(step==0) fields.add(step,raw,"initial",current);
        for(int p=0;p<2;++p) {
            const auto& m=s.pair[static_cast<size_t>(p)];
            if(m.edge_min<=0.0) throw std::runtime_error("ROI edge contains negative F");
            if(have_previous) {
                auto event=[&](int e,bool hit,const char* name) {
                    if(event_found[p][e]||!hit) return;
                    const auto io0=std::chrono::steady_clock::now();
                    fields.add(previous_step, previous_step*dt,
                               std::string("event_")+PAIR_NAMES[p]+"_"+name+"_before",previous_phi);
                    fields.add(step,raw,std::string("event_")+PAIR_NAMES[p]+"_"+name+"_after",current);
                    export_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-io0).count();
                    event_found[p][e]=true;
                };
                event(0,previous_central[p]<0.0 && m.central_f>=0.0,"central_contact");
                event(1,previous_area[p]>targets[p] && m.area<=targets[p],"geometric_T99");
                event(2,previous_area[p]>0.0 && m.area<=0.0,"geometric_zero_area");
            }
        previous_area[p]=m.area; previous_central[p]=m.central_f;
        }
        auto it=sched.fields.find(step); if(it!=sched.fields.end()) for(const auto& tag:it->second) fields.add(step,raw,tag,current);
        write_observation(summary,s); previous_phi=std::move(current); previous_step=step; have_previous=true;
        observer_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-ob0).count(); ++observer_count;
    };
    process_scalar(0);
    if(preflight_only) { fields.close(); summary.close(); haloout.close(); return; }
    for(int step=1;step<=terminal;++step) {
        if((step&31)==1 && (stop_requested || elapsed()>=max_seconds)) {
            const auto t0=std::chrono::steady_clock::now(); Sample partial=inspect(phase,ids,step-1,dt,factor,elapsed());
            auto current=capture(phase,ids); fields.add(step-1,(step-1)*dt,"partial_stop",current); write_observation(summary,partial);
            fields.close(); summary.close(); haloout.close(); std::ofstream status(summary_path+".status",std::ios::trunc);
            status<<"PARTIAL\tstep="<<(step-1)<<"\treason="<<(stop_requested?"external_guard":"soft_deadline")<<"\n"; status.close(); return;
        }
        const auto k0=std::chrono::steady_clock::now(); step_once(phase,ip,obstacle,boundary,drive,reg,ids,regularized,dt,halo,track_limiter);
        stepper_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-k0).count();
        write_halo_row(haloout,step,halo);
        if(sched.scalar.count(step)) process_scalar(step);
        else { auto it=sched.fields.find(step); if(it!=sched.fields.end()) { auto x=capture(phase,ids); save_fixed(step,step*dt,x); } }
    }
    fields.close(); summary.close(); haloout.close(); if(!summary||!haloout) throw std::runtime_error("scalar/halo finalization failed");
    struct rusage usage{}; if(getrusage(RUSAGE_SELF,&usage)!=0) throw std::runtime_error("getrusage failed");
    std::ofstream timing(summary_path+".timing",std::ios::trunc);
    timing<<std::setprecision(17)<<"schema\topenphase-full-route-timing-v001\nroute\t"<<route<<"\n"
          <<"updates\t"<<terminal<<"\nraw_dt_W2\t"<<dt<<"\nphysical_factor_a\t"<<factor<<"\n"
          <<"initialization_wall_s\t"<<init_wall<<"\nstepper_wall_s\t"<<stepper_s<<"\nobserver_inclusive_wall_s\t"<<observer_s<<"\nfield_export_subset_wall_s\t"<<export_s<<"\n"
          <<"observer_rows\t"<<observer_count<<"\nroute_wall_s\t"<<elapsed()<<"\nprocess_cpu_s\t"<<double(std::clock()-cpu_begin)/CLOCKS_PER_SEC<<"\npeak_rss_bytes\t"<<usage.ru_maxrss<<"\n"
          <<"halo_chi\t"<<chi<<"\nhalo_radius\t"<<radius<<"\nhalo_steps\t"<<halo.total.steps<<"\n"
          <<"halo_source_cells_sum\t"<<halo.total.source_cells_sum<<"\nhalo_cells_sum\t"<<halo.total.halo_cells_sum<<"\n"
          <<"halo_capillary_entries_total\t"<<halo.total.capillary_entries_total<<"\n"
          <<"halo_capillary_entries_in_halo\t"<<halo.total.capillary_entries_in_halo<<"\nhalo_scaled_entries\t"<<halo.total.scaled_entries<<"\n"
          <<"halo_pair01_entries\t"<<halo.total.pair_entries[0]<<"\nhalo_pair02_entries\t"<<halo.total.pair_entries[1]<<"\nhalo_pair12_entries\t"<<halo.total.pair_entries[2]<<"\n"
          <<"halo_nonzero_value2\t"<<halo.total.nonzero_value2<<"\n"
          <<"halo_limiter_candidate_nodes\t"<<halo.total.limiter_candidate_nodes<<"\nhalo_limiter_changed_nodes\t"<<halo.total.limiter_changed_nodes<<"\nhalo_limiter_changed_entries\t"<<halo.total.limiter_changed_entries<<"\n"
          <<"halo_max_limiter_delta\t"<<halo.total.max_limiter_delta<<"\n"
          <<"halo_max_source_product\t"<<halo.total.max_source_product<<"\nhalo_max_H\t"<<halo.total.max_H<<"\nhalo_max_multiplier\t"<<halo.total.max_multiplier<<"\n";
    timing.close();
}
} // namespace

int main(int argc,char** argv) {
    std::signal(SIGUSR1,request_stop);
    if(argc!=17) { std::cerr<<"usage: openphase_h24_full route opi input.raw summary.tsv fields.bin tags.tsv schedule.tsv factor t99short t99long max_seconds preflight_only chi radius track_limiter halo.tsv\n"; return 2; }
    try {
        self_test_halo_math();
        run_full(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],number(argv[8],"factor"),
            number(argv[9],"short T99"),number(argv[10],"long T99"),number(argv[11],"maximum seconds"),
            integer(argv[12],"preflight flag")!=0,number(argv[13],"halo chi"),integer(argv[14],"halo radius"),
            integer(argv[15],"track limiter flag")!=0,argv[16]);
    }
    catch(const std::exception& e) { std::cerr<<"OPENPHASE_H24_FULL_FAILURE\t"<<e.what()<<'\n'; return 1; }
    return 0;
}
