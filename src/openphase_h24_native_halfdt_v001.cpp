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
void step_once(PhaseField& phase,InterfaceProperties& ip,DoubleObstacle& obstacle,
               BoundaryConditions& boundary,DrivingForce& drive,InterfaceRegularization& reg,
               const size_t ids[3],bool regularized,double dt) {
    ip.Set(phase,boundary);
    if(regularized) { obstacle.CalculatePhaseFieldIncrements(phase,ip,reg); reg.Average(phase,boundary); reg.MergePhaseFieldIncrements(phase,ip); }
    else obstacle.CalculatePhaseFieldIncrements(phase,ip);
    drive.Clear();
    for(int j=0;j<NY;++j) for(int i=0;i<NX;++i) {
        const NodePF& n=phase.Fields(i,j,0); if(!n.wide_interface()) continue;
        const double u0=n.get_value(ids[0]),u1=n.get_value(ids[1]),u2=n.get_value(ids[2]);
        if(u0>0.0&&u1>0.0) drive.Force(i,j,0).set_all(ids[0],ids[1],30.0,0.0,0.0,std::sqrt(u0*u1));
        if(u0>0.0&&u2>0.0) drive.Force(i,j,0).set_all(ids[0],ids[2],30.0,0.0,0.0,std::sqrt(u0*u2));
    }
    drive.SetBoundaryConditions(boundary); drive.MergePhaseFieldIncrements(phase,ip);
    phase.NormalizeIncrements(boundary,dt); phase.MergeIncrements(boundary,dt);
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

void run_full(const std::string& route,const std::string& opi_path,const std::string& input_path,
              const std::string& summary_path,const std::string& fields_path,const std::string& tags_path,
              const std::string& schedule_path,double factor,double t99_short,double t99_long,double max_seconds,
              bool preflight_only) {
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
    if(!std::isfinite(factor) || factor<=0.0) throw std::runtime_error("physical clock factor must be positive");
    load_2d(phase,ids,input_path); phase.FinalizeInitialization(boundary); set_initial_stable_volumes(phase,ids); phase.Finalize(boundary);
    const double init_wall=elapsed();
    const Schedule sched=read_schedule(schedule_path); const int terminal=*sched.scalar.rbegin();
    if(preflight_only) {
        if(terminal!=0) throw std::runtime_error("preflight schedule must contain only step zero");
    } else if(terminal!=rtc.nSteps) throw std::runtime_error("OPI nSteps differs from nearest-step terminal");
    std::ofstream summary(summary_path,std::ios::trunc); if(!summary) throw std::runtime_error("cannot create scalar output"); write_scalar_header(summary);
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
    if(preflight_only) { fields.close(); summary.close(); return; }
    for(int step=1;step<=terminal;++step) {
        if((step&31)==1 && (stop_requested || elapsed()>=max_seconds)) {
            const auto t0=std::chrono::steady_clock::now(); Sample partial=inspect(phase,ids,step-1,dt,factor,elapsed());
            auto current=capture(phase,ids); fields.add(step-1,(step-1)*dt,"partial_stop",current); write_observation(summary,partial);
            fields.close(); summary.close(); std::ofstream status(summary_path+".status",std::ios::trunc);
            status<<"PARTIAL\tstep="<<(step-1)<<"\treason="<<(stop_requested?"external_guard":"soft_deadline")<<"\n"; status.close(); return;
        }
        const auto k0=std::chrono::steady_clock::now(); step_once(phase,ip,obstacle,boundary,drive,reg,ids,regularized,dt);
        stepper_s+=std::chrono::duration<double>(std::chrono::steady_clock::now()-k0).count();
        if(sched.scalar.count(step)) process_scalar(step);
        else { auto it=sched.fields.find(step); if(it!=sched.fields.end()) { auto x=capture(phase,ids); save_fixed(step,step*dt,x); } }
    }
    fields.close(); summary.close(); if(!summary) throw std::runtime_error("scalar finalization failed");
    struct rusage usage{}; if(getrusage(RUSAGE_SELF,&usage)!=0) throw std::runtime_error("getrusage failed");
    std::ofstream timing(summary_path+".timing",std::ios::trunc);
    timing<<std::setprecision(17)<<"schema\topenphase-full-route-timing-v001\nroute\t"<<route<<"\n"
          <<"updates\t"<<terminal<<"\nraw_dt_W2\t"<<dt<<"\nphysical_factor_a\t"<<factor<<"\n"
          <<"initialization_wall_s\t"<<init_wall<<"\nstepper_wall_s\t"<<stepper_s<<"\nobserver_inclusive_wall_s\t"<<observer_s<<"\nfield_export_subset_wall_s\t"<<export_s<<"\n"
          <<"observer_rows\t"<<observer_count<<"\nroute_wall_s\t"<<elapsed()<<"\nprocess_cpu_s\t"<<double(std::clock()-cpu_begin)/CLOCKS_PER_SEC<<"\npeak_rss_bytes\t"<<usage.ru_maxrss<<"\n";
    timing.close();
}
} // namespace

int main(int argc,char** argv) {
    std::signal(SIGUSR1,request_stop);
    if(argc!=13) { std::cerr<<"usage: openphase_h24_full route opi input.raw summary.tsv fields.bin tags.tsv schedule.tsv factor t99short t99long max_seconds preflight_only\n"; return 2; }
    try { run_full(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],number(argv[8],"factor"),number(argv[9],"short T99"),number(argv[10],"long T99"),number(argv[11],"maximum seconds"),integer(argv[12],"preflight flag")!=0); }
    catch(const std::exception& e) { std::cerr<<"OPENPHASE_H24_FULL_FAILURE\t"<<e.what()<<'\n'; return 1; }
    return 0;
}
