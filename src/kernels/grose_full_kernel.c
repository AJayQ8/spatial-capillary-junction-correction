#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* Full-grid, source-based Grose--Asimow path for three order parameters and
   two thermodynamic phase fields (alpha=phi0, beta=phi1+phi2).

   The capillary and bulk force uses the direct variational sums in manuscript
   Eqs. 17--18 and 32--34. Eq. 43 is applied to coexisting pairs local to the
   five-point support. Each node performs the explicit Appendix-H local event
   recursion with frozen macro-step neighbors. The phase-profile correction is
   a separate function, applied after the force update using a shared
   pre-correction snapshot, Appendix G normal extension, and Eqs. 37--42.

   This source is a new implementation; it is not the authors' code. */

static inline size_t ix(int y, int x, int nx) {
    return (size_t)y * (size_t)nx + (size_t)x;
}

static inline int clampi(int v, int n) {
    return v < 0 ? 0 : (v >= n ? n - 1 : v);
}

static inline double fdc2(double h, double eta) {
    const double z = M_PI * h / eta;
    if (fabs(z) < 1.0e-4) return 1.0 + z*z/12.0 + z*z*z*z/240.0;
    return z*z / (2.0 * (1.0 - cos(z)));
}

static void local_laplacian(const double *phi, int ny, int nx, int y, int x,
                            double h, double fdc, double lap[3],
                            double nb[4][3], double center[3]) {
    const size_t n = (size_t)ny * (size_t)nx;
    const int ym = y > 0 ? y - 1 : y;
    const int yp = y + 1 < ny ? y + 1 : y;
    const int xm = x > 0 ? x - 1 : x;
    const int xp = x + 1 < nx ? x + 1 : x;
    const int ys[4] = {ym, yp, y, y};
    const int xs[4] = {x, x, xm, xp};
    const double ih2 = 1.0 / (h * h);
    for (int p = 0; p < 3; ++p) {
        const size_t base = (size_t)p * n;
        center[p] = phi[base + ix(y, x, nx)];
        double sum = 0.0;
        for (int d = 0; d < 4; ++d) {
            nb[d][p] = phi[base + ix(ys[d], xs[d], nx)];
            sum += nb[d][p];
        }
        lap[p] = fdc * (sum - 4.0 * center[p]) * ih2;
    }
}

static inline void local_laplacian_from_frozen_neighbors(
    const double center[3], const double nb[4][3], double h, double fdc,
    double lap[3]) {
    const double ih2 = 1.0 / (h * h);
    for (int p = 0; p < 3; ++p)
        lap[p] = fdc * (nb[0][p] + nb[1][p] + nb[2][p] + nb[3][p]
                        - 4.0 * center[p]) * ih2;
}

/* Direct local equations. eligible is the current recursion participation
   mask after applying the positive five-point-neighborhood rule and retaining
   the per-node eliminated set. */
static void direct_rhs(const double p[3], const double lap[3],
                       const int eligible[3], const double sigma[9],
                       const double mobility[9], const double mu0[9],
                       double h, double eta, double amplitude, double rhs[3],
                       double trace[32]) {
    double mu_max = 0.0, lambda_max = 0.0;
    int nu = eligible[0] + eligible[1] + eligible[2];
    int pairs = 0;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            const int ij = 3*i + j;
            if (!eligible[i] || !eligible[j]) continue;
            ++pairs;
            const double a = fabs(mu0[ij]);
            if (a > mu_max) mu_max = a;
            if (sigma[ij] > 0.0) {
                const double l = h * a / sigma[ij];
                if (l > lambda_max) lambda_max = l;
            }
        }
    }
    double sa[9] = {0.0};
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            const int ij = 3*i + j, ji = 3*j + i;
            if (!eligible[i] || !eligible[j]) continue;
            const double ratio = mu_max > 0.0 ? fabs(mu0[ij]) / mu_max : 0.0;
            const double mult = mu_max > 0.0
                ? 1.0 + amplitude * lambda_max * (1.0 - ratio) : 1.0;
            sa[ij] = sa[ji] = sigma[ij] * mult;
        }
    }
    double weighted_mu = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            if (eligible[i] && eligible[j]) {
                const double q = p[i] * p[j] * mu0[3*i+j] * mu0[3*i+j];
                weighted_mu += q;
            }
    const double k2 = (M_PI / eta) * (M_PI / eta);
    double rpair[9] = {0.0};
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            const int ij = 3*i + j, ji = 3*j + i;
            if (!eligible[i] || !eligible[j]) continue;
            double mix_i = 0.0, mix_j = 0.0, grad_i = 0.0, grad_j = 0.0;
            for (int k = 0; k < 3; ++k) {
                if (k != i) { mix_i += sa[3*i+k] * p[k]; grad_i += sa[3*i+k] * lap[k]; }
                if (k != j) { mix_j += sa[3*j+k] * p[k]; grad_j += sa[3*j+k] * lap[k]; }
            }
            const double cap = (k2 * (mix_i - mix_j) + grad_i - grad_j) / (double)nu;
            const double gij = weighted_mu > 0.0
                ? (M_PI / eta) * p[i] * p[j] * fabs(mu0[ij]) / sqrt(weighted_mu)
                : 0.0;
            const double bulk = gij * mu0[ij];
            const double transfer = -mobility[ij] * (cap + bulk);
            rpair[ij] = transfer;
            rpair[ji] = -transfer;
            if (trace) {
                const int t = (i*3+j)*3;
                trace[t] = cap;
                trace[t+1] = bulk;
                trace[t+2] = transfer;
            }
        }
    }
    rhs[0] = rhs[1] = rhs[2] = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) rhs[i] += rpair[3*i+j];
    if (trace) {
        trace[27] = (double)nu;
        trace[28] = (double)pairs;
        trace[29] = mu_max;
        trace[30] = lambda_max;
        trace[31] = weighted_mu;
    }
}

int grose_rhs_field_v001(const double *phi, double *rhs, double *lap_out,
                         uint8_t *eligible_out, int ny, int nx,
                         double h, double eta, double amplitude,
                         const double *sigma, const double *mobility,
                         const double *mu0) {
    if (!phi || !rhs || !lap_out || !eligible_out || !sigma || !mobility || !mu0 ||
        ny < 2 || nx < 2 || !(h > 0.0) || !(eta > 0.0) || amplitude < 0.0) return 1;
    const size_t n = (size_t)ny * (size_t)nx;
    const double fdc = fdc2(h, eta);
    for (int y = 0; y < ny; ++y) for (int x = 0; x < nx; ++x) {
        double p[3], nb[4][3], lap[3], localrhs[3];
        local_laplacian(phi, ny, nx, y, x, h, fdc, lap, nb, p);
        int eligible[3];
        for (int i = 0; i < 3; ++i) {
            eligible[i] = p[i] > 0.0 || nb[0][i] > 0.0 || nb[1][i] > 0.0 ||
                          nb[2][i] > 0.0 || nb[3][i] > 0.0;
            eligible_out[(size_t)i*n + ix(y,x,nx)] = (uint8_t)eligible[i];
            lap_out[(size_t)i*n + ix(y,x,nx)] = lap[i];
        }
        direct_rhs(p, lap, eligible, sigma, mobility, mu0, h, eta, amplitude, localrhs, NULL);
        for (int i = 0; i < 3; ++i) rhs[(size_t)i*n + ix(y,x,nx)] = localrhs[i];
    }
    return 0;
}

int grose_step_v001(const double *phi, double *out, int ny, int nx,
                    double h, double eta, double dt, double amplitude,
                    const double *sigma, const double *mobility,
                    const double *mu0, double *stats) {
    if (!phi || !out || !sigma || !mobility || !mu0 || ny < 2 || nx < 2 ||
        !(h > 0.0) || !(eta > 0.0) || dt < 0.0 || amplitude < 0.0) return 1;
    const size_t n = (size_t)ny * (size_t)nx;
    if (dt == 0.0) { memcpy(out, phi, (size_t)3*n*sizeof(double)); if (stats) memset(stats,0,8*sizeof(double)); return 0; }
    if (stats) memset(stats, 0, 8*sizeof(double));
    const double fdc = fdc2(h, eta);
    for (int y = 0; y < ny; ++y) for (int x = 0; x < nx; ++x) {
        double p[3], nb[4][3], lap[3];
        local_laplacian(phi, ny, nx, y, x, h, fdc, lap, nb, p);
        int eliminated[3] = {0,0,0};
        double rem = dt, elapsed = 0.0;
        int substeps = 0;
        for (; substeps < 4 && rem > 0.0; ++substeps) {
            int eligible[3], full[3];
            for (int i = 0; i < 3; ++i) {
                full[i] = p[i] > 0.0 || nb[0][i] > 0.0 || nb[1][i] > 0.0 ||
                          nb[2][i] > 0.0 || nb[3][i] > 0.0;
                eligible[i] = full[i] && !eliminated[i];
            }
            const int neligible = eligible[0] + eligible[1] + eligible[2];
            if (neligible == 0) return 2;
            /* Appendix-H substeps update the center but retain macro-step
               neighbors; therefore the center-dependent stencil is rebuilt. */
            local_laplacian_from_frozen_neighbors(p, nb, h, fdc, lap);
            double rhs[3];
            direct_rhs(p, lap, eligible, sigma, mobility, mu0, h, eta, amplitude, rhs, NULL);
            double tmin = INFINITY;
            double tcross[3] = {INFINITY,INFINITY,INFINITY};
            for (int i = 0; i < 3; ++i) if (eligible[i] && rhs[i] < 0.0) {
                tcross[i] = p[i] == 0.0 ? 0.0 : p[i] / (-rhs[i]);
                if (tcross[i] < tmin) tmin = tcross[i];
            }
            const double scale = fmax(dt, fmax(rem, isfinite(tmin) ? fabs(tmin) : 0.0));
            const double tol = 16.0 * DBL_EPSILON * scale;
            if (!isfinite(tmin) || tmin > rem + tol) {
                for (int i = 0; i < 3; ++i) p[i] += rem * rhs[i];
                elapsed += rem; rem = 0.0;
            } else {
                double tau = fmin(rem, fmax(0.0, tmin));
                for (int i = 0; i < 3; ++i) p[i] += tau * rhs[i];
                int removed = 0;
                for (int i = 0; i < 3; ++i) if (fabs(tcross[i] - tmin) <= tol) {
                    p[i] = 0.0; eliminated[i] = 1; removed++;
                    if (stats) stats[0] += 1.0;
                }
                if (!removed) return 3;
                if (stats) stats[1] += 1.0;
                elapsed += tau;
                rem = fmax(0.0, rem - tau);
                double sum = p[0] + p[1] + p[2];
                const double resid = 1.0 - sum;
                int survivor = -1;
                for (int i = 0; i < 3; ++i) if (!eliminated[i] &&
                    (survivor < 0 || p[i] > p[survivor])) survivor = i;
                if (survivor < 0) return 4;
                p[survivor] += resid;
                if (stats && fabs(resid) > stats[4]) stats[4] = fabs(resid);
                if (rem <= tol) rem = 0.0;
            }
            if (stats && substeps + 1 > stats[2]) stats[2] = (double)(substeps + 1);
        }
        if (rem > 0.0) return 5;
        double sum = p[0] + p[1] + p[2];
        double minp = fmin(p[0], fmin(p[1], p[2]));
        if (minp < -2e-12 || !isfinite(sum)) return 6;
        if (minp < 0.0) {
            const double cleanup = -minp;
            for (int i = 0; i < 3; ++i) if (p[i] < 0.0) p[i] = 0.0;
            const double residual = 1.0 - (p[0]+p[1]+p[2]);
            int largest = p[1] > p[0] ? 1 : 0; if (p[2] > p[largest]) largest = 2;
            p[largest] += residual;
            if (stats && cleanup > stats[3]) stats[3] = cleanup;
        } else {
            const double residual = 1.0 - sum;
            int largest = p[1] > p[0] ? 1 : 0; if (p[2] > p[largest]) largest = 2;
            p[largest] += residual;
        }
        if (stats) {
            const double e = fabs((p[0]+p[1]+p[2]) - 1.0);
            if (e > stats[5]) stats[5] = e;
        }
        for (int i = 0; i < 3; ++i) out[(size_t)i*n + ix(y,x,nx)] = p[i];
        (void)elapsed;
    }
    return 0;
}

static inline double bilinear(const double *a, int ny, int nx,
                              double y, double x) {
    if (x < 0.0) x = 0.0; if (x > (double)(nx-1)) x = (double)(nx-1);
    if (y < 0.0) y = 0.0; if (y > (double)(ny-1)) y = (double)(ny-1);
    int x0 = (int)floor(x), y0 = (int)floor(y);
    int x1 = x0 + 1 < nx ? x0 + 1 : x0;
    int y1 = y0 + 1 < ny ? y0 + 1 : y0;
    double fx = x - x0, fy = y - y0;
    return (1.0-fy)*((1.0-fx)*a[ix(y0,x0,nx)] + fx*a[ix(y0,x1,nx)]) +
           fy*((1.0-fx)*a[ix(y1,x0,nx)] + fx*a[ix(y1,x1,nx)]);
}

static inline void bilinear_add(double *a, int ny, int nx,
                                double y, double x, double value) {
    if (x < 0.0) x = 0.0; if (x > (double)(nx-1)) x = (double)(nx-1);
    if (y < 0.0) y = 0.0; if (y > (double)(ny-1)) y = (double)(ny-1);
    int x0 = (int)floor(x), y0 = (int)floor(y);
    int x1 = x0 + 1 < nx ? x0 + 1 : x0;
    int y1 = y0 + 1 < ny ? y0 + 1 : y0;
    double fx = x - x0, fy = y - y0;
    a[ix(y0,x0,nx)] += (1.0-fy)*(1.0-fx)*value;
    a[ix(y0,x1,nx)] += (1.0-fy)*fx*value;
    a[ix(y1,x0,nx)] += fy*(1.0-fx)*value;
    a[ix(y1,x1,nx)] += fy*fx*value;
}

static inline int nearest(int q, int n) { return clampi(q, n); }

/* Output debug record per requested target and group has 20 doubles:
   [group,icf,extension_valid,anchor_found,target_psi,anchor_x,anchor_y,
    anchor_psi,signed_distance_cells,theta_anchor,theta_raw,theta_clipped,
    candidate_psi,delta_psi,beta1_fraction_at_anchor,beta2_fraction_at_anchor,
    corrected_psi,anchor_was_interpolated,fraction_sample_x,fraction_sample_y].
   The profile anchor is an exact ICF=2 grid node; beta fractions are sampled
   bilinearly at the first point on the ray whose nearest grid node is ICF=2. */
int grose_profile_correct_v001(
    const double *phi, double *out, int ny, int nx, double h, double eta,
    double *psi, uint8_t *icf, double *normal_x, double *normal_y,
    double *ext_x, double *ext_y, double *ext_weight, double *psi_corrected,
    double *anchor_x, double *anchor_y, double *fraction_x, double *fraction_y,
    uint8_t *anchor_valid,
    const int *watch_yx, int nwatch, double *debug, double *stats) {
    if (!phi || !out || !psi || !icf || !normal_x || !normal_y || !ext_x || !ext_y ||
        !ext_weight || !psi_corrected || !anchor_x || !anchor_y ||
        !fraction_x || !fraction_y || !anchor_valid ||
        ny < 2 || nx < 2 || !(h > 0.0) || !(eta > 0.0) || nwatch < 0) return 1;
    const size_t n = (size_t)ny * (size_t)nx;
    const double lcells = eta / h;
    const double coslim = cos(M_PI / lcells);
    const double high = 0.5 + 0.5*coslim;
    const double low = 1.0 - high;
    const double marker = 1.0 / sqrt(2.0); /* h/sqrt(2), in grid units */
    if (stats) memset(stats, 0, 16*sizeof(double));
    if (debug && nwatch > 0) for (int i = 0; i < nwatch*2*20; ++i) debug[i] = NAN;
    for (size_t q = 0; q < n; ++q) {
        psi[q] = phi[q];
        psi[n+q] = phi[n+q] + phi[2*n+q];
    }
    memcpy(out, phi, 3*n*sizeof(double));

    for (int g = 0; g < 2; ++g) {
        double *pg = psi + (size_t)g*n;
        uint8_t *mg = icf + (size_t)g*n;
        double *nxg = normal_x + (size_t)g*n;
        double *nyg = normal_y + (size_t)g*n;
        double *exg = ext_x + (size_t)g*n;
        double *eyg = ext_y + (size_t)g*n;
        double *ewg = ext_weight + (size_t)g*n;
        double *pcg = psi_corrected + (size_t)g*n;
        double *axg = anchor_x + (size_t)g*n;
        double *ayg = anchor_y + (size_t)g*n;
        double *fxg = fraction_x + (size_t)g*n;
        double *fyg = fraction_y + (size_t)g*n;
        uint8_t *avg = anchor_valid + (size_t)g*n;
        memset(exg,0,n*sizeof(double)); memset(eyg,0,n*sizeof(double));
        memset(ewg,0,n*sizeof(double)); memset(avg,0,n*sizeof(uint8_t));
        for (size_t q=0; q<n; ++q) { fxg[q]=NAN; fyg[q]=NAN; }
        for (int y = 0; y < ny; ++y) for (int x = 0; x < nx; ++x) {
            const size_t q = ix(y,x,nx);
            const double v = pg[q];
            if (v >= low && v <= high) mg[q] = 2;
            else if (v < low || v > high) mg[q] = 1;
            else mg[q] = 0;
            if (stats) { if (mg[q] == 1) stats[g] += 1.; else if (mg[q] == 2) stats[2+g] += 1.; }
            if (mg[q] != 2) { nxg[q] = nyg[q] = 0.0; continue; }
            const int ym = y > 0 ? y-1 : y, yp = y+1 < ny ? y+1 : y;
            const int xm = x > 0 ? x-1 : x, xp = x+1 < nx ? x+1 : x;
            const double gx = (pg[ix(y,xp,nx)] - pg[ix(y,xm,nx)])/(2.0*h);
            const double gy = (pg[ix(yp,x,nx)] - pg[ix(ym,x,nx)])/(2.0*h);
            const double norm = hypot(gx,gy);
            if (norm > 0.0 && isfinite(norm)) { nxg[q] = gx/norm; nyg[q] = gy/norm; }
            else { nxg[q] = nyg[q] = 0.0; }
        }
        /* Appendix G: a discrete “nearest to bulk” ICF=2 source is an
           interface-interior node adjacent to an ICF=1 edge node. Project
           outward from the interface: toward nil for psi<1/2 and toward bulk
           for psi>=1/2. Deposit three h/sqrt(2) marker steps, accumulating
           source normals with bilinear weights. */
        for (int y = 0; y < ny; ++y) for (int x = 0; x < nx; ++x) {
            const size_t q = ix(y,x,nx);
            if (mg[q] != 2 || (nxg[q] == 0.0 && nyg[q] == 0.0)) continue;
            int frontier = 0;
            for (int oy=-1; oy<=1 && !frontier; ++oy) for (int ox=-1; ox<=1; ++ox) {
                if (ox==0 && oy==0) continue;
                const int yy=y+oy, xx=x+ox;
                if (yy>=0 && yy<ny && xx>=0 && xx<nx && mg[ix(yy,xx,nx)]==1) {
                    frontier=1; break;
                }
            }
            if (!frontier) continue;
            if (stats) stats[12+g] += 1.0;
            const int outward = pg[q] < 0.5 ? -1 : 1;
            for (int k = 1; k <= 3; ++k) {
                const double xx = x + (double)outward * k * marker * nxg[q];
                const double yy = y + (double)outward * k * marker * nyg[q];
                /* Accumulate the source normal and bilinear weights separately. */
                double wx0 = fmin(fmax(xx,0.0),(double)(nx-1));
                double wy0 = fmin(fmax(yy,0.0),(double)(ny-1));
                int x0=(int)floor(wx0),y0=(int)floor(wy0),x1=x0+1<nx?x0+1:x0,y1=y0+1<ny?y0+1:y0;
                double fx=wx0-x0,fy=wy0-y0;
                const int tx[4]={x0,x1,x0,x1}, ty[4]={y0,y0,y1,y1};
                const double ww[4]={(1-fx)*(1-fy),fx*(1-fy),(1-fx)*fy,fx*fy};
                for(int z=0;z<4;++z){
                    const size_t t=ix(ty[z],tx[z],nx);
                    exg[t]+=ww[z]*nxg[q]; eyg[t]+=ww[z]*nyg[q]; ewg[t]+=ww[z];
                }
            }
        }
        for (size_t q = 0; q < n; ++q) {
            if (ewg[q] > 0.0) {
                exg[q] /= ewg[q]; eyg[q] /= ewg[q];
                const double norm = hypot(exg[q],eyg[q]);
                if (norm > 0.0) { exg[q] /= norm; eyg[q] /= norm; }
                else ewg[q] = 0.0;
            }
        }
        memcpy(pcg,pg,n*sizeof(double));
        for (int y = 0; y < ny; ++y) for (int x = 0; x < nx; ++x) {
            const size_t q=ix(y,x,nx);
            if (mg[q] != 1 || ewg[q] <= 0.0) continue;
            const double v=pg[q], nxv=exg[q], nyv=eyg[q];
            const double dir = v < 0.5 ? 1.0 : -1.0;
            double found=0.0, ax=(double)x, ay=(double)y, ap=0.0, dist=0.0;
            /* March toward the first interface-interior node along the extended
               normal. Samples are bilinear; the ICF test is at nearest node. */
            for (int k=1;k<=12;++k) {
                const double t=0.25*(double)k;
                const double sx=(double)x+dir*t*nxv;
                const double sy=(double)y+dir*t*nyv;
                const int xx=nearest((int)llround(sx),nx), yy=nearest((int)llround(sy),ny);
                if (mg[ix(yy,xx,nx)]==2) {
                    found=1.0;
                    /* Eq.39 uses the first interior NODE value. The ray hit
                       remains separately available for Eq.42's bilinear
                       grain-fraction estimate. */
                    ax=(double)xx; ay=(double)yy;
                    dist=-((ax-(double)x)*nxv+(ay-(double)y)*nyv);
                    ap=pg[ix(yy,xx,nx)];
                    fxg[q]=sx; fyg[q]=sy;
                    break;
                }
            }
            if (found==0.0) { if(stats) stats[6+g]+=1.; continue; }
            const double theta0=acos(fmin(1.0,fmax(-1.0,1.0-2.0*ap)));
            const double theta_raw=theta0+M_PI*dist/lcells;
            const double theta=fmin(M_PI,fmax(0.0,theta_raw));
            const double candidate=0.5-0.5*cos(theta);
            pcg[q]=candidate; axg[q]=ax; ayg[q]=ay; avg[q]=1;
            if(stats){ stats[4+g]+=1.; const double dd=fabs(candidate-v); if(dd>stats[8+g])stats[8+g]=dd; }
            for(int w=0;w<nwatch;++w){
                const int wy=watch_yx[2*w],wx=watch_yx[2*w+1];
                if(wy!=y||wx!=x||!debug)continue;
                double *d=debug+((size_t)w*2+(size_t)g)*20;
                d[0]=(double)g;d[1]=(double)mg[q];d[2]=1.;d[3]=found;d[4]=v;
                d[5]=ax;d[6]=ay;d[7]=ap;d[8]=dist;d[9]=theta0;d[10]=theta_raw;
                d[11]=theta;d[12]=candidate;d[13]=candidate-v;d[14]=NAN;d[15]=NAN;
                if (g == 1) {
                    double f1 = bilinear(phi+n,ny,nx,fyg[q],fxg[q]);
                    double f2 = bilinear(phi+2*n,ny,nx,fyg[q],fxg[q]);
                    const double fsum = f1 + f2;
                    if (fsum > 0.0) { d[14] = f1/fsum; d[15] = f2/fsum; }
                }
                d[16]=candidate;d[17]=0.;d[18]=fxg[q];d[19]=fyg[q];
            }
        }
    }

    /* Eq. 40 phase-field normalization. Preserve original psi if both candidate
       values vanish only through a malformed input/roundoff, and expose that. */
    for (size_t q=0;q<n;++q) {
        double a=psi_corrected[q],b=psi_corrected[n+q],s=a+b;
        if (!(s>0.0) || !isfinite(s)) { a=psi[q];b=psi[n+q];s=a+b; if(stats)stats[10]+=1.; }
        a/=s;b/=s;
        const double beta_old=psi[n+q];
        double f1=0.0,f2=0.0;
        if(beta_old>0.0){ f1=phi[n+q]/beta_old;f2=phi[2*n+q]/beta_old; }
        else if(b>0.0 && anchor_valid[n+q]) {
            f1=bilinear(phi+n,ny,nx,fraction_y[n+q],fraction_x[n+q]);
            f2=bilinear(phi+2*n,ny,nx,fraction_y[n+q],fraction_x[n+q]);
            const double fs=f1+f2;
            if(fs>0.0){f1/=fs;f2/=fs;}else {b=0.0;a=1.0;if(stats)stats[11]+=1.;}
        } else if(b>0.0) { b=0.0;a=1.0;if(stats)stats[11]+=1.; }
        double fsum=f1+f2;
        if(fsum>0.0){
            if(f1==0.0) f2=1.0;
            else if(f2==0.0) f1=1.0;
            else { f1/=fsum; f2/=fsum; }
        } else if(b>0.0) { b=0.0; a=1.0; if(stats)stats[11]+=1.; }
        out[q]=a;out[n+q]=b*f1;out[2*n+q]=b*f2;
        const double residual=1.0-(out[q]+out[n+q]+out[2*n+q]);
        int largest=0;
        if(out[n+q]>out[(size_t)largest*n+q])largest=1;
        if(out[2*n+q]>out[(size_t)largest*n+q])largest=2;
        out[(size_t)largest*n+q]+=residual;
    }
    if(debug && nwatch>0) for(int w=0;w<nwatch;++w){
        const int y=watch_yx[2*w],x=watch_yx[2*w+1];
        if(y<0||y>=ny||x<0||x>=nx)continue;
        for (int g=0; g<2; ++g) {
            double *d=debug+((size_t)w*2+(size_t)g)*20;
            if (!isfinite(d[0])) continue;
            d[16]=psi_corrected[(size_t)g*n+ix(y,x,nx)];
        }
    }
    return 0;
}
