#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

/* Fused, single-thread implementation of the same cell-centred B2 Euler step
   used by b2_rhs.py + project_simplex_roundoff_safe().  This is a performance
   path only; the NumPy routines remain the reference implementation. */
int b2_step(const double *phi, double *out, double *raw_rhs, double *lap,
            uint8_t *active, const double *gamma, const double *delta_f,
            const double *mobility, double *flux, int ny, int nx, double dx,
            double eta, double dt, double epsilon_factor) {
    if (!phi || !out || !raw_rhs || !lap || !active || !gamma || !delta_f ||
        !mobility || !flux || ny < 2 || nx < 2 || dx <= 0.0 || eta <= 0.0 ||
        dt < 0.0 || epsilon_factor < 0.0) return 1;

    const size_t n = (size_t)ny * (size_t)nx;
    const double inv_dx2 = 1.0 / (dx * dx);
    const double c2 = (M_PI / eta) * (M_PI / eta);
    const double b_prefactor = M_PI / eta;

    /* Edge-copy ghosts give zero normal derivative at the outer cell faces. */
    for (int p = 0; p < 3; ++p) {
        const size_t base = (size_t)p * n;
        for (int y = 0; y < ny; ++y) {
            const int ym = (y > 0) ? y - 1 : y;
            const int yp = (y + 1 < ny) ? y + 1 : y;
            for (int x = 0; x < nx; ++x) {
                const int xm = (x > 0) ? x - 1 : x;
                const int xp = (x + 1 < nx) ? x + 1 : x;
                const size_t q = base + (size_t)y * nx + x;
                const double center = phi[q];
                const double d2x = (phi[base + (size_t)y * nx + xp] - 2.0 * center +
                                    phi[base + (size_t)y * nx + xm]) * inv_dx2;
                const double d2y = (phi[base + (size_t)yp * nx + x] - 2.0 * center +
                                    phi[base + (size_t)ym * nx + x]) * inv_dx2;
                lap[q] = d2x + d2y;
                active[q] = (uint8_t)((center > 0.0) || (lap[q] != 0.0));
            }
        }
    }

    /* One unordered-pair pass preserves the antisymmetric Eq. (27) sum. */
    for (size_t q = 0; q < n; ++q) {
        int n_active = (int)active[q] + (int)active[n + q] + (int)active[2 * n + q];
        if (n_active <= 0) return 2;
        const double inv_n = 1.0 / (double)n_active;
        const double p0 = phi[q], pa = phi[n + q], pb = phi[2 * n + q];
        const double l0 = lap[q], la = lap[n + q], lb = lap[2 * n + q];
        double r0 = 0.0, ra = 0.0, rb = 0.0;

        if (active[q] && active[n + q]) {
            const double kappa = inv_n * (la - l0 + c2 * (pa - p0));
            const double junction = (gamma[1 * 3 + 2] - gamma[0 * 3 + 2]) *
                                    (lb + c2 * pb) * inv_n;
            const double bulk = b_prefactor * (p0 + pa) * sqrt(fmax(p0 * pa, 0.0));
            const double pair = mobility[0 * 3 + 1] *
                (bulk * delta_f[0 * 3 + 1] - gamma[0 * 3 + 1] * kappa + junction);
            r0 += pair; ra -= pair;
        }
        if (active[q] && active[2 * n + q]) {
            const double kappa = inv_n * (lb - l0 + c2 * (pb - p0));
            const double junction = (gamma[2 * 3 + 1] - gamma[0 * 3 + 1]) *
                                    (la + c2 * pa) * inv_n;
            const double bulk = b_prefactor * (p0 + pb) * sqrt(fmax(p0 * pb, 0.0));
            const double pair = mobility[0 * 3 + 2] *
                (bulk * delta_f[0 * 3 + 2] - gamma[0 * 3 + 2] * kappa + junction);
            r0 += pair; rb -= pair;
        }
        if (active[n + q] && active[2 * n + q]) {
            const double kappa = inv_n * (lb - la + c2 * (pb - pa));
            const double junction = (gamma[2 * 3 + 0] - gamma[1 * 3 + 0]) *
                                    (l0 + c2 * p0) * inv_n;
            const double bulk = b_prefactor * (pa + pb) * sqrt(fmax(pa * pb, 0.0));
            const double pair = mobility[1 * 3 + 2] *
                (bulk * delta_f[1 * 3 + 2] - gamma[1 * 3 + 2] * kappa + junction);
            ra += pair; rb -= pair;
        }
        raw_rhs[q] = r0;
        raw_rhs[n + q] = ra;
        raw_rhs[2 * n + q] = rb;
    }

    /* Roundoff-safe ternary Euclidean simplex projection.  The active-rank
       guard follows the NumPy KKT criterion; inactive support is exact zero. */
    const double guard_scale = epsilon_factor * DBL_EPSILON;
    double raw_phase0_sum = 0.0;
    double projected_phase0_delta_sum = 0.0;
    for (size_t q = 0; q < n; ++q) {
        const double trial[3] = {
            phi[q] + dt * raw_rhs[q],
            phi[n + q] + dt * raw_rhs[n + q],
            phi[2 * n + q] + dt * raw_rhs[2 * n + q]
        };
        double s1 = trial[0], s2 = trial[1], s3 = trial[2];
        if (s1 < s2) { double t = s1; s1 = s2; s2 = t; }
        if (s2 < s3) { double t = s2; s2 = s3; s3 = t; }
        if (s1 < s2) { double t = s1; s1 = s2; s2 = t; }
        const double scale = fmax(1.0, fmax(fabs(trial[0]), fmax(fabs(trial[1]), fabs(trial[2]))));
        const double tol = guard_scale * scale;
        const double margin2 = s2 - (s1 + s2 - 1.0) / 2.0;
        const double margin3 = s3 - (s1 + s2 + s3 - 1.0) / 3.0;
        double theta;
        if (margin3 > tol) theta = (s1 + s2 + s3 - 1.0) / 3.0;
        else if (margin2 > tol) theta = (s1 + s2 - 1.0) / 2.0;
        else theta = s1 - 1.0;

        double projected[3];
        double sum = 0.0;
        for (int p = 0; p < 3; ++p) {
            const double value = trial[p] - theta;
            projected[p] = (value > tol) ? value : 0.0;
            sum += projected[p];
        }
        if (!(sum > 0.0) || !isfinite(sum)) return 3;
        for (int p = 0; p < 3; ++p) projected[p] /= sum;
        int largest = 0;
        if (projected[1] > projected[largest]) largest = 1;
        if (projected[2] > projected[largest]) largest = 2;
        projected[largest] += 1.0 - (projected[0] + projected[1] + projected[2]);
        out[q] = projected[0];
        out[n + q] = projected[1];
        out[2 * n + q] = projected[2];
        raw_phase0_sum += raw_rhs[q];
        projected_phase0_delta_sum += projected[0] - phi[q];
    }
    flux[0] = raw_phase0_sum * dx * dx;
    flux[1] = (dt > 0.0) ? projected_phase0_delta_sum * dx * dx / dt : 0.0;
    flux[2] = flux[1] - flux[0];
    return 0;
}
