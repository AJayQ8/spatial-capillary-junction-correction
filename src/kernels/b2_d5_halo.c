#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

/* B2 pair-law kernel with an externally supplied, width-halo capillary
   multiplier.  The force remains (1 + lambda * max_filter(w0))*C + B. */

static inline int clamp_index(int i, int n) {
    return (i < 0) ? 0 : ((i >= n) ? n - 1 : i);
}

/* Exact centered max filter of radius `radius` with edge-copy extension.
   Monotone deque gives O(length) work per line, independent of filter width. */
static void max_filter_line(const double *src, double *dst, int length,
                            int stride, int radius, int *deque) {
    int head = 0;
    int tail = 0;
    for (int t = -radius; t < length + radius; ++t) {
        const int index = clamp_index(t, length);
        const double value = src[(size_t)index * (size_t)stride];
        while (tail > head) {
            const int previous_t = deque[tail - 1];
            const int previous_index = clamp_index(previous_t, length);
            if (src[(size_t)previous_index * (size_t)stride] > value) break;
            --tail;
        }
        deque[tail++] = t;

        const int first_allowed = t - 2 * radius;
        while (tail > head && deque[head] < first_allowed) ++head;

        const int center = t - radius;
        if (center >= 0 && center < length) {
            const int best_index = clamp_index(deque[head], length);
            dst[(size_t)center * (size_t)stride] =
                src[(size_t)best_index * (size_t)stride];
        }
    }
}

static void build_halo(const double *phi, double *weight0, double *horizontal,
                       double *halo, int *deque, int ny, int nx, int radius,
                       double *summary, double dx) {
    const size_t n = (size_t)ny * (size_t)nx;
    for (size_t q = 0; q < n; ++q)
        weight0[q] = 27.0 * phi[q] * phi[n + q] * phi[2 * n + q];

    for (int y = 0; y < ny; ++y) {
        const size_t row = (size_t)y * (size_t)nx;
        max_filter_line(weight0 + row, horizontal + row, nx, 1,
                        radius, deque);
    }
    for (int x = 0; x < nx; ++x)
        max_filter_line(horizontal + x, halo + x, ny, nx,
                        radius, deque);

    double maximum = 0.0;
    double sum = 0.0;
    double positive = 0.0;
    for (size_t q = 0; q < n; ++q) {
        const double v = halo[q];
        if (v > maximum) maximum = v;
        if (v > 0.0) positive += 1.0;
        sum += v;
    }
    summary[0] = maximum;
    summary[1] = positive;
    summary[2] = sum * dx * dx;
}

int b2_capillary_halo_step(
    const double *phi, double *out, double *raw_rhs, double *lap,
    uint8_t *active, const double *gamma, const double *delta_f,
    const double *mobility, double *flux, double *halo_summary,
    double *weight0, double *horizontal, double *halo, int *deque,
    int ny, int nx, double dx, double eta, double dt,
    double epsilon_factor, double lambda, int radius) {
    if (!phi || !out || !raw_rhs || !lap || !active || !gamma || !delta_f ||
        !mobility || !flux || !halo_summary || !weight0 || !horizontal ||
        !halo || !deque || ny < 2 || nx < 2 || dx <= 0.0 || eta <= 0.0 ||
        dt < 0.0 || epsilon_factor < 0.0 || lambda < 0.0 ||
        !isfinite(lambda) || radius < 0 || radius > ny || radius > nx)
        return 1;

    const size_t n = (size_t)ny * (size_t)nx;
    const double inv_dx2 = 1.0 / (dx * dx);
    const double c2 = (M_PI / eta) * (M_PI / eta);
    const double b_prefactor = M_PI / eta;

    build_halo(phi, weight0, horizontal, halo, deque, ny, nx, radius,
               halo_summary, dx);

    /* Edge-copy ghosts give zero normal derivative at outer cell faces. */
    for (int p = 0; p < 3; ++p) {
        const size_t base = (size_t)p * n;
        for (int y = 0; y < ny; ++y) {
            const int ym = (y > 0) ? y - 1 : y;
            const int yp = (y + 1 < ny) ? y + 1 : y;
            for (int x = 0; x < nx; ++x) {
                const int xm = (x > 0) ? x - 1 : x;
                const int xp = (x + 1 < nx) ? x + 1 : x;
                const size_t q = base + (size_t)y * (size_t)nx + (size_t)x;
                const double center = phi[q];
                const double d2x = (phi[base + (size_t)y * (size_t)nx + (size_t)xp]
                                    - 2.0 * center
                                    + phi[base + (size_t)y * (size_t)nx + (size_t)xm]) * inv_dx2;
                const double d2y = (phi[base + (size_t)yp * (size_t)nx + (size_t)x]
                                    - 2.0 * center
                                    + phi[base + (size_t)ym * (size_t)nx + (size_t)x]) * inv_dx2;
                lap[q] = d2x + d2y;
                active[q] = (uint8_t)((center > 0.0) || (lap[q] != 0.0));
            }
        }
    }

    /* One unordered-pair pass preserves the antisymmetric Eq. (27) sum. */
    for (size_t q = 0; q < n; ++q) {
        const int n_active = (int)active[q] + (int)active[n + q] +
                             (int)active[2 * n + q];
        if (n_active <= 0) return 2;
        const double inv_n = 1.0 / (double)n_active;
        const double p0 = phi[q], pa = phi[n + q], pb = phi[2 * n + q];
        const double l0 = lap[q], la = lap[n + q], lb = lap[2 * n + q];
        double r0 = 0.0, ra = 0.0, rb = 0.0;
        double b0 = 0.0, ba = 0.0, bb = 0.0;

        if (active[q] && active[n + q]) {
            const double kappa = inv_n * (la - l0 + c2 * (pa - p0));
            const double junction = (gamma[5] - gamma[2]) * (lb + c2 * pb) * inv_n;
            const double bulk = b_prefactor * (p0 + pa) * sqrt(fmax(p0 * pa, 0.0));
            const double pair = mobility[1] *
                (bulk * delta_f[1] - gamma[1] * kappa + junction);
            r0 += pair; ra -= pair;
            const double bulk_pair = mobility[1] * bulk * delta_f[1];
            b0 += bulk_pair; ba -= bulk_pair;
        }
        if (active[q] && active[2 * n + q]) {
            const double kappa = inv_n * (lb - l0 + c2 * (pb - p0));
            const double junction = (gamma[7] - gamma[1]) * (la + c2 * pa) * inv_n;
            const double bulk = b_prefactor * (p0 + pb) * sqrt(fmax(p0 * pb, 0.0));
            const double pair = mobility[2] *
                (bulk * delta_f[2] - gamma[2] * kappa + junction);
            r0 += pair; rb -= pair;
            const double bulk_pair = mobility[2] * bulk * delta_f[2];
            b0 += bulk_pair; bb -= bulk_pair;
        }
        if (active[n + q] && active[2 * n + q]) {
            const double kappa = inv_n * (lb - la + c2 * (pb - pa));
            const double junction = (gamma[6] - gamma[3]) * (l0 + c2 * p0) * inv_n;
            const double bulk = b_prefactor * (pa + pb) * sqrt(fmax(pa * pb, 0.0));
            const double pair = mobility[5] *
                (bulk * delta_f[5] - gamma[5] * kappa + junction);
            ra += pair; rb -= pair;
            const double bulk_pair = mobility[5] * bulk * delta_f[5];
            ba += bulk_pair; bb -= bulk_pair;
        }

        const double extra = lambda * halo[q];
        raw_rhs[q] = r0 + extra * (r0 - b0);
        raw_rhs[n + q] = ra + extra * (ra - ba);
        raw_rhs[2 * n + q] = rb + extra * (rb - bb);
    }

    /* Roundoff-safe ternary Euclidean simplex projection. */
    const double guard_scale = epsilon_factor * DBL_EPSILON;
    double raw_phase0_sum = 0.0;
    double projected_phase0_delta_sum = 0.0;
    for (size_t q = 0; q < n; ++q) {
        const double trial[3] = {phi[q] + dt * raw_rhs[q],
                                 phi[n + q] + dt * raw_rhs[n + q],
                                 phi[2 * n + q] + dt * raw_rhs[2 * n + q]};
        double s1 = trial[0], s2 = trial[1], s3 = trial[2];
        if (s1 < s2) { double t = s1; s1 = s2; s2 = t; }
        if (s2 < s3) { double t = s2; s2 = s3; s3 = t; }
        if (s1 < s2) { double t = s1; s1 = s2; s2 = t; }
        const double scale = fmax(1.0, fmax(fabs(trial[0]),
                                   fmax(fabs(trial[1]), fabs(trial[2]))));
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
