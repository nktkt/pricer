// pricer/optimize.hpp — a small Levenberg–Marquardt least-squares solver.
//
// Minimizes sum_i residual_i(p)^2 over parameters p. The Jacobian is computed by
// central finite differences, so any residual function works without supplying
// derivatives. This is the engine behind model calibration (e.g. fitting a vol
// model to market quotes).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace pricer::opt {

struct LMResult {
    std::vector<double> params;  // best parameters found
    double cost;                 // final sum of squared residuals
    int iterations;
    bool converged;
};

namespace detail {
// Solve the n×n system A x = b by Gaussian elimination with partial pivoting.
inline std::vector<double> solve_linear(std::vector<std::vector<double>> A, std::vector<double> b) {
    const int n = static_cast<int>(b.size());
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int r = col + 1; r < n; ++r)
            if (std::fabs(A[r][col]) > std::fabs(A[piv][col])) piv = r;
        std::swap(A[piv], A[col]);
        std::swap(b[piv], b[col]);
        const double d = A[col][col];
        if (std::fabs(d) < 1e-300) continue;  // singular column: leave as is
        for (int r = col + 1; r < n; ++r) {
            const double f = A[r][col] / d;
            for (int c = col; c < n; ++c) A[r][c] -= f * A[col][c];
            b[r] -= f * b[col];
        }
    }
    std::vector<double> x(n, 0.0);
    for (int r = n - 1; r >= 0; --r) {
        double s = b[r];
        for (int c = r + 1; c < n; ++c) s -= A[r][c] * x[c];
        x[r] = (std::fabs(A[r][r]) < 1e-300) ? 0.0 : s / A[r][r];
    }
    return x;
}

inline double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}
}  // namespace detail

// `residual(p, out)` must fill `out` (length m) with the residuals at p.
template <class Residual>
LMResult levenberg_marquardt(Residual residual, std::vector<double> p, int m,
                             int max_iter = 200, double tol = 1e-12) {
    const int n = static_cast<int>(p.size());
    std::vector<double> r(m);
    residual(p, r);
    double cost = detail::dot(r, r);
    double lambda = 1e-3;
    bool converged = false;
    int it = 0;

    for (; it < max_iter; ++it) {
        // Central-difference Jacobian J (m×n).
        std::vector<std::vector<double>> J(m, std::vector<double>(n));
        for (int j = 0; j < n; ++j) {
            const double h = 1e-6 * (std::fabs(p[j]) + 1e-6);
            std::vector<double> pp = p, pm = p, rp(m), rm(m);
            pp[j] += h; pm[j] -= h;
            residual(pp, rp);
            residual(pm, rm);
            for (int i = 0; i < m; ++i) J[i][j] = (rp[i] - rm[i]) / (2 * h);
        }

        // Normal equations: A = JᵀJ, g = Jᵀr.
        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        std::vector<double> g(n, 0.0);
        for (int i = 0; i < m; ++i)
            for (int a = 0; a < n; ++a) {
                g[a] += J[i][a] * r[i];
                for (int b = 0; b < n; ++b) A[a][b] += J[i][a] * J[i][b];
            }
        double gnorm = 0.0;
        for (double v : g) gnorm = std::max(gnorm, std::fabs(v));
        if (gnorm < tol) { converged = true; break; }

        // Damped step, increasing lambda until the cost decreases.
        bool stepped = false;
        for (int tries = 0; tries < 40; ++tries) {
            std::vector<std::vector<double>> Ad = A;
            for (int d = 0; d < n; ++d) Ad[d][d] += lambda * (A[d][d] + 1e-12);  // Marquardt scaling
            std::vector<double> neg_g(n);
            for (int d = 0; d < n; ++d) neg_g[d] = -g[d];
            const std::vector<double> step = detail::solve_linear(Ad, neg_g);

            std::vector<double> pnew(n), rnew(m);
            for (int d = 0; d < n; ++d) pnew[d] = p[d] + step[d];
            residual(pnew, rnew);
            const double cnew = detail::dot(rnew, rnew);

            if (cnew < cost) {
                const double old = cost;
                p = std::move(pnew); r = std::move(rnew); cost = cnew;
                lambda = std::max(lambda * 0.3, 1e-12);
                stepped = true;
                if (old - cnew < tol * (1.0 + old)) converged = true;
                break;
            }
            lambda *= 3.0;
            if (lambda > 1e12) break;
        }
        if (!stepped || converged) break;
    }
    return {std::move(p), cost, it, converged};
}

}  // namespace pricer::opt
