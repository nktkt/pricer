// pricer/american.hpp — American (early-exercise) option pricing.
//
// Unlike a European option, an American option may be exercised at any time up to
// expiry, so its value is an *optimal stopping* problem. Two complementary
// methods are provided:
//   * a Cox–Ross–Rubinstein binomial tree (`binomial_price`) — a fast, accurate
//     lattice that does backward induction with an early-exercise test at every
//     node; the European mode cross-checks against Black–Scholes;
//   * Longstaff–Schwartz least-squares Monte Carlo (`lsm_american`) — the
//     regression-based MC method, which estimates the continuation value by
//     regressing discounted future cash flows on a polynomial of the spot and
//     exercises when immediate payoff beats it. It reuses the counter-based RNG
//     and cross-checks against the tree.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "pricer/black_scholes.hpp"  // OptionType
#include "pricer/rng.hpp"            // cb_normal

namespace pricer {

// Intrinsic payoff of a vanilla option at spot S.
inline double vanilla_payoff(OptionType type, double S, double K) {
    return type == OptionType::Call ? std::max(S - K, 0.0) : std::max(K - S, 0.0);
}

// Cox–Ross–Rubinstein binomial price. `american=false` gives the European price
// (a backward induction with no early exercise), which converges to Black–Scholes
// as `steps` grows; `american=true` adds the early-exercise test at each node.
inline double binomial_price(OptionType type, double S, double K, double r, double sigma, double T,
                             int steps, bool american) {
    if (steps < 1) steps = 1;
    const double dt = T / steps;
    const double u = std::exp(sigma * std::sqrt(dt));
    const double d = 1.0 / u;
    const double disc = std::exp(-r * dt);
    const double p = (std::exp(r * dt) - d) / (u - d);  // risk-neutral up-probability

    // Terminal layer: payoff at each of the steps+1 final nodes (i = #down moves).
    std::vector<double> v(static_cast<std::size_t>(steps) + 1);
    for (int i = 0; i <= steps; ++i)
        v[i] = vanilla_payoff(type, S * std::pow(u, steps - i) * std::pow(d, i), K);

    // Backward induction.
    for (int step = steps - 1; step >= 0; --step) {
        for (int i = 0; i <= step; ++i) {
            double val = disc * (p * v[i] + (1.0 - p) * v[i + 1]);
            if (american) {
                const double spot = S * std::pow(u, step - i) * std::pow(d, i);
                val = std::max(val, vanilla_payoff(type, spot, K));
            }
            v[i] = val;
        }
    }
    return v[0];
}

// Greeks from a single CRR tree pass. The lattice already holds option values at
// neighbouring spot/time nodes, so delta, gamma and theta come "for free" by
// finite differences over those nodes (the standard Hull construction; with CRR
// u·d = 1 the middle step-2 node sits at the original spot, giving theta over
// 2·dt). Vega and rho — sensitivities to model inputs, not lattice coordinates —
// are central differences from re-pricing. Works for European and American.
inline Greeks binomial_greeks(OptionType type, double S, double K, double r, double sigma, double T,
                              int steps, bool american) {
    if (steps < 3) steps = 3;  // need the step-2 layer for gamma/theta
    const double dt = T / steps;
    const double u = std::exp(sigma * std::sqrt(dt));
    const double d = 1.0 / u;
    const double disc = std::exp(-r * dt);
    const double p = (std::exp(r * dt) - d) / (u - d);

    std::vector<double> v(static_cast<std::size_t>(steps) + 1);
    for (int i = 0; i <= steps; ++i)
        v[i] = vanilla_payoff(type, S * std::pow(u, steps - i) * std::pow(d, i), K);

    double v00 = 0, v10 = 0, v11 = 0, v20 = 0, v21 = 0, v22 = 0;  // captured early layers
    for (int step = steps - 1; step >= 0; --step) {
        for (int i = 0; i <= step; ++i) {
            double val = disc * (p * v[i] + (1.0 - p) * v[i + 1]);
            if (american)
                val = std::max(val, vanilla_payoff(type, S * std::pow(u, step - i) * std::pow(d, i), K));
            v[i] = val;
        }
        if (step == 2) { v20 = v[0]; v21 = v[1]; v22 = v[2]; }
        if (step == 1) { v10 = v[0]; v11 = v[1]; }
    }
    v00 = v[0];

    Greeks g{};
    g.price = v00;
    g.delta = (v10 - v11) / (S * u - S * d);
    const double su2 = S * u * u, sd2 = S * d * d;
    const double delta_up = (v20 - v21) / (su2 - S);
    const double delta_dn = (v21 - v22) / (S - sd2);
    g.gamma = (delta_up - delta_dn) / (0.5 * (su2 - sd2));
    g.theta = (v21 - v00) / (2.0 * dt);  // v21 is at spot S, two steps (2·dt) later

    // Vega (per 1.0 of vol) and rho (per 1.0 of rate) by central differences.
    const double hs = 0.01, hr = 1e-4;
    g.vega = (binomial_price(type, S, K, r, sigma + hs, T, steps, american) -
              binomial_price(type, S, K, r, sigma - hs, T, steps, american)) / (2.0 * hs);
    g.rho = (binomial_price(type, S, K, r + hr, sigma, T, steps, american) -
             binomial_price(type, S, K, r - hr, sigma, T, steps, american)) / (2.0 * hr);
    return g;
}

namespace detail {

// Solve a 3x3 system A x = b in place (Gaussian elimination, partial pivoting).
// Returns false if (near-)singular.
inline bool solve3x3(double A[3][3], double b[3], double x[3]) {
    for (int col = 0; col < 3; ++col) {
        int piv = col;
        for (int r2 = col + 1; r2 < 3; ++r2)
            if (std::fabs(A[r2][col]) > std::fabs(A[piv][col])) piv = r2;
        if (std::fabs(A[piv][col]) < 1e-12) return false;
        if (piv != col) {
            for (int c = 0; c < 3; ++c) std::swap(A[piv][c], A[col][c]);
            std::swap(b[piv], b[col]);
        }
        for (int r2 = col + 1; r2 < 3; ++r2) {
            const double f = A[r2][col] / A[col][col];
            for (int c = col; c < 3; ++c) A[r2][c] -= f * A[col][c];
            b[r2] -= f * b[col];
        }
    }
    for (int row = 2; row >= 0; --row) {
        double s = b[row];
        for (int c = row + 1; c < 3; ++c) s -= A[row][c] * x[c];
        x[row] = s / A[row][row];
    }
    return true;
}

}  // namespace detail

// Longstaff–Schwartz least-squares Monte Carlo for an American option. Simulates
// GBM paths with the counter-based RNG, then walks backward estimating the
// continuation value by regressing discounted future cash flows on {1, x, x^2}
// (x = S/K, normalized for conditioning) over the in-the-money paths, exercising
// where the immediate payoff exceeds the estimate. Reproducible for a given seed.
inline double lsm_american(OptionType type, double S, double K, double r, double sigma, double T,
                           long n_paths, int n_steps, std::uint64_t seed = 12345) {
    if (n_steps < 1) n_steps = 1;
    if (n_paths < 1) n_paths = 1;
    const double dt = T / n_steps;
    const double drift = (r - 0.5 * sigma * sigma) * dt;
    const double vol = sigma * std::sqrt(dt);
    const double disc = std::exp(-r * dt);  // one-step discount factor

    // Simulate and store the full paths (n_paths x (n_steps+1)).
    const std::size_t W = static_cast<std::size_t>(n_steps) + 1;
    std::vector<double> path(static_cast<std::size_t>(n_paths) * W);
    for (long pth = 0; pth < n_paths; ++pth) {
        double s = S;
        path[static_cast<std::size_t>(pth) * W] = s;
        for (int t = 1; t <= n_steps; ++t) {
            const std::uint64_t ctr = static_cast<std::uint64_t>(pth) * n_steps + (t - 1);
            s *= std::exp(drift + vol * cb_normal(seed, ctr));
            path[static_cast<std::size_t>(pth) * W + t] = s;
        }
    }

    // Cash flow per path and the step at which it is realized; start at expiry.
    std::vector<double> cashflow(static_cast<std::size_t>(n_paths));
    std::vector<int> exercise_step(static_cast<std::size_t>(n_paths), n_steps);
    for (long pth = 0; pth < n_paths; ++pth)
        cashflow[pth] = vanilla_payoff(type, path[static_cast<std::size_t>(pth) * W + n_steps], K);

    // Backward induction over the interior exercise dates.
    for (int t = n_steps - 1; t >= 1; --t) {
        double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        double rhs[3] = {0, 0, 0};
        bool any = false;
        for (long pth = 0; pth < n_paths; ++pth) {
            const double s = path[static_cast<std::size_t>(pth) * W + t];
            if (vanilla_payoff(type, s, K) <= 0.0) continue;  // regress on ITM paths only
            any = true;
            const double x = s / K;
            const double basis[3] = {1.0, x, x * x};
            // Discount the realized cash flow back from its exercise step to t.
            const double y = cashflow[pth] * std::pow(disc, exercise_step[pth] - t);
            for (int a = 0; a < 3; ++a) {
                for (int c = 0; c < 3; ++c) A[a][c] += basis[a] * basis[c];
                rhs[a] += basis[a] * y;
            }
        }
        double coef[3];
        if (!any || !detail::solve3x3(A, rhs, coef)) continue;  // keep holding this step
        for (long pth = 0; pth < n_paths; ++pth) {
            const double s = path[static_cast<std::size_t>(pth) * W + t];
            const double ex = vanilla_payoff(type, s, K);
            if (ex <= 0.0) continue;
            const double x = s / K;
            const double cont = coef[0] + coef[1] * x + coef[2] * x * x;
            if (ex > cont) {  // exercising now beats the estimated continuation value
                cashflow[pth] = ex;
                exercise_step[pth] = t;
            }
        }
    }

    // Value today: average of each path's cash flow discounted to t = 0.
    double sum = 0.0;
    for (long pth = 0; pth < n_paths; ++pth)
        sum += cashflow[pth] * std::pow(disc, exercise_step[pth]);
    return sum / static_cast<double>(n_paths);
}

}  // namespace pricer
