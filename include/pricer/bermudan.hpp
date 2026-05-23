// pricer/bermudan.hpp — Bermudan option pricing by Longstaff–Schwartz LSM.
//
// A Bermudan option may be exercised early, but only on a fixed, finite set of
// dates t_1 < t_2 < … < t_m = T (e.g. a swaption callable each quarter). It sits
// between a European option (exercise only at expiry) and an American option
// (exercise at any time): adding exercise dates can only raise the value, and as
// the dates fill in densely the price converges to the American one.
//
// This reuses the Longstaff–Schwartz machinery from american.hpp (regress the
// discounted continuation value on {1, S/K, (S/K)^2} over the in-the-money paths
// and exercise when the immediate payoff beats it) but over an *arbitrary*,
// possibly irregular, schedule. Because GBM is Markov, the underlying is jumped
// exactly from one exercise date to the next — no fine time grid is needed.
#pragma once
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "pricer/american.hpp"  // vanilla_payoff, detail::solve3x3
#include "pricer/black_scholes.hpp"  // OptionType
#include "pricer/rng.hpp"            // cb_normal

namespace pricer {

// Build m equally-spaced exercise dates T/m, 2T/m, …, T (a convenience for the
// common case; bermudan_lsm itself takes any sorted schedule).
inline std::vector<double> equally_spaced_dates(double T, int m) {
    if (m < 1) m = 1;
    std::vector<double> t(static_cast<std::size_t>(m));
    for (int i = 1; i <= m; ++i) t[i - 1] = T * i / m;
    return t;
}

// Longstaff–Schwartz price of a Bermudan option exercisable on the dates in
// `exercise_times` (sorted ascending, in years, the last one being expiry).
// Paths are simulated exactly from date to date with the counter-based RNG, so
// the result is reproducible for a given seed. With a single date this is the
// European price; with many equally-spaced dates it approaches the American one.
inline double bermudan_lsm(OptionType type, double S, double K, double r, double sigma,
                           const std::vector<double>& exercise_times, long n_paths,
                           std::uint64_t seed = 12345, double q = 0.0) {
    const std::size_t m = exercise_times.size();
    if (m == 0) throw std::invalid_argument("bermudan_lsm: need >= 1 exercise date");
    if (n_paths < 1) n_paths = 1;

    // Per-interval GBM jump parameters: drift_i and vol_i over [t_{i-1}, t_i].
    std::vector<double> drift(m), vol(m);
    double prev = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
        const double dt = exercise_times[i] - prev;
        if (dt <= 0.0) throw std::invalid_argument("bermudan_lsm: exercise_times must increase from 0");
        drift[i] = (r - q - 0.5 * sigma * sigma) * dt;
        vol[i] = sigma * std::sqrt(dt);
        prev = exercise_times[i];
    }

    // Simulate the underlying at each exercise date (n_paths x m), exact jumps.
    std::vector<double> path(static_cast<std::size_t>(n_paths) * m);
    for (long p = 0; p < n_paths; ++p) {
        double s = S;
        for (std::size_t i = 0; i < m; ++i) {
            const std::uint64_t ctr = static_cast<std::uint64_t>(p) * m + i;
            s *= std::exp(drift[i] + vol[i] * cb_normal(seed, ctr));
            path[static_cast<std::size_t>(p) * m + i] = s;
        }
    }

    // Cash flow per path and the date (in years) at which it is realized; start at
    // expiry (the last date).
    std::vector<double> cashflow(static_cast<std::size_t>(n_paths));
    std::vector<double> exercise_t(static_cast<std::size_t>(n_paths), exercise_times[m - 1]);
    for (long p = 0; p < n_paths; ++p)
        cashflow[p] = vanilla_payoff(type, path[static_cast<std::size_t>(p) * m + (m - 1)], K);

    // Backward induction over the interior exercise dates.
    for (std::size_t idx = m - 1; idx-- > 0;) {  // idx = m-2 … 0
        const double ti = exercise_times[idx];
        double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        double rhs[3] = {0, 0, 0};
        bool any = false;
        for (long p = 0; p < n_paths; ++p) {
            const double s = path[static_cast<std::size_t>(p) * m + idx];
            if (vanilla_payoff(type, s, K) <= 0.0) continue;  // regress on ITM paths only
            any = true;
            const double x = s / K;
            const double basis[3] = {1.0, x, x * x};
            const double y = cashflow[p] * std::exp(-r * (exercise_t[p] - ti));  // discount to t_i
            for (int a = 0; a < 3; ++a) {
                for (int c = 0; c < 3; ++c) A[a][c] += basis[a] * basis[c];
                rhs[a] += basis[a] * y;
            }
        }
        double coef[3];
        if (!any || !detail::solve3x3(A, rhs, coef)) continue;  // keep holding at this date
        for (long p = 0; p < n_paths; ++p) {
            const double s = path[static_cast<std::size_t>(p) * m + idx];
            const double ex = vanilla_payoff(type, s, K);
            if (ex <= 0.0) continue;
            const double x = s / K;
            const double cont = coef[0] + coef[1] * x + coef[2] * x * x;
            if (ex > cont) {  // exercising now beats the estimated continuation value
                cashflow[p] = ex;
                exercise_t[p] = ti;
            }
        }
    }

    // Value today: average of each path's cash flow discounted to t = 0.
    double sum = 0.0;
    for (long p = 0; p < n_paths; ++p) sum += cashflow[p] * std::exp(-r * exercise_t[p]);
    return sum / static_cast<double>(n_paths);
}

}  // namespace pricer
