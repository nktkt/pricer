// variance_reduction.cpp — reach a target accuracy with far fewer paths.
//
// A variance-reduction technique is judged by its *error spread* over many runs,
// not a single lucky/unlucky seed. So this demo estimates the RMSE (root-mean-
// square error vs. the analytic price) of each estimator across many seeds, at a
// fixed payoff-evaluation budget per run. Quasi-Monte Carlo is deterministic, so
// it has a single error.
//
// The "variance-reduction factor" = (RMSE_plain / RMSE_method)^2 is how many
// times more plain-MC paths you would need to match that method's accuracy.
#include "pricer/black_scholes.hpp"
#include "pricer/monte_carlo.hpp"
#include "pricer/qmc.hpp"
#include "pricer/variance_reduction.hpp"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <initializer_list>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    auto call = [K](double ST) { return ST > K ? ST - K : 0.0; };
    const double bs = black_scholes_call(S, K, r, sigma, T);
    const int seeds = 64;  // independent runs used to estimate RMSE

    std::printf("analytic (Black-Scholes) = %.6f   (RMSE over %d seeds)\n\n", bs, seeds);
    std::printf("%10s | %12s | %12s | %12s | %12s\n",
                "N", "plain", "antithetic", "control", "qmc(1 run)");
    std::printf("-----------|--------------|--------------|--------------|--------------\n");

    for (long N : {10000L, 100000L, 1000000L}) {
        double se_plain = 0, se_anti = 0, se_ctrl = 0;
        for (int k = 0; k < seeds; ++k) {
            const std::uint64_t sd = 1000 + k;
            const double ep = mc::price_terminal(call, S, r, sigma, T, N, sd) - bs;
            const double ea = mc::price_terminal_antithetic(call, S, r, sigma, T, N / 2, sd) - bs;
            const double ec = mc::price_terminal_control(call, S, r, sigma, T, N, sd) - bs;
            se_plain += ep * ep; se_anti += ea * ea; se_ctrl += ec * ec;
        }
        const double rp = std::sqrt(se_plain / seeds);
        const double ra = std::sqrt(se_anti / seeds);
        const double rc = std::sqrt(se_ctrl / seeds);
        const double qe = std::fabs(mc::price_terminal_qmc(call, S, r, sigma, T, N) - bs);
        std::printf("%10ld | %12.6f | %12.6f | %12.6f | %12.6f\n", N, rp, ra, rc, qe);
        if (N == 1000000L)
            std::printf("\nvariance-reduction factor at N=1e6 (plain/method)^2: "
                        "antithetic %.1fx, control %.1fx, qmc %.0fx\n",
                        (rp / ra) * (rp / ra), (rp / rc) * (rp / rc), (rp / qe) * (rp / qe));
    }

    std::printf("\nSame budget, far less error: control variate and QMC reach the analytic\n"
                "price with much smaller RMSE, i.e. equivalent accuracy from fewer paths.\n");
    return 0;
}
