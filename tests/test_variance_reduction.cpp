// Tests for variance-reduction estimators and quasi-Monte Carlo.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/qmc.hpp"
#include "pricer/variance_reduction.hpp"

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    auto call = [K](double ST) { return ST > K ? ST - K : 0.0; };
    const double bs = black_scholes_call(S, K, r, sigma, T);

    const double anti = mc::price_terminal_antithetic(call, S, r, sigma, T, 1'000'000, 7);
    const double ctrl = mc::price_terminal_control(call, S, r, sigma, T, 1'000'000, 7);
    const double qmc  = mc::price_terminal_qmc(call, S, r, sigma, T, 1'000'000);

    check::approx("antithetic vs BS",      anti, bs, 0.02);
    check::approx("control variate vs BS", ctrl, bs, 0.01);
    check::approx("QMC vs BS",             qmc,  bs, 0.005);

    // Inverse-normal sanity: inv_norm_cdf is the inverse of the CDF.
    check::approx("inv_norm_cdf(0.5)",   inv_norm_cdf(0.5),   0.0,       1e-9);
    check::approx("inv_norm_cdf(0.975)", inv_norm_cdf(0.975), 1.959964,  1e-4);

    return check::report("variance_reduction");
}
