// pricer/greeks_mc.hpp — Monte Carlo Greeks.
//
// Two approaches:
//   * Bump-and-revalue with COMMON RANDOM NUMBERS (CRN): reprice with the same
//     seed at S±h (or sigma±h) so the random draws cancel and only the
//     sensitivity remains — far lower variance than independent bumps. Generic:
//     works for any payoff callable.
//   * Pathwise (infinitesimal-perturbation) derivative: differentiate the payoff
//     along each path. This is exactly what adjoint AD computes; here it is shown
//     in closed form for the vanilla call (delta = disc * E[1_{S_T>K} * S_T/S]).
#pragma once
#include <cmath>
#include <cstdint>
#include <random>
#include "pricer/monte_carlo.hpp"

namespace pricer::mc {

// Central-difference delta with common random numbers.
template <class Payoff>
double mc_delta_crn(Payoff payoff, double S, double r, double sigma, double T, long n,
                    double h = 1e-2, std::uint64_t seed = 12345) {
    const double up = price_terminal(payoff, S + h, r, sigma, T, n, seed);
    const double dn = price_terminal(payoff, S - h, r, sigma, T, n, seed);
    return (up - dn) / (2 * h);
}

// Central-difference vega with common random numbers (per 1.0 of vol).
template <class Payoff>
double mc_vega_crn(Payoff payoff, double S, double r, double sigma, double T, long n,
                   double h = 1e-4, std::uint64_t seed = 12345) {
    const double up = price_terminal(payoff, S, r, sigma + h, T, n, seed);
    const double dn = price_terminal(payoff, S, r, sigma - h, T, n, seed);
    return (up - dn) / (2 * h);
}

// Pathwise delta of a European call (no bumping; one pass).
inline double mc_call_delta_pathwise(double S, double K, double r, double sigma, double T,
                                     long n, std::uint64_t seed = 12345) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);
    const double drift = (r - 0.5 * sigma * sigma) * T, vol = sigma * std::sqrt(T);
    double sum = 0.0;
    for (long i = 0; i < n; ++i) {
        const double ST = S * std::exp(drift + vol * Z(rng));
        if (ST > K) sum += ST / S;  // d/dS max(ST-K,0) = 1_{ST>K} * ST/S
    }
    return std::exp(-r * T) * (sum / static_cast<double>(n));
}

}  // namespace pricer::mc
