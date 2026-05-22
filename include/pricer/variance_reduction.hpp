// pricer/variance_reduction.hpp — same accuracy from fewer paths.
//
// Two classic estimators that cut Monte Carlo variance without bias:
//   * Antithetic variates: evaluate each draw Z and its mirror -Z, cancelling
//     the odd part of the error.
//   * Control variate: subtract a correlated quantity with a known mean. Here
//     the control is the discounted terminal price disc*S_T, whose risk-neutral
//     expectation is exactly S; the optimal coefficient beta is estimated from
//     the sample in a single pass.
#pragma once
#include <cmath>
#include <cstdint>
#include <random>

namespace pricer::mc {

// Antithetic-variates estimator over `n_pairs` mirrored draws (2*n_pairs paths).
template <class Payoff>
double price_terminal_antithetic(Payoff payoff, double S, double r, double sigma, double T,
                                 long n_pairs, std::uint64_t seed = 12345) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);
    const double drift = (r - 0.5 * sigma * sigma) * T, vol = sigma * std::sqrt(T);
    double sum = 0.0;
    for (long i = 0; i < n_pairs; ++i) {
        const double z = Z(rng);
        const double up = S * std::exp(drift + vol * z);
        const double dn = S * std::exp(drift - vol * z);
        sum += 0.5 * (payoff(up) + payoff(dn));
    }
    return std::exp(-r * T) * (sum / static_cast<double>(n_pairs));
}

// Control-variate estimator using disc*S_T (mean S). Returns the discounted,
// variance-reduced price. Beta is estimated from the same n samples.
template <class Payoff>
double price_terminal_control(Payoff payoff, double S, double r, double sigma, double T,
                              long n, std::uint64_t seed = 12345) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);
    const double drift = (r - 0.5 * sigma * sigma) * T, vol = sigma * std::sqrt(T);
    const double disc = std::exp(-r * T);

    // Y = discounted payoff, X = discounted terminal price (E[X] = S).
    double sumY = 0, sumX = 0, sumXX = 0, sumXY = 0;
    for (long i = 0; i < n; ++i) {
        const double ST = S * std::exp(drift + vol * Z(rng));
        const double Y = disc * payoff(ST);
        const double X = disc * ST;
        sumY += Y; sumX += X; sumXX += X * X; sumXY += X * Y;
    }
    const double dn = static_cast<double>(n);
    const double covXY = sumXY - sumX * sumY / dn;
    const double varX = sumXX - sumX * sumX / dn;
    const double beta = (varX > 0) ? covXY / varX : 0.0;
    // mean(Y) - beta * (mean(X) - E[X]),  E[X] = S
    return sumY / dn - beta * (sumX / dn - S);
}

}  // namespace pricer::mc
