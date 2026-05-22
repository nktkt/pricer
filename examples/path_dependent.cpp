// path_dependent.cpp — price several path-dependent products from formulas alone.
//
// The Monte Carlo engine simulates a full price path and exposes four aggregates
// to the payoff DSL as variables:
//   ST    terminal price
//   avg   arithmetic average over the observation dates
//   Smax  running maximum
//   Smin  running minimum
// plus the strike K. With those, whole families of exotics are one-liners:
//
//   max(ST - K, 0)                  vanilla call (sanity check vs Black–Scholes)
//   max(avg - K, 0)                 arithmetic Asian call (cheaper: averaging cuts vol)
//   max(ST - K, 0) * (Smax < 130)   up-and-out barrier call (knocked out above 130)
//   ST - Smin                       floating-strike lookback call
//   (ST > K) * 10                   cash-or-nothing digital
//
// Nothing about these instruments is compiled into C++: each formula is turned
// into native code at runtime by pricer::PayoffJit.
#include "pricer/black_scholes.hpp"
#include "pricer/payoff_jit.hpp"
#include <cstdio>
#include <random>
#include <vector>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long n_paths = 2'000'000;
    const int  n_steps = 250;  // ~daily observations

    const double dt = T / n_steps;
    const double drift = (r - 0.5 * sigma * sigma) * dt;
    const double vol = sigma * std::sqrt(dt);
    const double disc = std::exp(-r * T);

    // Variable layout shared by every compiled payoff.
    const std::vector<std::string> vars = {"ST", "avg", "Smax", "Smin", "K"};

    struct Product { const char* name; const char* formula; };
    const Product products[] = {
        {"vanilla call",        "max(ST - K, 0)"},
        {"Asian call (avg)",    "max(avg - K, 0)"},
        {"up-and-out @130",     "max(ST - K, 0) * (Smax < 130)"},
        {"lookback (float K)",  "ST - Smin"},
        {"digital cash-or-not", "(ST > K) * 10"},
    };

    PayoffJit jit;
    std::vector<PayoffJit::Fn> fns;
    for (const auto& p : products) fns.push_back(jit.compile(p.formula, vars));

    // One shared set of paths prices all products simultaneously.
    std::mt19937_64 rng(2024);
    std::normal_distribution<double> Z(0.0, 1.0);
    std::vector<double> sums(std::size(products), 0.0);

    for (long p = 0; p < n_paths; ++p) {
        double price = S, sum = 0.0, smax = S, smin = S;
        for (int s = 0; s < n_steps; ++s) {
            price *= std::exp(drift + vol * Z(rng));
            sum += price;
            if (price > smax) smax = price;
            if (price < smin) smin = price;
        }
        double v[5] = {price, sum / n_steps, smax, smin, K};  // ST, avg, Smax, Smin, K
        for (size_t i = 0; i < fns.size(); ++i) sums[i] += fns[i](v);
    }

    std::printf("S=%.0f K=%.0f r=%.2f sigma=%.2f T=%.0f  paths=%ld steps=%d\n\n",
                S, K, r, sigma, T, n_paths, n_steps);
    std::printf("%-22s | %-32s | %10s\n", "product", "formula", "price");
    std::printf("-----------------------|----------------------------------|-----------\n");
    for (size_t i = 0; i < std::size(products); ++i)
        std::printf("%-22s | %-32s | %10.6f\n",
                    products[i].name, products[i].formula, disc * sums[i] / n_paths);

    std::printf("\nreference: Black-Scholes vanilla call = %.6f\n",
                black_scholes_call(S, K, r, sigma, T));
    std::printf("\nEach price above came from a formula string compiled to native code.\n");
    return 0;
}
