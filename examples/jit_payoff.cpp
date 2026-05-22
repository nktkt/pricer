// jit_payoff.cpp — compile a terminal payoff formula to native code at runtime.
//
// Thin demo over pricer::PayoffJit: parse a formula in the variables {ST, K},
// JIT-compile it, print the generated LLVM IR, then price it by Monte Carlo and
// (for the default call) compare with the analytic Black–Scholes value.
//
// Change the formula and you price a different instrument — no C++ recompile:
//   ./jit_payoff                              # call:  max(ST - K, 0)
//   ./jit_payoff "max(K - ST, 0)"            # put
//   ./jit_payoff "max(ST-K,0)+max(K-ST,0)"   # straddle
//   ./jit_payoff "(ST > K) * 10"             # cash-or-nothing digital
#include "pricer/black_scholes.hpp"
#include "pricer/payoff_jit.hpp"
#include <chrono>
#include <cstdio>
#include <random>

using namespace pricer;

int main(int argc, char** argv) {
    const std::string expr = (argc > 1) ? argv[1] : "max(ST - K, 0)";

    PayoffJit jit;
    PayoffJit::Fn payoff;
    try {
        payoff = jit.compile(expr, {"ST", "K"});  // v[0]=ST, v[1]=K
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }

    std::printf("formula : \"%s\"\n\n--- generated LLVM IR ---\n%s-------------------------\n\n",
                expr.c_str(), jit.last_ir().c_str());

    auto call = [&](double ST, double K) { double v[2] = {ST, K}; return payoff(v); };
    std::printf("JIT-compiled. payoff(ST=120, K=100) = %.4f\n", call(120.0, 100.0));

    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long n = 10'000'000;
    std::mt19937_64 rng(777);
    std::normal_distribution<double> Z(0.0, 1.0);
    const double drift = (r - 0.5 * sigma * sigma) * T, vol = sigma * std::sqrt(T);

    const auto t0 = std::chrono::high_resolution_clock::now();
    double sum = 0.0;
    for (long i = 0; i < n; ++i) {
        const double ST = S * std::exp(drift + vol * Z(rng));
        sum += call(ST, K);
    }
    const double price = std::exp(-r * T) * (sum / n);
    const auto t1 = std::chrono::high_resolution_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::printf("\nMonte Carlo price (%ld paths) = %.6f  (%.1f ms)\n", n, price, ms);
    if (expr == "max(ST - K, 0)")
        std::printf("analytic (BS call)            = %.6f  <- should match\n",
                    black_scholes_call(S, K, r, sigma, T));
    return 0;
}
