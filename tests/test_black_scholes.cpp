// Tests for the closed-form Black–Scholes pricer.
#include "check.hpp"
#include "pricer/black_scholes.hpp"

using namespace pricer;

int main() {
    // Textbook reference: S=K=100, r=5%, sigma=20%, T=1y.
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;

    const double call = black_scholes_call(S, K, r, sigma, T);
    const double put  = black_scholes_put(S, K, r, sigma, T);

    check::approx("call value", call, 10.4506, 1e-3);
    check::approx("put value",  put,   5.5735, 1e-3);

    // Put–call parity: C - P = S - K*e^{-rT}
    const double parity_lhs = call - put;
    const double parity_rhs = S - K * std::exp(-r * T);
    check::approx("put-call parity", parity_lhs, parity_rhs, 1e-9);

    // Deep in-the-money call ~ S - K*e^{-rT}; deep OTM call ~ 0.
    check::approx("deep ITM call", black_scholes_call(200, 100, r, sigma, T),
                  200 - 100 * std::exp(-r * T), 1e-2);
    check::is_true("deep OTM call ~ 0",
                   black_scholes_call(10, 100, r, sigma, T) < 1e-6);

    return check::report("black_scholes");
}
