// portfolio_aad.cpp — full book risk for a multi-name book in ONE reverse sweep.
//
// Reverse-mode adjoint AD computes the gradient of one scalar (the book value)
// with respect to every input in a single backward pass, at a cost independent
// of the number of inputs. So a book of options on several underlyings yields
// its value and every position's delta and vega from one tape sweep — the
// efficiency that makes adjoint AD the production standard for large books. Here
// we cross-check the one-pass risk against the closed-form per-position Greeks.
#include "pricer/black_scholes.hpp"
#include "pricer/portfolio.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

using namespace pricer;

int main() {
    const std::vector<Position> book = {
        {OptionType::Call, 100.0, 100.0, 0.05, 0.20, 1.0, 10.0},  // long 10 ATM calls
        {OptionType::Put, 90.0, 95.0, 0.05, 0.25, 0.5, -5.0},     // short 5 puts
        {OptionType::Call, 120.0, 110.0, 0.03, 0.18, 2.0, 7.0},   // long 7 ITM calls
    };

    const BookGreeks g = book_greeks_aad(book);

    std::printf("Multi-name book: full first-order risk from ONE reverse-mode AAD sweep\n\n");
    std::printf("Book value = %.4f\n\n", g.value);
    std::printf("%-4s  %-5s  %8s  %6s  %14s  %14s\n", "pos", "type", "S", "qty",
                "delta (1-pass)", "vega (1-pass)");
    std::printf("--------------------------------------------------------------------\n");
    for (std::size_t k = 0; k < book.size(); ++k) {
        const Position& p = book[k];
        std::printf("%-4zu  %-5s  %8.1f  %6.0f  %14.5f  %14.5f\n", k,
                    p.type == OptionType::Call ? "call" : "put", p.S, p.qty, g.delta[k], g.vega[k]);
    }

    // Cross-check: the one-pass gradient equals the closed-form per-position
    // Greeks scaled by quantity.
    double max_err = 0.0;
    for (std::size_t k = 0; k < book.size(); ++k) {
        const Position& p = book[k];
        const Greeks gr = black_scholes_greeks(p.type, p.S, p.K, p.r, p.sigma, p.T);
        const double de = std::fabs(g.delta[k] - p.qty * gr.delta);
        const double ve = std::fabs(g.vega[k] - p.qty * gr.vega);
        if (de > max_err) max_err = de;
        if (ve > max_err) max_err = ve;
    }
    std::printf("\nMax error vs. closed-form delta/vega = %.2e (all from one tape sweep)\n", max_err);
    return 0;
}
