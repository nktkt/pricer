// Tests for book-level Greeks in one reverse-mode AAD sweep (portfolio.hpp).
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/portfolio.hpp"

#include <cstddef>
#include <vector>

using namespace pricer;

int main() {
    // A small multi-name book: different underlyings, strikes, vols, maturities,
    // and signed quantities (long and short).
    const std::vector<Position> book = {
        {OptionType::Call, 100.0, 100.0, 0.05, 0.20, 1.0, 10.0},
        {OptionType::Put, 90.0, 95.0, 0.05, 0.25, 0.5, -5.0},
        {OptionType::Call, 120.0, 110.0, 0.03, 0.18, 2.0, 7.0},
    };

    const BookGreeks g = book_greeks_aad(book);

    // Book value = Σ qty · BS price.
    double value = 0.0;
    for (const Position& p : book)
        value += p.qty * black_scholes_price(p.type, p.S, p.K, p.r, p.sigma, p.T);
    check::approx("book value", g.value, value, 1e-9);

    // The single reverse sweep reproduces each position's analytic delta/vega
    // (scaled by quantity) — full book risk in one pass.
    for (std::size_t k = 0; k < book.size(); ++k) {
        const Position& p = book[k];
        const Greeks gr = black_scholes_greeks(p.type, p.S, p.K, p.r, p.sigma, p.T);
        check::approx("one-pass delta matches analytic", g.delta[k], p.qty * gr.delta, 1e-7);
        check::approx("one-pass vega matches analytic", g.vega[k], p.qty * gr.vega, 1e-7);
    }

    return check::report("portfolio");
}
