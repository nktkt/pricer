// pricer/portfolio.hpp — book-level risk in one reverse-mode AAD sweep.
//
// A trading book holds many positions on many underlyings. Reverse-mode adjoint
// AD computes the gradient of a single scalar output with respect to *all* inputs
// in one backward pass, at a cost independent of the number of inputs — so the
// whole book's value and every position's delta and vega come out of one sweep
// over one tape. This is the property that makes adjoint AD the standard for
// production risk on large books.
#pragma once
#include <cstddef>
#include <vector>

#include "pricer/adjoint.hpp"        // Tape, Var, make_var
#include "pricer/black_scholes.hpp"  // OptionType
#include "pricer/greeks_ad.hpp"      // bs_price_ad<Num>

namespace pricer {

// A vanilla option position in a multi-name book.
struct Position {
    OptionType type;
    double S, K, r, sigma, T;
    double qty = 1.0;
};

// Book value and per-position first-order risk, all from ONE backward sweep.
struct BookGreeks {
    double value = 0.0;
    std::vector<double> delta;  // dValue/dS_k    (already scaled by qty_k)
    std::vector<double> vega;   // dValue/dsigma_k
};

// Price the whole book V = Σ qty_k · BS(S_k, K_k, r_k, sigma_k, T_k) on a single
// tape, then one reverse sweep yields the book's sensitivity to every position's
// spot and volatility simultaneously.
inline BookGreeks book_greeks_aad(const std::vector<Position>& book) {
    Tape tape;
    std::vector<Var> Sv, sigv;
    Sv.reserve(book.size());
    sigv.reserve(book.size());

    Var value = make_var(tape, 0.0);
    for (const Position& p : book) {
        const Var S = make_var(tape, p.S);
        const Var K = make_var(tape, p.K);
        const Var r = make_var(tape, p.r);
        const Var sg = make_var(tape, p.sigma);
        const Var T = make_var(tape, p.T);
        value = value + p.qty * bs_price_ad<Var>(p.type, S, K, r, sg, T);
        Sv.push_back(S);
        sigv.push_back(sg);
    }

    const std::vector<double> adj = tape.grad(value.idx);
    BookGreeks g;
    g.value = value.v;
    g.delta.resize(book.size());
    g.vega.resize(book.size());
    for (std::size_t k = 0; k < book.size(); ++k) {
        g.delta[k] = adj[Sv[k].idx];
        g.vega[k] = adj[sigv[k].idx];
    }
    return g;
}

}  // namespace pricer
