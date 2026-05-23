// pricer/adjoint.hpp — reverse-mode automatic differentiation (AAD).
//
// A tape records each operation's parents and local partial derivatives during
// the forward pass; one backward sweep then propagates the output adjoint to
// every input. For a single output (the price) and many inputs (S, sigma, r, T),
// this yields all first-order Greeks in ONE reverse pass — the efficiency that
// makes adjoint AD the standard for production risk.
//
// The active type `Var` overloads the same operators/functions as the dual
// numbers, so the generic `bs_price_ad` in greeks_ad.hpp is reused verbatim.
#pragma once
#include <cmath>
#include <vector>

#include "pricer/black_scholes.hpp"  // OptionType, Greeks
#include "pricer/greeks_ad.hpp"      // bs_price_ad<Num>, black_scholes_greeks_ad (for gamma)

namespace pricer {

class Tape {
public:
    struct Node { int p0, p1; double w0, w1; };  // up to two parents + local partials

    int leaf() { nodes_.push_back({-1, -1, 0.0, 0.0}); return last(); }
    int unary(int a, double wa) { nodes_.push_back({a, -1, wa, 0.0}); return last(); }
    int binary(int a, double wa, int b, double wb) { nodes_.push_back({a, b, wa, wb}); return last(); }

    // Reverse sweep from output `out`; returns the adjoint of every node.
    std::vector<double> grad(int out) const {
        std::vector<double> adj(nodes_.size(), 0.0);
        if (out >= 0) adj[out] = 1.0;
        for (int i = static_cast<int>(nodes_.size()) - 1; i >= 0; --i) {
            const Node& n = nodes_[i];
            if (n.p0 >= 0) adj[n.p0] += adj[i] * n.w0;
            if (n.p1 >= 0) adj[n.p1] += adj[i] * n.w1;
        }
        return adj;
    }

private:
    int last() const { return static_cast<int>(nodes_.size()) - 1; }
    std::vector<Node> nodes_;
};

// Active variable: a value plus its node index on a tape.
struct Var {
    double v;
    int idx;
    Tape* t;
};

inline Var make_var(Tape& t, double x) { return {x, t.leaf(), &t}; }

// --- operators (record value + local partials) ---
inline Var operator+(const Var& a, const Var& b) { return {a.v + b.v, a.t->binary(a.idx, 1.0, b.idx, 1.0), a.t}; }
inline Var operator-(const Var& a, const Var& b) { return {a.v - b.v, a.t->binary(a.idx, 1.0, b.idx, -1.0), a.t}; }
inline Var operator*(const Var& a, const Var& b) { return {a.v * b.v, a.t->binary(a.idx, b.v, b.idx, a.v), a.t}; }
inline Var operator/(const Var& a, const Var& b) {
    const double inv = 1.0 / b.v;
    return {a.v * inv, a.t->binary(a.idx, inv, b.idx, -a.v * inv * inv), a.t};
}
inline Var operator-(const Var& a) { return {-a.v, a.t->unary(a.idx, -1.0), a.t}; }
inline Var operator*(double s, const Var& a) { return {s * a.v, a.t->unary(a.idx, s), a.t}; }
inline Var operator*(const Var& a, double s) { return s * a; }
inline Var operator+(const Var& a, double s) { return {a.v + s, a.t->unary(a.idx, 1.0), a.t}; }
inline Var operator+(double s, const Var& a) { return a + s; }
inline Var operator-(const Var& a, double s) { return {a.v - s, a.t->unary(a.idx, 1.0), a.t}; }
inline Var operator-(double s, const Var& a) { return {s - a.v, a.t->unary(a.idx, -1.0), a.t}; }

// --- math functions ---
inline Var exp(const Var& a) { const double e = std::exp(a.v); return {e, a.t->unary(a.idx, e), a.t}; }
inline Var log(const Var& a) { return {std::log(a.v), a.t->unary(a.idx, 1.0 / a.v), a.t}; }
inline Var sqrt(const Var& a) { const double s = std::sqrt(a.v); return {s, a.t->unary(a.idx, 0.5 / s), a.t}; }
inline Var erfc(const Var& a) {
    const double val = std::erfc(a.v);
    const double d = -1.1283791670955126 * std::exp(-a.v * a.v);  // d/dx erfc = -2/√π e^{-x²}
    return {val, a.t->unary(a.idx, d), a.t};
}

// Greeks via reverse-mode AAD: delta/vega/rho/theta from one backward sweep.
// (Gamma, a second derivative, is supplied by forward-mode AD.)
inline Greeks black_scholes_greeks_aad(OptionType type, double S, double K, double r, double sigma,
                                       double T, double q = 0.0) {
    Tape tape;
    const Var Sv = make_var(tape, S), Kv = make_var(tape, K), rv = make_var(tape, r),
              sigv = make_var(tape, sigma), Tv = make_var(tape, T), qv = make_var(tape, q);
    const Var price = bs_price_ad<Var>(type, Sv, Kv, rv, sigv, Tv, qv);
    const std::vector<double> adj = tape.grad(price.idx);

    Greeks g{};
    g.price = price.v;
    g.delta = adj[Sv.idx];
    g.vega = adj[sigv.idx];
    g.rho = adj[rv.idx];
    g.theta = -adj[Tv.idx];  // time decay
    g.gamma = black_scholes_greeks_ad(type, S, K, r, sigma, T, q).gamma;  // 2nd order via forward AD
    return g;
}

}  // namespace pricer
