// pricer/greeks_ad.hpp — exact Black–Scholes Greeks by forward-mode AD.
//
// The Black–Scholes price is written once as a generic function of any numeric
// type; evaluating it in dual numbers yields exact derivatives (no bumping):
//   * one Dual<4> pass gives delta, vega, rho and dPrice/dT (theta) together;
//   * a nested Dual<1,Dual<1>> pass in S gives gamma (a second derivative).
#pragma once
#include <cmath>

#include "pricer/black_scholes.hpp"  // OptionType, Greeks
#include "pricer/dual.hpp"

namespace pricer {

// Standard-normal CDF generic over the numeric type (uses erfc).
template <class Num>
Num norm_cdf_ad(const Num& x) {
    using std::erfc;
    return 0.5 * erfc((-1.0 / std::sqrt(2.0)) * x);
}

// Black–Scholes price generic over the numeric type Num (double or Dual<...>),
// with a continuous dividend yield q (q = 0 is no-dividend Black–Scholes).
template <class Num>
Num bs_price_ad(OptionType type, const Num& S, const Num& K, const Num& r, const Num& sigma,
                const Num& T, const Num& q) {
    using std::sqrt;
    const Num sqrtT = sqrt(T);
    const Num d1 = (log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const Num d2 = d1 - sigma * sqrtT;
    const Num disc = exp((-1.0) * r * T);
    const Num discq = exp((-1.0) * q * T);
    if (type == OptionType::Call)
        return S * discq * norm_cdf_ad(d1) - K * disc * norm_cdf_ad(d2);
    return K * disc * norm_cdf_ad((-1.0) * d2) - S * discq * norm_cdf_ad((-1.0) * d1);
}

// Greeks via automatic differentiation; matches the closed-form values exactly.
// q is carried as a constant (delta/vega/rho/theta reflect it; the dividend
// sensitivity itself is not reported).
inline Greeks black_scholes_greeks_ad(OptionType type, double S, double K, double r, double sigma,
                                      double T, double q = 0.0) {
    // First order: seed S, sigma, r, T (K and q constant) in one Dual<4> evaluation.
    using D = Dual<4, double>;
    const D price = bs_price_ad<D>(type, ad_var<4>(S, 0), ad_const<4>(K), ad_var<4>(r, 2),
                                   ad_var<4>(sigma, 1), ad_var<4>(T, 3), ad_const<4>(q));
    Greeks g{};
    g.price = price.v;
    g.delta = price.d[0];
    g.vega = price.d[1];
    g.rho = price.d[2];
    g.theta = -price.d[3];  // time decay = -∂Price/∂T

    // Second order in S for gamma: nested dual, S active at both levels.
    using D1 = Dual<1, double>;
    using DD = Dual<1, D1>;
    auto dd_const = [](double x) { DD c; c.v = D1(x); c.d[0] = D1(0.0); return c; };
    DD Snn;                       // value carries ∂/∂S, outer carries ∂²/∂S²
    Snn.v = D1(S); Snn.v.d[0] = 1.0;
    Snn.d[0] = D1(1.0);          // (its own d[0] stays 0)
    const DD pr = bs_price_ad<DD>(type, Snn, dd_const(K), dd_const(r), dd_const(sigma),
                                  dd_const(T), dd_const(q));
    g.gamma = pr.d[0].d[0];
    return g;
}

}  // namespace pricer
