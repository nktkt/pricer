// pricer/local_vol.hpp — Dupire local volatility from a call-price surface.
//
// Dupire's formula recovers the local volatility consistent with a surface of
// European call prices C(K, T) (zero dividends):
//                 ∂C/∂T + r·K·∂C/∂K
//   σ_loc²(K,T) = ----------------------
//                    ½·K²·∂²C/∂K²
// The derivatives are taken by central finite differences of any price function
// C(K, T). For a flat implied-vol surface this returns that constant vol.
#pragma once
#include <algorithm>
#include <cmath>

namespace pricer {

// `price(K, T)` returns the European call price at strike K, expiry T.
template <class PriceFn>
double dupire_local_vol(PriceFn price, double K, double T, double r, double hK = 1.0,
                        double hT = 1e-3) {
    const double dC_dT = (price(K, T + hT) - price(K, T - hT)) / (2 * hT);
    const double dC_dK = (price(K + hK, T) - price(K - hK, T)) / (2 * hK);
    const double d2C_dK2 = (price(K + hK, T) - 2 * price(K, T) + price(K - hK, T)) / (hK * hK);

    const double numerator = dC_dT + r * K * dC_dK;
    const double denominator = 0.5 * K * K * d2C_dK2;
    if (denominator <= 0.0) return 0.0;
    return std::sqrt(std::max(numerator / denominator, 0.0));
}

}  // namespace pricer
