// pricer/qmc.hpp — quasi-Monte Carlo (low-discrepancy) terminal pricing.
//
// Instead of pseudo-random points, QMC uses a deterministic low-discrepancy
// sequence that fills space more evenly, mapped to normals through the inverse
// CDF. For low-dimensional problems the error decays roughly like 1/N instead of
// 1/sqrt(N), so far fewer points reach a target accuracy.
#pragma once
#include <cmath>
#include <cstdint>

namespace pricer {

// Inverse standard-normal CDF — Acklam's rational approximation (~1e-9 abs err).
inline double inv_norm_cdf(double p) {
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
                               -2.759285104469687e+02, 1.383577518672690e+02,
                               -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
                               -1.556989798598866e+02, 6.680131188771972e+01,
                               -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                               -2.400758277161838e+00, -2.549732539343734e+00,
                               4.374664141464968e+00,  2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
                               2.445134137142996e+00, 3.754408661907416e+00};
    const double plow = 0.02425, phigh = 1.0 - plow;
    if (p < plow) {
        const double q = std::sqrt(-2 * std::log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    } else if (p <= phigh) {
        const double q = p - 0.5, rr = q * q;
        return (((((a[0] * rr + a[1]) * rr + a[2]) * rr + a[3]) * rr + a[4]) * rr + a[5]) * q /
               (((((b[0] * rr + b[1]) * rr + b[2]) * rr + b[3]) * rr + b[4]) * rr + 1);
    } else {
        const double q = std::sqrt(-2 * std::log(1 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
}

namespace mc {

// Van der Corput radical inverse in base 2 — the 1-D low-discrepancy sequence.
inline double radical_inverse_base2(std::uint64_t i) {
    double f = 0.0, base = 0.5;
    while (i) { f += (i & 1u) * base; i >>= 1; base *= 0.5; }
    return f;
}

// 1-D quasi-Monte Carlo: low-discrepancy points -> inverse-CDF normals -> price.
template <class Payoff>
double price_terminal_qmc(Payoff payoff, double S, double r, double sigma, double T, long n) {
    const double drift = (r - 0.5 * sigma * sigma) * T, vol = sigma * std::sqrt(T);
    double sum = 0.0;
    for (long i = 1; i <= n; ++i) {  // start at 1 so the point is never 0
        const double u = radical_inverse_base2(static_cast<std::uint64_t>(i));
        const double z = inv_norm_cdf(u);
        const double ST = S * std::exp(drift + vol * z);
        sum += payoff(ST);
    }
    return std::exp(-r * T) * (sum / static_cast<double>(n));
}

}  // namespace mc
}  // namespace pricer
