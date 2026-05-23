// Tests for rainbow (two-asset best-of / worst-of) options (rainbow.hpp). The
// bivariate normal CDF is checked against its analytic special cases; the Stulz
// formulas are pinned by the call-on-max + call-on-min == two vanilla calls
// parity and cross-checked against a correlated two-asset Monte Carlo.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/rainbow.hpp"

using namespace pricer;

int main() {
    // === Bivariate normal CDF special cases ===
    check::approx("M(0,0;0) == 1/4", bivariate_normal_cdf(0, 0, 0), 0.25, 1e-6);
    // M(0,0;rho) = 1/4 + asin(rho)/(2*pi); for rho=0.5 this is 1/3.
    check::approx("M(0,0;0.5) == 1/3", bivariate_normal_cdf(0, 0, 0.5), 1.0 / 3.0, 1e-6);
    check::approx("M(0,0;-0.5)", bivariate_normal_cdf(0, 0, -0.5), 0.25 - std::asin(0.5) / (2 * M_PI),
                  1e-6);
    // Independence: M(a,b;0) = N(a) N(b).
    check::approx("M(a,b;0) == N(a)N(b)", bivariate_normal_cdf(0.3, -0.7, 0.0),
                  norm_cdf(0.3) * norm_cdf(-0.7), 1e-6);
    // Perfect correlation: M(a,b;1) = N(min(a,b)); M(a,b;-1) = max(0, N(a)+N(b)-1).
    check::approx("M(a,b;1) == N(min)", bivariate_normal_cdf(0.4, 0.9, 0.999999),
                  norm_cdf(0.4), 1e-4);
    check::approx("M(a,b;-1) == max(0,N(a)+N(b)-1)", bivariate_normal_cdf(0.4, 0.9, -0.999999),
                  std::max(0.0, norm_cdf(0.4) + norm_cdf(0.9) - 1.0), 1e-4);
    // Marginal: M(a, +inf; rho) = N(a).
    check::approx("M(a,+inf;rho) == N(a)", bivariate_normal_cdf(0.3, 10.0, 0.5), norm_cdf(0.3), 1e-6);

    // === Rainbow options ===
    const double S1 = 100, S2 = 95, K = 100, r = 0.05, sig1 = 0.20, sig2 = 0.30, rho = 0.4, T = 1.0;

    const double cmax = rainbow_price(OptionType::Call, RainbowType::Max, S1, S2, K, r, sig1, sig2, rho, T);
    const double cmin = rainbow_price(OptionType::Call, RainbowType::Min, S1, S2, K, r, sig1, sig2, rho, T);
    const double c1 = black_scholes_call(S1, K, r, sig1, T);
    const double c2 = black_scholes_call(S2, K, r, sig2, T);

    // Parity that pins both Stulz call formulas: max-call + min-call == two vanillas.
    check::approx("Cmax + Cmin == c1 + c2", cmax + cmin, c1 + c2, 1e-9);

    // A call on the best-of dominates either single call; the worst-of is dominated.
    check::is_true("Cmax >= max(c1,c2)", cmax >= std::max(c1, c2) - 1e-9);
    check::is_true("Cmin <= min(c1,c2)", cmin <= std::min(c1, c2) + 1e-9);
    check::is_true("Cmax >= Cmin", cmax >= cmin);

    // As rho -> 1 with identical assets, max and min coincide with the single call.
    // (The gap shrinks like sqrt(1-rho), so rho must be very close to 1.)
    const double cmax_hi = rainbow_price(OptionType::Call, RainbowType::Max, 100, 100, K, r, 0.2, 0.2, 0.99999, T);
    const double cmin_hi = rainbow_price(OptionType::Call, RainbowType::Min, 100, 100, K, r, 0.2, 0.2, 0.99999, T);
    const double c_one = black_scholes_call(100, K, r, 0.2, T);
    check::approx("rho->1: Cmax ~ vanilla", cmax_hi, c_one, 0.05);
    check::approx("rho->1: Cmin ~ vanilla", cmin_hi, c_one, 0.05);

    // Dispersion: lower correlation widens the gap (max up, min down).
    const double cmax_lo = rainbow_price(OptionType::Call, RainbowType::Max, S1, S2, K, r, sig1, sig2, -0.5, T);
    const double cmax_hr = rainbow_price(OptionType::Call, RainbowType::Max, S1, S2, K, r, sig1, sig2, 0.9, T);
    check::is_true("Cmax falls with correlation", cmax_lo > cmax_hr);
    const double cmin_lo = rainbow_price(OptionType::Call, RainbowType::Min, S1, S2, K, r, sig1, sig2, -0.5, T);
    const double cmin_hr = rainbow_price(OptionType::Call, RainbowType::Min, S1, S2, K, r, sig1, sig2, 0.9, T);
    check::is_true("Cmin rises with correlation", cmin_hr > cmin_lo);

    // === Monte Carlo cross-checks the closed form (call & put, max & min) ===
    const long n = 800'000;
    check::approx("MC call-max", rainbow_price_mc(OptionType::Call, RainbowType::Max, S1, S2, K, r, sig1, sig2, rho, T, n), cmax, 0.15);
    check::approx("MC call-min", rainbow_price_mc(OptionType::Call, RainbowType::Min, S1, S2, K, r, sig1, sig2, rho, T, n), cmin, 0.15);
    const double pmax = rainbow_price(OptionType::Put, RainbowType::Max, S1, S2, K, r, sig1, sig2, rho, T);
    const double pmin = rainbow_price(OptionType::Put, RainbowType::Min, S1, S2, K, r, sig1, sig2, rho, T);
    check::approx("MC put-max", rainbow_price_mc(OptionType::Put, RainbowType::Max, S1, S2, K, r, sig1, sig2, rho, T, n), pmax, 0.15);
    check::approx("MC put-min", rainbow_price_mc(OptionType::Put, RainbowType::Min, S1, S2, K, r, sig1, sig2, rho, T, n), pmin, 0.15);

    // Reproducible for a fixed seed.
    check::is_true("rainbow MC reproducible",
                   rainbow_price_mc(OptionType::Call, RainbowType::Max, S1, S2, K, r, sig1, sig2, rho, T, 100000) ==
                       rainbow_price_mc(OptionType::Call, RainbowType::Max, S1, S2, K, r, sig1, sig2, rho, T, 100000));

    return check::report("rainbow");
}
