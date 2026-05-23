// Tests for multi-asset options (basket.hpp): basket options on correlated GBM
// and spread/exchange options. Monte Carlo engines are cross-checked against the
// geometric-basket closed form and Margrabe's exact exchange-option formula.
#include "check.hpp"
#include "pricer/basket.hpp"
#include "pricer/black_scholes.hpp"

using namespace pricer;

int main() {
    const double r = 0.05, T = 1.0;

    // === Basket ===
    // One asset, weight 1: the geometric basket reduces to Black–Scholes.
    check::approx("geo basket n=1 == BS call",
                  geometric_basket_price(OptionType::Call, {100.0}, {1.0}, 100.0, r, {0.20},
                                         {{1.0}}, T),
                  black_scholes_call(100.0, 100.0, r, 0.20, T), 1e-9);

    // Three correlated assets, equal weights.
    const std::vector<double> S{100, 95, 105}, w{1.0 / 3, 1.0 / 3, 1.0 / 3}, sig{0.20, 0.25, 0.30};
    const std::vector<std::vector<double>> corr{
        {1.0, 0.5, 0.3}, {0.5, 1.0, 0.4}, {0.3, 0.4, 1.0}};
    const double K = 100.0;

    // Geometric basket: Monte Carlo matches the exact closed form (this also
    // validates the Cholesky-correlated simulation, since V uses the full matrix).
    const double geo_cf = geometric_basket_price(OptionType::Call, S, w, K, r, sig, corr, T);
    const double geo_mc =
        basket_price_mc(OptionType::Call, AverageType::Geometric, S, w, K, r, sig, corr, T, 400'000);
    check::approx("geo basket MC vs closed form", geo_mc, geo_cf, 0.10);

    // The arithmetic basket (the real product) is worth at least the geometric one.
    const double ari_mc =
        basket_price_mc(OptionType::Call, AverageType::Arithmetic, S, w, K, r, sig, corr, T, 400'000);
    check::is_true("arithmetic basket >= geometric", ari_mc > geo_mc);

    // A basket is a sum, so higher correlation raises its variance and its value.
    const std::vector<std::vector<double>> lo{{1.0, 0.1, 0.1}, {0.1, 1.0, 0.1}, {0.1, 0.1, 1.0}};
    const std::vector<std::vector<double>> hi{{1.0, 0.9, 0.9}, {0.9, 1.0, 0.9}, {0.9, 0.9, 1.0}};
    const double basket_lo =
        basket_price_mc(OptionType::Call, AverageType::Arithmetic, S, w, K, r, sig, lo, T, 400'000);
    const double basket_hi =
        basket_price_mc(OptionType::Call, AverageType::Arithmetic, S, w, K, r, sig, hi, T, 400'000);
    check::is_true("basket value rises with correlation", basket_hi > basket_lo);

    // Reproducible for a fixed seed.
    const double geo_mc2 =
        basket_price_mc(OptionType::Call, AverageType::Geometric, S, w, K, r, sig, corr, T, 400'000);
    check::is_true("basket MC reproducible", geo_mc == geo_mc2);

    // === Spread / exchange ===
    const double S1 = 100, S2 = 100, sig1 = 0.20, sig2 = 0.30, rho = 0.4;

    // Kirk reduces to Margrabe exactly at K = 0 (an exchange option).
    const double marg = margrabe_exchange_price(S1, S2, sig1, sig2, rho, T);
    check::approx("Kirk(K=0) == Margrabe",
                  spread_kirk_price(OptionType::Call, S1, S2, 0.0, r, sig1, sig2, rho, T), marg,
                  1e-10);

    // Monte Carlo of the exchange payoff matches Margrabe's exact price.
    const double exch_mc = spread_price_mc(OptionType::Call, S1, S2, 0.0, r, sig1, sig2, rho, T, 400'000);
    check::approx("exchange MC vs Margrabe", exch_mc, marg, 0.10);

    // Spread with a non-zero strike: Kirk's approximation matches MC closely.
    const double kirk = spread_kirk_price(OptionType::Call, S1, S2, 5.0, r, sig1, sig2, rho, T);
    const double spread_mc = spread_price_mc(OptionType::Call, S1, S2, 5.0, r, sig1, sig2, rho, T, 400'000);
    check::approx("spread Kirk vs MC (K=5)", kirk, spread_mc, 0.15);

    // A spread narrows as the legs co-move, so higher correlation lowers the value.
    const double spread_lorho = spread_kirk_price(OptionType::Call, S1, S2, 5.0, r, sig1, sig2, 0.1, T);
    const double spread_hirho = spread_kirk_price(OptionType::Call, S1, S2, 5.0, r, sig1, sig2, 0.9, T);
    check::is_true("spread value falls with correlation", spread_hirho < spread_lorho);

    return check::report("basket");
}
