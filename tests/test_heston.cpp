// Tests for the Heston model: Black–Scholes limit, parity, and calibration.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/heston.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

using namespace pricer;

int main() {
    const double S = 100, r = 0.05, T = 1.0;

    // 1) Black–Scholes limit: with tiny vol-of-vol and v0 = theta, Heston ≈ BS(sqrt(v0)).
    {
        HestonParams p{1.0, 0.04, 0.01, 0.0, 0.04};  // sqrt(0.04) = 0.20 vol
        for (double K : {80.0, 100.0, 120.0}) {
            char nm[32]; std::snprintf(nm, sizeof nm, "BS-limit K=%.0f", K);
            check::approx(nm, heston_call(p, S, K, r, T), black_scholes_call(S, K, r, 0.20, T), 2e-3);
        }
    }

    // 2) Put–call parity: C - P = S - K e^{-rT}.
    {
        HestonParams p{2.0, 0.04, 0.30, -0.60, 0.04};
        const double c = heston_price(OptionType::Call, p, S, 100, r, T);
        const double put = heston_price(OptionType::Put, p, S, 100, r, T);
        check::approx("parity", c - put, S - 100 * std::exp(-r * T), 1e-8);
    }

    // 3) Calibration: recover market prices generated from known Heston params.
    {
        const HestonParams truth{2.0, 0.04, 0.30, -0.60, 0.04};
        std::vector<double> Ks = {90, 100, 110, 90, 100, 110};
        std::vector<double> Ts = {0.5, 0.5, 0.5, 1.0, 1.0, 1.0};
        std::vector<double> prices;
        for (std::size_t i = 0; i < Ks.size(); ++i)
            prices.push_back(heston_call(truth, S, Ks[i], r, Ts[i]));

        const HestonParams fit =
            calibrate_heston(S, r, Ks, Ts, prices, {1.5, 0.05, 0.40, -0.40, 0.05});

        double max_err = 0.0;
        for (std::size_t i = 0; i < Ks.size(); ++i)
            max_err = std::max(max_err, std::fabs(heston_call(fit, S, Ks[i], r, Ts[i]) - prices[i]));
        check::is_true("calibration reproduces quotes (<1e-2)", max_err < 1e-2);
        std::printf("  [info] heston calibration max price error = %.2e\n", max_err);
    }

    return check::report("heston");
}
