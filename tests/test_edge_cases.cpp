// Robustness tests: boundary inputs and error handling across the library.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/csv.hpp"
#include "pricer/curve.hpp"
#include "pricer/implied_vol.hpp"
#include "pricer/risk.hpp"
#include "pricer/smile.hpp"
#include <stdexcept>
#include <vector>

using namespace pricer;

template <class F>
static bool throws(F f) {
    try { f(); return false; } catch (const std::exception&) { return true; }
}

int main() {
    // VaR/ES on an empty sample must not crash (returns zeros).
    {
        const RiskMeasures m = var_es({});
        check::approx("empty VaR", m.var, 0.0, 0.0);
        check::approx("empty ES", m.es, 0.0, 0.0);
    }

    // Discount curve input validation.
    check::is_true("curve rejects size mismatch",
                   throws([] { DiscountCurve c({1.0, 2.0}, {0.03}); }));
    check::is_true("curve rejects non-increasing times",
                   throws([] { DiscountCurve c({2.0, 1.0}, {0.03, 0.04}); }));
    {
        DiscountCurve c({1.0, 2.0}, {0.03, 0.04});
        check::is_true("forward_rate rejects t2<=t1", throws([&] { c.forward_rate(2.0, 1.0); }));
    }

    // Implied vol round-trips even for deep in/out-of-the-money quotes.
    {
        const double S = 100, r = 0.05, T = 1.0;
        for (double K : {60.0, 150.0}) {
            const double price = black_scholes_call(S, K, r, 0.35, T);
            check::approx("deep iv round-trip", implied_vol(OptionType::Call, price, S, K, r, T),
                          0.35, 1e-5);
        }
    }

    // Smile calibration needs at least three quotes.
    check::is_true("smile rejects <3 quotes",
                   throws([] { calibrate_smile(100.0, {90.0, 110.0}, {0.2, 0.2}); }));

    // CSV reader reports a missing file rather than silently returning empty.
    check::is_true("read_csv rejects missing file",
                   throws([] { read_csv("/no/such/pricer_file_xyz.csv"); }));

    return check::report("edge_cases");
}
