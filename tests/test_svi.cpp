// Tests for SVI smile calibration via Levenberg–Marquardt.
#include "check.hpp"
#include "pricer/svi.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

using namespace pricer;

int main() {
    const double S = 100.0, T = 1.0;

    // Generate market vols from a known SVI smile, then calibrate to them.
    const SVIParams truth{0.04, 0.10, -0.40, 0.00, 0.20};
    std::vector<double> Ks = {70, 80, 90, 100, 110, 120, 140};
    std::vector<double> vols;
    for (double K : Ks) vols.push_back(truth.vol(K, S, T));

    const SVIParams fit = calibrate_svi(S, T, Ks, vols);

    // The fit must reproduce the quotes within a tight tolerance.
    double max_err = 0.0;
    for (std::size_t i = 0; i < Ks.size(); ++i) {
        const double e = std::fabs(fit.vol(Ks[i], S, T) - vols[i]);
        if (e > max_err) max_err = e;
        char nm[32];
        std::snprintf(nm, sizeof nm, "reprice K=%.0f", Ks[i]);
        check::approx(nm, fit.vol(Ks[i], S, T), vols[i], 1e-4);
    }
    check::is_true("max reproduction error tiny", max_err < 1e-4);

    // ATM total variance is recovered.
    check::approx("ATM total var", fit.total_variance(0.0), truth.total_variance(0.0), 1e-5);

    return check::report("svi");
}
