// Tests for quadratic smile calibration.
#include "check.hpp"
#include "pricer/smile.hpp"
#include <cmath>
#include <vector>

using namespace pricer;

int main() {
    const double S = 100.0;

    // Generate vols from a known smile, then recover its parameters.
    const double a = 0.20, b = -0.05, c = 0.30;
    std::vector<double> Ks = {70, 85, 100, 115, 130, 150};
    std::vector<double> vols;
    for (double K : Ks) {
        const double k = std::log(K / S);
        vols.push_back(a + b * k + c * k * k);
    }

    const Smile fit = calibrate_smile(S, Ks, vols);
    check::approx("recover a", fit.a, a, 1e-9);
    check::approx("recover b", fit.b, b, 1e-9);
    check::approx("recover c", fit.c, c, 1e-9);

    // Reproduces the quotes within tolerance (the calibration exit criterion).
    for (size_t i = 0; i < Ks.size(); ++i) {
        char nm[32];
        std::snprintf(nm, sizeof nm, "reprice K=%.0f", Ks[i]);
        check::approx(nm, fit.vol(Ks[i]), vols[i], 1e-9);
    }

    return check::report("smile");
}
