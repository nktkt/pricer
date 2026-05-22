// Tests for the Levenberg–Marquardt least-squares solver.
#include "check.hpp"
#include "pricer/optimize.hpp"
#include <cmath>
#include <vector>

using namespace pricer;

int main() {
    // Nonlinear fit: recover (p0, p1) of y = p0 * exp(p1 * x) from clean data.
    const double p0_true = 2.0, p1_true = 0.5;
    std::vector<double> xs = {0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0};
    std::vector<double> ys;
    for (double x : xs) ys.push_back(p0_true * std::exp(p1_true * x));

    auto residual = [&](const std::vector<double>& p, std::vector<double>& out) {
        for (std::size_t i = 0; i < xs.size(); ++i)
            out[i] = p[0] * std::exp(p[1] * xs[i]) - ys[i];
    };

    const opt::LMResult res =
        opt::levenberg_marquardt(residual, {1.0, 1.0}, static_cast<int>(xs.size()));

    check::is_true("converged", res.converged);
    check::approx("recover p0", res.params[0], p0_true, 1e-5);
    check::approx("recover p1", res.params[1], p1_true, 1e-5);
    check::is_true("cost ~ 0", res.cost < 1e-12);

    return check::report("optimize");
}
