// Tests for the discount curve.
#include "check.hpp"
#include "pricer/curve.hpp"
#include <cmath>

using namespace pricer;

int main() {
    DiscountCurve c({0.5, 1.0, 2.0, 5.0}, {0.030, 0.035, 0.040, 0.045});

    // At pillars, df = exp(-z*t) and zero_rate returns the pillar rate.
    check::approx("zero at pillar", c.zero_rate(1.0), 0.035, 1e-12);
    check::approx("df at pillar", c.df(2.0), std::exp(-0.040 * 2.0), 1e-12);

    // Linear interpolation halfway between 1.0 and 2.0.
    check::approx("zero interp", c.zero_rate(1.5), 0.0375, 1e-12);

    // Flat extrapolation past the ends.
    check::approx("flat short end", c.zero_rate(0.1), 0.030, 1e-12);
    check::approx("flat long end", c.zero_rate(10.0), 0.045, 1e-12);

    // df(0) = 1 and discounting is monotone decreasing for positive rates.
    check::approx("df at 0", c.df(0.0), 1.0, 1e-12);
    check::is_true("df decreasing", c.df(1.0) > c.df(2.0) && c.df(2.0) > c.df(5.0));

    // Forward rate consistency: df(t1)*exp(-f*(t2-t1)) == df(t2).
    const double f = c.forward_rate(1.0, 2.0);
    check::approx("forward consistency", c.df(1.0) * std::exp(-f * 1.0), c.df(2.0), 1e-12);

    return check::report("curve");
}
