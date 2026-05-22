// Tests for Gauss–Legendre quadrature.
#include "check.hpp"
#include "pricer/quadrature.hpp"
#include <cmath>

using namespace pricer;

static double integrate(int n, double a, double b, double (*f)(double)) {
    auto gl = gauss_legendre(n, a, b);
    double s = 0.0;
    for (std::size_t i = 0; i < gl.first.size(); ++i) s += gl.second[i] * f(gl.first[i]);
    return s;
}

int main() {
    // Polynomials up to degree 2n-1 are integrated exactly.
    check::approx("∫_0^1 x^3", integrate(4, 0, 1, [](double x){ return x * x * x; }), 0.25, 1e-12);
    check::approx("∫_0^2 x^5", integrate(6, 0, 2, [](double x){ return std::pow(x, 5); }),
                  64.0 / 6.0, 1e-9);

    // Smooth non-polynomial: ∫_{-1}^1 e^x = e - 1/e.
    check::approx("∫ exp", integrate(16, -1, 1, [](double x){ return std::exp(x); }),
                  std::exp(1.0) - std::exp(-1.0), 1e-12);

    // Weights sum to the interval length.
    auto gl = gauss_legendre(10, 2.0, 7.0);
    double wsum = 0.0;
    for (double w : gl.second) wsum += w;
    check::approx("weights sum to b-a", wsum, 5.0, 1e-12);

    return check::report("quadrature");
}
