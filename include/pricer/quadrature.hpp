// pricer/quadrature.hpp — Gauss–Legendre numerical integration.
//
// Computes nodes and weights for n-point Gauss–Legendre quadrature on [a, b]
// (roots of the Legendre polynomial found by Newton's method). Used to integrate
// the Heston characteristic-function inversion, but it is a general-purpose tool.
#pragma once
#include <cmath>
#include <utility>
#include <vector>

namespace pricer {

// Returns (nodes, weights) for n-point Gauss–Legendre quadrature on [a, b]:
//   integral_a^b f ≈ sum_i weights[i] * f(nodes[i]).
inline std::pair<std::vector<double>, std::vector<double>> gauss_legendre(int n, double a,
                                                                          double b) {
    std::vector<double> x(n), w(n);
    const double pi = 3.14159265358979323846;
    for (int i = 0; i < n; ++i) {
        double xi = std::cos(pi * (i + 0.75) / (n + 0.5));  // initial guess for the i-th root
        for (int it = 0; it < 100; ++it) {                  // Newton refinement on P_n
            double p0 = 1.0, p1 = xi;
            for (int k = 2; k <= n; ++k) {
                const double p2 = ((2 * k - 1) * xi * p1 - (k - 1) * p0) / k;
                p0 = p1; p1 = p2;
            }
            const double dp = n * (xi * p1 - p0) / (xi * xi - 1.0);  // P_n'(xi)
            const double dx = -p1 / dp;
            xi += dx;
            if (std::fabs(dx) < 1e-15) break;
        }
        // Recompute P_n'(xi) for the weight.
        double p0 = 1.0, p1 = xi;
        for (int k = 2; k <= n; ++k) {
            const double p2 = ((2 * k - 1) * xi * p1 - (k - 1) * p0) / k;
            p0 = p1; p1 = p2;
        }
        const double dp = n * (xi * p1 - p0) / (xi * xi - 1.0);
        x[i] = 0.5 * (a + b) + 0.5 * (b - a) * xi;
        w[i] = (b - a) / ((1.0 - xi * xi) * dp * dp);
    }
    return {std::move(x), std::move(w)};
}

}  // namespace pricer
