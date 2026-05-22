// pricer/normal.hpp — standard normal distribution helpers.
#pragma once
#include <cmath>

namespace pricer {

// Probability density function of the standard normal: n(x) = e^{-x^2/2} / sqrt(2*pi)
inline double norm_pdf(double x) {
    constexpr double inv_sqrt_2pi = 0.3989422804014327;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

// Cumulative distribution function of the standard normal, via erfc.
inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

}  // namespace pricer
