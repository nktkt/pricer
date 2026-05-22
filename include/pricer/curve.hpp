// pricer/curve.hpp — a discount (yield) curve.
//
// Built from pillar times and continuously-compounded zero rates. The zero rate
// is linearly interpolated in time (flat beyond the ends); discount factors and
// forward rates follow from it. This is the term-structure input that pricing
// and risk ultimately discount against.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace pricer {

class DiscountCurve {
public:
    // `times` must be strictly increasing and positive; `zeros` are the
    // continuously-compounded zero rates at those times.
    DiscountCurve(std::vector<double> times, std::vector<double> zeros)
        : t_(std::move(times)), z_(std::move(zeros)) {
        if (t_.size() != z_.size() || t_.empty())
            throw std::invalid_argument("DiscountCurve: size mismatch or empty");
        for (std::size_t i = 1; i < t_.size(); ++i)
            if (!(t_[i] > t_[i - 1]))
                throw std::invalid_argument("DiscountCurve: times must be increasing");
    }

    // Continuously-compounded zero rate at time t (linear interp, flat ends).
    double zero_rate(double t) const {
        if (t <= t_.front()) return z_.front();
        if (t >= t_.back()) return z_.back();
        const std::size_t i = static_cast<std::size_t>(
            std::upper_bound(t_.begin(), t_.end(), t) - t_.begin());
        const double w = (t - t_[i - 1]) / (t_[i] - t_[i - 1]);
        return z_[i - 1] + w * (z_[i] - z_[i - 1]);
    }

    // Discount factor P(0,t) = exp(-z(t) * t).
    double df(double t) const { return std::exp(-zero_rate(t) * t); }

    // Continuously-compounded forward rate over [t1, t2].
    double forward_rate(double t1, double t2) const {
        if (t2 <= t1) throw std::invalid_argument("forward_rate: need t2 > t1");
        return (zero_rate(t2) * t2 - zero_rate(t1) * t1) / (t2 - t1);
    }

private:
    std::vector<double> t_, z_;
};

}  // namespace pricer
