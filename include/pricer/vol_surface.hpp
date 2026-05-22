// pricer/vol_surface.hpp — an implied-volatility surface from market prices.
//
// Given a grid of expiries, strikes and (call) market prices, back out the
// Black–Scholes implied vol at each node, then interpolate bilinearly in
// (strike, expiry). Querying the surface gives an implied vol — and hence a
// price — at any (K, T) inside the grid (clamped at the edges).
#pragma once
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "pricer/black_scholes.hpp"
#include "pricer/implied_vol.hpp"

namespace pricer {

class VolSurface {
public:
    // expiries (ascending) × strikes (ascending); prices[i][j] is the market
    // call price for expiry i, strike j. S = spot, r = rate.
    VolSurface(double S, double r, std::vector<double> expiries, std::vector<double> strikes,
               const std::vector<std::vector<double>>& call_prices)
        : S_(S), r_(r), T_(std::move(expiries)), K_(std::move(strikes)) {
        if (T_.empty() || K_.empty() || call_prices.size() != T_.size())
            throw std::invalid_argument("VolSurface: grid size mismatch");
        iv_.assign(T_.size(), std::vector<double>(K_.size(), 0.0));
        for (std::size_t i = 0; i < T_.size(); ++i) {
            if (call_prices[i].size() != K_.size())
                throw std::invalid_argument("VolSurface: row size mismatch");
            for (std::size_t j = 0; j < K_.size(); ++j)
                iv_[i][j] = implied_vol(OptionType::Call, call_prices[i][j], S_, K_[j], r_, T_[i]);
        }
    }

    // Implied vol at (K, T) by bilinear interpolation (clamped to the grid).
    double vol(double K, double T) const {
        std::size_t i0, i1; double wt;
        std::size_t j0, j1; double wk;
        bracket(T_, T, i0, i1, wt);
        bracket(K_, K, j0, j1, wk);
        const double v0 = iv_[i0][j0] + wk * (iv_[i0][j1] - iv_[i0][j0]);
        const double v1 = iv_[i1][j0] + wk * (iv_[i1][j1] - iv_[i1][j0]);
        return v0 + wt * (v1 - v0);
    }

    // Black–Scholes call price using the interpolated vol.
    double call_price(double K, double T) const {
        return black_scholes_call(S_, K, r_, vol(K, T), T);
    }

private:
    // Find indices lo/hi in ascending `xs` bracketing x (clamped) and weight w.
    static void bracket(const std::vector<double>& xs, double x, std::size_t& lo,
                        std::size_t& hi, double& w) {
        if (x <= xs.front() || xs.size() == 1) { lo = hi = 0; w = 0.0; return; }
        if (x >= xs.back()) { lo = hi = xs.size() - 1; w = 0.0; return; }
        std::size_t i = 1;
        while (i < xs.size() && xs[i] < x) ++i;
        lo = i - 1; hi = i;
        w = (x - xs[lo]) / (xs[hi] - xs[lo]);
    }

    double S_, r_;
    std::vector<double> T_, K_;
    std::vector<std::vector<double>> iv_;
};

}  // namespace pricer
