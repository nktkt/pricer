// Tests for the implied-volatility surface.
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/vol_surface.hpp"
#include <cstdio>
#include <vector>

using namespace pricer;

int main() {
    const double S = 100, r = 0.05;
    std::vector<double> Ts = {0.25, 1.0, 2.0};
    std::vector<double> Ks = {80, 100, 120};

    // Case 1: a flat surface. Prices generated at constant vol 0.20 must imply 0.20
    // at every node, and interpolated points must too.
    {
        std::vector<std::vector<double>> prices(Ts.size(), std::vector<double>(Ks.size()));
        for (size_t i = 0; i < Ts.size(); ++i)
            for (size_t j = 0; j < Ks.size(); ++j)
                prices[i][j] = black_scholes_call(S, Ks[j], r, 0.20, Ts[i]);
        VolSurface surf(S, r, Ts, Ks, prices);
        check::approx("flat node vol", surf.vol(100, 1.0), 0.20, 1e-6);
        check::approx("flat interp vol", surf.vol(110, 1.5), 0.20, 1e-6);
        check::approx("flat clamp vol", surf.vol(60, 0.1), 0.20, 1e-6);
        check::approx("reprice via surface", surf.call_price(100, 1.0),
                      black_scholes_call(S, 100, r, 0.20, 1.0), 1e-6);
    }

    // Case 2: a strike-dependent smile. Vol varies by strike; nodes must recover it.
    {
        auto node_vol = [](double K) { return 0.20 + 0.001 * (100.0 - K); };  // skew
        std::vector<std::vector<double>> prices(Ts.size(), std::vector<double>(Ks.size()));
        for (size_t i = 0; i < Ts.size(); ++i)
            for (size_t j = 0; j < Ks.size(); ++j)
                prices[i][j] = black_scholes_call(S, Ks[j], r, node_vol(Ks[j]), Ts[i]);
        VolSurface surf(S, r, Ts, Ks, prices);
        check::approx("skew node 80",  surf.vol(80, 1.0),  node_vol(80),  1e-6);
        check::approx("skew node 120", surf.vol(120, 1.0), node_vol(120), 1e-6);
    }

    return check::report("vol_surface");
}
