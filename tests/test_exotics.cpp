// Tests for path-dependent exotic options (exotics.hpp): Asian, barrier, lookback.
// Each Monte Carlo engine is cross-checked against a closed form that is exact
// (geometric Asian) or exact under continuous monitoring (barrier, lookback).
#include "check.hpp"
#include "pricer/black_scholes.hpp"
#include "pricer/exotics.hpp"

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;

    // === Asian ===
    // Geometric-average Asian closed form with n=1 reduces to Black–Scholes.
    check::approx("geo Asian n=1 == BS call",
                  geometric_asian_price(OptionType::Call, S, K, r, sigma, T, 1),
                  black_scholes_call(S, K, r, sigma, T), 1e-9);

    // Monte Carlo geometric Asian matches its exact closed form (same monitoring).
    const int n_steps = 50;
    const double geo_cf = geometric_asian_price(OptionType::Call, S, K, r, sigma, T, n_steps);
    const double geo_mc =
        asian_price_mc(OptionType::Call, AverageType::Geometric, S, K, r, sigma, T, n_steps, 300'000);
    check::approx("geo Asian MC vs closed form", geo_mc, geo_cf, 0.05);

    // Arithmetic average >= geometric (AM–GM), so the arithmetic Asian call is dearer;
    // both are cheaper than the vanilla European (averaging dampens the terminal vol).
    const double ari_mc =
        asian_price_mc(OptionType::Call, AverageType::Arithmetic, S, K, r, sigma, T, n_steps, 300'000);
    check::is_true("arithmetic Asian >= geometric", ari_mc > geo_mc);
    check::is_true("Asian call < vanilla European call", ari_mc < black_scholes_call(S, K, r, sigma, T));

    // MC is reproducible for a fixed seed.
    const double geo_mc2 =
        asian_price_mc(OptionType::Call, AverageType::Geometric, S, K, r, sigma, T, n_steps, 300'000);
    check::is_true("Asian MC reproducible", geo_mc == geo_mc2);

    // === Barrier ===
    // Knock-in + knock-out = vanilla (parity), for every up/down × call/put case.
    const double bs_call = black_scholes_call(S, K, r, sigma, T);
    const double bs_put = black_scholes_put(S, K, r, sigma, T);
    const double Bup = 130.0, Bdn = 80.0;
    check::approx("UpIn + UpOut call == vanilla",
                  barrier_price(OptionType::Call, BarrierType::UpIn, S, K, Bup, r, sigma, T) +
                      barrier_price(OptionType::Call, BarrierType::UpOut, S, K, Bup, r, sigma, T),
                  bs_call, 1e-9);
    check::approx("DownIn + DownOut call == vanilla",
                  barrier_price(OptionType::Call, BarrierType::DownIn, S, K, Bdn, r, sigma, T) +
                      barrier_price(OptionType::Call, BarrierType::DownOut, S, K, Bdn, r, sigma, T),
                  bs_call, 1e-9);
    check::approx("UpIn + UpOut put == vanilla",
                  barrier_price(OptionType::Put, BarrierType::UpIn, S, K, Bup, r, sigma, T) +
                      barrier_price(OptionType::Put, BarrierType::UpOut, S, K, Bup, r, sigma, T),
                  bs_put, 1e-9);
    check::approx("DownIn + DownOut put == vanilla",
                  barrier_price(OptionType::Put, BarrierType::DownIn, S, K, Bdn, r, sigma, T) +
                      barrier_price(OptionType::Put, BarrierType::DownOut, S, K, Bdn, r, sigma, T),
                  bs_put, 1e-9);

    // A far-away barrier never knocks: up-and-out call ~ vanilla; down-and-out call ~ vanilla.
    check::approx("UpOut call, far barrier ~ vanilla",
                  barrier_price(OptionType::Call, BarrierType::UpOut, S, K, 1.0e6, r, sigma, T),
                  bs_call, 1e-6);
    check::approx("DownOut call, far barrier ~ vanilla",
                  barrier_price(OptionType::Call, BarrierType::DownOut, S, K, 1.0e-6, r, sigma, T),
                  bs_call, 1e-6);

    // MC barrier (with BGK continuity correction) matches the continuous closed form.
    const int bsteps = 250;
    const long bpaths = 400'000;
    const double do_cf = barrier_price(OptionType::Call, BarrierType::DownOut, S, K, 90.0, r, sigma, T);
    const double do_mc = barrier_price_mc(OptionType::Call, BarrierType::DownOut, S, K, 90.0, r, sigma,
                                          T, bsteps, bpaths);
    check::approx("DownOut call MC vs analytic", do_mc, do_cf, 0.15);

    const double uo_cf = barrier_price(OptionType::Call, BarrierType::UpOut, S, K, 130.0, r, sigma, T);
    const double uo_mc = barrier_price_mc(OptionType::Call, BarrierType::UpOut, S, K, 130.0, r, sigma,
                                          T, bsteps, bpaths);
    check::approx("UpOut call MC vs analytic", uo_mc, uo_cf, 0.15);

    // A knock-out is never worth more than the vanilla.
    check::is_true("knock-out <= vanilla", do_cf <= bs_call + 1e-9 && uo_cf <= bs_call + 1e-9);

    // === Lookback (floating strike) ===
    const double lb_call_cf = lookback_floating_price(OptionType::Call, S, r, sigma, T);
    const double lb_put_cf = lookback_floating_price(OptionType::Put, S, r, sigma, T);
    const double lb_call_mc =
        lookback_floating_price_mc(OptionType::Call, S, r, sigma, T, bsteps, bpaths);
    const double lb_put_mc =
        lookback_floating_price_mc(OptionType::Put, S, r, sigma, T, bsteps, bpaths);
    check::approx("lookback call MC vs closed form", lb_call_mc, lb_call_cf, 0.30);
    check::approx("lookback put MC vs closed form", lb_put_mc, lb_put_cf, 0.30);

    // A floating-strike lookback call (buy at the lowest) dominates the vanilla call.
    check::is_true("lookback call > vanilla call", lb_call_cf > bs_call);
    check::is_true("lookback payoffs positive", lb_call_cf > 0.0 && lb_put_cf > 0.0);

    return check::report("exotics");
}
