// risk_demo.cpp — 1-day VaR & Expected Shortfall of an option portfolio.
// Revalue a small multi-position book under many one-day Monte Carlo moves of a
// single underlying, build the P&L distribution, then read off VaR and ES.
#include "pricer/black_scholes.hpp"
#include "pricer/risk.hpp"
#include <vector>
#include <random>
#include <cmath>
#include <cstdio>

using namespace pricer;

// One option line in the book: contract spec plus signed quantity (long > 0).
struct Position {
    OptionType type;
    double strike;
    double T;     // time to expiry today (years)
    double qty;   // signed number of contracts
};

// Value the whole book at spot `S` with maturities shifted by `dt` (0 = today).
static double portfolio_value(const std::vector<Position>& book, double S,
                              double r, double sigma, double dt) {
    double v = 0.0;
    for (const Position& p : book)
        v += p.qty * black_scholes_price(p.type, S, p.strike, r, sigma, p.T - dt);
    return v;
}

int main() {
    // Market and underlying.
    const double S0 = 100.0, r = 0.05, sigma = 0.20;

    // Portfolio: a few positions on the same underlying.
    const std::vector<Position> book = {
        {OptionType::Call, 100.0, 1.0,  10.0},  // long 10 calls,  K=100
        {OptionType::Call, 110.0, 1.0,  -5.0},  // short 5 calls,  K=110
        {OptionType::Put,   90.0, 1.0,   8.0},  // long 8 puts,    K=90
    };

    // Current book value (no maturity shift).
    const double V0 = portfolio_value(book, S0, r, sigma, 0.0);

    // Net portfolio delta — qualitative read on directional exposure.
    double net_delta = 0.0;
    for (const Position& p : book)
        net_delta += p.qty * black_scholes_greeks(p.type, S0, p.strike, r, sigma, p.T).delta;

    // One-day risk horizon and Monte Carlo setup.
    const double dt = 1.0 / 252.0;          // one trading day
    const long n_scen = 500'000;            // number of scenarios
    std::mt19937_64 rng(20240521ULL);       // fixed seed for reproducibility
    std::normal_distribution<double> Z(0.0, 1.0);

    // Real-world-ish one-day GBM step shared by every scenario.
    const double drift = (r - 0.5 * sigma * sigma) * dt;
    const double vol   = sigma * std::sqrt(dt);

    // Simulate one-day P&L: revalue with the same options at maturity T - dt.
    std::vector<double> pnl;
    pnl.reserve(static_cast<std::size_t>(n_scen));
    for (long i = 0; i < n_scen; ++i) {
        const double S1 = S0 * std::exp(drift + vol * Z(rng));
        const double V1 = portfolio_value(book, S1, r, sigma, dt);
        pnl.push_back(V1 - V0);
    }

    // Risk measures at two confidence levels.
    const RiskMeasures r99 = var_es(pnl, 0.99);
    const RiskMeasures r95 = var_es(pnl, 0.95);

    // Report.
    std::printf("Current portfolio value V0 = %.4f\n", V0);
    std::printf("Net portfolio delta        = %+.4f (%s underlying)\n",
                net_delta, net_delta >= 0.0 ? "long" : "short");
    std::printf("Scenarios                  = %ld (1-day horizon, dt = 1/252)\n\n", n_scen);

    std::printf("%10s | %12s | %12s\n", "confidence", "VaR", "ES");
    std::printf("-----------|--------------|-------------\n");
    std::printf("%9.0f%% | %12.4f | %12.4f\n", 95.0, r95.var, r95.es);
    std::printf("%9.0f%% | %12.4f | %12.4f\n", 99.0, r99.var, r99.es);

    std::printf("\nVaR/ES are reported as positive one-day loss amounts; "
                "ES >= VaR and the 99%% loss exceeds the 95%% loss.\n");
    return 0;
}
