// Tests for VaR / ES against the analytic normal values.
#include "check.hpp"
#include "pricer/risk.hpp"
#include <random>
#include <vector>

using namespace pricer;

int main() {
    // Standard-normal P&L: loss ~ N(0,1). Analytic 99% VaR = z_0.99 = 2.3263,
    // ES = phi(z)/(1-c) = 2.6652.
    std::mt19937_64 rng(2024);
    std::normal_distribution<double> N(0.0, 1.0);
    std::vector<double> pnl(2'000'000);
    for (auto& x : pnl) x = N(rng);

    const RiskMeasures m = var_es(pnl, 0.99);
    check::approx("VaR 99% vs normal", m.var, 2.3263, 2e-2);
    check::approx("ES 99% vs normal",  m.es,  2.6652, 2e-2);

    // ES must be at least as large as VaR (tail mean >= quantile).
    check::is_true("ES >= VaR", m.es >= m.var);

    return check::report("risk");
}
