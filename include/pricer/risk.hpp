// pricer/risk.hpp — portfolio risk from a P&L sample.
//
// Value-at-Risk (VaR) and Expected Shortfall (ES, a.k.a. CVaR) at a confidence
// level, computed from a sample of profit-and-loss outcomes. Both are reported
// as positive loss numbers: VaR is the loss quantile, ES the average loss in the
// tail beyond it. Feed in P&L scenarios from any pricing/simulation routine.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace pricer {

struct RiskMeasures {
    double var;  // loss not exceeded with `confidence` probability
    double es;   // mean loss given the loss is worse than VaR
};

// Historical/simulation VaR & ES from P&L outcomes (profit positive, loss
// negative). `confidence` e.g. 0.99 for a 99% measure.
inline RiskMeasures var_es(std::vector<double> pnl, double confidence = 0.99) {
    RiskMeasures m{0.0, 0.0};
    if (pnl.empty()) return m;

    std::vector<double> loss(pnl.size());
    for (std::size_t i = 0; i < pnl.size(); ++i) loss[i] = -pnl[i];
    std::sort(loss.begin(), loss.end());  // ascending: tail losses at the end

    const std::size_t n = loss.size();
    std::size_t idx = static_cast<std::size_t>(std::ceil(confidence * n));
    if (idx >= n) idx = n - 1;            // index of the VaR quantile
    m.var = loss[idx];

    double sum = 0.0;
    std::size_t cnt = 0;
    for (std::size_t i = idx; i < n; ++i) { sum += loss[i]; ++cnt; }
    m.es = cnt ? sum / static_cast<double>(cnt) : m.var;
    return m;
}

}  // namespace pricer
