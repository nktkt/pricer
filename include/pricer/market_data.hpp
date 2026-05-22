// pricer/market_data.hpp — load market data and persist results (CSV adapters).
//
// File-based adapters for getting quotes and a discount curve into the engine,
// and for writing priced results back out. CSV keeps it dependency-free; the
// same interfaces could be backed by a database or a live feed later.
#pragma once
#include <string>
#include <vector>

#include "pricer/black_scholes.hpp"  // OptionType
#include "pricer/csv.hpp"
#include "pricer/curve.hpp"

namespace pricer {

// One market option quote.
struct OptionQuote {
    OptionType type;
    double strike;
    double expiry;
    double price;
};

inline OptionType parse_option_type(const std::string& s) {
    return (s == "put" || s == "Put" || s == "PUT" || s == "P" || s == "p") ? OptionType::Put
                                                                            : OptionType::Call;
}

// Load option quotes from a CSV with header: type,strike,expiry,price
inline std::vector<OptionQuote> load_option_quotes(const std::string& path) {
    const auto rows = read_csv(path);
    std::vector<OptionQuote> quotes;
    for (size_t i = 1; i < rows.size(); ++i) {  // row 0 is the header
        const auto& r = rows[i];
        if (r.size() < 4) continue;
        quotes.push_back({parse_option_type(r[0]), std::stod(r[1]), std::stod(r[2]), std::stod(r[3])});
    }
    return quotes;
}

// Load a discount curve from a CSV with header: time,zero_rate
inline DiscountCurve load_curve(const std::string& path) {
    const auto rows = read_csv(path);
    std::vector<double> t, z;
    for (size_t i = 1; i < rows.size(); ++i) {
        if (rows[i].size() < 2) continue;
        t.push_back(std::stod(rows[i][0]));
        z.push_back(std::stod(rows[i][1]));
    }
    return DiscountCurve(std::move(t), std::move(z));
}

// A priced result row, for persistence.
struct PricedResult {
    double strike, expiry, price, delta, gamma, vega, theta, rho;
};

inline void save_results(const std::string& path, const std::vector<PricedResult>& results) {
    std::vector<std::vector<std::string>> rows;
    rows.push_back({"strike", "expiry", "price", "delta", "gamma", "vega", "theta", "rho"});
    for (const auto& r : results)
        rows.push_back({std::to_string(r.strike), std::to_string(r.expiry), std::to_string(r.price),
                        std::to_string(r.delta), std::to_string(r.gamma), std::to_string(r.vega),
                        std::to_string(r.theta), std::to_string(r.rho)});
    write_csv(path, rows);
}

inline std::vector<PricedResult> load_results(const std::string& path) {
    const auto rows = read_csv(path);
    std::vector<PricedResult> out;
    for (size_t i = 1; i < rows.size(); ++i) {
        const auto& r = rows[i];
        if (r.size() < 8) continue;
        out.push_back({std::stod(r[0]), std::stod(r[1]), std::stod(r[2]), std::stod(r[3]),
                       std::stod(r[4]), std::stod(r[5]), std::stod(r[6]), std::stod(r[7])});
    }
    return out;
}

}  // namespace pricer
