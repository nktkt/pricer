// market_data_demo.cpp — end-to-end market-data workflow.
// Write sample option quotes to CSV, load them back, recover each quote's
// implied volatility, compute its Greeks, and round-trip the priced results
// through save_results / load_results — all on temp files that are cleaned up.
#include "pricer/market_data.hpp"
#include "pricer/implied_vol.hpp"
#include "pricer/black_scholes.hpp"

#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

using namespace pricer;

int main() {
    const double S = 100, r = 0.05;

    // Temp file paths, so we don't litter the source tree.
    namespace fs = std::filesystem;
    const fs::path quotes_path = fs::temp_directory_path() / "pricer_quotes.csv";
    const fs::path results_path = fs::temp_directory_path() / "pricer_results.csv";

    // 1) Build a sample quotes CSV. Prices are computed at 20% vol with the
    //    Black–Scholes model so the implied vols we recover later are realistic.
    const double true_vol = 0.20;
    std::vector<std::vector<std::string>> rows;
    rows.push_back({"type", "strike", "expiry", "price"});
    struct SampleQuote { OptionType type; const char* name; double K; double T; };
    const std::vector<SampleQuote> samples = {
        {OptionType::Call, "call", 90.0, 0.5},
        {OptionType::Put, "put", 100.0, 0.5},
        {OptionType::Call, "call", 110.0, 1.0},
        {OptionType::Put, "put", 90.0, 1.0},
        {OptionType::Call, "call", 100.0, 1.0},
    };
    for (const auto& q : samples) {
        const double price = black_scholes_greeks(q.type, S, q.K, r, true_vol, q.T).price;
        rows.push_back({q.name, std::to_string(q.K), std::to_string(q.T), std::to_string(price)});
    }
    write_csv(quotes_path.string(), rows);
    std::printf("Wrote %zu sample quotes to %s\n\n", samples.size(), quotes_path.string().c_str());

    // 2) Load the quotes back from disk.
    const std::vector<OptionQuote> quotes = load_option_quotes(quotes_path.string());
    std::printf("Loaded quotes (S = %.2f, r = %.2f):\n", S, r);
    std::printf("%6s | %8s | %8s | %12s\n", "type", "strike", "expiry", "price");
    std::printf("-------|----------|----------|-------------\n");
    for (const auto& q : quotes) {
        std::printf("%6s | %8.2f | %8.2f | %12.6f\n",
                    q.type == OptionType::Call ? "call" : "put",
                    q.strike, q.expiry, q.price);
    }
    std::printf("\n");

    // 3) Recover implied vol per quote and compute Greeks from it.
    std::vector<PricedResult> results;
    std::printf("Recovered implied volatilities:\n");
    std::printf("%6s | %8s | %8s | %12s\n", "type", "strike", "expiry", "implied_vol");
    std::printf("-------|----------|----------|-------------\n");
    for (const auto& q : quotes) {
        const double iv = implied_vol(q.type, q.price, S, q.strike, r, q.expiry);
        const Greeks g = black_scholes_greeks(q.type, S, q.strike, r, iv, q.expiry);
        results.push_back({q.strike, q.expiry, g.price, g.delta, g.gamma, g.vega, g.theta, g.rho});
        std::printf("%6s | %8.2f | %8.2f | %12.6f\n",
                    q.type == OptionType::Call ? "call" : "put",
                    q.strike, q.expiry, iv);
    }
    std::printf("\n");

    // 4) Persist the priced results, then load them back to show the round-trip.
    save_results(results_path.string(), results);
    std::printf("Saved %zu priced results to %s\n\n", results.size(), results_path.string().c_str());

    const std::vector<PricedResult> loaded = load_results(results_path.string());
    std::printf("Round-tripped priced results:\n");
    std::printf("%8s | %8s | %10s | %8s | %8s | %8s | %10s | %10s\n",
                "strike", "expiry", "price", "delta", "gamma", "vega", "theta", "rho");
    std::printf("---------|----------|------------|----------|----------|----------|------------|------------\n");
    for (const auto& res : loaded) {
        std::printf("%8.2f | %8.2f | %10.6f | %8.4f | %8.4f | %8.4f | %10.4f | %10.4f\n",
                    res.strike, res.expiry, res.price, res.delta, res.gamma,
                    res.vega, res.theta, res.rho);
    }
    std::printf("\n");

    // 5) Clean up the temp files.
    fs::remove(quotes_path);
    fs::remove(results_path);
    std::printf("Cleaned up temp files. Recovered implied vols should all be ~%.2f.\n", true_vol);

    return 0;
}
