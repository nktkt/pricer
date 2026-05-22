// Tests for the CSV market-data adapters and result persistence.
#include "check.hpp"
#include "pricer/market_data.hpp"
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

using namespace pricer;

int main() {
    const auto dir = std::filesystem::temp_directory_path();
    const std::string quotes_path = (dir / "pricer_quotes.csv").string();
    const std::string curve_path = (dir / "pricer_curve.csv").string();
    const std::string results_path = (dir / "pricer_results.csv").string();

    // --- option quotes round-trip ---
    {
        std::ofstream f(quotes_path);
        f << "type,strike,expiry,price\n"
          << "call,100,1.0,10.45\n"
          << "put,90,0.5,2.31\n";
    }
    const auto quotes = load_option_quotes(quotes_path);
    check::is_true("loaded 2 quotes", quotes.size() == 2);
    check::is_true("quote0 is call", quotes[0].type == OptionType::Call);
    check::approx("quote0 strike", quotes[0].strike, 100.0, 1e-12);
    check::is_true("quote1 is put", quotes[1].type == OptionType::Put);
    check::approx("quote1 price", quotes[1].price, 2.31, 1e-12);

    // --- curve from CSV ---
    {
        std::ofstream f(curve_path);
        f << "time,zero\n0.5,0.03\n1.0,0.035\n2.0,0.04\n";
    }
    const DiscountCurve curve = load_curve(curve_path);
    check::approx("curve df at pillar", curve.df(1.0), std::exp(-0.035 * 1.0), 1e-12);

    // --- results persistence round-trip ---
    {
        std::vector<PricedResult> out = {{100, 1.0, 10.450584, 0.6368, 0.01876, 37.524, -6.414, 53.23},
                                         {110, 1.0, 5.5, -0.4, 0.018, 36.0, -5.0, 40.0}};
        save_results(results_path, out);
        const auto back = load_results(results_path);
        check::is_true("results round-trip count", back.size() == out.size());
        check::approx("result0 price", back[0].price, 10.450584, 1e-6);
        check::approx("result1 delta", back[1].delta, -0.4, 1e-6);
    }

    std::filesystem::remove(quotes_path);
    std::filesystem::remove(curve_path);
    std::filesystem::remove(results_path);
    return check::report("market_data");
}
