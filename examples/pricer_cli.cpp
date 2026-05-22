// pricer_cli.cpp — a small command-line front-end for the header-only `pricer`
// library. Supports closed-form Black-Scholes pricing + Greeks, implied
// volatility, and parallel Monte Carlo pricing. Manual argv parsing only,
// no external dependencies.

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>

#include "pricer/black_scholes.hpp"
#include "pricer/implied_vol.hpp"
#include "pricer/parallel.hpp"

namespace {

// Print usage to stdout.
void print_usage() {
    std::printf(
        "pricer_cli - command-line interface for the pricer library\n"
        "\n"
        "Usage:\n"
        "  pricer_cli price --type {call|put} --S <n> --K <n> --r <n> --sigma <n> --T <n>\n"
        "      Print the Black-Scholes price and all Greeks.\n"
        "\n"
        "  pricer_cli iv    --type {call|put} --price <n> --S <n> --K <n> --r <n> --T <n>\n"
        "      Print the implied volatility.\n"
        "\n"
        "  pricer_cli mc    --type {call|put} --S <n> --K <n> --r <n> --sigma <n> --T <n>\n"
        "                   [--paths <n>] [--threads <n>]\n"
        "      Print the Monte Carlo (parallel) price.\n"
        "      Defaults: --paths 1000000, --threads 0 (= all cores).\n"
        "\n"
        "  pricer_cli --help\n"
        "      Print this message.\n");
}

// Collect every \"--flag value\" pair into a map. Flags may appear in any order.
// Throws std::runtime_error if a flag is missing its value.
std::map<std::string, std::string> parse_flags(int argc, char** argv, int start) {
    std::map<std::string, std::string> flags;
    for (int i = start; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--", 0) == 0) {
            std::string key = arg.substr(2);
            if (i + 1 >= argc)
                throw std::runtime_error("missing value for flag --" + key);
            flags[key] = argv[++i];
        } else {
            throw std::runtime_error("unexpected argument: " + arg);
        }
    }
    return flags;
}

// Look up a required flag; throw if absent.
const std::string& require(const std::map<std::string, std::string>& flags,
                           const std::string& key) {
    auto it = flags.find(key);
    if (it == flags.end())
        throw std::runtime_error("missing required flag --" + key);
    return it->second;
}

// Parse a required numeric flag via std::stod, with a clear error message.
double require_num(const std::map<std::string, std::string>& flags,
                   const std::string& key) {
    const std::string& s = require(flags, key);
    try {
        std::size_t pos = 0;
        double v = std::stod(s, &pos);
        if (pos != s.size())
            throw std::invalid_argument("trailing characters");
        return v;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid numeric value for --" + key + ": " + s);
    }
}

// Parse an optional numeric flag; return `def` when absent.
double optional_num(const std::map<std::string, std::string>& flags,
                    const std::string& key, double def) {
    if (flags.find(key) == flags.end()) return def;
    return require_num(flags, key);
}

// Lower-case helper for case-insensitive flag values.
std::string to_lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// Parse the --type flag (case-insensitive). Throws on an unknown value.
pricer::OptionType parse_type(const std::map<std::string, std::string>& flags) {
    std::string t = to_lower(require(flags, "type"));
    if (t == "call") return pricer::OptionType::Call;
    if (t == "put")  return pricer::OptionType::Put;
    throw std::runtime_error("invalid --type (expected call|put): " + t);
}

// --- Subcommands. Each returns the process exit code. ---

int cmd_price(const std::map<std::string, std::string>& flags) {
    pricer::OptionType type = parse_type(flags);
    double S = require_num(flags, "S");
    double K = require_num(flags, "K");
    double r = require_num(flags, "r");
    double sigma = require_num(flags, "sigma");
    double T = require_num(flags, "T");

    pricer::Greeks g = pricer::black_scholes_greeks(type, S, K, r, sigma, T);
    // First line is machine-checkable: exactly "price: <value>" with 6 decimals.
    std::printf("price: %.6f\n", g.price);
    std::printf("delta: %.6f\n", g.delta);
    std::printf("gamma: %.6f\n", g.gamma);
    std::printf("vega: %.6f\n", g.vega);
    std::printf("theta: %.6f\n", g.theta);
    std::printf("rho: %.6f\n", g.rho);
    return 0;
}

int cmd_iv(const std::map<std::string, std::string>& flags) {
    pricer::OptionType type = parse_type(flags);
    double price = require_num(flags, "price");
    double S = require_num(flags, "S");
    double K = require_num(flags, "K");
    double r = require_num(flags, "r");
    double T = require_num(flags, "T");

    double iv = pricer::implied_vol(type, price, S, K, r, T);
    std::printf("implied_vol: %.6f\n", iv);
    return 0;
}

int cmd_mc(const std::map<std::string, std::string>& flags) {
    pricer::OptionType type = parse_type(flags);
    double S = require_num(flags, "S");
    double K = require_num(flags, "K");
    double r = require_num(flags, "r");
    double sigma = require_num(flags, "sigma");
    double T = require_num(flags, "T");
    long n_paths = static_cast<long>(optional_num(flags, "paths", 1000000.0));
    unsigned n_threads = static_cast<unsigned>(optional_num(flags, "threads", 0.0));

    double mc_price;
    if (type == pricer::OptionType::Call) {
        // Call payoff: max(ST - K, 0).
        auto payoff = [K](double ST) { return ST > K ? ST - K : 0.0; };
        mc_price = pricer::mc::price_terminal_parallel(payoff, S, r, sigma, T,
                                                       n_paths, 12345, n_threads);
    } else {
        // Put payoff: max(K - ST, 0).
        auto payoff = [K](double ST) { return ST < K ? K - ST : 0.0; };
        mc_price = pricer::mc::price_terminal_parallel(payoff, S, r, sigma, T,
                                                       n_paths, 12345, n_threads);
    }
    std::printf("mc_price: %.6f\n", mc_price);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    // No args, or an explicit help request -> print usage and succeed.
    if (argc < 2) {
        print_usage();
        return 0;
    }
    std::string cmd = argv[1];
    if (cmd == "--help" || cmd == "-h" || cmd == "help") {
        print_usage();
        return 0;
    }

    try {
        std::map<std::string, std::string> flags = parse_flags(argc, argv, 2);
        if (cmd == "price") return cmd_price(flags);
        if (cmd == "iv")    return cmd_iv(flags);
        if (cmd == "mc")    return cmd_mc(flags);
        // Unknown command is a genuine error.
        throw std::runtime_error("unknown command: " + cmd);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        std::fprintf(stderr, "Run 'pricer_cli --help' for usage.\n");
        return 1;
    }
}
