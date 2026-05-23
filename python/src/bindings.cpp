// bindings.cpp — Python bindings for the pricer C++ core (via pybind11).
//
// Exposes pricing, Greeks, implied volatility, Monte Carlo (serial / parallel /
// quasi-MC) and portfolio risk so a quant can price a book from Python without
// touching C++. Monte Carlo payoffs are built in C++ (call/put), so no Python
// callback is invoked per path.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // std::vector <-> list

#include <cstdint>
#include <string>

#include "pricer/black_scholes.hpp"
#include "pricer/implied_vol.hpp"
#include "pricer/monte_carlo.hpp"
#include "pricer/parallel.hpp"
#include "pricer/qmc.hpp"
#include "pricer/risk.hpp"

namespace py = pybind11;
using namespace pricer;

// Vanilla call/put payoff as a C++ lambda (kept off the Python side for speed).
static auto payoff_for(OptionType t, double K) {
    return [t, K](double ST) {
        return (t == OptionType::Call) ? (ST > K ? ST - K : 0.0) : (ST < K ? K - ST : 0.0);
    };
}

PYBIND11_MODULE(_pricer, m) {
    m.doc() = "pricer — option pricing & risk engine (C++ core via pybind11)";
    m.attr("__version__") = "0.2.0";

    py::enum_<OptionType>(m, "OptionType")
        .value("Call", OptionType::Call)
        .value("Put", OptionType::Put);

    py::class_<Greeks>(m, "Greeks")
        .def_readonly("price", &Greeks::price)
        .def_readonly("delta", &Greeks::delta)
        .def_readonly("gamma", &Greeks::gamma)
        .def_readonly("vega", &Greeks::vega)
        .def_readonly("theta", &Greeks::theta)
        .def_readonly("rho", &Greeks::rho)
        .def("__repr__", [](const Greeks& g) {
            return "Greeks(price=" + std::to_string(g.price) + ", delta=" + std::to_string(g.delta) +
                   ", gamma=" + std::to_string(g.gamma) + ", vega=" + std::to_string(g.vega) +
                   ", theta=" + std::to_string(g.theta) + ", rho=" + std::to_string(g.rho) + ")";
        });

    py::class_<RiskMeasures>(m, "RiskMeasures")
        .def_readonly("var", &RiskMeasures::var)
        .def_readonly("es", &RiskMeasures::es)
        .def("__repr__", [](const RiskMeasures& r) {
            return "RiskMeasures(var=" + std::to_string(r.var) + ", es=" + std::to_string(r.es) + ")";
        });

    // --- closed form ---
    m.def("black_scholes_price", &black_scholes_price, py::arg("type"), py::arg("S"),
          py::arg("K"), py::arg("r"), py::arg("sigma"), py::arg("T"));
    m.def("black_scholes_call", &black_scholes_call, py::arg("S"), py::arg("K"), py::arg("r"),
          py::arg("sigma"), py::arg("T"));
    m.def("black_scholes_put", &black_scholes_put, py::arg("S"), py::arg("K"), py::arg("r"),
          py::arg("sigma"), py::arg("T"));
    m.def("black_scholes_greeks", &black_scholes_greeks, py::arg("type"), py::arg("S"),
          py::arg("K"), py::arg("r"), py::arg("sigma"), py::arg("T"));
    m.def("implied_vol", &implied_vol, py::arg("type"), py::arg("price"), py::arg("S"),
          py::arg("K"), py::arg("r"), py::arg("T"), py::arg("tol") = 1e-8,
          py::arg("max_iter") = 100);

    // --- Monte Carlo (vanilla call/put) ---
    m.def("mc_price",
          [](OptionType t, double S, double K, double r, double sig, double T, long n,
             std::uint64_t seed) {
              return mc::price_terminal(payoff_for(t, K), S, r, sig, T, n, seed);
          },
          py::arg("type"), py::arg("S"), py::arg("K"), py::arg("r"), py::arg("sigma"),
          py::arg("T"), py::arg("n_paths") = 1000000, py::arg("seed") = 12345);

    m.def("mc_price_parallel",
          [](OptionType t, double S, double K, double r, double sig, double T, long n,
             std::uint64_t seed, unsigned threads) {
              return mc::price_terminal_parallel(payoff_for(t, K), S, r, sig, T, n, seed, threads);
          },
          py::arg("type"), py::arg("S"), py::arg("K"), py::arg("r"), py::arg("sigma"),
          py::arg("T"), py::arg("n_paths") = 1000000, py::arg("seed") = 12345,
          py::arg("threads") = 0, py::call_guard<py::gil_scoped_release>());

    m.def("qmc_price",
          [](OptionType t, double S, double K, double r, double sig, double T, long n) {
              return mc::price_terminal_qmc(payoff_for(t, K), S, r, sig, T, n);
          },
          py::arg("type"), py::arg("S"), py::arg("K"), py::arg("r"), py::arg("sigma"),
          py::arg("T"), py::arg("n") = 1000000);

    // --- risk ---
    m.def("var_es", &var_es, py::arg("pnl"), py::arg("confidence") = 0.99);
}
