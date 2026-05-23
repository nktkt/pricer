// bindings.cpp — Python bindings for the pricer C++ core (via pybind11).
//
// Exposes pricing, Greeks, implied volatility, Monte Carlo (serial / parallel /
// quasi-MC) and portfolio risk so a quant can price a book from Python without
// touching C++. Monte Carlo payoffs are built in C++ (call/put), so no Python
// callback is invoked per path.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // std::vector <-> list

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "pricer/black_scholes.hpp"
#include "pricer/curve.hpp"
#include "pricer/implied_vol.hpp"
#include "pricer/monte_carlo.hpp"
#include "pricer/parallel.hpp"
#include "pricer/parallel_simd.hpp"
#include "pricer/portfolio.hpp"
#include "pricer/qmc.hpp"
#include "pricer/risk.hpp"
#include "pricer/simd_mc.hpp"
#include "pricer/xva.hpp"

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
    m.attr("__version__") = "0.6.0";

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

    // --- portfolio (book-level AAD risk) ---
    py::class_<Position>(m, "Position")
        .def(py::init([](OptionType type, double S, double K, double r, double sigma, double T,
                         double qty) { return Position{type, S, K, r, sigma, T, qty}; }),
             py::arg("type"), py::arg("S"), py::arg("K"), py::arg("r"), py::arg("sigma"),
             py::arg("T"), py::arg("qty") = 1.0)
        .def_readwrite("type", &Position::type)
        .def_readwrite("S", &Position::S)
        .def_readwrite("K", &Position::K)
        .def_readwrite("r", &Position::r)
        .def_readwrite("sigma", &Position::sigma)
        .def_readwrite("T", &Position::T)
        .def_readwrite("qty", &Position::qty);

    py::class_<BookGreeks>(m, "BookGreeks")
        .def_readonly("value", &BookGreeks::value)
        .def_readonly("delta", &BookGreeks::delta)  // per-position dValue/dS
        .def_readonly("vega", &BookGreeks::vega)     // per-position dValue/dsigma
        .def("__repr__", [](const BookGreeks& g) {
            return "BookGreeks(value=" + std::to_string(g.value) + ", n_positions=" +
                   std::to_string(g.delta.size()) + ")";
        });

    // --- term structure & credit (for xVA) ---
    py::class_<DiscountCurve>(m, "DiscountCurve")
        .def(py::init<std::vector<double>, std::vector<double>>(), py::arg("times"),
             py::arg("zeros"))
        .def("df", &DiscountCurve::df, py::arg("t"))
        .def("zero_rate", &DiscountCurve::zero_rate, py::arg("t"))
        .def("forward_rate", &DiscountCurve::forward_rate, py::arg("t1"), py::arg("t2"));

    py::class_<SurvivalCurve>(m, "SurvivalCurve")
        .def(py::init([](double hazard) { return SurvivalCurve{hazard}; }), py::arg("hazard"))
        .def_static("from_spread", &SurvivalCurve::from_spread, py::arg("spread"),
                    py::arg("recovery"))
        .def_readwrite("hazard", &SurvivalCurve::hazard)
        .def("survival", &SurvivalCurve::survival, py::arg("t"))
        .def("default_prob", &SurvivalCurve::default_prob, py::arg("a"), py::arg("b"));

    py::class_<ExposureProfile>(m, "ExposureProfile")
        .def_readonly("times", &ExposureProfile::times)
        .def_readonly("epe", &ExposureProfile::epe)
        .def_readonly("ene", &ExposureProfile::ene);

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

    // --- SIMD / multicore Monte Carlo (counter-based path generation) ---
    m.def("mc_price_simd",
          [](OptionType t, double S, double K, double r, double sig, double T, long n,
             std::uint64_t seed) {
              return mc::price_terminal_cb_simd(payoff_for(t, K), S, r, sig, T, n, seed);
          },
          py::arg("type"), py::arg("S"), py::arg("K"), py::arg("r"), py::arg("sigma"),
          py::arg("T"), py::arg("n_paths") = 1000000, py::arg("seed") = 12345,
          py::call_guard<py::gil_scoped_release>());

    m.def("mc_price_parallel_simd",
          [](OptionType t, double S, double K, double r, double sig, double T, long n,
             std::uint64_t seed, unsigned threads) {
              return mc::price_terminal_cb_parallel_simd(payoff_for(t, K), S, r, sig, T, n, seed,
                                                         threads);
          },
          py::arg("type"), py::arg("S"), py::arg("K"), py::arg("r"), py::arg("sigma"),
          py::arg("T"), py::arg("n_paths") = 1000000, py::arg("seed") = 12345,
          py::arg("threads") = 0, py::call_guard<py::gil_scoped_release>());

    // --- book risk: all positions' delta/vega from one reverse-mode AAD sweep ---
    m.def("book_greeks_aad", &book_greeks_aad, py::arg("book"));

    // --- risk ---
    m.def("var_es", &var_es, py::arg("pnl"), py::arg("confidence") = 0.99);

    // --- xVA: exposure simulation + CVA / DVA / BCVA ---
    m.def("european_exposure_profile",
          [](OptionType type, double S, double K, double r, double sigma, double T,
             std::vector<double> grid, long n_paths, std::uint64_t seed) {
              // Mark-to-market of a vanilla European option, computed in C++ (no
              // per-path Python callback): Black–Scholes value while alive, payoff at expiry.
              auto mtm = [type, K, r, sigma, T](double t, double St) {
                  if (T - t > 1e-8) return black_scholes_price(type, St, K, r, sigma, T - t);
                  return type == OptionType::Call ? std::max(St - K, 0.0) : std::max(K - St, 0.0);
              };
              return exposure_profile_gbm(mtm, S, r, sigma, grid, n_paths, seed);
          },
          py::arg("type"), py::arg("S"), py::arg("K"), py::arg("r"), py::arg("sigma"), py::arg("T"),
          py::arg("grid"), py::arg("n_paths") = 200000, py::arg("seed") = 12345,
          py::call_guard<py::gil_scoped_release>());

    m.def("cva", &cva, py::arg("exposure"), py::arg("counterparty"), py::arg("disc"),
          py::arg("recovery"));
    m.def("dva", &dva, py::arg("exposure"), py::arg("own"), py::arg("disc"),
          py::arg("recovery_own"));
    m.def("bcva", &bcva, py::arg("exposure"), py::arg("counterparty"), py::arg("own"),
          py::arg("disc"), py::arg("recovery_cp"), py::arg("recovery_own"));
}
