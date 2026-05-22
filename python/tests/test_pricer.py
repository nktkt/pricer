"""Smoke tests for the pricer Python bindings.

Run with the package importable, e.g. after a CMake build:
    PYTHONPATH=python python python/tests/test_pricer.py
or after `pip install .`:
    python python/tests/test_pricer.py
"""
import pricer
from pricer import OptionType


def approx(a, b, tol):
    assert abs(a - b) <= tol, f"{a} != {b} (tol {tol})"


def main():
    S, K, r, sigma, T = 100.0, 100.0, 0.05, 0.20, 1.0

    # Closed form
    call = pricer.black_scholes_call(S, K, r, sigma, T)
    put = pricer.black_scholes_put(S, K, r, sigma, T)
    approx(call, 10.450584, 1e-3)
    approx(call - put, S - K * pow(2.718281828459045, -r * T), 1e-6)  # put-call parity

    # Greeks
    g = pricer.black_scholes_greeks(OptionType.Call, S, K, r, sigma, T)
    approx(g.delta, 0.636831, 1e-4)
    approx(g.price, call, 1e-9)

    # Implied vol round-trip
    iv = pricer.implied_vol(OptionType.Call, call, S, K, r, T)
    approx(iv, sigma, 1e-5)

    # Monte Carlo (serial, parallel, quasi)
    approx(pricer.mc_price(OptionType.Call, S, K, r, sigma, T, n_paths=2_000_000), call, 0.05)
    approx(pricer.mc_price_parallel(OptionType.Call, S, K, r, sigma, T, n_paths=2_000_000, threads=4),
           call, 0.05)
    approx(pricer.qmc_price(OptionType.Call, S, K, r, sigma, T, n=1_000_000), call, 0.01)

    # Risk: VaR/ES on a tiny P&L sample
    rm = pricer.var_es([1.0, -2.0, -3.0, 0.5, -5.0, -1.0], 0.8)
    assert rm.es >= rm.var > 0.0

    print("all python tests passed; pricer", pricer.__version__)


if __name__ == "__main__":
    main()
