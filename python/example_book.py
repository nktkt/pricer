"""Price and risk a small option book from Python using the `pricer` package.

This script reads like a quant pricing a tiny book on a single underlying.
Everything below is driven from Python (no C++) via the `pricer` public API:

  * closed-form Black-Scholes prices and Greeks,
  * a Monte Carlo / Quasi-Monte Carlo cross-check of one price,
  * an implied-vol round-trip sanity check,
  * a 1-day 99% historical-style VaR/ES via a pure-Python Monte Carlo
    revaluation of the whole book.

Dependencies: standard library (math, random) plus `pricer`. No numpy.

Run:
  .venv/bin/python python/example_book.py
"""

import math
import random

import pricer
from pricer import OptionType


# ---------------------------------------------------------------------------
# Market and model parameters (single underlying)
# ---------------------------------------------------------------------------
S0 = 100.0       # spot price of the underlying today
R = 0.05         # continuously compounded risk-free rate
SIGMA = 0.20     # flat Black-Scholes volatility
SEED = 12345     # fixed seed -> reproducible Monte Carlo

# The book: a list of (option type, strike K, maturity T, signed quantity).
# Positive qty = long, negative qty = short.
POSITIONS = [
    (OptionType.Call, 100.0, 1.0, +10),
    (OptionType.Call, 110.0, 1.0, -5),
    (OptionType.Put, 90.0, 1.0, +8),
]


def type_name(option_type):
    """Return a short human-readable label for an OptionType."""
    return "Call" if option_type == OptionType.Call else "Put"


def book_value(spot, time_shift=0.0):
    """Mark-to-market value of the whole book.

    `time_shift` shortens every option's remaining maturity (used to
    revalue the book one day forward). Maturities are floored at a tiny
    positive number so we never ask for a non-positive time to expiry.
    """
    total = 0.0
    for option_type, strike, maturity, qty in POSITIONS:
        ttm = max(maturity - time_shift, 1e-8)
        unit = pricer.black_scholes_price(option_type, spot, strike, R, SIGMA, ttm)
        total += qty * unit
    return total


def main():
    # -------------------------------------------------------------------
    # 1) + 2) Per-position table, total book value and total delta
    # -------------------------------------------------------------------
    print("=" * 72)
    print("OPTION BOOK  (S0=%.2f, r=%.4f, sigma=%.4f)" % (S0, R, SIGMA))
    print("=" * 72)
    header = "%-5s %8s %5s %12s %14s %12s" % (
        "Type", "Strike", "Qty", "UnitPrice", "PosValue", "PosDelta",
    )
    print(header)
    print("-" * 72)

    total_value = 0.0
    total_delta = 0.0
    for option_type, strike, maturity, qty in POSITIONS:
        unit = pricer.black_scholes_price(option_type, S0, strike, R, SIGMA, maturity)
        greeks = pricer.black_scholes_greeks(option_type, S0, strike, R, SIGMA, maturity)
        pos_value = qty * unit
        pos_delta = greeks.delta * qty
        total_value += pos_value
        total_delta += pos_delta
        print("%-5s %8.2f %5d %12.6f %14.6f %12.6f" % (
            type_name(option_type), strike, qty, unit, pos_value, pos_delta,
        ))

    print("-" * 72)
    print("Total book value : %14.6f" % total_value)
    print("Total book delta : %14.6f" % total_delta)
    print()

    # -------------------------------------------------------------------
    # 3) Cross-check one option's price: closed form vs MC vs QMC
    # -------------------------------------------------------------------
    chk_type, chk_K, chk_T, _ = POSITIONS[0]
    cf = pricer.black_scholes_price(chk_type, S0, chk_K, R, SIGMA, chk_T)
    mc = pricer.mc_price(chk_type, S0, chk_K, R, SIGMA, chk_T,
                         n_paths=1000000, seed=SEED)
    qmc = pricer.qmc_price(chk_type, S0, chk_K, R, SIGMA, chk_T, n=1000000)
    print("Price cross-check for %s K=%.2f T=%.2f:" % (
        type_name(chk_type), chk_K, chk_T))
    print("  Closed form : %.6f" % cf)
    print("  Monte Carlo : %.6f  (diff %+.6f)" % (mc, mc - cf))
    print("  QMC         : %.6f  (diff %+.6f)" % (qmc, qmc - cf))
    print()

    # -------------------------------------------------------------------
    # 4) Implied-vol round-trip: price -> implied vol -> should recover SIGMA
    # -------------------------------------------------------------------
    iv = pricer.implied_vol(chk_type, cf, S0, chk_K, R, chk_T)
    print("Implied-vol round-trip for %s K=%.2f T=%.2f:" % (
        type_name(chk_type), chk_K, chk_T))
    print("  Input sigma     : %.6f" % SIGMA)
    print("  Recovered sigma : %.6f  (error %+.2e)" % (iv, iv - SIGMA))
    print()

    # -------------------------------------------------------------------
    # 5) 1-day 99% VaR / ES via a pure-Python Monte Carlo book revaluation
    # -------------------------------------------------------------------
    n_sims = 200000
    dt = 1.0 / 252.0
    drift = (R - 0.5 * SIGMA * SIGMA) * dt
    diffusion = SIGMA * math.sqrt(dt)

    v0 = book_value(S0, time_shift=0.0)
    rng = random.Random(SEED)

    pnl = []
    for _ in range(n_sims):
        z = rng.gauss(0.0, 1.0)                       # standard normal shock
        s1 = S0 * math.exp(drift + diffusion * z)     # one-day-ahead spot
        v1 = book_value(s1, time_shift=dt)            # revalue book one day on
        pnl.append(v1 - v0)                           # realised P&L

    risk = pricer.var_es(pnl, 0.99)
    print("1-day 99%% VaR / ES  (%d Monte Carlo paths, dt=1/252):" % n_sims)
    print("  Book value V0 : %14.6f" % v0)
    print("  VaR (99%%)     : %14.6f" % risk.var)
    print("  ES  (99%%)     : %14.6f" % risk.es)
    print("=" * 72)


if __name__ == "__main__":
    main()
