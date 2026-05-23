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

    # Continuous dividend yield q: dividend-adjusted put-call parity, and q lowers the call.
    e = 2.718281828459045
    cq = pricer.black_scholes_call(S, K, r, sigma, T, q=0.03)
    pq = pricer.black_scholes_put(S, K, r, sigma, T, q=0.03)
    approx(cq - pq, S * e ** (-0.03 * T) - K * e ** (-r * T), 1e-6)
    assert cq < call, "a dividend yield should lower the call value"
    approx(pricer.implied_vol(OptionType.Call, cq, S, K, r, T, q=0.03), sigma, 1e-5)

    # Monte Carlo (serial, parallel, quasi, SIMD, multicore+SIMD)
    approx(pricer.mc_price(OptionType.Call, S, K, r, sigma, T, n_paths=2_000_000), call, 0.05)
    approx(pricer.mc_price_parallel(OptionType.Call, S, K, r, sigma, T, n_paths=2_000_000, threads=4),
           call, 0.05)
    approx(pricer.qmc_price(OptionType.Call, S, K, r, sigma, T, n=1_000_000), call, 0.01)
    approx(pricer.mc_price_simd(OptionType.Call, S, K, r, sigma, T, n_paths=2_000_000), call, 0.05)
    # Multicore+SIMD is deterministic in the thread count (bit-identical).
    p_a = pricer.mc_price_parallel_simd(OptionType.Call, S, K, r, sigma, T, n_paths=2_000_000, threads=2)
    p_b = pricer.mc_price_parallel_simd(OptionType.Call, S, K, r, sigma, T, n_paths=2_000_000, threads=8)
    approx(p_a, call, 0.05)
    assert p_a == p_b, "parallel+SIMD must be bit-identical across thread counts"

    # Risk: VaR/ES on a tiny P&L sample
    rm = pricer.var_es([1.0, -2.0, -3.0, 0.5, -5.0, -1.0], 0.8)
    assert rm.es >= rm.var > 0.0

    # Book-level Greeks from one reverse-mode AAD sweep.
    book = [
        pricer.Position(OptionType.Call, 100.0, 100.0, 0.05, 0.20, 1.0, 10.0),
        pricer.Position(OptionType.Put, 90.0, 95.0, 0.05, 0.25, 0.5, -5.0),
    ]
    bg = pricer.book_greeks_aad(book)
    book_value = sum(
        p.qty * pricer.black_scholes_price(p.type, p.S, p.K, p.r, p.sigma, p.T) for p in book
    )
    approx(bg.value, book_value, 1e-9)
    g0 = pricer.black_scholes_greeks(OptionType.Call, 100.0, 100.0, 0.05, 0.20, 1.0)
    approx(bg.delta[0], 10.0 * g0.delta, 1e-7)  # one-pass delta == qty * analytic delta
    approx(bg.vega[0], 10.0 * g0.vega, 1e-7)

    # xVA: exposure profile -> CVA / DVA for a long European call.
    grid = [i / 12.0 for i in range(1, 13)]
    ep = pricer.european_exposure_profile(OptionType.Call, S, K, 0.03, sigma, T, grid, n_paths=200_000)
    assert len(ep.times) == len(grid) and max(ep.ene) <= 1e-12  # long option: no negative exposure
    disc = pricer.DiscountCurve([0.5, 1.0, 2.0], [0.03, 0.03, 0.03])
    cp = pricer.SurvivalCurve.from_spread(0.02, 0.4)
    assert pricer.cva(ep, cp, disc, 0.4) > 0.0
    approx(pricer.dva(ep, cp, disc, 0.4), 0.0, 1e-12)  # DVA == 0 for a long option

    # Path-dependent exotics: each Monte Carlo engine agrees with its closed form.
    from pricer import AverageType, BarrierType

    # Asian: geometric MC matches the exact discrete-geometric closed form; n=1 == BS.
    approx(pricer.geometric_asian_price(OptionType.Call, S, K, r, sigma, T, 1), call, 1e-9)
    geo = pricer.geometric_asian_price(OptionType.Call, S, K, r, sigma, T, 50)
    geo_mc = pricer.asian_price_mc(OptionType.Call, AverageType.Geometric, S, K, r, sigma, T,
                                   n_steps=50, n_paths=300_000)
    approx(geo_mc, geo, 0.05)

    # Barrier: knock-in + knock-out == vanilla, and MC (BGK-corrected) matches the analytic.
    ui = pricer.barrier_price(OptionType.Call, BarrierType.UpIn, S, K, 130.0, r, sigma, T)
    uo = pricer.barrier_price(OptionType.Call, BarrierType.UpOut, S, K, 130.0, r, sigma, T)
    approx(ui + uo, call, 1e-9)
    uo_mc = pricer.barrier_price_mc(OptionType.Call, BarrierType.UpOut, S, K, 130.0, r, sigma, T,
                                    n_steps=250, n_paths=400_000)
    approx(uo_mc, uo, 0.15)

    # Lookback (floating strike): MC matches the closed form and dominates the vanilla call.
    lb = pricer.lookback_floating_price(OptionType.Call, S, r, sigma, T)
    lb_mc = pricer.lookback_floating_price_mc(OptionType.Call, S, r, sigma, T,
                                              n_steps=250, n_paths=400_000)
    approx(lb_mc, lb, 0.30)
    assert lb > call, "a floating-strike lookback call dominates the vanilla call"

    # SABR: beta=1, nu=0 is a flat lognormal smile; calibration recovers the params.
    F, Tf = 100.0, 1.0
    flat = pricer.SabrParams(0.25, 1.0, 0.0, 0.0)
    approx(pricer.sabr_implied_vol(F, 80.0, Tf, flat), 0.25, 1e-10)
    truth = pricer.SabrParams(0.22, 0.5, -0.35, 0.45)
    ks = [70.0, 85.0, 100.0, 115.0, 130.0, 150.0]
    mkt = [pricer.sabr_implied_vol(F, k, Tf, truth) for k in ks]
    fit = pricer.calibrate_sabr(F, Tf, ks, mkt, beta=0.5)
    approx(fit.params.alpha, 0.22, 1e-3)
    approx(fit.params.rho, -0.35, 1e-3)
    approx(fit.params.nu, 0.45, 1e-3)
    # Forward put-call parity through the SABR/Black-76 pricer.
    sc = pricer.sabr_black_price(OptionType.Call, F, 100.0, Tf, 1.0, truth)
    sp = pricer.sabr_black_price(OptionType.Put, F, 100.0, Tf, 1.0, truth)
    approx(sc - sp, F - 100.0, 1e-9)
    # The SABR SDE Monte Carlo cross-checks the Hagan/Black-76 price (moderate tol).
    approx(pricer.sabr_price_mc(OptionType.Call, F, 100.0, Tf, 1.0, truth, n_paths=400_000),
           sc, 0.25)

    # Bermudan: one date == European; many dates approach the American value.
    eur_put = pricer.black_scholes_put(S, K, r, sigma, T)
    b1 = pricer.bermudan_lsm(OptionType.Put, S, K, r, sigma, pricer.equally_spaced_dates(T, 1),
                             n_paths=200_000)
    approx(b1, eur_put, 0.06)
    b50 = pricer.bermudan_lsm(OptionType.Put, S, K, r, sigma, pricer.equally_spaced_dates(T, 50),
                              n_paths=200_000)
    assert b50 > eur_put, "a Bermudan put with many dates is worth more than the European"
    # An irregular schedule is accepted too.
    b_irr = pricer.bermudan_lsm(OptionType.Put, S, K, r, sigma, [0.25, 0.5, 0.75, 1.0],
                                n_paths=200_000)
    assert eur_put < b_irr < b50 + 0.1

    # Multi-asset: geometric basket MC matches the closed form; n=1 == Black–Scholes.
    from pricer import AverageType as AT
    Sv, wv, sigv = [100.0, 95.0, 105.0], [1 / 3, 1 / 3, 1 / 3], [0.20, 0.25, 0.30]
    corr = [[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]]
    approx(pricer.geometric_basket_price(OptionType.Call, [100.0], [1.0], 100.0, r, [0.20],
                                         [[1.0]], T), call, 1e-9)
    gcf = pricer.geometric_basket_price(OptionType.Call, Sv, wv, 100.0, r, sigv, corr, T)
    gmc = pricer.basket_price_mc(OptionType.Call, AT.Geometric, Sv, wv, 100.0, r, sigv, corr, T,
                                 n_paths=400_000)
    approx(gmc, gcf, 0.10)

    # Spread / exchange: Kirk reduces to Margrabe at K=0, and MC matches Margrabe.
    marg = pricer.margrabe_exchange_price(100.0, 100.0, 0.20, 0.30, 0.4, T)
    approx(pricer.spread_kirk_price(OptionType.Call, 100.0, 100.0, 0.0, r, 0.20, 0.30, 0.4, T),
           marg, 1e-10)
    approx(pricer.spread_price_mc(OptionType.Call, 100.0, 100.0, 0.0, r, 0.20, 0.30, 0.4, T,
                                  n_paths=400_000), marg, 0.10)

    print("all python tests passed; pricer", pricer.__version__)


if __name__ == "__main__":
    main()
