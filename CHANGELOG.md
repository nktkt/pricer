# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/) and the project aims for
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.15.0] — 2026-05-23

### Added
- **Bachelier (normal) model** (`bachelier.hpp`): arithmetic Brownian motion on
  the forward (`F_T = F + σ_n·√T·Z`), so negative forwards and strikes are fine —
  the convention for interest-rate options and spreads, where Black–Scholes'
  lognormal assumption breaks down. The volatility `σ_n` is an absolute "normal"
  vol and the rate enters only through the discount factor `df`.
  - **Closed form** (`bachelier_price`): the exact call/put price under the normal
    model; the ATM price reduces to `df·σ_n·√T·φ(0)` and tends to the discounted
    intrinsic as `σ_n → 0`.
  - **Greeks** (`bachelier_greeks`): analytic delta, gamma and vega, with theta the
    time decay at a fixed `df` (rho is left 0 since the rate enters only via `df`).
  - **Normal implied vol** (`bachelier_implied_vol`): recovers the absolute (normal)
    vol from a price by a safeguarded Newton iteration seeded from the exact ATM
    inversion.
  - **Monte Carlo** (`bachelier_price_mc`): simulates the arithmetic forward (which
    can go negative) and cross-checks the closed form.
  - Anchored by the **validating relationships**: put–call parity
    (`C − P = df·(F − K)`), the ATM formula, the zero-vol intrinsic limit, a
    negative-rate floorlet, the Greeks vs. finite differences, an implied-vol
    round-trip, and the Monte Carlo cross-check.
  - Exposed through the Python bindings (`bachelier_price`, `bachelier_greeks`,
    `bachelier_implied_vol`, `bachelier_price_mc`). Covered by `test_bachelier`
    and demonstrated in `examples/bachelier_demo.cpp`.

## [0.14.0] — 2026-05-23

### Added
- **Rainbow (two-asset best-of / worst-of) options** (`rainbow.hpp`): options on
  the maximum or minimum of two correlated underlyings, with the Stulz (1982)
  closed form cross-checked against a Monte Carlo engine.
  - **Bivariate normal CDF** (`bivariate_normal_cdf`): the standard bivariate
    normal `M(a, b; ρ)` via Drezner's (1978) reduction to a 2-D Gauss quadrature,
    checked against its analytic special cases.
  - **Stulz closed form** (`rainbow_price`, with the `RainbowType` `Max` / `Min`
    selector): exact best-of (max) and worst-of (min) call formulas; the puts come
    from put–call parity using the discounted forward on the max/min (a vanilla
    leg plus/minus a Margrabe exchange option). Anchored by the **validating
    parity** that pins both call formulas at once — a call on the max plus a call
    on the min equals two vanilla calls (`C_max(K) + C_min(K) = c(S1, K) +
    c(S2, K)`) — with best-of dominating either single call, worst-of dominated by
    it, and the expected correlation effects (best-of falls / worst-of rises with
    correlation).
  - **Monte Carlo** (`rainbow_price_mc`): a correlated two-asset GBM engine that
    cross-checks the closed form for call/put on max/min.
  - Exposed through the Python bindings (`RainbowType`, `bivariate_normal_cdf`,
    `rainbow_price`, `rainbow_price_mc`); all carry the continuous dividend yields
    `q1` / `q2`. Covered by `test_rainbow` and demonstrated in
    `examples/rainbow_options.cpp`.

## [0.13.0] — 2026-05-23

### Added
- **Digital (binary) options** (`digital.hpp`): cash-or-nothing and
  asset-or-nothing options, each with a closed form and a counter-based Monte
  Carlo engine (`digital_price`, `digital_price_mc`) that cross-check.
  - **Cash-or-nothing** pays a fixed `cash` if in the money
    (`cash·e^{-rT}·N(±d2)`); **asset-or-nothing** pays the terminal `S_T` if in
    the money (`S·e^{-qT}·N(±d1)`).
  - Anchored by the **exact vanilla decomposition** (a vanilla call equals the
    asset-or-nothing call minus `K` cash-or-nothing calls with `cash = 1`), so the
    digitals reprice Black–Scholes exactly, plus **cash parity**
    (call + put = `e^{-rT}`), **asset parity** (call + put = `S·e^{-qT}`), and the
    **digital = −∂C/∂K** identity (a cash-or-nothing call is the negative strike
    derivative of the vanilla call).
  - Exposed through the Python bindings (`DigitalType`, `digital_price`,
    `digital_price_mc`); both carry the continuous dividend yield `q`. Covered by
    `test_digital` and demonstrated in `examples/digital_options.cpp`.

## [0.12.0] — 2026-05-23

### Added
- **Multi-asset options** (`basket.hpp`): options on several correlated
  underlyings, with correlated Geometric Brownian Motion driven by a Cholesky
  factor of the correlation matrix and the counter-based RNG.
  - **Basket** (weighted average of the assets): an exact lognormal closed form
    for the geometric-average basket (`geometric_basket_price` — the geometric
    average of correlated GBMs is itself lognormal) plus a Monte Carlo engine for
    arithmetic *and* geometric averaging (`basket_price_mc`). The geometric basket
    MC cross-checks the closed form, the arithmetic basket sits above the geometric
    one (AM ≥ GM), and the basket value rises with correlation.
  - **Spread / exchange** (two assets): Margrabe's exact exchange-option formula
    for `max(S1 − S2, 0)` (`margrabe_exchange_price`), Kirk's approximation for the
    spread `max(S1 − S2 − K, 0)` (`spread_kirk_price`, which reduces to Margrabe
    exactly at `K = 0`), and a correlated two-asset Monte Carlo engine
    (`spread_price_mc`). Margrabe matches Kirk(K=0) and the MC, Kirk tracks the MC
    for non-zero strikes, and the spread value falls as correlation rises.
  - Exposed through the Python bindings (`geometric_basket_price`,
    `basket_price_mc`, `margrabe_exchange_price`, `spread_kirk_price`,
    `spread_price_mc`); all carry the continuous dividend yield `q`. Covered by
    `test_basket` and demonstrated in `examples/basket_options.cpp`.

## [0.11.0] — 2026-05-23

### Added
- **Bermudan option pricing** (`bermudan.hpp`): `bermudan_lsm` prices an option
  exercisable on a finite, possibly irregular, schedule of dates by
  Longstaff–Schwartz least-squares Monte Carlo. It reuses the LSM regression
  ({1, S/K, (S/K)²} over the in-the-money paths) from `american.hpp` but over an
  arbitrary schedule, jumping the underlying exactly from one exercise date to the
  next (GBM is Markov, so no fine time grid is needed). A Bermudan sits between the
  European and American values: one exercise date reproduces the European price,
  many equally-spaced dates approach the American (binomial) price, and adding
  dates only raises the value — all cross-checked in `test_bermudan` and shown in
  `examples/bermudan_option.cpp`. A helper `equally_spaced_dates` builds the common
  uniform schedule, and the dividend yield `q` is supported. Exposed through the
  Python bindings (`bermudan_lsm`, `equally_spaced_dates`).

## [0.10.0] — 2026-05-23

### Added
- **SABR stochastic-volatility model** (`sabr.hpp`): the industry-standard model
  for interest-rate and FX smiles.
  - **Hagan 2002 implied-volatility approximation** (`sabr_implied_vol`): the
    closed-form Black (lognormal) vol for a strike on a forward, with the
    at-the-money branch as the continuous K → F limit; reduces to a flat smile
    `σ = α` when `β = 1, ν = 0`.
  - **Option pricing** (`sabr_black_price`): Black-76 on the forward using the
    SABR implied vol.
  - **Calibration** (`calibrate_sabr`): fits `(α, ρ, ν)` at a fixed `β` to a smile
    of market vols via the shared Levenberg–Marquardt solver, seeding `α` from the
    ATM vol (`α ≈ σ_ATM·F^(1-β)`) so it is robust across `β` scales — a model
    round-trip recovers the parameters to ~1e-12.
  - **Monte Carlo** (`sabr_price_mc`): simulates the SABR SDE directly (the vol
    process is lognormal so it is stepped exactly; the forward uses a
    full-truncation Euler step) with the counter-based RNG, cross-checking the
    asymptotic Hagan price.
  - Exposed through the Python bindings (`SabrParams`, `SabrFit`,
    `sabr_implied_vol`, `sabr_black_price`, `calibrate_sabr`, `sabr_price_mc`);
    `test_sabr` and `examples/sabr_smile.cpp` cover smile shape, calibration and
    the MC cross-check.

## [0.9.0] — 2026-05-23

### Added
- **Path-dependent exotic options** (`exotics.hpp`): Asian, barrier and lookback
  options, each priced two ways that cross-check.
  - **Asian** (average-price, fixed strike): an exact closed form for the
    discretely-monitored geometric average (the geometric average of GBM is
    lognormal) plus a Monte Carlo engine for arithmetic *and* geometric averaging.
  - **Barrier** (single barrier, continuous monitoring): the Reiner–Rubinstein
    closed form for all eight up/down × in/out × call/put cases (knock-in plus the
    matching knock-out equals the vanilla — parity), plus a Monte Carlo engine.
  - **Lookback** (floating strike): the Conze–Viswanathan / Goldman–Sosin–Gatto
    closed form plus a Monte Carlo engine.
  - The barrier and lookback Monte Carlo engines apply the
    **Broadie–Glasserman–Kou continuity correction** (shifting the monitored
    barrier / running extreme by `exp(±0.5826·σ·√dt)`) so discrete stepping
    converges to the continuous-monitoring closed form; cross-checked in
    `test_exotics` and demonstrated in `examples/exotic_options.cpp`.
  - Exposed through the Python bindings (`geometric_asian_price`, `asian_price_mc`,
    `barrier_price`, `barrier_price_mc`, `lookback_floating_price`,
    `lookback_floating_price_mc`, with `AverageType` / `BarrierType` enums). All
    carry the continuous dividend yield `q`.

## [0.8.0] — 2026-05-23

### Added
- **Continuous dividend yield (`q`)** (Merton extension) across the closed-form
  pricer and Greeks (closed-form, forward-AD and adjoint-AD), implied vol, the
  American binomial tree / LSM, and the Python and REST surfaces — backward
  compatible (`q` defaults to 0). Verified by dividend-adjusted put–call parity
  and the early-exercise premium an American call gains when `q` is high
  (`test_dividends`).
- **American option Greeks** (`american.hpp` `binomial_greeks`): delta, gamma and
  theta read from a single CRR tree pass (finite differences over adjacent
  lattice nodes), with vega and rho by central-difference re-pricing. The European
  mode matches the closed-form Black–Scholes Greeks; `examples/american_option.cpp`
  prints the American put's Greeks.

## [0.7.0] — 2026-05-23

### Added
- **American (early-exercise) option pricing** (`american.hpp`): a Cox–Ross–Rubinstein
  binomial tree and Longstaff–Schwartz least-squares Monte Carlo, cross-checked
  against each other and Black–Scholes (`examples/american_option.cpp`,
  `test_american`). Captures the early-exercise premium an American put carries
  over its European counterpart.

## [0.6.0] — 2026-05-23

### Added
- **Kubernetes manifests** (`k8s/`): Deployment, Service and a CPU
  HorizontalPodAutoscaler for the REST service — non-root pods, read-only root
  filesystem, `/health` probes, and Prometheus scrape annotations for `/metrics`
  — with a `k8s/README.md`.
- **Project docs**: `ARCHITECTURE.md` (how the layers, the one-AST/two-backend
  DSL, the determinism model and the surfaces fit together) and `CONTRIBUTING.md`
  (build/test, conventions, adding a feature, and the release process), linked
  from the README.

## [0.5.0] — 2026-05-23

### Added
- **Container deployment**: a multi-stage `Dockerfile` (slim, non-root runtime
  carrying only the server binary) and a `docker-compose.yml`, plus a `docker` CI
  job that builds the image and checks the container serves `/health`, `/price`
  and `/metrics`.
- **DSL parser robustness tests** (`test_payoff_fuzz`): property and fuzz tests
  asserting `parse()` never crashes on arbitrary input (20k random + adversarial
  cases), that syntax vs. semantic errors are rejected at the right layer, that
  valid formulas round-trip through the pretty-printer, and that deep nesting is
  bounded.

### Changed
- **Hardened the payoff parser** (an untrusted-input surface, since the JIT
  compiles user formulas): recursion is now depth-bounded, so adversarially
  nested input gets a clean parse error instead of a stack overflow, and a
  malformed number literal raises the parser's own error type rather than leaking
  `std::invalid_argument` / `std::out_of_range` from `std::stod`.

## [0.4.0] — 2026-05-23

### Added
- **Server observability**: the REST service exposes a Prometheus `/metrics`
  endpoint (request, response-class, per-endpoint latency, async-job and uptime
  counters) and writes a structured JSON access-log line to stderr for every
  request; the smoke test now covers `/metrics`.
- **Python bindings for the v0.3.0 features**: SIMD and multicore+SIMD Monte
  Carlo (`mc_price_simd`, `mc_price_parallel_simd`), book-level Greeks in one AAD
  sweep (`Position`, `BookGreeks`, `book_greeks_aad`), and xVA (`DiscountCurve`,
  `SurvivalCurve`, `ExposureProfile`, `european_exposure_profile`, `cva`, `dva`,
  `bcva`). `python/example_book.py` now shows one-pass book Greeks and CVA.

## [0.3.0] — 2026-05-23

### Added
- **SIMD path generation** (`simd.hpp`, `simd_mc.hpp`): a portable SIMD layer
  built on the GCC/Clang vector extensions (AVX2 on x86, NEON on ARM) with
  vectorized `exp`/`log`/`sqrt` and inverse-normal CDF, plus a vectorized
  counter-based Monte Carlo engine (`mc::price_terminal_cb_simd`) that generates
  W=4 paths per step. ~2.2× over the stateful baseline / ~1.3× over the scalar
  counter-based engine, with results matching to floating-point tolerance.
- **Multicore + SIMD engine** (`parallel_simd.hpp`):
  `mc::price_terminal_cb_parallel_simd` stacks SIMD path generation across all
  cores with a deterministic, thread-count-independent (bit-identical) result.
  `examples/scale_benchmark.cpp` shows the full speedup ladder — ~12.8× over the
  Phase-1 baseline on a 10-core machine, meeting the Phase-3 ≥10× CPU target.
- **xVA** (`xva.hpp`): a GBM exposure-simulation scenario engine plus CVA, DVA
  and bilateral CVA against a hazard-rate survival curve and a discount curve
  (`examples/xva_demo.cpp`).
- **Book-level Greeks in one AAD sweep** (`portfolio.hpp`): `book_greeks_aad`
  returns a multi-name book's value and every position's delta/vega from a single
  reverse-mode pass (`examples/portfolio_aad.cpp`), completing the Phase-4 exit
  criterion (full Greeks for a multi-asset book in one pass).

## [0.2.0] — 2026-05-23

### Added
- **Payoff DSL typed AST** (`payoff_ast.hpp`): tokenizer + parser + AST +
  tree-walking interpreter (no LLVM); the JIT now walks this shared AST.
- **Counter-based RNG** (`rng.hpp`) for reproducible, parallel/SIMD-friendly path
  generation (~1.9× faster than `std::mt19937_64`).
- **Distributed Monte Carlo** (`distributed.hpp`) with bit-identical results
  regardless of how blocks are sharded, plus a multi-process driver.
- **Reverse-mode AAD** (`adjoint.hpp`) and forward-mode AD (`dual.hpp`) for exact
  Greeks; **Dupire local volatility** (`local_vol.hpp`).
- **Calibration & vol surface**: Levenberg–Marquardt solver, SVI and Heston
  models (pricing + calibration), plus a discount curve, implied-vol surface,
  and quadratic smile fit.
- **REST pricing service** (`server/`, opt-in `-DPRICER_BUILD_SERVER=ON`): a
  dependency-free HTTP server (POSIX sockets) with `/price`, `/impliedvol`,
  `/mc`, and an async Monte Carlo job API (`/submit`, `/job`); a `server` CI job
  builds it and smoke-tests the endpoints on Linux/macOS.
- **CSV market-data adapters and result persistence** (`market_data.hpp`).

### Changed
- First-party targets now build with `-Wall -Wextra -Werror`; added boundary /
  error-handling tests.

## [0.1.0] — 2026-05-21

Initial public release, built up in phases.

### Added
- **Core library** (`include/pricer/`, header-only): standard-normal helpers,
  Black–Scholes pricing & Greeks, a generic terminal-value Monte Carlo engine.
- **Payoff DSL / LLVM JIT** (`payoff_jit.hpp`): compile a payoff formula string to
  native code at runtime; functions `exp/log/sqrt/abs/pow`, comparisons, `if`;
  path-aware variables; a compiled-kernel cache; and vectorized `<W x double>`
  (SIMD) codegen.
- **Performance & scale**: deterministic multithreaded Monte Carlo
  (`parallel.hpp`), variance reduction (`variance_reduction.hpp`: antithetic &
  control variates), and quasi-Monte Carlo (`qmc.hpp`).
- **Risk & calibration**: implied-volatility solver (`implied_vol.hpp`), Monte
  Carlo Greeks (`greeks_mc.hpp`: bump+CRN and pathwise), and Value-at-Risk /
  Expected Shortfall (`risk.hpp`).
- **Python bindings** (`pip install .`) via pybind11 + scikit-build-core:
  pricing, Greeks, implied volatility, Monte Carlo (serial/parallel/quasi-MC)
  and VaR/ES — usable from Python with no C++.
- **Command-line tool** `pricer_cli`: `price` (price + Greeks), `iv` (implied
  vol), and `mc` (parallel Monte Carlo) sub-commands.
- Examples for every topic, a CTest suite, a CMake build, and CI (C++ and
  Python) on Linux and macOS.

[Unreleased]: https://github.com/nktkt/pricer/compare/v0.15.0...HEAD
[0.15.0]: https://github.com/nktkt/pricer/compare/v0.14.0...v0.15.0
[0.14.0]: https://github.com/nktkt/pricer/compare/v0.13.0...v0.14.0
[0.13.0]: https://github.com/nktkt/pricer/compare/v0.12.0...v0.13.0
[0.12.0]: https://github.com/nktkt/pricer/compare/v0.11.0...v0.12.0
[0.11.0]: https://github.com/nktkt/pricer/compare/v0.10.0...v0.11.0
[0.10.0]: https://github.com/nktkt/pricer/compare/v0.9.0...v0.10.0
[0.9.0]: https://github.com/nktkt/pricer/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/nktkt/pricer/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/nktkt/pricer/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/nktkt/pricer/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/nktkt/pricer/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/nktkt/pricer/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/nktkt/pricer/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/nktkt/pricer/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/nktkt/pricer/releases/tag/v0.1.0
