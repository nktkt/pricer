# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/) and the project aims for
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
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

[Unreleased]: https://github.com/nktkt/pricer/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/nktkt/pricer/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/nktkt/pricer/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/nktkt/pricer/releases/tag/v0.1.0
