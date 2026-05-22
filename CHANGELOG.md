# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/) and the project aims for
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- **REST pricing service** (`server/`, opt-in `-DPRICER_BUILD_SERVER=ON`): a
  dependency-free HTTP server (POSIX sockets) with `/price`, `/impliedvol`,
  `/mc`, and an async Monte Carlo job API (`/submit`, `/job`); a `server` CI job
  builds it and smoke-tests the endpoints on Linux/macOS.
- Vol-surface/curve/calibration stack: discount curve, implied-vol surface,
  quadratic & SVI smile calibration, Heston pricing + calibration (Levenberg–
  Marquardt), and exact Greeks by forward-mode automatic differentiation.

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

[Unreleased]: https://github.com/nktkt/pricer/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/nktkt/pricer/releases/tag/v0.1.0
