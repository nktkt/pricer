# Architecture

How `pricer` is put together. For *what* it does see [`README.md`](README.md);
for *where it is going* see [`ROADMAP.md`](ROADMAP.md); for the per-symbol API see
the [Doxygen reference](https://nktkt.github.io/pricer/).

## Overview

```
        ┌──────────────────────────────────────────────────────────┐
        │                     User-facing surfaces                  │
        │   CLI (pricer_cli) · Python (pybind11) · REST server      │
        │                       · Docker image                      │
        └───────────────────────────┬──────────────────────────────┘
                                     │
        ┌───────────────────────────▼──────────────────────────────┐
        │                  Core library (header-only)               │
        │                                                            │
        │  Payoff DSL:  text → AST ──┬── tree-walking interpreter    │
        │                            └── LLVM IR → JIT (native)      │
        │                                                            │
        │  Pricing:     Black–Scholes · Monte Carlo · Heston ·       │
        │               American (binomial / LSM)                    │
        │  Risk:        Greeks (closed-form / AD / AAD) · VaR/ES ·   │
        │               xVA (CVA/DVA) · book-level one-pass Greeks   │
        │  Calibration: curves · vol surface · SVI · Heston · Dupire │
        └───────────────────────────┬──────────────────────────────┘
                                     │
        ┌───────────────────────────▼──────────────────────────────┐
        │                     Compute / scaling                     │
        │   scalar · SIMD (vector ext.) · multicore · multi-process │
        │            counter-based RNG → reproducible everywhere     │
        └──────────────────────────────────────────────────────────┘
```

The core is **header-only** (a CMake `INTERFACE` target): there is nothing to
link, and every layer is just `#include "pricer/<x>.hpp"`. Heavy or
platform-specific capabilities — the LLVM JIT, the Python bindings, the REST
server — are **opt-in** behind CMake options, so the base library has zero
third-party dependencies.

## Design principles

- **One IR, many backends.** A payoff is parsed once into a typed AST; an
  interpreter and an LLVM JIT consume the *same* AST. Pricing methods are
  templated on the payoff/number type, so the same code runs scalar, vectorized,
  or under automatic differentiation.
- **Reproducible numerics.** Randomness comes from a *counter-based* RNG: draw
  `i` is a pure function of `(seed, i)`, with no carried state. Combined with
  fixed work partitioning and fixed-order reduction, results are bit-identical
  regardless of how the work is split across SIMD lanes, threads, or processes.
- **Pay only for what you use.** JIT / Python / server are opt-in; the core
  pulls in nothing heavy.
- **Correctness is checked.** Every method is cross-checked against an analytic
  solution or a trusted reference, with a test per topic.

## The payoff DSL: one AST, two backends

This is the seed idea — describe an instrument as a formula, get native code.

```
"max(ST - K, 0)"
      │  tokenize + recursive-descent parse        (payoff_ast.hpp)
      ▼
   typed AST  ──────────────┬─────────────────────────────┐
                            │                              │
        eval(node, env)     │       Codegen walks the AST  │
   tree-walking interpreter │       → LLVM IR → ORC JIT    │
       (payoff_ast.hpp)     │          (payoff_jit.hpp)    │
                            ▼                              ▼
                      double result                native function ptr
```

- `payoff_ast.hpp` has **no LLVM dependency**: it tokenizes, parses, pretty-prints
  and interprets. It is the shared front end and is always available (and CI-tested
  without LLVM).
- `payoff_jit.hpp` reuses that AST and walks it to emit LLVM IR, JIT-compiling a
  scalar `double(const double*)` kernel or a vectorized `<W x double>` batch
  kernel, with a compiled-kernel cache keyed by `(width, variables, formula)`.
- Both backends validate the same way: the parser rejects **syntax** errors
  (with a depth-bounded recursion guard against adversarial nesting), while
  unknown variables / functions / arities are **semantic** errors raised by the
  interpreter (`eval`) and the JIT codegen.

## Core library layers

All under `include/pricer/`:

| Layer | Headers | Notes |
|-------|---------|-------|
| Math primitives | `normal.hpp` | standard-normal pdf/cdf |
| Closed-form pricing | `black_scholes.hpp` | price + Greeks struct |
| Monte Carlo | `monte_carlo.hpp` | generic terminal-value engine, templated on payoff |
| Early exercise | `american.hpp` | American options: CRR binomial tree + Longstaff–Schwartz LSM |
| Variance reduction | `variance_reduction.hpp`, `qmc.hpp` | antithetic / control variate; Sobol-style QMC + inverse-normal CDF |
| RNG & scaling | `rng.hpp`, `parallel.hpp`, `simd.hpp`, `simd_mc.hpp`, `parallel_simd.hpp`, `distributed.hpp` | counter-based RNG; deterministic thread pool; portable SIMD layer + vectorized MC; multicore×SIMD; block-sharded distributed MC |
| Risk | `greeks_mc.hpp`, `risk.hpp`, `xva.hpp` | MC Greeks (bump+CRN / pathwise); VaR/ES; exposure simulation + CVA/DVA/BCVA |
| Automatic differentiation | `dual.hpp`, `greeks_ad.hpp`, `adjoint.hpp`, `portfolio.hpp` | forward-mode duals; AD Greeks; reverse-mode AAD tape; book-level Greeks in one sweep |
| Curves & calibration | `curve.hpp`, `vol_surface.hpp`, `smile.hpp`, `optimize.hpp`, `svi.hpp`, `quadrature.hpp`, `heston.hpp`, `local_vol.hpp` | discount curve; vol surface; LM least-squares; SVI / Heston / Dupire |
| Market data | `csv.hpp`, `market_data.hpp` | CSV adapters + result persistence |
| Payoff DSL | `payoff_ast.hpp`, `payoff_jit.hpp` | front end + two backends (above) |

### Templated numerics → reuse across backends

The Black–Scholes price is written once as `bs_price_ad<Num>` over an abstract
number type. Instantiating `Num = double` gives the ordinary price; `Num = Dual`
gives forward-mode Greeks; `Num = Var` records onto the AAD tape. The same source
yields prices and exact Greeks with no duplication — and `portfolio.hpp` builds a
whole book on one tape so its value and *every* position's delta and vega come
out of a single reverse sweep.

### Determinism & scaling

`rng.hpp`'s `cb_normal(seed, i)` is stateless, so path generation is
embarrassingly parallel and order-independent. The scaling engines build on it:

- `parallel.hpp` / `parallel_simd.hpp` split work into a **fixed** number of
  blocks (independent of the thread count), sum each block, and reduce in block
  order → bit-identical for any number of threads.
- `simd_mc.hpp` generates `W` paths per step through the `simd.hpp` vector layer
  (vectorized RNG, inverse-normal CDF and `exp`).
- `distributed.hpp` shards blocks across processes with a canonical-order global
  reduction → identical price for any worker count.

## User-facing surfaces

- **CLI** — `examples/pricer_cli.cpp`: `price` / `iv` / `mc` subcommands.
- **Python** — `python/` (pybind11 + scikit-build-core, `pip install .`): exposes
  pricing, Greeks, implied vol, Monte Carlo (serial / parallel / QMC / SIMD /
  multicore+SIMD), VaR/ES, book-level AAD Greeks, and xVA.
- **REST server** — `server/pricer_server.cpp`: a dependency-free HTTP server on
  POSIX sockets with `/price`, `/impliedvol`, `/mc`, an async job API
  (`/submit`, `/job`), a Prometheus `/metrics` endpoint and a structured JSON
  access log.
- **Docker** — a multi-stage `Dockerfile` (slim, non-root runtime carrying only
  the server binary) and `docker-compose.yml`.

## Build, test & CI

- **Build**: CMake ≥ 3.16, C++17. The core is an `INTERFACE` target; options
  `PRICER_BUILD_EXAMPLES`, `PRICER_BUILD_TESTS`, `PRICER_ENABLE_JIT`,
  `PRICER_BUILD_PYTHON`, `PRICER_BUILD_SERVER` gate the rest. First-party targets
  use `-Wall -Wextra -Werror` via the `pricer_warnings` interface target.
- **Tests**: dependency-free executables using `tests/check.hpp`, registered with
  CTest (30 locally with JIT on, 29 in CI with JIT off).
- **CI** (`.github/workflows/`): a C++ build/test matrix and `python`, `server`
  and `docker` jobs on Linux + macOS, plus a `docs` workflow that publishes the
  Doxygen reference to GitHub Pages. JIT is disabled in CI (`-DPRICER_ENABLE_JIT=OFF`)
  because the hosted runners' LLVM cannot host a JIT; the JIT path is exercised
  locally.

## Repository map

```
include/pricer/   header-only core library (see the table above)
examples/         runnable demos, one per topic (also the examples gallery)
tests/            CTest suite (check.hpp + one test_*.cpp per topic)
python/           pybind11 bindings + package + tests
server/           REST service + smoke test
.github/workflows ci.yml (build/python/server/docker) + docs.yml
Dockerfile        multi-stage container build for the server
CMakeLists.txt    top-level build (core INTERFACE lib + options)
```
