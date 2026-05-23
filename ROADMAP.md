# Roadmap — `pricer`

A long-range plan to grow this from a set of teaching samples into a
**scalable, high-performance pricing & risk engine**.

The north star: *let a user describe a financial instrument as a formula and get
fast, correct prices and risk — on a laptop or across a cluster — without
rewriting C++.* The LLVM JIT payoff compiler is the seed of that idea.

---

## Vision

```
            ┌─────────────────────────────────────────────┐
            │            User-facing surfaces             │
            │   CLI   ·   Python API   ·   REST / gRPC     │
            └───────────────┬─────────────────────────────┘
                            │
            ┌───────────────▼─────────────────────────────┐
            │                 Core engine                  │
            │  Payoff DSL → AST → LLVM IR → JIT (native)   │
            │  Pricing methods: Analytic · MC · PDE        │
            │  Risk: Greeks (bump / AAD), VaR, xVA         │
            └───────────────┬─────────────────────────────┘
                            │
            ┌───────────────▼─────────────────────────────┐
            │              Compute backends                │
            │   Scalar · SIMD · Multicore · GPU · Cluster  │
            └─────────────────────────────────────────────┘
```

Design principles:

- **One IR, many backends.** Express a payoff once; run it scalar, vectorized,
  on GPU, or distributed.
- **Correctness is a feature.** Every method is cross-checked against an
  analytic solution or a trusted reference (e.g. QuantLib).
- **Pay only for what you use.** The library core has zero heavy dependencies;
  JIT/GPU/cluster are opt-in.
- **Reproducible numerics.** Seeded RNG, deterministic aggregation, pinned
  results in regression tests.

---

## Phases

### Phase 0 — Foundation ✅ (done)
The current educational samples: Black–Scholes, Monte Carlo, convergence,
parallelism, Greeks, a barrier option, and the LLVM JIT payoff compiler.

### Phase 1 — Core library (0–3 months)
Turn scripts into a reusable library.
- [ ] CMake build; split `include/` (public API) and `src/`
- [ ] `Instrument` / `Model` / `PricingEngine` abstractions
- [ ] Pluggable RNG interface (Mersenne, PCG, Philox)
- [ ] Unit tests (GoogleTest/Catch2) + GitHub Actions CI (Linux + macOS)
- [ ] Numerical regression suite with pinned tolerances
- **Exit criteria:** `find_package(pricer)` works; CI green; vanilla options
  priced via a stable public API.

### Phase 2 — Payoff DSL & JIT engine (3–6 months)
Make the JIT compiler the real core.
- [x] Formal grammar + tokenizer/parser → typed AST
      *(payoff_ast.hpp: tokenizer + parser + typed AST + tree-walking interpreter; JIT now walks the AST)*
- [x] Functions: `exp log sqrt pow abs`, comparisons, `if/else`
- [x] Path-aware payoffs (time-indexed `S[t]`, multi-asset `S1, S2, …`)
- [x] Vectorized codegen (LLVM vector types → SIMD kernels)
- [x] Compiled-kernel cache keyed by formula hash
- **Exit criteria:** barrier/asian/basket payoffs expressed in the DSL match
  hand-written C++ to floating-point tolerance.

### Phase 3 — Performance & scale-up (6–9 months)
Single-node speed.
- [x] Thread-pool MC with deterministic reduction
- [x] SIMD path generation (AVX2/AVX-512/NEON via LLVM)
      *(simd.hpp: GCC/Clang vector extensions → AVX2 on x86, NEON on ARM; vectorized RNG, inverse-normal CDF and exp generate W=4 paths per step — ~2.2× over the stateful baseline, ~1.3× over the scalar counter-based engine)*
- [ ] GPU backend (LLVM NVPTX / SPIR-V) sharing the same payoff IR
- [x] Variance reduction: antithetic, control variates, Sobol QMC
- **Exit criteria:** ≥10× over the Phase 1 baseline on CPU; GPU path for the
  vanilla MC engine.
      *(≥10× CPU half DONE: `parallel_simd.hpp` stacks SIMD × multicore into one
      deterministic engine — `examples/scale_benchmark.cpp` measures ~12.8× over
      the Phase-1 stateful-mt19937 baseline on a 10-core machine. GPU path still
      open, so the phase stays open.)*

### Phase 4 — Risk & calibration (9–12 months)
From price to risk.
- [x] Greeks via adjoint algorithmic differentiation (AAD)
      *(tape-based reverse-mode AAD: first-order Greeks in one backward sweep; gamma via forward AD; 2nd-order adjoint is future)*
- [x] Yield-curve and volatility-surface construction
- [x] Model calibration (Heston, local vol) with least-squares solvers
      *(LM solver + SVI smile + Heston char-function calibration + Dupire local vol)*
- [x] Portfolio risk: VaR/ES, scenario engine, basic xVA
      *(VaR/ES `risk.hpp`; `xva.hpp` — GBM exposure-simulation scenario engine + CVA/DVA/BCVA against a hazard-rate survival curve)*
- **Exit criteria:** full Greeks for a multi-asset book in one pass; calibration
  reproduces market quotes within tolerance. ✓
      *(`portfolio.hpp` `book_greeks_aad` returns a multi-name book's value and every position's delta/vega from ONE reverse-mode sweep — matches closed-form to ~1e-13; Heston/SVI calibration reproduces quotes within tolerance.)*

### Phase 5 — Productization (12–18 months)
Make it usable by others.
- [x] Python bindings (pybind11/nanobind) + `pip install pricer`
- [x] CLI and REST/gRPC service with a job API
      *(pricer_cli + a POSIX-socket REST service with an async job API; gRPC not done)*
- [x] Market-data adapters and result persistence
      *(CSV quote/curve loaders + CSV result persistence; extensible to feeds/DBs)*
- [x] Versioned releases (SemVer), docs site, examples gallery
      *(v0.1.0 release + CHANGELOG; Doxygen docs site live at nktkt.github.io/pricer; examples/ dir as the gallery)*
- **Exit criteria:** a quant prices a book from Python without touching C++.  ✓

### Phase 6 — Distributed scale-out (18–24+ months)
Beyond one machine.
- [ ] Distributed MC (MPI / Ray / Kubernetes) with sharded path generation
      *(block-sharded MC + a multi-process driver landed; cross-node MPI/Ray/k8s transport pending)*
- [ ] Work scheduling, fault tolerance, deterministic global aggregation
      *(deterministic global aggregation done — bit-identical for any worker count; scheduling/fault-tolerance pending)*
- [ ] Observability: metrics, tracing, cost/perf dashboards
      *(REST server exposes Prometheus `/metrics` — request/latency/job/uptime counters — plus a structured JSON access log; distributed tracing and dashboards pending)*
- [ ] Cloud-native deployment (containers, autoscaling)
      *(multi-stage `Dockerfile` + `docker-compose.yml` ship the REST service as a slim non-root container, built and endpoint-checked by a `docker` CI job; orchestration/autoscaling pending)*
- **Exit criteria:** linear scaling of a large MC job across N nodes with
  reproducible results.

---

## Cross-cutting concerns

| Area | Commitment |
|------|------------|
| **Testing** | Unit + numerical regression + property-based tests on every PR |
| **Benchmarking** | Tracked microbenchmarks; perf regressions block merge |
| **Accuracy** | Reference checks vs. analytic solutions and QuantLib |
| **Docs** | Doxygen API + tutorials; every public symbol documented |
| **Security** | Dependency scanning; no untrusted code in JIT without sandboxing |
| **Governance** | Conventional Commits, semantic versioning, CHANGELOG |

---

## Success metrics

- **Speed:** vanilla MC throughput (paths/sec) per core, and total on GPU/cluster
- **Accuracy:** max relative error vs. reference across the instrument suite
- **Productivity:** lines of user code to price a new exotic (target: a one-line
  formula)
- **Scale:** parallel efficiency (%) at 1 → N cores → N nodes
- **Adoption:** PyPI installs, reproducible example notebooks

---

## Near-term next steps

1. Introduce CMake + CI and a first unit test around `black_scholes_call`.
2. Extract a `pricer::mc::price(payoff, model, n_paths)` API from the samples.
3. Grow the DSL with `exp/log/sqrt` and add an Asian-option example.
4. Add a benchmark harness so every later optimization is measured, not guessed.

> Timelines are directional, not commitments — they describe sequencing and
> dependencies, and will shift with scope and contributors.
