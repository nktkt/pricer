# pricer

Option pricing in modern C++ — from the Black–Scholes formula and Monte Carlo
simulation to a **runtime payoff compiler built on LLVM JIT**.

Each program is small, self-contained, and heavily commented. Together they walk
from the basics of pricing to the kind of just-in-time numerical engines used in
quantitative finance.

This is the seed of a larger goal: a **scalable pricing & risk engine** where you
describe an instrument as a formula and get fast, correct results — on a laptop
or across a cluster. See [`ROADMAP.md`](ROADMAP.md) for the long-range plan.

## Project layout

```
include/pricer/   header-only core library
  normal.hpp        standard-normal pdf / cdf
  black_scholes.hpp closed-form pricing + Greeks
  monte_carlo.hpp   generic terminal-value MC engine
  payoff_jit.hpp    payoff-formula → LLVM IR → native-code compiler (needs LLVM)
examples/         runnable demos built on the library
tests/            CTest suite (dependency-free; DSL test needs LLVM)
```

| Example | Topic | Highlight |
|---------|-------|-----------|
| `black_scholes_demo` | Black–Scholes vs. Monte Carlo | Two independent methods agree to ~0.05% |
| `convergence` | Accuracy vs. speed | Error shrinks like `1/sqrt(N)` (×100 paths → ~1/10 error) |
| `parallel_mc` | Multithreaded Monte Carlo | ~8× faster across 10 logical cores |
| `benchmark` | MC throughput harness | ~275 Mpaths/s; baseline for tracking optimizations |
| `greeks` | Risk sensitivities (Greeks) | Closed-form vs. finite-difference cross-check |
| `barrier_option` | Path-dependent product | Up-and-out barrier call via stepped MC |
| `jit_payoff` | **LLVM JIT** payoff compiler | Parses a formula string → LLVM IR → native code at runtime |
| `path_dependent` | Exotics from formulas | Asian / barrier / lookback / digital, each a one-line formula |

## Requirements

- CMake 3.16+ and a C++17 compiler (`clang++` or `g++`)
- **LLVM** (optional, only for the `jit_payoff` / `path_dependent` examples and the
  DSL test). On macOS: `brew install llvm`

## Build, test & run

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure   # run the test suite
./build/examples/black_scholes_demo          # run a demo
```

The JIT examples are built automatically when CMake finds LLVM (point it there
with `-DLLVM_DIR=$(llvm-config --cmakedir)` if needed); otherwise they are
skipped and the rest still builds. CI builds and tests on Linux and macOS.

## The LLVM JIT highlight

`jit_payoff` takes a payoff **formula as a string**, turns it into LLVM IR,
JIT-compiles it to native code at runtime, and calls it from a Monte Carlo loop.
Change the formula and you price a different instrument — without recompiling C++.

```sh
./build/examples/jit_payoff                            # call:  max(ST - K, 0)  (matches Black–Scholes)
./build/examples/jit_payoff "max(K - ST, 0)"          # put
./build/examples/jit_payoff "max(ST-K,0)+max(K-ST,0)" # straddle
./build/examples/jit_payoff "(ST > K) * 10"           # cash-or-nothing digital
```

For `max(ST - K, 0)` (variables read from the `double* v` argument) it generates:

```llvm
define double @payoff_0(ptr %v) {
entry:
  %0 = getelementptr inbounds double, ptr %v, i64 0
  %1 = load double, ptr %0, align 8        ; ST
  %2 = getelementptr inbounds double, ptr %v, i64 1
  %3 = load double, ptr %2, align 8        ; K
  %4 = fsub double %1, %3
  %5 = fcmp ogt double %4, 0.000000e+00
  %6 = select i1 %5, double %4, double 0.000000e+00
  ret double %6
}
```

`path_dependent` goes further: the engine simulates whole paths and exposes
`ST`, `avg`, `Smax`, `Smin` (plus `K`) to the DSL, so exotics become one-liners —
e.g. an arithmetic Asian call is `max(avg - K, 0)` and an up-and-out barrier is
`max(ST - K, 0) * (Smax < 130)`.

Supported grammar: numbers, variables (whichever names you bind), operators
`+ - * /` and unary minus, comparisons `< > <= >= == !=` (yielding `1.0`/`0.0`),
parentheses, and functions `exp log sqrt abs` (1 arg), `max min pow` (2 args),
`if(cond, a, b)` (3 args).

> Note: the DSL parser uses exceptions for error reporting, so the `jit_payoff`
> target is compiled with `-fexceptions` (LLVM itself often ships built without).

## Disclaimer

Educational code. Not investment advice and not production-grade financial
software.

## License

MIT
