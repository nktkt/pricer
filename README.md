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
examples/         runnable demos built on the library
tests/            CTest suite (dependency-free)
```

| Example | Topic | Highlight |
|---------|-------|-----------|
| `black_scholes_demo` | Black–Scholes vs. Monte Carlo | Two independent methods agree to ~0.05% |
| `convergence` | Accuracy vs. speed | Error shrinks like `1/sqrt(N)` (×100 paths → ~1/10 error) |
| `parallel_mc` | Multithreaded Monte Carlo | ~8× faster across 10 logical cores |
| `greeks` | Risk sensitivities (Greeks) | Closed-form vs. finite-difference cross-check |
| `barrier_option` | Path-dependent product | Up-and-out barrier call via stepped MC |
| `jit_payoff` | **LLVM JIT** payoff compiler | Parses a formula string → LLVM IR → native code at runtime |

## Requirements

- CMake 3.16+ and a C++17 compiler (`clang++` or `g++`)
- **LLVM** (optional, only for the `jit_payoff` example). On macOS: `brew install llvm`

## Build, test & run

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure   # run the test suite
./build/examples/black_scholes_demo          # run a demo
```

The `jit_payoff` example is built automatically when CMake finds LLVM (point it
there with `-DLLVM_DIR=$(llvm-config --cmakedir)` if needed); otherwise it is
skipped and the rest still builds. CI builds and tests on Linux and macOS.

## The LLVM JIT highlight

`jit_payoff` takes a payoff **formula as a string**, turns it into LLVM IR,
JIT-compiles it to native code at runtime, and calls it from a Monte Carlo loop.
Change the formula and you price a different instrument — without recompiling C++.

```sh
./build/examples/jit_payoff                            # call:  max(ST - K, 0)  (matches Black–Scholes)
./build/examples/jit_payoff "max(K - ST, 0)"          # put
./build/examples/jit_payoff "max(ST-K,0)+max(K-ST,0)" # straddle
```

For `max(ST - K, 0)` it generates and compiles:

```llvm
define double @payoff(double %ST, double %K) {
entry:
  %0 = fsub double %ST, %K
  %1 = fcmp ogt double %0, 0.000000e+00
  %2 = select i1 %1, double %0, double 0.000000e+00
  ret double %2
}
```

Supported grammar: numbers, variables `ST` / `K`, operators `+ - * /`, unary
minus, parentheses, and `max(a, b)` / `min(a, b)`.

> Note: the DSL parser uses exceptions for error reporting, so the `jit_payoff`
> target is compiled with `-fexceptions` (LLVM itself often ships built without).

## Disclaimer

Educational code. Not investment advice and not production-grade financial
software.

## License

MIT
