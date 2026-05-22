# cpp-llvm-option-pricing

Option pricing experiments in modern C++ — from the Black–Scholes formula and
Monte Carlo simulation to a **runtime payoff compiler built on LLVM JIT**.

Each program is small, self-contained, and heavily commented. Together they walk
from the basics of pricing to the kind of just-in-time numerical engines used in
quantitative finance.

## What's inside

| File | Topic | Highlight |
|------|-------|-----------|
| `option_pricing.cpp` | Black–Scholes vs. Monte Carlo | Two independent methods agree to ~0.05% |
| `convergence.cpp` | Accuracy vs. speed | Error shrinks like `1/sqrt(N)` (×100 paths → ~1/10 error) |
| `parallel_mc.cpp` | Multithreaded Monte Carlo | ~8× faster across 10 logical cores |
| `greeks.cpp` | Risk sensitivities (Greeks) | Closed-form vs. finite-difference cross-check |
| `barrier_option.cpp` | Path-dependent product | Up-and-out barrier call via stepped MC |
| `jit_payoff.cpp` | **LLVM JIT** payoff compiler | Parses a formula string → LLVM IR → native code at runtime |
| `bs_common.hpp` | Shared helpers | `norm_cdf`, `norm_pdf`, `black_scholes_call` |

## Requirements

- A C++17 compiler (`clang++` or `g++`)
- **LLVM 18+** (only for `jit_payoff`). On macOS: `brew install llvm`

The `Makefile` assumes Homebrew LLVM at `/opt/homebrew/opt/llvm`. Adjust the
`LLVM` variable if yours lives elsewhere.

## Build & run

```sh
make            # build everything
make run        # build, then run every sample in order
make clean      # remove binaries
```

Or build a single sample:

```sh
clang++ -std=c++17 -O2 -pthread convergence.cpp -o convergence && ./convergence
```

## The LLVM JIT highlight

`jit_payoff` takes a payoff **formula as a string**, turns it into LLVM IR,
JIT-compiles it to native code at runtime, and calls it from a Monte Carlo loop.
Change the formula and you price a different instrument — without recompiling C++.

```sh
./jit_payoff                              # call:     max(ST - K, 0)   (matches Black–Scholes)
./jit_payoff "max(K - ST, 0)"            # put
./jit_payoff "max(ST-K,0)+max(K-ST,0)"   # straddle
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

> Note: `llvm-config --cxxflags` emits `-fno-exceptions`; the build appends
> `-fexceptions` afterwards so the parser can throw on malformed input.

## Disclaimer

Educational code. Not investment advice and not production-grade financial
software.

## License

MIT
