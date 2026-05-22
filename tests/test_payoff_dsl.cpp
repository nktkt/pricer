// Tests for the LLVM-JIT payoff DSL (pricer::PayoffJit).
// Compiles formulas to native code and checks the returned values.
#include "check.hpp"
#include "pricer/payoff_jit.hpp"
#include <vector>

using namespace pricer;

int main() {
    PayoffJit jit;
    const std::vector<std::string> stk = {"ST", "K"};

    // Arithmetic and the call payoff.
    auto call = jit.compile("max(ST - K, 0)", stk);
    { double v[] = {120, 100}; check::approx("call ITM", call(v), 20.0, 1e-12); }
    { double v[] = {80, 100};  check::approx("call OTM", call(v), 0.0,  1e-12); }

    // Math intrinsics.
    auto idf = jit.compile("exp(log(ST))", {"ST"});
    { double v[] = {50}; check::approx("exp(log(x))", idf(v), 50.0, 1e-9); }
    auto rt = jit.compile("sqrt(ST)", {"ST"});
    { double v[] = {16}; check::approx("sqrt(16)", rt(v), 4.0, 1e-12); }
    auto pw = jit.compile("pow(ST, 3)", {"ST"});
    { double v[] = {2}; check::approx("pow(2,3)", pw(v), 8.0, 1e-12); }
    auto ab = jit.compile("abs(ST - K)", stk);
    { double v[] = {80, 100}; check::approx("abs(-20)", ab(v), 20.0, 1e-12); }

    // Comparisons yield 1.0 / 0.0 and compose arithmetically (digital).
    auto digi = jit.compile("(ST > K) * 10", stk);
    { double v[] = {120, 100}; check::approx("digital ITM", digi(v), 10.0, 1e-12); }
    { double v[] = {80, 100};  check::approx("digital OTM", digi(v), 0.0,  1e-12); }

    // Conditional builtin.
    auto cond = jit.compile("if(ST > K, ST - K, 0)", stk);
    { double v[] = {130, 100}; check::approx("if true",  cond(v), 30.0, 1e-12); }
    { double v[] = {90, 100};  check::approx("if false", cond(v), 0.0,  1e-12); }

    // Malformed input must throw rather than misbehave.
    bool threw = false;
    try { jit.compile("max(ST - K)", stk); } catch (const std::exception&) { threw = true; }
    check::is_true("bad arity throws", threw);

    // Compiled-kernel cache: recompiling a formula returns the same pointer and
    // does not trigger a second JIT compile.
    const unsigned before = jit.compiles();
    auto c1 = jit.compile("ST + K * 2", stk);
    auto c2 = jit.compile("ST + K * 2", stk);
    check::is_true("cache returns same fn", c1 == c2);
    check::is_true("cache avoids recompile", jit.compiles() - before == 1);

    // Vectorized (batch) kernel must match the scalar result lane-for-lane.
    const unsigned Wd = 4;
    auto sc = jit.compile("max(ST - K, 0)", stk);
    auto bt = jit.compile_batch("max(ST - K, 0)", stk, Wd);
    double v[8] = {120, 80, 100, 150,   100, 100, 100, 100};  // [ST x4][K x4]
    double out[4];
    bt(v, out);
    for (unsigned l = 0; l < Wd; ++l) {
        double sv[2] = {v[l], v[Wd + l]};
        char nm[32]; std::snprintf(nm, sizeof nm, "batch lane %u", l);
        check::approx(nm, out[l], sc(sv), 1e-12);
    }

    return check::report("payoff_dsl");
}
