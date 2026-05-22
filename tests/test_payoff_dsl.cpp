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

    return check::report("payoff_dsl");
}
