// Tests for the payoff DSL tokenizer/parser/AST and interpreter (no LLVM).
#include "check.hpp"
#include "pricer/payoff_ast.hpp"
#include <map>
#include <string>

using namespace pricer;

static double ev(const std::string& f, std::map<std::string, double> env) {
    return ast::eval(ast::parse(f), env);
}

int main() {
    const std::map<std::string, double> mkt = {{"ST", 120}, {"K", 100}};

    // Arithmetic, precedence, parentheses, unary minus.
    check::approx("precedence", ev("1 + 2 * 3", {}), 7.0, 1e-12);
    check::approx("parens", ev("(1 + 2) * 3", {}), 9.0, 1e-12);
    check::approx("unary minus", ev("-ST + K", mkt), -20.0, 1e-12);

    // Variables and the call payoff.
    check::approx("call ITM", ev("max(ST - K, 0)", mkt), 20.0, 1e-12);
    check::approx("call OTM", ev("max(ST - K, 0)", {{"ST", 80}, {"K", 100}}), 0.0, 1e-12);

    // Math functions.
    check::approx("exp(log)", ev("exp(log(ST))", {{"ST", 50}}), 50.0, 1e-9);
    check::approx("sqrt", ev("sqrt(ST)", {{"ST", 16}}), 4.0, 1e-12);
    check::approx("pow", ev("pow(ST, 3)", {{"ST", 2}}), 8.0, 1e-12);
    check::approx("abs", ev("abs(ST - K)", {{"ST", 80}, {"K", 100}}), 20.0, 1e-12);

    // Comparisons yield 1/0 and compose; conditional builtin.
    check::approx("digital", ev("(ST > K) * 10", mkt), 10.0, 1e-12);
    check::approx("if true", ev("if(ST > K, ST - K, 0)", mkt), 20.0, 1e-12);
    check::approx("if false", ev("if(ST > K, ST - K, 0)", {{"ST", 90}, {"K", 100}}), 0.0, 1e-12);

    // Pretty-printer round-trips structurally.
    check::is_true("to_string", ast::to_string(*ast::parse("max(ST-K,0)")) == "max((ST - K), 0.000000)");

    // Malformed inputs must throw.
    auto throws = [](const std::string& f) {
        try { ast::parse(f); return false; } catch (const std::exception&) { return true; }
    };
    check::is_true("missing paren throws", throws("max(ST - K"));
    check::is_true("bad token throws", throws("ST $ K"));
    bool bad_arity = false;
    try { ev("max(ST)", mkt); } catch (const std::exception&) { bad_arity = true; }
    check::is_true("bad arity throws (eval)", bad_arity);

    return check::report("payoff_ast");
}
