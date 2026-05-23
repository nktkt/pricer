// Robustness tests for the payoff DSL parser (payoff_ast.hpp).
//
// The JIT compiles user-supplied payoff formulas, so the parser is an untrusted
// input surface. Its contract: for ANY input string, parse() either returns a
// valid AST or throws std::exception — it must never crash, hang, or leak some
// other exception type. These tests assert that contract on random and
// adversarial inputs, and check that valid inputs round-trip and evaluate.
#include "check.hpp"
#include "pricer/payoff_ast.hpp"

#include <cstddef>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace pricer::ast;

namespace {

const std::map<std::string, double> kEnv{{"S", 100.0}, {"K", 90.0}};

// Parse (and, if it parses, also pretty-print and evaluate) the input. Returns
// true if the call completed "cleanly" — either a valid result or a thrown
// std::exception. Returns false only if a non-std::exception escaped.
bool handled_cleanly(const std::string& src) {
    try {
        NodePtr ast = parse(src);
        (void)to_string(*ast);
        try {
            (void)eval(*ast, kEnv);  // may throw (unknown var, etc.) — that is fine
        } catch (const std::exception&) {
        }
        return true;
    } catch (const std::exception&) {
        return true;  // clean rejection
    } catch (...) {
        return false;  // a non-std::exception leaked: contract violation
    }
}

// A small recursive generator of *valid* formulas (numbers, vars, binary ops,
// unary minus, function calls, parens) up to a bounded depth.
std::string gen_valid(std::mt19937& rng, int depth) {
    std::uniform_int_distribution<int> leaf_or_node(0, depth > 0 ? 6 : 1);
    switch (leaf_or_node(rng)) {
        case 0:
            return std::to_string(std::uniform_int_distribution<int>(0, 200)(rng));
        case 1:
            return (rng() & 1u) ? "S" : "K";
        case 2: {
            static const char* ops[] = {"+", "-", "*", "/", "<", ">", "<=", ">=", "==", "!="};
            const char* op = ops[std::uniform_int_distribution<int>(0, 9)(rng)];
            return "(" + gen_valid(rng, depth - 1) + op + gen_valid(rng, depth - 1) + ")";
        }
        case 3:
            return "(-" + gen_valid(rng, depth - 1) + ")";
        case 4: {
            static const char* f1[] = {"exp", "log", "sqrt", "abs"};
            return std::string(f1[std::uniform_int_distribution<int>(0, 3)(rng)]) + "(" +
                   gen_valid(rng, depth - 1) + ")";
        }
        case 5: {
            static const char* f2[] = {"max", "min", "pow"};
            return std::string(f2[std::uniform_int_distribution<int>(0, 2)(rng)]) + "(" +
                   gen_valid(rng, depth - 1) + "," + gen_valid(rng, depth - 1) + ")";
        }
        default:
            return "if(" + gen_valid(rng, depth - 1) + "," + gen_valid(rng, depth - 1) + "," +
                   gen_valid(rng, depth - 1) + ")";
    }
}

}  // namespace

int main() {
    std::mt19937 rng(20260523);

    // 1) Random byte soup from the DSL alphabet must never crash the parser.
    const std::string alphabet = "0123456789.+-*/()<>=!,abcSKxy_ exp log sqrt max min pow if";
    bool soup_ok = true;
    for (int t = 0; t < 20000; ++t) {
        std::string s;
        const int len = std::uniform_int_distribution<int>(0, 40)(rng);
        for (int k = 0; k < len; ++k)
            s += alphabet[std::uniform_int_distribution<size_t>(0, alphabet.size() - 1)(rng)];
        if (!handled_cleanly(s)) { soup_ok = false; break; }
    }
    check::is_true("random fuzz: no uncaught crash (20k inputs)", soup_ok);

    // 2) Syntactically malformed inputs must be rejected by parse() with a
    //    std::exception. (Semantic errors — bad arity, unknown names — are a
    //    separate layer, checked in step 2b.)
    const std::vector<std::string> bad = {
        "",   "(",    ")",      "1+",  "*2",     "max(", "(1",  "1 2",
        ",",  "<=",   "1<",     "a b", "1++2",   "()",   ".",   "/",
        "(()", "max(1,", "S K", "1)",  "(1+)",
    };
    bool bad_ok = true;
    for (const std::string& b : bad) {
        bool threw = false;
        try {
            parse(b);
        } catch (const std::exception&) {
            threw = true;
        } catch (...) {
            bad_ok = false;  // wrong exception type
        }
        if (!threw) { std::printf("  (not rejected: %s)\n", b.c_str()); bad_ok = false; }
    }
    check::is_true("syntax errors throw std::exception", bad_ok);

    // 2b) Syntactically valid but semantically invalid: parse() accepts them
    //     (they are well-formed expressions), but eval() must reject them with a
    //     std::exception — unknown variable, unknown function, or wrong arity.
    const std::vector<std::string> semantic_bad = {
        "exp", "unknownvar", "max()", "if(1,2)", "max(1,2,3)", "foo(1)", "sqrt(1,2)",
    };
    bool sem_ok = true;
    for (const std::string& s : semantic_bad) {
        bool threw = false;
        try {
            (void)eval(*parse(s), kEnv);
        } catch (const std::exception&) {
            threw = true;
        } catch (...) {
            sem_ok = false;
        }
        if (!threw) { std::printf("  (semantic error not caught: %s)\n", s.c_str()); sem_ok = false; }
    }
    check::is_true("semantic errors throw at eval", sem_ok);

    // 3) Adversarial deep nesting must throw cleanly (not overflow the stack).
    std::string deep_parens(5000, '('), deep_unary(5000, '-');
    deep_parens += "S";
    deep_parens += std::string(5000, ')');
    deep_unary += "S";
    bool deep_ok = true;
    for (const std::string& d : {deep_parens, deep_unary}) {
        try { parse(d); deep_ok = false; }       // should not succeed
        catch (const std::exception&) {}          // bounded-depth rejection: good
        catch (...) { deep_ok = false; }
    }
    check::is_true("deep nesting rejected without crashing", deep_ok);

    // 4) A long run of digits (std::stod out-of-range) is a clean lexer error.
    check::is_true("huge number literal rejected cleanly", handled_cleanly(std::string(500, '9')));

    // 5) Property: generated valid formulas parse, evaluate without throwing, and
    //    pretty-print idempotently (parse -> print -> parse -> print is stable).
    bool valid_ok = true, roundtrip_ok = true;
    for (int t = 0; t < 5000; ++t) {
        const std::string f = gen_valid(rng, 5);
        try {
            const NodePtr a = parse(f);
            (void)eval(*a, kEnv);  // all vars are S/K and all arities valid -> no throw
            const std::string p1 = to_string(*a);
            const std::string p2 = to_string(*parse(p1));
            if (p1 != p2) { roundtrip_ok = false; break; }
        } catch (const std::exception&) {
            valid_ok = false;
            std::printf("  (valid formula failed: %s)\n", f.c_str());
            break;
        }
    }
    check::is_true("generated valid formulas parse & eval", valid_ok);
    check::is_true("pretty-print round-trips (parse->print stable)", roundtrip_ok);

    // 6) A moderately deep (but in-bounds) nesting still parses fine.
    {
        std::string ok_deep(50, '('), tail(50, ')');
        ok_deep += "S";
        ok_deep += tail;
        bool ok = true;
        try { parse(ok_deep); } catch (...) { ok = false; }
        check::is_true("in-bounds nesting (50 deep) still parses", ok);
    }

    return check::report("payoff_fuzz");
}
