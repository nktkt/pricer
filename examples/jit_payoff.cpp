// jit_payoff.cpp — compile a payoff *formula* to native code at runtime with LLVM.
//
// What it does:
//   1. Parse a payoff formula string,
//   2. Emit an LLVM IR function  double payoff(double ST, double K),
//   3. JIT-compile it to native code with LLVM ORC,
//   4. Take the function pointer and call it from a Monte Carlo loop.
//
// Change the formula and you price a different instrument — without recompiling
// the C++. This is the seed of the project's "payoff DSL" (see ROADMAP Phase 2).
//
// Grammar: numbers, variables ST / K, operators + - * /, unary minus,
//          parentheses, max(a,b), min(a,b).
// Examples: "max(ST - K, 0)"            (call)
//           "max(K - ST, 0)"            (put)
//           "max(ST-K,0)+max(K-ST,0)"   (straddle)

#include "pricer/black_scholes.hpp"
#include <cstdio>
#include <string>
#include <map>
#include <memory>
#include <random>
#include <chrono>
#include <stdexcept>

#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Verifier.h"
#include "llvm/Support/Error.h"
#include "llvm/Support/TargetSelect.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/ExecutionEngine/Orc/LLJIT.h"
#include "llvm/ExecutionEngine/Orc/ThreadSafeModule.h"

using namespace llvm;
using namespace pricer;

// ---- recursive-descent compiler from a formula string to LLVM IR ----
class PayoffCompiler {
public:
    PayoffCompiler(IRBuilder<>& b, LLVMContext& c, std::map<std::string, Value*>& vars)
        : B(b), Ctx(c), Vars(vars) {}

    Value* compile(const std::string& src) {
        S = src; pos = 0;
        Value* v = parseExpr();
        skipWs();
        if (pos != S.size()) throw std::runtime_error("trailing characters: " + S.substr(pos));
        return v;
    }

private:
    IRBuilder<>& B;
    LLVMContext& Ctx;
    std::map<std::string, Value*>& Vars;
    std::string S;
    size_t pos = 0;

    void skipWs() { while (pos < S.size() && std::isspace((unsigned char)S[pos])) pos++; }
    char peek() { skipWs(); return pos < S.size() ? S[pos] : '\0'; }
    Constant* dbl(double x) { return ConstantFP::get(Ctx, APFloat(x)); }

    Value* parseExpr() {                              // expr := term (('+'|'-') term)*
        Value* lhs = parseTerm();
        while (true) {
            char c = peek();
            if (c == '+') { pos++; lhs = B.CreateFAdd(lhs, parseTerm()); }
            else if (c == '-') { pos++; lhs = B.CreateFSub(lhs, parseTerm()); }
            else return lhs;
        }
    }
    Value* parseTerm() {                              // term := factor (('*'|'/') factor)*
        Value* lhs = parseFactor();
        while (true) {
            char c = peek();
            if (c == '*') { pos++; lhs = B.CreateFMul(lhs, parseFactor()); }
            else if (c == '/') { pos++; lhs = B.CreateFDiv(lhs, parseFactor()); }
            else return lhs;
        }
    }
    Value* parseFactor() {                            // factor := '-' factor | primary
        if (peek() == '-') { pos++; return B.CreateFNeg(parseFactor()); }
        return parsePrimary();
    }
    Value* parsePrimary() {                           // primary := num | '(' expr ')' | ident[ '(' args ')' ]
        char c = peek();
        if (c == '(') {
            pos++; Value* v = parseExpr();
            if (peek() != ')') throw std::runtime_error("missing ')'");
            pos++; return v;
        }
        if (std::isdigit((unsigned char)c) || c == '.') return parseNumber();
        if (std::isalpha((unsigned char)c)) return parseIdent();
        throw std::runtime_error(std::string("unexpected character: ") + c);
    }
    Value* parseNumber() {
        skipWs();
        size_t start = pos;
        while (pos < S.size() && (std::isdigit((unsigned char)S[pos]) || S[pos] == '.')) pos++;
        return dbl(std::stod(S.substr(start, pos - start)));
    }
    Value* parseIdent() {
        skipWs();
        size_t start = pos;
        while (pos < S.size() && std::isalnum((unsigned char)S[pos])) pos++;
        std::string name = S.substr(start, pos - start);

        if (peek() == '(') {                          // function call: max / min
            pos++;
            Value* a = parseExpr();
            if (peek() != ',') throw std::runtime_error(name + " takes 2 arguments");
            pos++;
            Value* b = parseExpr();
            if (peek() != ')') throw std::runtime_error("missing ')'");
            pos++;
            if (name == "max") return B.CreateSelect(B.CreateFCmpOGT(a, b), a, b);
            if (name == "min") return B.CreateSelect(B.CreateFCmpOLT(a, b), a, b);
            throw std::runtime_error("unknown function: " + name);
        }
        auto it = Vars.find(name);                    // variable: ST / K
        if (it == Vars.end()) throw std::runtime_error("unknown variable: " + name);
        return it->second;
    }
};

// Monte Carlo using a JIT-compiled payoff(ST, K).
static double mc_with_jit(double (*payoff)(double, double),
                          double S, double K, double r, double sigma, double T, long n) {
    std::mt19937_64 rng(777);
    std::normal_distribution<double> Z(0.0, 1.0);
    const double drift = (r - 0.5 * sigma * sigma) * T, vol = sigma * std::sqrt(T);
    double sum = 0.0;
    for (long i = 0; i < n; ++i) {
        const double ST = S * std::exp(drift + vol * Z(rng));
        sum += payoff(ST, K);                         // call the runtime-compiled native fn
    }
    return std::exp(-r * T) * (sum / n);
}

int main(int argc, char** argv) {
    const std::string expr = (argc > 1) ? argv[1] : "max(ST - K, 0)";

    InitializeNativeTarget();
    InitializeNativeTargetAsmPrinter();
    ExitOnError ExitOnErr;

    auto jit = ExitOnErr(orc::LLJITBuilder().create());

    // Build IR: define double @payoff(double %ST, double %K)
    auto ctx = std::make_unique<LLVMContext>();
    auto mod = std::make_unique<Module>("payoff_module", *ctx);
    Type* dty = Type::getDoubleTy(*ctx);
    FunctionType* fty = FunctionType::get(dty, {dty, dty}, false);
    Function* fn = Function::Create(fty, Function::ExternalLinkage, "payoff", mod.get());
    fn->getArg(0)->setName("ST");
    fn->getArg(1)->setName("K");

    IRBuilder<> builder(BasicBlock::Create(*ctx, "entry", fn));
    std::map<std::string, Value*> vars{{"ST", fn->getArg(0)}, {"K", fn->getArg(1)}};

    try {
        PayoffCompiler compiler(builder, *ctx, vars);
        builder.CreateRet(compiler.compile(expr));
    } catch (const std::exception& e) {
        std::fprintf(stderr, "formula error: %s\n", e.what());
        return 1;
    }
    if (verifyFunction(*fn, &errs())) { std::fprintf(stderr, "IR verification failed\n"); return 1; }

    std::printf("formula : \"%s\"\n", expr.c_str());
    std::printf("\n--- generated LLVM IR ---\n");
    mod->print(outs(), nullptr);
    std::printf("-------------------------\n\n");

    ExitOnErr(jit->addIRModule(orc::ThreadSafeModule(std::move(mod), std::move(ctx))));
    auto sym = ExitOnErr(jit->lookup("payoff"));
    auto payoff = sym.toPtr<double (*)(double, double)>();

    std::printf("JIT-compiled. payoff(ST=120, K=100) = %.4f\n", payoff(120.0, 100.0));

    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long n = 10'000'000;
    const auto t0 = std::chrono::high_resolution_clock::now();
    const double price = mc_with_jit(payoff, S, K, r, sigma, T, n);
    const auto t1 = std::chrono::high_resolution_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::printf("\nMonte Carlo price (%ld paths) = %.6f  (%.1f ms)\n", n, price, ms);
    if (expr == "max(ST - K, 0)")
        std::printf("analytic (BS call)            = %.6f  <- should match\n",
                    black_scholes_call(S, K, r, sigma, T));
    return 0;
}

/*  Built automatically by CMake when LLVM is found (see examples/CMakeLists.txt).
    Manual build:
      clang++ -std=c++17 -O2 -Iinclude examples/jit_payoff.cpp \
        $(llvm-config --cxxflags --ldflags --libs core orcjit native --system-libs) \
        -fexceptions -o jit_payoff
*/
