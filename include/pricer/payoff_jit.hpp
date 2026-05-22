// pricer/payoff_jit.hpp — compile a payoff *formula* to native code with LLVM.
//
// `PayoffJit::compile(formula, varNames)` parses a formula string, emits an LLVM
// IR function `double payoff(const double* v)` that reads each named variable
// from `v[index]`, JIT-compiles it, and returns a raw function pointer. The
// caller fills `v` per Monte Carlo path, so the same engine prices any
// instrument expressible in the grammar — including path-dependent ones, by
// exposing path aggregates (terminal, average, max, min) as variables.
//
// Grammar (precedence low → high):
//   expr   := cmp
//   cmp    := add (('<'|'>'|'<='|'>='|'=='|'!=') add)?   // yields 1.0 / 0.0
//   add    := mul (('+'|'-') mul)*
//   mul    := factor (('*'|'/') factor)*
//   factor := '-' factor | primary
//   primary:= number | '(' expr ')' | ident | ident '(' args ')'
// Functions: exp log sqrt abs (1 arg), max min pow (2 args), if(cond,a,b) (3 args).
//
// This header requires LLVM; only targets that link LLVM should include it.
#pragma once
#include <cctype>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "llvm/ExecutionEngine/Orc/LLJIT.h"
#include "llvm/ExecutionEngine/Orc/ThreadSafeModule.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Intrinsics.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Verifier.h"
#include "llvm/Support/Error.h"
#include "llvm/Support/TargetSelect.h"
#include "llvm/Support/raw_ostream.h"

namespace pricer {

namespace detail {

// Recursive-descent translator: formula text -> LLVM IR values, reading
// variables from a `double*` function argument.
class PayoffParser {
public:
    PayoffParser(llvm::IRBuilder<>& b, llvm::LLVMContext& ctx, llvm::Value* args,
                 const std::map<std::string, unsigned>& vars)
        : B(b), Ctx(ctx), Args(args), Vars(vars) {}

    llvm::Value* run(const std::string& src) {
        S = src;
        pos = 0;
        llvm::Value* v = parseExpr();
        skipWs();
        if (pos != S.size()) err("trailing characters: '" + S.substr(pos) + "'");
        return v;
    }

private:
    llvm::IRBuilder<>& B;
    llvm::LLVMContext& Ctx;
    llvm::Value* Args;
    const std::map<std::string, unsigned>& Vars;
    std::string S;
    size_t pos = 0;

    llvm::Type* dty() { return llvm::Type::getDoubleTy(Ctx); }
    llvm::Constant* dbl(double x) { return llvm::ConstantFP::get(Ctx, llvm::APFloat(x)); }
    [[noreturn]] void err(const std::string& m) { throw std::runtime_error("payoff formula: " + m); }

    void skipWs() { while (pos < S.size() && std::isspace((unsigned char)S[pos])) pos++; }
    char peek() { skipWs(); return pos < S.size() ? S[pos] : '\0'; }

    llvm::Value* parseExpr() { return parseCmp(); }

    llvm::Value* parseCmp() {
        llvm::Value* lhs = parseAdd();
        skipWs();
        if (pos < S.size()) {
            const char c = S[pos];
            llvm::CmpInst::Predicate pred;
            int len = 0;
            const bool eq = (pos + 1 < S.size() && S[pos + 1] == '=');
            if (c == '<') { pred = eq ? llvm::CmpInst::FCMP_OLE : llvm::CmpInst::FCMP_OLT; len = eq ? 2 : 1; }
            else if (c == '>') { pred = eq ? llvm::CmpInst::FCMP_OGE : llvm::CmpInst::FCMP_OGT; len = eq ? 2 : 1; }
            else if (c == '=' && eq) { pred = llvm::CmpInst::FCMP_OEQ; len = 2; }
            else if (c == '!' && eq) { pred = llvm::CmpInst::FCMP_ONE; len = 2; }
            if (len) {
                pos += len;
                llvm::Value* rhs = parseAdd();
                // Map the boolean result to 1.0 / 0.0 so it composes arithmetically.
                return B.CreateUIToFP(B.CreateFCmp(pred, lhs, rhs), dty());
            }
        }
        return lhs;
    }

    llvm::Value* parseAdd() {
        llvm::Value* lhs = parseMul();
        while (true) {
            const char c = peek();
            if (c == '+') { pos++; lhs = B.CreateFAdd(lhs, parseMul()); }
            else if (c == '-') { pos++; lhs = B.CreateFSub(lhs, parseMul()); }
            else return lhs;
        }
    }

    llvm::Value* parseMul() {
        llvm::Value* lhs = parseFactor();
        while (true) {
            const char c = peek();
            if (c == '*') { pos++; lhs = B.CreateFMul(lhs, parseFactor()); }
            else if (c == '/') { pos++; lhs = B.CreateFDiv(lhs, parseFactor()); }
            else return lhs;
        }
    }

    llvm::Value* parseFactor() {
        if (peek() == '-') { pos++; return B.CreateFNeg(parseFactor()); }
        return parsePrimary();
    }

    llvm::Value* parsePrimary() {
        const char c = peek();
        if (c == '(') {
            pos++;
            llvm::Value* v = parseExpr();
            if (peek() != ')') err("missing ')'");
            pos++;
            return v;
        }
        if (std::isdigit((unsigned char)c) || c == '.') return parseNumber();
        if (std::isalpha((unsigned char)c) || c == '_') return parseIdent();
        err(std::string("unexpected character '") + c + "'");
    }

    llvm::Value* parseNumber() {
        skipWs();
        const size_t start = pos;
        while (pos < S.size() && (std::isdigit((unsigned char)S[pos]) || S[pos] == '.')) pos++;
        return dbl(std::stod(S.substr(start, pos - start)));
    }

    llvm::Value* parseIdent() {
        skipWs();
        const size_t start = pos;
        while (pos < S.size() && (std::isalnum((unsigned char)S[pos]) || S[pos] == '_')) pos++;
        const std::string name = S.substr(start, pos - start);

        if (peek() == '(') {  // function call
            pos++;
            std::vector<llvm::Value*> a;
            if (peek() != ')') {
                a.push_back(parseExpr());
                while (peek() == ',') { pos++; a.push_back(parseExpr()); }
            }
            if (peek() != ')') err("missing ')' in call to " + name);
            pos++;
            return emitCall(name, a);
        }

        auto it = Vars.find(name);  // variable: load from v[index]
        if (it == Vars.end()) err("unknown variable: " + name);
        llvm::Value* gep = B.CreateConstInBoundsGEP1_64(dty(), Args, it->second);
        return B.CreateLoad(dty(), gep);
    }

    llvm::Value* emitCall(const std::string& name, const std::vector<llvm::Value*>& a) {
        const size_t n = a.size();
        if (name == "max" && n == 2) return B.CreateSelect(B.CreateFCmpOGT(a[0], a[1]), a[0], a[1]);
        if (name == "min" && n == 2) return B.CreateSelect(B.CreateFCmpOLT(a[0], a[1]), a[0], a[1]);
        if (name == "pow" && n == 2) return B.CreateBinaryIntrinsic(llvm::Intrinsic::pow, a[0], a[1]);
        if (name == "exp" && n == 1) return B.CreateUnaryIntrinsic(llvm::Intrinsic::exp, a[0]);
        if (name == "log" && n == 1) return B.CreateUnaryIntrinsic(llvm::Intrinsic::log, a[0]);
        if (name == "sqrt" && n == 1) return B.CreateUnaryIntrinsic(llvm::Intrinsic::sqrt, a[0]);
        if (name == "abs" && n == 1) return B.CreateUnaryIntrinsic(llvm::Intrinsic::fabs, a[0]);
        if (name == "if" && n == 3)
            return B.CreateSelect(B.CreateFCmpONE(a[0], dbl(0.0)), a[1], a[2]);
        err("unknown function or wrong arity: " + name + "/" + std::to_string(n));
    }
};

}  // namespace detail

// JIT engine that turns payoff formulas into callable native functions.
class PayoffJit {
public:
    using Fn = double (*)(const double*);  // payoff(v) where v holds the variables

    PayoffJit() {
        llvm::InitializeNativeTarget();
        llvm::InitializeNativeTargetAsmPrinter();
        auto jit = llvm::orc::LLJITBuilder().create();
        if (!jit) { llvm::consumeError(jit.takeError()); throw std::runtime_error("LLJIT creation failed"); }
        jit_ = std::move(*jit);
    }

    // Compile `formula`, binding the given variable names to indices 0..n-1.
    // The returned pointer reads those variables from its `const double*` argument.
    Fn compile(const std::string& formula, const std::vector<std::string>& varNames) {
        auto ctx = std::make_unique<llvm::LLVMContext>();
        auto mod = std::make_unique<llvm::Module>("payoff_module", *ctx);
        llvm::Type* dty = llvm::Type::getDoubleTy(*ctx);
        llvm::Type* pty = llvm::PointerType::getUnqual(*ctx);
        llvm::FunctionType* fty = llvm::FunctionType::get(dty, {pty}, false);

        const std::string name = "payoff_" + std::to_string(counter_++);
        llvm::Function* fn =
            llvm::Function::Create(fty, llvm::Function::ExternalLinkage, name, mod.get());
        fn->getArg(0)->setName("v");

        llvm::IRBuilder<> builder(llvm::BasicBlock::Create(*ctx, "entry", fn));
        std::map<std::string, unsigned> vars;
        for (unsigned i = 0; i < varNames.size(); ++i) vars[varNames[i]] = i;

        detail::PayoffParser parser(builder, *ctx, fn->getArg(0), vars);
        builder.CreateRet(parser.run(formula));

        { llvm::raw_string_ostream os(last_ir_); last_ir_.clear(); mod->print(os, nullptr); }
        if (llvm::verifyFunction(*fn)) throw std::runtime_error("IR verification failed");

        if (auto e = jit_->addIRModule(llvm::orc::ThreadSafeModule(std::move(mod), std::move(ctx)))) {
            llvm::consumeError(std::move(e));
            throw std::runtime_error("addIRModule failed");
        }
        auto sym = jit_->lookup(name);
        if (!sym) { llvm::consumeError(sym.takeError()); throw std::runtime_error("symbol lookup failed"); }
        return sym->toPtr<Fn>();
    }

    // LLVM IR text produced by the most recent compile() (for inspection/printing).
    const std::string& last_ir() const { return last_ir_; }

private:
    std::unique_ptr<llvm::orc::LLJIT> jit_;
    unsigned counter_ = 0;
    std::string last_ir_;
};

}  // namespace pricer
