// pricer/payoff_jit.hpp — JIT-compile a payoff formula to native code with LLVM.
//
// The formula is tokenized and parsed into a typed AST (pricer/payoff_ast.hpp);
// this header then walks that AST to emit LLVM IR and JIT-compiles it. One
// grammar, two backends: the AST interpreter (no LLVM) and this code generator.
//
// `compile(formula, varNames)` yields `double payoff(const double* v)` reading
// each named variable from v[index]. `compile_batch(formula, varNames, W)`
// yields `void payoff_v(const double* v, double* out)` over `<W x double>` IR
// (structure-of-arrays input). Compiled kernels are cached by (width, vars,
// formula). Requires LLVM; only targets that link LLVM should include it.
#pragma once
#include <cstdint>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "pricer/payoff_ast.hpp"

#include "llvm/ExecutionEngine/Orc/LLJIT.h"
#include "llvm/ExecutionEngine/Orc/ThreadSafeModule.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Intrinsics.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Verifier.h"
#include "llvm/Support/Alignment.h"
#include "llvm/Support/Error.h"
#include "llvm/Support/TargetSelect.h"
#include "llvm/Support/raw_ostream.h"

namespace pricer {

namespace detail {

// Walk a payoff AST and emit LLVM IR. Works for scalar (Width==1, ValTy=double)
// and vector (ValTy = <Width x double>) codegen — the IRBuilder ops are type-
// polymorphic, so only the leaves (constants, variable loads) differ.
class Codegen {
public:
    Codegen(llvm::IRBuilder<>& b, llvm::Value* args, const std::map<std::string, unsigned>& vars,
            llvm::Type* valTy, llvm::Type* eltTy, unsigned width)
        : B(b), Args(args), Vars(vars), ValTy(valTy), EltTy(eltTy), Width(width) {}

    llvm::Value* gen(const ast::Node& n) {
        switch (n.type) {
            case ast::NodeType::Number: return llvm::ConstantFP::get(ValTy, n.number);  // splats if vector
            case ast::NodeType::Var: return loadVar(n.name);
            case ast::NodeType::Unary: return B.CreateFNeg(gen(*n.kids[0]));
            case ast::NodeType::Binary: return genBinary(n);
            case ast::NodeType::Call: return genCall(n);
        }
        err("bad AST node");
    }

private:
    [[noreturn]] void err(const std::string& m) { throw std::runtime_error("payoff codegen: " + m); }

    llvm::Value* loadVar(const std::string& name) {
        auto it = Vars.find(name);
        if (it == Vars.end()) err("unknown variable: " + name);
        llvm::Value* ptr =
            B.CreateConstInBoundsGEP1_64(EltTy, Args, static_cast<uint64_t>(it->second) * Width);
        return B.CreateAlignedLoad(ValTy, ptr, llvm::Align(8));
    }

    llvm::Value* genBinary(const ast::Node& n) {
        llvm::Value* a = gen(*n.kids[0]);
        llvm::Value* b = gen(*n.kids[1]);
        const std::string& o = n.op;
        if (o == "+") return B.CreateFAdd(a, b);
        if (o == "-") return B.CreateFSub(a, b);
        if (o == "*") return B.CreateFMul(a, b);
        if (o == "/") return B.CreateFDiv(a, b);
        llvm::CmpInst::Predicate p;
        if (o == "<") p = llvm::CmpInst::FCMP_OLT;
        else if (o == ">") p = llvm::CmpInst::FCMP_OGT;
        else if (o == "<=") p = llvm::CmpInst::FCMP_OLE;
        else if (o == ">=") p = llvm::CmpInst::FCMP_OGE;
        else if (o == "==") p = llvm::CmpInst::FCMP_OEQ;
        else if (o == "!=") p = llvm::CmpInst::FCMP_ONE;
        else err("bad operator " + o);
        // Map the boolean to 1.0 / 0.0 so it composes arithmetically.
        return B.CreateUIToFP(B.CreateFCmp(p, a, b), ValTy);
    }

    llvm::Value* genCall(const ast::Node& n) {
        std::vector<llvm::Value*> a;
        a.reserve(n.kids.size());
        for (const auto& k : n.kids) a.push_back(gen(*k));
        const std::string& f = n.name;
        const size_t k = a.size();
        if (f == "max" && k == 2) return B.CreateSelect(B.CreateFCmpOGT(a[0], a[1]), a[0], a[1]);
        if (f == "min" && k == 2) return B.CreateSelect(B.CreateFCmpOLT(a[0], a[1]), a[0], a[1]);
        if (f == "pow" && k == 2) return B.CreateBinaryIntrinsic(llvm::Intrinsic::pow, a[0], a[1]);
        if (f == "exp" && k == 1) return B.CreateUnaryIntrinsic(llvm::Intrinsic::exp, a[0]);
        if (f == "log" && k == 1) return B.CreateUnaryIntrinsic(llvm::Intrinsic::log, a[0]);
        if (f == "sqrt" && k == 1) return B.CreateUnaryIntrinsic(llvm::Intrinsic::sqrt, a[0]);
        if (f == "abs" && k == 1) return B.CreateUnaryIntrinsic(llvm::Intrinsic::fabs, a[0]);
        if (f == "if" && k == 3)
            return B.CreateSelect(B.CreateFCmpONE(a[0], llvm::ConstantFP::get(ValTy, 0.0)), a[1], a[2]);
        err("unknown function or wrong arity: " + f + "/" + std::to_string(k));
    }

    llvm::IRBuilder<>& B;
    llvm::Value* Args;
    const std::map<std::string, unsigned>& Vars;
    llvm::Type* ValTy;   // double, or <Width x double>
    llvm::Type* EltTy;   // always double (the storage element)
    unsigned Width;
};

}  // namespace detail

// JIT engine that turns payoff formulas into callable native functions.
class PayoffJit {
public:
    using Fn = double (*)(const double*);              // scalar:  payoff(v)
    using BatchFn = void (*)(const double*, double*);  // vector:  payoff_v(v, out)

    PayoffJit() {
        llvm::InitializeNativeTarget();
        llvm::InitializeNativeTargetAsmPrinter();
        auto jit = llvm::orc::LLJITBuilder().create();
        if (!jit) { llvm::consumeError(jit.takeError()); throw std::runtime_error("LLJIT creation failed"); }
        jit_ = std::move(*jit);
    }

    // Compile `formula` (scalar). Variable names map to indices 0..n-1, read from
    // the `const double*` argument. Cached by formula + variable list.
    Fn compile(const std::string& formula, const std::vector<std::string>& varNames) {
        return reinterpret_cast<Fn>(static_cast<uintptr_t>(build(formula, varNames, /*width=*/1)));
    }

    // Compile a vectorized kernel processing `width` paths per call (SoA layout).
    BatchFn compile_batch(const std::string& formula, const std::vector<std::string>& varNames,
                          unsigned width) {
        if (width < 2) throw std::runtime_error("compile_batch requires width >= 2");
        return reinterpret_cast<BatchFn>(static_cast<uintptr_t>(build(formula, varNames, width)));
    }

    const std::string& last_ir() const { return last_ir_; }  // IR of the most recent build()
    unsigned compiles() const { return compiles_; }          // actual JIT compiles (cache misses)

private:
    // Build (or fetch from cache) a kernel; returns its native address.
    std::uint64_t build(const std::string& formula, const std::vector<std::string>& varNames,
                        unsigned width) {
        std::string key = std::to_string(width);
        for (const auto& v : varNames) key += "|" + v;
        key += "#" + formula;
        if (auto it = cache_.find(key); it != cache_.end()) return it->second;

        const ast::NodePtr root = ast::parse(formula);  // tokenize + parse -> typed AST

        auto ctx = std::make_unique<llvm::LLVMContext>();
        auto mod = std::make_unique<llvm::Module>("payoff_module", *ctx);
        llvm::Type* eltTy = llvm::Type::getDoubleTy(*ctx);
        llvm::Type* valTy = (width <= 1)
            ? eltTy
            : static_cast<llvm::Type*>(llvm::FixedVectorType::get(eltTy, width));
        llvm::Type* pty = llvm::PointerType::getUnqual(*ctx);

        const bool scalar = (width <= 1);
        llvm::FunctionType* fty = scalar
            ? llvm::FunctionType::get(eltTy, {pty}, false)
            : llvm::FunctionType::get(llvm::Type::getVoidTy(*ctx), {pty, pty}, false);

        const std::string name = (scalar ? "payoff_" : "payoff_v_") + std::to_string(counter_++);
        llvm::Function* fn =
            llvm::Function::Create(fty, llvm::Function::ExternalLinkage, name, mod.get());
        fn->getArg(0)->setName("v");
        if (!scalar) fn->getArg(1)->setName("out");

        llvm::IRBuilder<> builder(llvm::BasicBlock::Create(*ctx, "entry", fn));
        std::map<std::string, unsigned> vars;
        for (unsigned i = 0; i < varNames.size(); ++i) vars[varNames[i]] = i;

        detail::Codegen cg(builder, fn->getArg(0), vars, valTy, eltTy, width);
        llvm::Value* result = cg.gen(*root);
        if (scalar) {
            builder.CreateRet(result);
        } else {
            builder.CreateAlignedStore(result, fn->getArg(1), llvm::Align(8));
            builder.CreateRetVoid();
        }

        last_ir_.clear();
        { llvm::raw_string_ostream os(last_ir_); mod->print(os, nullptr); }
        if (llvm::verifyFunction(*fn)) throw std::runtime_error("IR verification failed");

        if (auto e = jit_->addIRModule(llvm::orc::ThreadSafeModule(std::move(mod), std::move(ctx)))) {
            llvm::consumeError(std::move(e));
            throw std::runtime_error("addIRModule failed");
        }
        auto sym = jit_->lookup(name);
        if (!sym) { llvm::consumeError(sym.takeError()); throw std::runtime_error("symbol lookup failed"); }

        ++compiles_;
        const std::uint64_t addr = sym->getValue();
        cache_.emplace(std::move(key), addr);
        return addr;
    }

    std::unique_ptr<llvm::orc::LLJIT> jit_;
    std::unordered_map<std::string, std::uint64_t> cache_;
    unsigned counter_ = 0;
    unsigned compiles_ = 0;
    std::string last_ir_;
};

}  // namespace pricer
