// jit_payoff.cpp  — テーマ5: LLVM JIT につなげる（最初の話題への回帰）
//
// やること:
//   1. コマンドラインで受け取ったペイオフ「数式の文字列」をパースし、
//   2. LLVM IR の関数  double payoff(double ST, double K)  を生成、
//   3. LLVM の JIT(ORC) で実行時にネイティブコードへコンパイル、
//   4. 関数ポインタを取り出し、モンテカルロのループから呼んで価格を計算する。
//
// これが「クオンツが入力した数式を、実行時に高速コードへコンパイルして大量計算する」
// という、第1回で話した LLVM の代表的な使い道のミニ実装。
//
// 使える文法:  数値, 変数 ST / K, 演算子 + - * / , 単項マイナス, 括弧, max(a,b), min(a,b)
// 例:  "max(ST - K, 0)"            … 通常コール
//      "max(K - ST, 0)"            … 通常プット
//      "max(ST-K,0)+max(K-ST,0)"   … ストラドル
//
// ビルド: 末尾のコメント参照（llvm-config を使う）

#include "bs_common.hpp"
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

// ---- 数式 → LLVM IR への再帰下降パーサ ----
// 文字列を読みながら、その場で IRBuilder を使って IR の値(Value*)を組み立てる。
class PayoffCompiler {
public:
    PayoffCompiler(IRBuilder<>& b, LLVMContext& c, std::map<std::string, Value*>& vars)
        : B(b), Ctx(c), Vars(vars) {}

    // 入力式をコンパイルし、結果(double)の Value* を返す
    Value* compile(const std::string& src) {
        S = src; pos = 0;
        Value* v = parseExpr();
        skipWs();
        if (pos != S.size()) throw std::runtime_error("余分な文字: " + S.substr(pos));
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

    // expr := term (('+'|'-') term)*
    Value* parseExpr() {
        Value* lhs = parseTerm();
        while (true) {
            char c = peek();
            if (c == '+') { pos++; lhs = B.CreateFAdd(lhs, parseTerm()); }
            else if (c == '-') { pos++; lhs = B.CreateFSub(lhs, parseTerm()); }
            else return lhs;
        }
    }
    // term := factor (('*'|'/') factor)*
    Value* parseTerm() {
        Value* lhs = parseFactor();
        while (true) {
            char c = peek();
            if (c == '*') { pos++; lhs = B.CreateFMul(lhs, parseFactor()); }
            else if (c == '/') { pos++; lhs = B.CreateFDiv(lhs, parseFactor()); }
            else return lhs;
        }
    }
    // factor := '-' factor | primary
    Value* parseFactor() {
        if (peek() == '-') { pos++; return B.CreateFNeg(parseFactor()); }
        return parsePrimary();
    }
    // primary := number | '(' expr ')' | ident | ident '(' args ')'
    Value* parsePrimary() {
        char c = peek();
        if (c == '(') {
            pos++; Value* v = parseExpr();
            if (peek() != ')') throw std::runtime_error("')' がない");
            pos++; return v;
        }
        if (std::isdigit((unsigned char)c) || c == '.') return parseNumber();
        if (std::isalpha((unsigned char)c)) return parseIdent();
        throw std::runtime_error(std::string("予期しない文字: ") + c);
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

        if (peek() == '(') {                 // 関数呼び出し max(a,b) / min(a,b)
            pos++;
            Value* a = parseExpr();
            if (peek() != ',') throw std::runtime_error(name + " は引数2つ");
            pos++;
            Value* b = parseExpr();
            if (peek() != ')') throw std::runtime_error("')' がない");
            pos++;
            if (name == "max") return B.CreateSelect(B.CreateFCmpOGT(a, b), a, b);
            if (name == "min") return B.CreateSelect(B.CreateFCmpOLT(a, b), a, b);
            throw std::runtime_error("未知の関数: " + name);
        }
        auto it = Vars.find(name);            // 変数 ST / K
        if (it == Vars.end()) throw std::runtime_error("未知の変数: " + name);
        return it->second;
    }
};

// JIT コンパイル済み payoff を使ったモンテカルロ
static double mc_with_jit(double (*payoff)(double, double),
                          double S, double K, double r, double sigma, double T, long n) {
    std::mt19937_64 rng(777);
    std::normal_distribution<double> Z(0.0, 1.0);
    double drift = (r - 0.5 * sigma * sigma) * T, vol = sigma * std::sqrt(T);
    double sum = 0.0;
    for (long i = 0; i < n; ++i) {
        double ST = S * std::exp(drift + vol * Z(rng));
        sum += payoff(ST, K);             // ← 実行時にコンパイルされたネイティブ関数を呼ぶ
    }
    return std::exp(-r * T) * (sum / n);
}

int main(int argc, char** argv) {
    std::string expr = (argc > 1) ? argv[1] : "max(ST - K, 0)";

    // --- LLVM 初期化 ---
    InitializeNativeTarget();
    InitializeNativeTargetAsmPrinter();
    ExitOnError ExitOnErr;

    auto jit = ExitOnErr(orc::LLJITBuilder().create());

    // --- IR 構築: define double @payoff(double %ST, double %K) ---
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
        std::fprintf(stderr, "数式エラー: %s\n", e.what());
        return 1;
    }
    if (verifyFunction(*fn, &errs())) { std::fprintf(stderr, "IR 検証失敗\n"); return 1; }

    std::printf("入力した数式 : \"%s\"\n", expr.c_str());
    std::printf("\n--- 生成された LLVM IR ---\n");
    mod->print(outs(), nullptr);          // 実際に生成された IR を表示
    std::printf("--------------------------\n\n");

    // --- JIT: モジュールを登録し、ネイティブコードにコンパイルして関数ポインタを得る ---
    ExitOnErr(jit->addIRModule(orc::ThreadSafeModule(std::move(mod), std::move(ctx))));
    auto sym = ExitOnErr(jit->lookup("payoff"));
    auto payoff = sym.toPtr<double (*)(double, double)>();

    // 動作確認: payoff(120, 100) を実行
    std::printf("実行時コンパイル完了。payoff(ST=120, K=100) = %.4f\n", payoff(120.0, 100.0));

    // --- JIT 関数を使ってモンテカルロ ---
    double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    long n = 10'000'000;
    auto t0 = std::chrono::high_resolution_clock::now();
    double price = mc_with_jit(payoff, S, K, r, sigma, T, n);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::printf("\nモンテカルロ価格(%ld本) = %.6f  (%.1f ms)\n", n, price, ms);
    if (expr == "max(ST - K, 0)")
        std::printf("解析解(BSコール)        = %.6f  ← 一致すれば成功\n",
                    black_scholes_call(S, K, r, sigma, T));
    return 0;
}

/*  ビルド例:
    clang++ -std=c++17 -O2 jit_payoff.cpp \
      $(/opt/homebrew/opt/llvm/bin/llvm-config --cxxflags --ldflags --libs core orcjit native --system-libs) \
      -o jit_payoff
    実行:   ./jit_payoff                       # 既定: コールのペイオフ
            ./jit_payoff "max(K - ST, 0)"      # プット
*/
