// payoff_interpret.cpp — the payoff DSL without LLVM.
//
// Tokenize and parse a formula into a typed AST, pretty-print it, then evaluate
// it with the tree-walking interpreter. This is the same grammar the JIT
// compiles, but here it runs in pure C++ — handy for testing, validation, or
// environments without LLVM. (For hot loops, JIT the AST instead; see jit_payoff.)
#include "pricer/payoff_ast.hpp"
#include <cstdio>
#include <map>
#include <string>

using namespace pricer;

int main(int argc, char** argv) {
    const std::string formula = (argc > 1) ? argv[1] : "max(ST - K, 0)";

    ast::NodePtr tree;
    try {
        tree = ast::parse(formula);  // tokenizer + parser -> typed AST
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }

    std::printf("formula : \"%s\"\n", formula.c_str());
    std::printf("AST     : %s\n\n", ast::to_string(*tree).c_str());

    // Evaluate the same instrument at a few market states.
    std::printf("%8s | %8s | %12s\n", "ST", "K", "payoff");
    std::printf("---------|----------|-------------\n");
    try {
        for (double ST : {80.0, 100.0, 120.0, 140.0}) {
            const std::map<std::string, double> env = {{"ST", ST}, {"K", 100.0}};
            std::printf("%8.2f | %8.2f | %12.6f\n", ST, 100.0, ast::eval(*tree, env));
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s (this demo only binds ST and K)\n", e.what());
        return 1;
    }
    std::printf("\nParsed once into an AST, then interpreted — no LLVM required.\n");
    return 0;
}
