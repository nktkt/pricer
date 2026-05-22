// pricer/payoff_ast.hpp — tokenizer, parser, and typed AST for the payoff DSL.
//
// Splits the payoff language into stages: tokenize → parse → a typed AST that a
// backend then consumes. This header has NO LLVM dependency: it can tokenize,
// parse, pretty-print, and *interpret* a formula in pure C++ (used by tests and
// the no-LLVM example). The JIT (payoff_jit.hpp) reuses this same AST and walks
// it to emit LLVM IR, so there is one grammar with two backends.
//
// Grammar (precedence low → high):
//   expr   := cmp
//   cmp    := add (('<'|'>'|'<='|'>='|'=='|'!=') add)?   // yields 1.0 / 0.0
//   add    := mul (('+'|'-') mul)*
//   mul    := factor (('*'|'/') factor)*
//   factor := '-' factor | primary
//   primary:= number | '(' expr ')' | ident | ident '(' args ')'
// Functions: exp log sqrt abs (1), max min pow (2), if(cond,a,b) (3).
#pragma once
#include <cctype>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace pricer::ast {

// ---- tokens ----
enum class Tok { Num, Ident, Op, LParen, RParen, Comma, End };
struct Token {
    Tok kind;
    double num = 0.0;
    std::string text;  // identifier name or operator spelling
};

inline std::vector<Token> tokenize(const std::string& s) {
    std::vector<Token> out;
    size_t i = 0;
    auto err = [](const std::string& m) { throw std::runtime_error("payoff lexer: " + m); };
    while (i < s.size()) {
        const char c = s[i];
        if (std::isspace((unsigned char)c)) { ++i; continue; }
        if (std::isdigit((unsigned char)c) || c == '.') {
            size_t j = i;
            while (j < s.size() && (std::isdigit((unsigned char)s[j]) || s[j] == '.')) ++j;
            out.push_back({Tok::Num, std::stod(s.substr(i, j - i)), {}});
            i = j;
        } else if (std::isalpha((unsigned char)c) || c == '_') {
            size_t j = i;
            while (j < s.size() && (std::isalnum((unsigned char)s[j]) || s[j] == '_')) ++j;
            out.push_back({Tok::Ident, 0.0, s.substr(i, j - i)});
            i = j;
        } else if (c == '(') { out.push_back({Tok::LParen, 0, "("}); ++i; }
        else if (c == ')') { out.push_back({Tok::RParen, 0, ")"}); ++i; }
        else if (c == ',') { out.push_back({Tok::Comma, 0, ","}); ++i; }
        else if (c == '+' || c == '-' || c == '*' || c == '/') {
            out.push_back({Tok::Op, 0, std::string(1, c)});
            ++i;
        } else if (c == '<' || c == '>' || c == '=' || c == '!') {
            const bool eq = (i + 1 < s.size() && s[i + 1] == '=');
            if ((c == '=' || c == '!') && !eq) err(std::string("stray '") + c + "'");
            out.push_back({Tok::Op, 0, eq ? std::string{c, '='} : std::string(1, c)});
            i += eq ? 2 : 1;
        } else {
            err(std::string("unexpected character '") + c + "'");
        }
    }
    out.push_back({Tok::End, 0, {}});
    return out;
}

// ---- AST ----
enum class NodeType { Number, Var, Unary, Binary, Call };
struct Node;
using NodePtr = std::shared_ptr<Node>;
struct Node {
    NodeType type;
    double number = 0.0;        // Number
    std::string name;           // Var name / Call function name
    std::string op;             // Unary ("neg") / Binary ("+","-","*","/","<","<=",...)
    std::vector<NodePtr> kids;  // Unary:1, Binary:2, Call:n
};

// ---- parser (recursive descent over the token stream) ----
class Parser {
public:
    explicit Parser(std::vector<Token> toks) : t_(std::move(toks)) {}
    NodePtr parse() {
        NodePtr n = expr();
        expect(Tok::End);
        return n;
    }

private:
    std::vector<Token> t_;
    size_t p_ = 0;

    const Token& peek() const { return t_[p_]; }
    Token next() { return t_[p_++]; }
    [[noreturn]] void err(const std::string& m) { throw std::runtime_error("payoff parser: " + m); }
    void expect(Tok k) { if (peek().kind != k) err("unexpected token"); }
    bool is_op(const char* s) const { return peek().kind == Tok::Op && peek().text == s; }

    static NodePtr make(NodeType t) { auto n = std::make_shared<Node>(); n->type = t; return n; }

    NodePtr expr() { return cmp(); }
    NodePtr cmp() {
        NodePtr a = add();
        if (peek().kind == Tok::Op) {
            const std::string& o = peek().text;
            if (o == "<" || o == ">" || o == "<=" || o == ">=" || o == "==" || o == "!=") {
                next();
                auto n = make(NodeType::Binary);
                n->op = o; n->kids = {a, add()};
                return n;
            }
        }
        return a;
    }
    NodePtr add() {
        NodePtr a = mul();
        while (is_op("+") || is_op("-")) {
            const std::string o = next().text;
            auto n = make(NodeType::Binary);
            n->op = o; n->kids = {a, mul()};
            a = n;
        }
        return a;
    }
    NodePtr mul() {
        NodePtr a = factor();
        while (is_op("*") || is_op("/")) {
            const std::string o = next().text;
            auto n = make(NodeType::Binary);
            n->op = o; n->kids = {a, factor()};
            a = n;
        }
        return a;
    }
    NodePtr factor() {
        if (is_op("-")) {
            next();
            auto n = make(NodeType::Unary);
            n->op = "neg"; n->kids = {factor()};
            return n;
        }
        return primary();
    }
    NodePtr primary() {
        const Token& tk = peek();
        if (tk.kind == Tok::Num) { auto n = make(NodeType::Number); n->number = next().num; return n; }
        if (tk.kind == Tok::LParen) {
            next();
            NodePtr v = expr();
            if (peek().kind != Tok::RParen) err("missing ')'");
            next();
            return v;
        }
        if (tk.kind == Tok::Ident) {
            const std::string name = next().text;
            if (peek().kind == Tok::LParen) {  // function call
                next();
                auto n = make(NodeType::Call);
                n->name = name;
                if (peek().kind != Tok::RParen) {
                    n->kids.push_back(expr());
                    while (peek().kind == Tok::Comma) { next(); n->kids.push_back(expr()); }
                }
                if (peek().kind != Tok::RParen) err("missing ')' in call to " + name);
                next();
                return n;
            }
            auto n = make(NodeType::Var);
            n->name = name;
            return n;
        }
        err("expected a value");
    }
};

inline NodePtr parse(const std::string& formula) { return Parser(tokenize(formula)).parse(); }

// ---- tree-walking interpreter (pure C++, no LLVM) ----
inline double eval(const Node& n, const std::map<std::string, double>& env) {
    auto e = [&](int k) { return eval(*n.kids[k], env); };
    switch (n.type) {
        case NodeType::Number: return n.number;
        case NodeType::Var: {
            auto it = env.find(n.name);
            if (it == env.end()) throw std::runtime_error("payoff eval: unknown variable " + n.name);
            return it->second;
        }
        case NodeType::Unary: return -e(0);
        case NodeType::Binary: {
            const double a = e(0), b = e(1);
            const std::string& o = n.op;
            if (o == "+") return a + b;
            if (o == "-") return a - b;
            if (o == "*") return a * b;
            if (o == "/") return a / b;
            if (o == "<") return a < b ? 1.0 : 0.0;
            if (o == ">") return a > b ? 1.0 : 0.0;
            if (o == "<=") return a <= b ? 1.0 : 0.0;
            if (o == ">=") return a >= b ? 1.0 : 0.0;
            if (o == "==") return a == b ? 1.0 : 0.0;
            if (o == "!=") return a != b ? 1.0 : 0.0;
            throw std::runtime_error("payoff eval: bad operator " + o);
        }
        case NodeType::Call: {
            const std::string& f = n.name;
            const size_t k = n.kids.size();
            if (f == "max" && k == 2) return std::fmax(e(0), e(1));
            if (f == "min" && k == 2) return std::fmin(e(0), e(1));
            if (f == "pow" && k == 2) return std::pow(e(0), e(1));
            if (f == "exp" && k == 1) return std::exp(e(0));
            if (f == "log" && k == 1) return std::log(e(0));
            if (f == "sqrt" && k == 1) return std::sqrt(e(0));
            if (f == "abs" && k == 1) return std::fabs(e(0));
            if (f == "if" && k == 3) return e(0) != 0.0 ? e(1) : e(2);
            throw std::runtime_error("payoff eval: unknown function/arity " + f);
        }
    }
    return 0.0;  // unreachable
}

inline double eval(const NodePtr& n, const std::map<std::string, double>& env) { return eval(*n, env); }

// ---- pretty-printer (for debugging / inspection) ----
inline std::string to_string(const Node& n) {
    switch (n.type) {
        case NodeType::Number: {
            std::string s = std::to_string(n.number);
            return s;
        }
        case NodeType::Var: return n.name;
        case NodeType::Unary: return "(-" + to_string(*n.kids[0]) + ")";
        case NodeType::Binary: return "(" + to_string(*n.kids[0]) + " " + n.op + " " + to_string(*n.kids[1]) + ")";
        case NodeType::Call: {
            std::string s = n.name + "(";
            for (size_t i = 0; i < n.kids.size(); ++i) { if (i) s += ", "; s += to_string(*n.kids[i]); }
            return s + ")";
        }
    }
    return {};
}

}  // namespace pricer::ast
