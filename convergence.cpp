// convergence.cpp  — テーマ1: 精度と速度の関係を見る
// モンテカルロのシナリオ数 N を 1万 → 1億 と増やし、
//   ・解析解（ブラック・ショールズ）との誤差
//   ・計算時間
// がどう変わるかを観察する。
// 理論上、モンテカルロの誤差は 1/sqrt(N) で減る（N を 100 倍にして誤差は 1/10）。
#include "bs_common.hpp"
#include <cstdio>
#include <random>
#include <chrono>

static double monte_carlo_call(double S, double K, double r, double sigma, double T,
                               long n_paths) {
    std::mt19937_64 rng(12345);
    std::normal_distribution<double> Z(0.0, 1.0);
    double drift = (r - 0.5 * sigma * sigma) * T;
    double vol   = sigma * std::sqrt(T);
    double sum = 0.0;
    for (long i = 0; i < n_paths; ++i) {
        double ST = S * std::exp(drift + vol * Z(rng));
        if (ST > K) sum += ST - K;
    }
    return std::exp(-r * T) * (sum / static_cast<double>(n_paths));
}

int main() {
    double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    double exact = black_scholes_call(S, K, r, sigma, T);
    std::printf("解析解(ブラック・ショールズ) = %.6f\n\n", exact);
    std::printf("%12s | %12s | %10s | %10s\n", "シナリオ数", "MC価格", "誤差", "時間(ms)");
    std::printf("-------------|--------------|------------|----------\n");

    for (long n = 10'000; n <= 100'000'000; n *= 10) {
        auto t0 = std::chrono::high_resolution_clock::now();
        double mc = monte_carlo_call(S, K, r, sigma, T, n);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::printf("%12ld | %12.6f | %10.6f | %10.1f\n", n, mc, mc - exact, ms);
    }
    std::printf("\nN を 100 倍にするたび、誤差はおよそ 1/10 に縮む（1/sqrt(N) の法則）。\n");
    return 0;
}
