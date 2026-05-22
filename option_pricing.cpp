// option_pricing.cpp
// ヨーロピアン・コールオプションの理論価格を 2 通りの方法で求めて比較する。
//   (1) ブラック・ショールズ式   … 数式を直接計算する「解析解」
//   (2) モンテカルロ法           … 乱数で満期株価を大量に生成して平均をとる近似
// 同じ前提なら両者はほぼ一致するはず、というのを確認するのが狙い。

#include <cmath>
#include <cstdio>
#include <random>
#include <chrono>

// 標準正規分布の累積分布関数 N(x)。
// erfc（相補誤差関数）を使って計算できる。
static double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// ---- (1) ブラック・ショールズ式によるコール価格 ----
// S     : 現在の株価
// K     : 権利行使価格（ストライク）
// r     : 無リスク金利（年率）
// sigma : ボラティリティ（年率の標準偏差）
// T     : 満期までの年数
static double black_scholes_call(double S, double K, double r, double sigma, double T) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}

// ---- (2) モンテカルロ法によるコール価格 ----
// 幾何ブラウン運動のもとでの満期株価
//   S_T = S * exp((r - 0.5*sigma^2)*T + sigma*sqrt(T)*Z),  Z ~ 標準正規分布
// を n_paths 回サンプリングし、ペイオフ max(S_T - K, 0) の平均を割り引く。
static double monte_carlo_call(double S, double K, double r, double sigma, double T,
                               long n_paths) {
    std::mt19937_64 rng(12345);                 // 再現性のため固定シード
    std::normal_distribution<double> Z(0.0, 1.0);

    double drift = (r - 0.5 * sigma * sigma) * T;
    double vol   = sigma * std::sqrt(T);

    double payoff_sum = 0.0;
    for (long i = 0; i < n_paths; ++i) {
        double ST = S * std::exp(drift + vol * Z(rng));
        double payoff = ST - K;
        if (payoff > 0.0) payoff_sum += payoff;  // max(ST - K, 0)
    }
    double mean_payoff = payoff_sum / static_cast<double>(n_paths);
    return std::exp(-r * T) * mean_payoff;       // 現在価値に割り引く
}

int main() {
    // 前提条件（パラメータ）
    double S = 100.0;     // 現在株価
    double K = 100.0;     // 権利行使価格
    double r = 0.05;      // 金利 5%
    double sigma = 0.20;  // ボラティリティ 20%
    double T = 1.0;       // 満期 1 年
    long n_paths = 10'000'000;  // モンテカルロのシナリオ数（1000万）

    // (1) 解析解
    double bs = black_scholes_call(S, K, r, sigma, T);

    // (2) モンテカルロ（実行時間も測る）
    auto t0 = std::chrono::high_resolution_clock::now();
    double mc = monte_carlo_call(S, K, r, sigma, T, n_paths);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::printf("前提: S=%.1f K=%.1f r=%.2f sigma=%.2f T=%.1f\n", S, K, r, sigma, T);
    std::printf("シナリオ数: %ld\n\n", n_paths);
    std::printf("(1) ブラック・ショールズ : %.6f\n", bs);
    std::printf("(2) モンテカルロ法       : %.6f\n", mc);
    std::printf("    誤差                 : %.6f (%.4f%%)\n", mc - bs, (mc - bs) / bs * 100.0);
    std::printf("    モンテカルロ計算時間 : %.1f ms\n", ms);
    return 0;
}
