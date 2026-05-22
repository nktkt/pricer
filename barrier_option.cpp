// barrier_option.cpp  — テーマ4: 複雑な商品に拡張（バリアオプション）
// 「アップ・アンド・アウト・コール」:
//   通常のコールだが、満期までに一度でも株価がバリア B 以上に達したら無効（ノックアウト）。
// 満期株価だけでなく「途中の経路」に依存する（パス依存）ため、
// 単純な公式では扱いづらく、モンテカルロで経路を刻んでシミュレーションするのが自然。
//
// 比較のため、同じ乱数で「通常のコール」も計算する。
// ノックアウトされる分だけ、バリア付きの方が必ず安く（価値が低く）なる。
#include "bs_common.hpp"
#include <cstdio>
#include <random>

int main() {
    double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    double B = 130.0;                 // バリア水準（これに触れたら無効）
    long n_paths = 2'000'000;         // 経路の本数
    int  n_steps = 250;               // 1年を 250 営業日に分割

    double dt = T / n_steps;
    double drift = (r - 0.5 * sigma * sigma) * dt;
    double vol   = sigma * std::sqrt(dt);
    double disc  = std::exp(-r * T);

    std::mt19937_64 rng(2024);
    std::normal_distribution<double> Z(0.0, 1.0);

    double sum_barrier = 0.0, sum_vanilla = 0.0;
    for (long p = 0; p < n_paths; ++p) {
        double price = S;
        bool knocked_out = false;
        for (int s = 0; s < n_steps; ++s) {        // 1日ずつ経路を進める
            price *= std::exp(drift + vol * Z(rng));
            if (price >= B) knocked_out = true;     // 一度でも触れたら無効
        }
        double payoff = (price > K) ? (price - K) : 0.0;
        sum_vanilla += payoff;                      // 通常コール
        if (!knocked_out) sum_barrier += payoff;    // バリア付きは生き残った経路のみ
    }

    double vanilla_mc = disc * sum_vanilla / n_paths;
    double barrier_mc = disc * sum_barrier / n_paths;
    double vanilla_bs = black_scholes_call(S, K, r, sigma, T);

    std::printf("前提: S=%.0f K=%.0f バリア B=%.0f 経路=%ld ステップ=%d\n\n",
                S, K, B, n_paths, n_steps);
    std::printf("通常コール (解析解BS)        : %.6f\n", vanilla_bs);
    std::printf("通常コール (モンテカルロ)    : %.6f  ← BS とほぼ一致\n", vanilla_mc);
    std::printf("バリア付き (アップ&アウト)   : %.6f  ← ノックアウト分だけ安い\n", barrier_mc);
    std::printf("差（ノックアウトの影響）     : %.6f\n", vanilla_mc - barrier_mc);
    std::printf("\nバリア付きは公式が複雑/存在しないことも多く、経路を刻むMCが活きる典型例。\n");
    return 0;
}
