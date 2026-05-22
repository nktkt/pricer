// parallel_mc.cpp  — テーマ2: 並列化で高速化
// 同じモンテカルロ計算を「1スレッド」と「全CPUコア」で実行し、速度差を見る。
// 各スレッドは別々の乱数シードを持ち、担当ぶんのシナリオを計算して、
// 最後にペイオフの合計を足し合わせる（並列リダクション）。
#include "bs_common.hpp"
#include <cstdio>
#include <random>
#include <chrono>
#include <thread>
#include <vector>

// シードと担当シナリオ数を受け取り、ペイオフ合計を返す（割引前）
static double payoff_sum(unsigned seed, long n, double S, double K,
                         double drift, double vol) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);
    double sum = 0.0;
    for (long i = 0; i < n; ++i) {
        double ST = S * std::exp(drift + vol * Z(rng));
        if (ST > K) sum += ST - K;
    }
    return sum;
}

int main() {
    double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    long total = 200'000'000;            // 合計2億シナリオ
    double drift = (r - 0.5 * sigma * sigma) * T;
    double vol   = sigma * std::sqrt(T);
    double disc  = std::exp(-r * T);
    double exact = black_scholes_call(S, K, r, sigma, T);

    unsigned ncores = std::thread::hardware_concurrency();
    std::printf("解析解 = %.6f / 合計シナリオ = %ld / 論理コア数 = %u\n\n",
                exact, total, ncores);

    // --- 1スレッド ---
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        double price = disc * payoff_sum(1, total, S, K, drift, vol) / total;
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::printf("1スレッド   : 価格=%.6f  時間=%8.1f ms\n", price, ms);
    }

    // --- 全コア並列 ---
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::thread> threads;
        std::vector<double> partial(ncores, 0.0);
        long per = total / ncores;
        for (unsigned t = 0; t < ncores; ++t) {
            long n = (t == ncores - 1) ? (total - per * (ncores - 1)) : per;
            threads.emplace_back([&, t, n] {
                partial[t] = payoff_sum(100 + t, n, S, K, drift, vol);
            });
        }
        for (auto& th : threads) th.join();
        double sum = 0.0;
        for (double p : partial) sum += p;
        double price = disc * sum / total;
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::printf("%u スレッド  : 価格=%.6f  時間=%8.1f ms\n", ncores, price, ms);
    }
    std::printf("\nコア数にほぼ比例して速くなる（メモリ帯域などで頭打ちはある）。\n");
    return 0;
}
