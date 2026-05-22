// greeks.cpp  — テーマ3: グリークス（リスク指標）を追加
// 「前提が少し動いたとき、オプション価格はどれだけ動くか」を表す感応度。
//   Delta : 株価 S が 1 動くと価格はいくら動くか
//   Gamma : Delta 自体が S でどれだけ動くか（曲がり具合）
//   Vega  : ボラティリティ sigma が 1%(0.01) 動くと価格はいくら動くか
//   Theta : 時間が 1 日 経つと価格はいくら減るか（時間価値の目減り）
//   Rho   : 金利 r が 1%(0.01) 動くと価格はいくら動くか
// 解析解（公式）と、数値微分（前提を少しずらして差をとる）を並べて検算する。
#include "bs_common.hpp"
#include <cstdio>

int main() {
    double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;

    double sqrtT = std::sqrt(T);
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    double d2 = d1 - sigma * sqrtT;

    // --- 解析解（公式）---
    double price = black_scholes_call(S, K, r, sigma, T);
    double delta = norm_cdf(d1);
    double gamma = norm_pdf(d1) / (S * sigma * sqrtT);
    double vega  = S * norm_pdf(d1) * sqrtT;                       // sigma 1.0 あたり
    double theta = -(S * norm_pdf(d1) * sigma) / (2 * sqrtT)
                   - r * K * std::exp(-r * T) * norm_cdf(d2);      // 1年あたり
    double rho   = K * T * std::exp(-r * T) * norm_cdf(d2);        // r 1.0 あたり

    // --- 数値微分（中心差分）で検算: f'(x) ≈ (f(x+h) - f(x-h)) / 2h ---
    auto C = [&](double S_, double K_, double r_, double sig_, double T_) {
        return black_scholes_call(S_, K_, r_, sig_, T_);
    };
    double hS = 0.01, hSig = 1e-4, hT = 1e-4, hr = 1e-4;
    double fd_delta = (C(S + hS, K, r, sigma, T) - C(S - hS, K, r, sigma, T)) / (2 * hS);
    double fd_gamma = (C(S + hS, K, r, sigma, T) - 2 * price + C(S - hS, K, r, sigma, T)) / (hS * hS);
    double fd_vega  = (C(S, K, r, sigma + hSig, T) - C(S, K, r, sigma - hSig, T)) / (2 * hSig);
    double fd_theta = -(C(S, K, r, sigma, T + hT) - C(S, K, r, sigma, T - hT)) / (2 * hT);
    double fd_rho   = (C(S, K, r + hr, sigma, T) - C(S, K, r - hr, sigma, T)) / (2 * hr);

    std::printf("コール価格 = %.6f\n\n", price);
    std::printf("%-7s | %12s | %12s | %s\n", "指標", "公式", "数値微分", "意味");
    std::printf("--------|--------------|--------------|------------------------------\n");
    std::printf("%-7s | %12.6f | %12.6f | %s\n", "Delta", delta, fd_delta, "株価+1 で価格 +Delta");
    std::printf("%-7s | %12.6f | %12.6f | %s\n", "Gamma", gamma, fd_gamma, "Delta の変化率");
    std::printf("%-7s | %12.6f | %12.6f | %s\n", "Vega",  vega,  fd_vega,  "vol+1.0 で価格 +Vega (1%なら/100)");
    std::printf("%-7s | %12.6f | %12.6f | %s\n", "Theta", theta, fd_theta, "1年あたりの時間価値の減少");
    std::printf("%-7s | %12.6f | %12.6f | %s\n", "Rho",   rho,   fd_rho,   "金利+1.0 で価格 +Rho");
    std::printf("\n公式と数値微分がほぼ一致 → 実装が正しいことの確認になる。\n");
    std::printf("実務での目安: Theta/365 = 1日あたり %.4f, Vega/100 = vol1%%で %.4f\n",
                theta / 365.0, vega / 100.0);
    return 0;
}
