// bs_common.hpp
// オプション価格計算の共通部品。各サンプルから #include して使う。
#pragma once
#include <cmath>

// 標準正規分布の確率密度関数 n(x) = exp(-x^2/2)/sqrt(2*pi)
inline double norm_pdf(double x) {
    static const double inv_sqrt_2pi = 0.3989422804014327;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

// 標準正規分布の累積分布関数 N(x)
inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// ブラック・ショールズ式によるヨーロピアン・コール価格（解析解）
//   S: 株価, K: 行使価格, r: 金利, sigma: ボラティリティ, T: 満期(年)
inline double black_scholes_call(double S, double K, double r, double sigma, double T) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}
