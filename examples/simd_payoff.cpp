// simd_payoff.cpp — vectorized payoff codegen (Phase 2 SIMD).
//
// PayoffJit can emit two kernels from the same formula:
//   scalar:  double payoff(const double* v)            — one path per call
//   batch:   void   payoff_v(const double* v, double*) — `W` paths per call,
//            using <W x double> IR (structure-of-arrays input).
//
// This benchmark isolates payoff *evaluation* (paths are pre-generated) and
// compares scalar vs. batch throughput. The batch kernel wins by amortizing the
// per-call overhead over W lanes and by issuing SIMD instructions for the
// arithmetic — both are properties of the generated IR, not hand-written code.
#include "pricer/payoff_jit.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace pricer;

int main() {
    const unsigned W = 8;            // lanes per batch call
    const long M = 4'000'000;        // number of paths (multiple of W)
    const int  R = 10;               // repeats, for stable timing
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const char* formula = "max(ST - K, 0) + max(K - ST, 0)";  // straddle (pure arithmetic)

    // Pre-generate terminal prices once; both kernels evaluate the same data.
    std::mt19937_64 rng(123);
    std::normal_distribution<double> Z(0.0, 1.0);
    const double drift = (r - 0.5 * sigma * sigma) * T, vol = sigma * std::sqrt(T);
    std::vector<double> ST(M);
    for (long i = 0; i < M; ++i) ST[i] = S * std::exp(drift + vol * Z(rng));

    // Scalar input layout: [ST, K] per evaluation.
    std::vector<double> sbuf(2 * M);
    for (long i = 0; i < M; ++i) { sbuf[2 * i] = ST[i]; sbuf[2 * i + 1] = K; }

    // Batch input layout (structure-of-arrays per chunk): [ST x W][K x W].
    const long chunks = M / W;
    std::vector<double> vbuf(static_cast<size_t>(chunks) * 2 * W);
    for (long c = 0; c < chunks; ++c) {
        for (unsigned l = 0; l < W; ++l) vbuf[c * 2 * W + l] = ST[c * W + l];
        for (unsigned l = 0; l < W; ++l) vbuf[c * 2 * W + W + l] = K;
    }

    PayoffJit jit;
    auto fs = jit.compile(formula, {"ST", "K"});
    auto fv = jit.compile_batch(formula, {"ST", "K"}, W);

    std::printf("formula : \"%s\"   (W=%u lanes)\n\n--- generated vector IR ---\n%s---------------------------\n\n",
                formula, W, jit.last_ir().c_str());

    // Scalar timing.
    double ssum = 0.0;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int rep = 0; rep < R; ++rep)
        for (long i = 0; i < M; ++i) ssum += fs(&sbuf[2 * i]);
    auto t1 = std::chrono::high_resolution_clock::now();
    const double sms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Batch timing.
    double vsum = 0.0;
    double out[64];  // >= W
    t0 = std::chrono::high_resolution_clock::now();
    for (int rep = 0; rep < R; ++rep)
        for (long c = 0; c < chunks; ++c) {
            fv(&vbuf[c * 2 * W], out);
            for (unsigned l = 0; l < W; ++l) vsum += out[l];
        }
    t1 = std::chrono::high_resolution_clock::now();
    const double vms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    const double evals = static_cast<double>(M) * R;
    std::printf("%-8s | %10s | %12s | %s\n", "kernel", "time(ms)", "Mevals/s", "checksum");
    std::printf("---------|------------|--------------|------------------\n");
    std::printf("%-8s | %10.1f | %12.1f | %.6f\n", "scalar", sms, evals / sms / 1e3, ssum);
    std::printf("%-8s | %10.1f | %12.1f | %.6f\n", "batch", vms, evals / vms / 1e3, vsum);
    std::printf("\nspeedup (scalar/batch time) = %.2fx\n", sms / vms);
    std::printf("checksums match: %s\n", (std::abs(ssum - vsum) < 1e-3) ? "yes" : "NO");
    return 0;
}
