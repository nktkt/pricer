// pricer/simd.hpp — a small, portable SIMD layer for path generation.
//
// Built on the GCC/Clang vector extensions (`__attribute__((vector_size))`),
// which the compiler lowers to AVX2/AVX-512 on x86-64 and to NEON on ARM — one
// source, many backends, no intrinsics and no third-party math library. We need
// vectorized transcendentals to make Monte Carlo path generation actually use
// the SIMD units (a plain loop over `std::exp` does not vectorize without a
// vector libm), so this header provides `v_exp`, `v_log`, `v_sqrt` and a
// branchless vectorized inverse-normal CDF.
//
// The width is fixed at 4 doubles (256-bit). The transcendentals are accurate to
// ~1e-11 relative — well inside Monte Carlo's 1/sqrt(N) sampling error — and the
// integer RNG mixing is bit-identical to the scalar path, so the SIMD and scalar
// engines agree to floating-point tolerance.
#pragma once
#include <cstdint>

#if defined(__clang__) || defined(__GNUC__)
#define PRICER_HAVE_SIMD 1
#else
#define PRICER_HAVE_SIMD 0
#endif

namespace pricer::simd {

#if PRICER_HAVE_SIMD

// 4-wide (256-bit) lanes. `vd` doubles, `vi`/`vu` signed/unsigned 64-bit ints.
inline constexpr int kWidth = 4;
typedef double   vd __attribute__((vector_size(sizeof(double) * 4)));
typedef int64_t  vi __attribute__((vector_size(sizeof(int64_t) * 4)));
typedef uint64_t vu __attribute__((vector_size(sizeof(uint64_t) * 4)));

// --- reinterpret (bit-cast) and broadcast helpers ---------------------------
inline vu as_u(vd x) { union { vd d; vu u; } c; c.d = x; return c.u; }
inline vd as_d(vu x) { union { vd d; vu u; } c; c.u = x; return c.d; }
inline vu as_u(vi x) { union { vi i; vu u; } c; c.i = x; return c.u; }
inline vi as_i(vu x) { union { vi i; vu u; } c; c.u = x; return c.i; }
inline vd as_d(vi x) { union { vd d; vi i; } c; c.i = x; return c.d; }

inline vd splat(double c) { return vd{c, c, c, c}; }
inline vu splat_u(uint64_t c) { return vu{c, c, c, c}; }
inline vi splat_i(int64_t c) { return vi{c, c, c, c}; }

// Comparison masks (all-ones / all-zeros lanes), as an unsigned bit pattern.
inline vu gt(vd a, vd b) { return as_u(static_cast<vi>(a > b)); }
inline vu lt(vd a, vd b) { return as_u(static_cast<vi>(a < b)); }
inline vu eq(vd a, vd b) { return as_u(static_cast<vi>(a == b)); }

// Lane-wise select: `mask ? a : b`, using the bit pattern of `mask`.
inline vd select(vu mask, vd a, vd b) {
    return as_d((mask & as_u(a)) | (~mask & as_u(b)));
}

// Horizontal sum of the four lanes.
inline double hsum(vd v) { return v[0] + v[1] + v[2] + v[3]; }

// --- transcendentals --------------------------------------------------------

// floor, via truncating int conversion with a correction for negatives.
inline vd v_floor(vd x) {
    vi t = __builtin_convertvector(x, vi);          // truncates toward zero
    vd tf = __builtin_convertvector(t, vd);
    return select(gt(tf, x), tf - splat(1.0), tf);  // step down where we rounded up
}

// e^x. Range-reduce x = n*ln2 + r with |r| <= ln2/2 (Cody–Waite split of ln2),
// approximate e^r by a degree-9 Taylor series, then scale by 2^n via the
// exponent bits. Inputs are clamped to a safe range first.
inline vd v_exp(vd x) {
    x = select(gt(x, splat(709.0)), splat(709.0), x);
    x = select(lt(x, splat(-708.0)), splat(-708.0), x);

    const vd log2e = splat(1.4426950408889634);
    const vd ln2_hi = splat(6.93145751953125e-1);    // ln2 split for an accurate r
    const vd ln2_lo = splat(1.42860682030941723212e-6);

    vd nf = v_floor(x * log2e + splat(0.5));          // nearest integer n
    vd r = x - nf * ln2_hi - nf * ln2_lo;

    vd p = splat(1.0 / 362880.0);                     // Horner: sum r^k / k!
    p = p * r + splat(1.0 / 40320.0);
    p = p * r + splat(1.0 / 5040.0);
    p = p * r + splat(1.0 / 720.0);
    p = p * r + splat(1.0 / 120.0);
    p = p * r + splat(1.0 / 24.0);
    p = p * r + splat(1.0 / 6.0);
    p = p * r + splat(0.5);
    p = p * r + splat(1.0);
    p = p * r + splat(1.0);

    vi n = __builtin_convertvector(nf, vi);
    vd pow2 = as_d((n + splat_i(1023)) << 52);        // 2^n from the exponent field
    return p * pow2;
}

// ln(x) for x > 0. Decompose x = m * 2^e with m in [sqrt(0.5), sqrt(2)), then
// ln(m) via the fast-converging atanh series ln(m) = 2*(s + s^3/3 + ...),
// s = (m-1)/(m+1) (|s| <= 0.172), and add e*ln2.
inline vd v_log(vd x) {
    vu xi = as_u(x);
    vd m = as_d((xi & splat_u(0x000fffffffffffffull)) | splat_u(0x3ff0000000000000ull));  // [1,2)
    vi e = as_i((xi >> 52) & splat_u(0x7ff)) - splat_i(1023);

    vu big = gt(m, splat(1.4142135623730951));        // pull m into [sqrt(0.5), sqrt(2))
    m = select(big, m * splat(0.5), m);
    e = e - as_i(big);                                // mask is -1 where true, so this adds 1
    vd ef = __builtin_convertvector(e, vd);

    vd s = (m - splat(1.0)) / (m + splat(1.0));
    vd s2 = s * s;
    vd poly = splat(1.0 / 11.0);
    poly = poly * s2 + splat(1.0 / 9.0);
    poly = poly * s2 + splat(1.0 / 7.0);
    poly = poly * s2 + splat(1.0 / 5.0);
    poly = poly * s2 + splat(1.0 / 3.0);
    poly = poly * s2 + splat(1.0);
    vd logm = splat(2.0) * s * poly;

    return ef * splat(0.6931471805599453) + logm;
}

// sqrt(x) for x >= 0. Bit-trick initial rsqrt estimate refined by three Newton
// steps, then sqrt = x * rsqrt; the x == 0 lane is forced to 0.
inline vd v_sqrt(vd x) {
    vd y = as_d(splat_u(0x5fe6eb50c7b537a9ull) - (as_u(x) >> 1));  // ~rsqrt(x)
    const vd half = splat(0.5), three_half = splat(1.5);
    y = y * (three_half - half * x * y * y);
    y = y * (three_half - half * x * y * y);
    y = y * (three_half - half * x * y * y);
    return select(eq(x, splat(0.0)), splat(0.0), x * y);
}

// Inverse standard-normal CDF, branchless — Acklam's rational approximation
// evaluated for all lanes and blended by region (matches the scalar inv_norm_cdf
// in qmc.hpp; the central region, ~95% of draws, is bit-identical).
inline vd inv_norm_cdf(vd p) {
    const vd a0 = splat(-3.969683028665376e+01), a1 = splat(2.209460984245205e+02),
             a2 = splat(-2.759285104469687e+02), a3 = splat(1.383577518672690e+02),
             a4 = splat(-3.066479806614716e+01), a5 = splat(2.506628277459239e+00);
    const vd b0 = splat(-5.447609879822406e+01), b1 = splat(1.615858368580409e+02),
             b2 = splat(-1.556989798598866e+02), b3 = splat(6.680131188771972e+01),
             b4 = splat(-1.328068155288572e+01);
    const vd c0 = splat(-7.784894002430293e-03), c1 = splat(-3.223964580411365e-01),
             c2 = splat(-2.400758277161838e+00), c3 = splat(-2.549732539343734e+00),
             c4 = splat(4.374664141464968e+00), c5 = splat(2.938163982698783e+00);
    const vd d0 = splat(7.784695709041462e-03), d1 = splat(3.224671290700398e-01),
             d2 = splat(2.445134137142996e+00), d3 = splat(3.754408661907416e+00);

    const vd one = splat(1.0);
    const vu lo = lt(p, splat(0.02425));
    const vu hi = gt(p, splat(0.97575));

    // central region
    vd q = p - splat(0.5), rr = q * q;
    vd mid = (((((a0 * rr + a1) * rr + a2) * rr + a3) * rr + a4) * rr + a5) * q /
             (((((b0 * rr + b1) * rr + b2) * rr + b3) * rr + b4) * rr + one);

    // lower tail
    vd ql = v_sqrt(splat(-2.0) * v_log(p));
    vd low = (((((c0 * ql + c1) * ql + c2) * ql + c3) * ql + c4) * ql + c5) /
             ((((d0 * ql + d1) * ql + d2) * ql + d3) * ql + one);

    // upper tail (negate the lower-tail form evaluated at 1-p)
    vd qh = v_sqrt(splat(-2.0) * v_log(one - p));
    vd high = -((((((c0 * qh + c1) * qh + c2) * qh + c3) * qh + c4) * qh + c5) /
                ((((d0 * qh + d1) * qh + d2) * qh + d3) * qh + one));

    vd r = select(lo, low, mid);
    return select(hi, high, r);
}

#endif  // PRICER_HAVE_SIMD

}  // namespace pricer::simd
