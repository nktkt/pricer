// pricer/dual.hpp — forward-mode automatic differentiation (dual numbers).
//
// `Dual<N, T>` carries a value plus an N-dimensional gradient, propagating exact
// derivatives through any computation built from the overloaded operators and
// math functions below. One evaluation in `Dual<N>` yields all N first-order
// partials at once. The scalar type `T` may itself be a `Dual`, so nesting
// (e.g. `Dual<1, Dual<1>>`) gives exact second derivatives.
//
// This is forward-mode AD; it is exact (no bumping). True adjoint/reverse-mode
// (AAD) is a separate, heavier construction.
#pragma once
#include <array>
#include <cmath>

namespace pricer {

template <int N, class T = double>
struct Dual {
    T v{};                 // value
    std::array<T, N> d{};  // gradient (partials)

    Dual() = default;
    Dual(const T& val) : v(val) { d.fill(T(0)); }  // constant (zero gradient)
};

// --- construction helpers ---
template <int N>
inline Dual<N, double> ad_var(double x, int i) {  // active variable: ∂/∂i = 1
    Dual<N, double> r;
    r.v = x;
    r.d[i] = 1.0;
    return r;
}
template <int N>
inline Dual<N, double> ad_const(double x) {  // constant
    Dual<N, double> r;
    r.v = x;
    return r;
}

// --- arithmetic (Dual ⊕ Dual) ---
template <int N, class T>
Dual<N, T> operator+(const Dual<N, T>& a, const Dual<N, T>& b) {
    Dual<N, T> r; r.v = a.v + b.v;
    for (int i = 0; i < N; ++i) r.d[i] = a.d[i] + b.d[i];
    return r;
}
template <int N, class T>
Dual<N, T> operator-(const Dual<N, T>& a, const Dual<N, T>& b) {
    Dual<N, T> r; r.v = a.v - b.v;
    for (int i = 0; i < N; ++i) r.d[i] = a.d[i] - b.d[i];
    return r;
}
template <int N, class T>
Dual<N, T> operator*(const Dual<N, T>& a, const Dual<N, T>& b) {
    Dual<N, T> r; r.v = a.v * b.v;
    for (int i = 0; i < N; ++i) r.d[i] = a.v * b.d[i] + b.v * a.d[i];
    return r;
}
template <int N, class T>
Dual<N, T> operator/(const Dual<N, T>& a, const Dual<N, T>& b) {
    Dual<N, T> r; r.v = a.v / b.v;
    for (int i = 0; i < N; ++i) r.d[i] = (a.d[i] * b.v - a.v * b.d[i]) / (b.v * b.v);
    return r;
}
template <int N, class T>
Dual<N, T> operator-(const Dual<N, T>& a) {
    Dual<N, T> r; r.v = -a.v;
    for (int i = 0; i < N; ++i) r.d[i] = -a.d[i];
    return r;
}

// --- arithmetic with a plain double (constant) ---
template <int N, class T>
Dual<N, T> operator*(double s, const Dual<N, T>& a) {
    Dual<N, T> r; r.v = T(s) * a.v;
    for (int i = 0; i < N; ++i) r.d[i] = T(s) * a.d[i];
    return r;
}
template <int N, class T>
Dual<N, T> operator*(const Dual<N, T>& a, double s) { return s * a; }
template <int N, class T>
Dual<N, T> operator+(const Dual<N, T>& a, double s) {
    Dual<N, T> r = a; r.v = a.v + T(s); return r;
}
template <int N, class T>
Dual<N, T> operator+(double s, const Dual<N, T>& a) { return a + s; }
template <int N, class T>
Dual<N, T> operator-(const Dual<N, T>& a, double s) {
    Dual<N, T> r = a; r.v = a.v - T(s); return r;
}
template <int N, class T>
Dual<N, T> operator-(double s, const Dual<N, T>& a) { return (-a) + s; }
template <int N, class T>
Dual<N, T> operator/(double s, const Dual<N, T>& a) {
    Dual<N, T> num; num.v = T(s);  // constant s
    return num / a;
}
template <int N, class T>
Dual<N, T> operator/(const Dual<N, T>& a, double s) { return a * (1.0 / s); }

// --- math functions (chain rule); `using std::f` dispatches scalar vs nested Dual ---
template <int N, class T>
Dual<N, T> exp(const Dual<N, T>& a) {
    using std::exp;
    const T e = exp(a.v);
    Dual<N, T> r; r.v = e;
    for (int i = 0; i < N; ++i) r.d[i] = e * a.d[i];
    return r;
}
template <int N, class T>
Dual<N, T> log(const Dual<N, T>& a) {
    using std::log;
    Dual<N, T> r; r.v = log(a.v);
    for (int i = 0; i < N; ++i) r.d[i] = a.d[i] / a.v;
    return r;
}
template <int N, class T>
Dual<N, T> sqrt(const Dual<N, T>& a) {
    using std::sqrt;
    const T s = sqrt(a.v);
    Dual<N, T> r; r.v = s;
    for (int i = 0; i < N; ++i) r.d[i] = a.d[i] / (T(2.0) * s);
    return r;
}
template <int N, class T>
Dual<N, T> erfc(const Dual<N, T>& a) {
    using std::erfc;
    using std::exp;
    const T val = erfc(a.v);
    const T factor = T(-1.1283791670955126) * exp(-(a.v * a.v));  // d/dx erfc = -2/√π e^{-x²}
    Dual<N, T> r; r.v = val;
    for (int i = 0; i < N; ++i) r.d[i] = factor * a.d[i];
    return r;
}

}  // namespace pricer
