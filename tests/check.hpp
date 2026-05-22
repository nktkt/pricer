// tests/check.hpp — a tiny, dependency-free assertion helper for CTest.
// Each test is an executable that returns non-zero on failure. This keeps the
// project buildable offline; it can be swapped for Catch2/GoogleTest later
// (see ROADMAP Phase 1).
#pragma once
#include <cmath>
#include <cstdio>

namespace check {

inline int& failures() { static int f = 0; return f; }

// Assert |got - want| <= tol.
inline void approx(const char* name, double got, double want, double tol) {
    const double err = std::fabs(got - want);
    if (err > tol) {
        std::printf("  FAIL %-28s got %.10g  want %.10g  (err %.3g > tol %.3g)\n",
                    name, got, want, err, tol);
        failures()++;
    } else {
        std::printf("  ok   %-28s %.10g\n", name, got);
    }
}

inline void is_true(const char* name, bool cond) {
    if (!cond) { std::printf("  FAIL %s\n", name); failures()++; }
    else       { std::printf("  ok   %s\n", name); }
}

// Call from main(): returns a process exit code (0 = all passed).
inline int report(const char* suite) {
    if (failures()) {
        std::printf("[%s] %d failure(s)\n", suite, failures());
        return 1;
    }
    std::printf("[%s] all passed\n", suite);
    return 0;
}

}  // namespace check
