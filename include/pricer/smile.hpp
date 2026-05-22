// pricer/smile.hpp — least-squares calibration of a volatility smile.
//
// Fits a quadratic in log-moneyness k = log(K / S):
//     vol(k) ≈ a + b·k + c·k²
// to a set of (strike, implied-vol) quotes by ordinary least squares. The fit is
// closed form (3×3 normal equations), so it is exact and dependency-free — a
// small but genuine calibration that reproduces market quotes within tolerance.
#pragma once
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pricer {

struct Smile {
    double a, b, c;  // vol(k) = a + b*k + c*k^2,  k = log(K/S)
    double S;        // reference spot for moneyness

    double vol(double K) const {
        const double k = std::log(K / S);
        return a + b * k + c * k * k;
    }
};

// Calibrate a quadratic smile to (strikes, vols) quotes around spot S.
inline Smile calibrate_smile(double S, const std::vector<double>& strikes,
                             const std::vector<double>& vols) {
    if (strikes.size() != vols.size() || strikes.size() < 3)
        throw std::invalid_argument("calibrate_smile: need >= 3 matching quotes");

    // Accumulate the normal-equation moments for basis [1, k, k^2].
    double n = 0, Sx = 0, Sx2 = 0, Sx3 = 0, Sx4 = 0, Sy = 0, Sxy = 0, Sx2y = 0;
    for (std::size_t i = 0; i < strikes.size(); ++i) {
        const double k = std::log(strikes[i] / S), y = vols[i];
        const double k2 = k * k;
        n += 1; Sx += k; Sx2 += k2; Sx3 += k2 * k; Sx4 += k2 * k2;
        Sy += y; Sxy += k * y; Sx2y += k2 * y;
    }

    double M[3][3] = {{n, Sx, Sx2}, {Sx, Sx2, Sx3}, {Sx2, Sx3, Sx4}};
    double v[3] = {Sy, Sxy, Sx2y};

    // Gaussian elimination with partial pivoting.
    for (int col = 0; col < 3; ++col) {
        int piv = col;
        for (int r = col + 1; r < 3; ++r)
            if (std::fabs(M[r][col]) > std::fabs(M[piv][col])) piv = r;
        if (piv != col) {
            for (int c = 0; c < 3; ++c) std::swap(M[piv][c], M[col][c]);
            std::swap(v[piv], v[col]);
        }
        if (std::fabs(M[col][col]) < 1e-300)
            throw std::runtime_error("calibrate_smile: singular system");
        for (int r = col + 1; r < 3; ++r) {
            const double f = M[r][col] / M[col][col];
            for (int c = col; c < 3; ++c) M[r][c] -= f * M[col][c];
            v[r] -= f * v[col];
        }
    }
    double x[3];
    for (int r = 2; r >= 0; --r) {
        double s = v[r];
        for (int c = r + 1; c < 3; ++c) s -= M[r][c] * x[c];
        x[r] = s / M[r][r];
    }
    return Smile{x[0], x[1], x[2], S};
}

}  // namespace pricer
