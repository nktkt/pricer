// bermudan_option.cpp — a Bermudan option as the bridge between European and
// American exercise.
//
// A Bermudan may be exercised early, but only on a finite schedule of dates. With
// one date it is European; as the schedule fills in it climbs toward the American
// value. This demo prices a Bermudan put on increasingly fine schedules with
// Longstaff–Schwartz LSM and shows the value rising from the European price
// (Black–Scholes) toward the American one (a binomial tree).
#include "pricer/american.hpp"
#include "pricer/bermudan.hpp"
#include "pricer/black_scholes.hpp"

#include <cstdio>

using namespace pricer;

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.20, T = 1.0;
    const long n_paths = 400'000;

    const double eur = black_scholes_put(S, K, r, sigma, T);
    const double am = binomial_price(OptionType::Put, S, K, r, sigma, T, 2000, true);

    std::printf("Bermudan put   S=%.0f K=%.0f r=%.2f sigma=%.2f T=%.1f   (%ld LSM paths)\n\n",
                S, K, r, sigma, T, n_paths);
    std::printf("  European (Black-Scholes)      : %.4f\n", eur);
    std::printf("  exercise dates -> Bermudan (LSM):\n");
    for (int m : {1, 2, 4, 12, 50, 250}) {
        const double v = bermudan_lsm(OptionType::Put, S, K, r, sigma, equally_spaced_dates(T, m), n_paths);
        std::printf("    %3d date(s)                 : %.4f\n", m, v);
    }
    std::printf("  American (binomial, 2000 steps): %.4f\n", am);

    // An irregular (non-uniform) exercise schedule is priced the same way.
    const std::vector<double> sched{0.25, 0.5, 0.75, 1.0};
    const double v = bermudan_lsm(OptionType::Put, S, K, r, sigma, sched, n_paths);
    std::printf("\n  quarterly schedule {0.25,0.5,0.75,1.0}: %.4f\n", v);
    return 0;
}
