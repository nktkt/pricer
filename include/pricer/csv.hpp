// pricer/csv.hpp — minimal CSV read/write (comma-separated, no embedded commas).
//
// A tiny, dependency-free helper for the market-data adapters and result
// persistence. Lines are split on commas; surrounding whitespace and trailing
// carriage returns are trimmed; blank lines are skipped.
#pragma once
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace pricer {

inline std::string csv_trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r')) ++a;
    while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r')) --b;
    return s.substr(a, b - a);
}

// Read a CSV file into rows of string fields.
inline std::vector<std::vector<std::string>> read_csv(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("read_csv: cannot open " + path);
    std::vector<std::vector<std::string>> rows;
    std::string line;
    while (std::getline(f, line)) {
        if (csv_trim(line).empty()) continue;
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) fields.push_back(csv_trim(cell));
        rows.push_back(std::move(fields));
    }
    return rows;
}

// Write rows of string fields to a CSV file.
inline void write_csv(const std::string& path, const std::vector<std::vector<std::string>>& rows) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("write_csv: cannot open " + path);
    for (const auto& row : rows) {
        for (size_t i = 0; i < row.size(); ++i) {
            if (i) f << ',';
            f << row[i];
        }
        f << '\n';
    }
}

}  // namespace pricer
