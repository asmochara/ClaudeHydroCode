#include "table2d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

Table2D::Table2D(const std::string& path, int nblocks) : path_(path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Table2D: cannot open " + path);

    std::vector<double> vals;
    std::string line;
    while (std::getline(in, line)) {
        auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        std::istringstream ss(line);
        double v;
        while (ss >> v) vals.push_back(v);
    }
    if (vals.size() < 2) throw std::runtime_error("Table2D: truncated file " + path);
    const size_t NR = static_cast<size_t>(vals[0]);
    const size_t NT = static_cast<size_t>(vals[1]);
    const size_t need = 2 + NR + NT + static_cast<size_t>(nblocks) * NR * NT;
    if (vals.size() != need)
        throw std::runtime_error("Table2D: expected " + std::to_string(need) +
                                 " numbers in " + path + " (" + std::to_string(nblocks) +
                                 " blocks), got " + std::to_string(vals.size()));

    size_t k = 2;
    lnRho_.resize(NR);
    for (size_t i = 0; i < NR; ++i) lnRho_[i] = std::log(vals[k++]);
    lnT_.resize(NT);
    for (size_t j = 0; j < NT; ++j) lnT_[j] = std::log(vals[k++]);
    for (int b = 0; b < nblocks; ++b) {
        blocks_.emplace_back(vals.begin() + k, vals.begin() + k + NR * NT);
        k += NR * NT;
    }

    for (size_t i = 1; i < NR; ++i)
        if (lnRho_[i] <= lnRho_[i - 1])
            throw std::runtime_error("Table2D: rho grid not ascending in " + path);
    for (size_t j = 1; j < NT; ++j)
        if (lnT_[j] <= lnT_[j - 1])
            throw std::runtime_error("Table2D: T grid not ascending in " + path);

    Tmin_ = std::exp(lnT_.front());
    Tmax_ = std::exp(lnT_.back());
}

namespace {
inline void locate(const std::vector<double>& g, double x, size_t& i, double& w) {
    if (x <= g.front()) { i = 0; w = 0.0; return; }
    if (x >= g.back())  { i = g.size() - 2; w = 1.0; return; }
    i = static_cast<size_t>(std::upper_bound(g.begin(), g.end(), x) - g.begin()) - 1;
    w = (x - g[i]) / (g[i + 1] - g[i]);
}
}  // namespace

double Table2D::interp(int block, double rho, double T) const {
    const auto& tab = blocks_[block];
    const size_t NT = lnT_.size();
    size_t i, j;
    double wi, wj;
    locate(lnRho_, std::log(std::max(rho, 1e-300)), i, wi);
    locate(lnT_,   std::log(std::max(T,   1e-300)), j, wj);
    const double f00 = tab[i * NT + j],       f01 = tab[i * NT + j + 1];
    const double f10 = tab[(i + 1) * NT + j], f11 = tab[(i + 1) * NT + j + 1];
    return (1 - wi) * ((1 - wj) * f00 + wj * f01) + wi * ((1 - wj) * f10 + wj * f11);
}

void Table2D::logBlock(int block) {
    for (auto& v : blocks_[block]) {
        if (v <= 0.0)
            throw std::runtime_error("Table2D: non-positive entry in log block of " + path_);
        v = std::log(v);
    }
}

double Table2D::interpLog(int block, double rho, double T) const {
    return std::exp(interp(block, rho, T));
}
