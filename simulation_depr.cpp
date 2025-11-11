#include "mc.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <chrono>
#include <random>
#include <cmath>
#include <map>
#include <array>
#include <algorithm>
#include <deque>
#include <iomanip>
#include <cstring>
#include <cctype>
#include <stack>
#include <unordered_map>
#include <functional>
#include <limits>

using std::array;
using std::vector;



constexpr double kB = 8.617333262e-11;
constexpr double th = 2.5e-8;
constexpr double MEV_TO_J = 1.602176634e-13;
constexpr double M_N = 674927498e-27;


size_t pickIndex(const std::vector<double>& cum, double u) {
    auto it = std::upper_bound(cum.begin(), cum.end(), u);
    if (it == cum.end()) return cum.empty() ? 0 : cum.size() - 1;  // clamp
    return size_t(it - cum.begin());
}

std::vector<double> linspace(double start, double end, int n) {
    std::vector<double> x;
    x.reserve(n);

    if (n == 1) {
        x.push_back(start);
        return x;
    }

    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; ++i) {
        x.push_back(start + i * step);
    }
    return x;
}

std::vector<double> logspace(double exp_start, double exp_end, int n) {
    std::vector<double> exps = linspace(exp_start, exp_end, n);
    std::vector<double> x;
    x.reserve(n);

    for (double e : exps) {
        x.push_back(std::pow(10.0, e));
    }
    return x;
}

double valueInterp(const std::vector<std::pair<double,double>>& data, double target) {
    const size_t n = data.size();
    if (n < 2) return 0.0;
    if (target < data.front().first || target > data.back().first) return 0.0;

    auto it = std::lower_bound(
        data.begin(), data.end(), target,
        [](const auto& p, double t){ return p.first < t; }
    );

    if (it == data.begin()) return it->second;
    if (it == data.end())   return 0.0;

    const auto& [E2, s2] = *it;
    const auto& [E1, s1] = *(it - 1);
    const double denom = (E2 - E1);
    if (denom == 0.0) return s1;
    return s1 + (s2 - s1) * (target - E1) / denom;
}

double interpMT(const std::map<int, std::vector<std::pair<double,double>>>& mt, int code, double E) {
    auto it = mt.find(code);
    return it == mt.end() ? 0.0 : valueInterp(it->second, E);
}

double randomVal(float min, float max) {
    static thread_local std::mt19937 gen(std::random_device{}()); // engine
    std::uniform_real_distribution<double> dist(min, max);        // [0,1)

    return dist(gen);
}

double materialWeight(const Material& m, double E, const std::vector<int>& mts_total) {
    double micro = 0.0;
    for (int mt : mts_total) micro += interpMT(m.mt, mt, E);
    return std::max(0.0, m.proportion * m.rho * micro);
}

void buildMaterialCum(const std::vector<Material>& mats, double E,
                               const std::vector<int>& mts_total,
                               std::vector<double>& cum, double& total) {
    cum.clear(); cum.reserve(mats.size());
    double run = 0.0;
    for (const auto& m : mats) { run += materialWeight(m, E, mts_total); cum.push_back(run); }
    total = run;
}

void buildReactionCum(const Material& m, double E,
                               const std::vector<int>& mts_sample,
                               std::vector<int>& labels,
                               std::vector<double>& cum, double& total) {
    labels.clear(); cum.clear();
    double run = 0.0;
    for (int mt : mts_sample) {
        double w = interpMT(m.mt, mt, E);
        if (w <= 0.0) continue;
        run += w;
        labels.push_back(mt);
        cum.push_back(run);
    }
    total = run;
}

int sampleMultiplicity(double nuBar) {
    const int n = (int)std::floor(nuBar);
    const double frac = nuBar - n;
    return n + (randomVal() < frac ? 1 : 0);
}

std::array<double,3> iso_dir() {
    const double u1 = randomVal(), u2 = randomVal();
    const double mu = 2.0*u1 - 1.0;
    const double phi = 2.0*M_PI*u2;
    const double s = std::sqrt(std::max(0.0, 1.0 - mu*mu));
    return { s*std::cos(phi), s*std::sin(phi), mu };
}

std::array<double,3> scale(const std::array<double,3>& a, double s){
    return {s*a[0], s*a[1], s*a[2]};
}

std::array<double,3> add(const std::array<double,3>& a, const std::array<double,3>& b){
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

double norm2(const std::array<double,3>& a){ return a[0]*a[0]+a[1]*a[1]+a[2]*a[2]; }

std::array<double,3> normalize(const std::array<double,3>& a){
    double n = std::sqrt(std::max(0.0, norm2(a)));
    return (n>0)? std::array<double,3>{a[0]/n, a[1]/n, a[2]/n} : std::array<double,3>{1,0,0};
}

double sampleMaxwellE(double T) {
    double xi1, xi2, xi3, xi4, R;
    do {
        xi1 = randomVal(); xi2 = randomVal();
        R = xi1*xi1 + xi2*xi2;
    } while (R > 1.0 || R == 0.0);
    xi3 = randomVal();
    xi4 = randomVal();
    // logs are negative; the minus sign makes E positive
    return -T * ( (xi1*xi1) * std::log(xi3)/R + std::log(xi4) );
}

void elasticScatter(double En, double A, double TK, double Efg, double& Eout)
{
    std::array<double,3> VL{0,0,0};
    const bool use_free_gas = (En < Efg && A <= 10.0);
    if (use_free_gas) {
        const double T_mev = kB * TK;
        const double Et = sampleMaxwellE(T_mev);
        const double vmag = std::sqrt(std::max(0.0, 2.0*Et / A));
        VL = scale(iso_dir(), vmag);
    }

    const std::array<double,3> n_hat_in = std::array<double,3>{0,0,1};
    const double vL_mag = std::sqrt(std::max(0.0, 2.0*En));
    const std::array<double,3> vL = scale(n_hat_in, vL_mag);

    const std::array<double,3> VCM = scale(add(vL, scale(VL, A)), 1.0/(1.0 + A));
    const std::array<double,3> vC  = add(vL, scale(VCM, -1.0));
    const double vC_mag = std::sqrt(std::max(0.0, norm2(vC)));
    const std::array<double,3> vCprime = scale(iso_dir(), vC_mag);
    const std::array<double,3> vLprime = add(vCprime, VCM);

    Eout = 0.5 * norm2(vLprime);
}

void recordCollision(Collisions& c, int collIdx, double E) {
    if (collIdx >= (int)c.num.size()) {
        c.num.resize(collIdx + 1, 0);
        c.sumEnergy.resize(collIdx + 1, 0.0);
    }
    c.num[collIdx] += 1;
    c.sumEnergy[collIdx] += E;
}

double massRatioA(const Material& m) {
    return (m.a > 0 ? double(m.a) : std::max(1, m.z));
}

double elasticEnergyStationary(double En, double A) {
    const double a = (A-1.0)/(A+1.0);
    const double alpha = a*a;
    const double u = randomVal();
    return (alpha + (1.0 - alpha)*u) * En;
}

double inelasticEnergyStationary(double En, double A, double Delta) {
    if (Delta <= 0.0) return elasticEnergyStationary(En, A);
    if (En <= Delta)  return 0.0;
    const double Ein = En - Delta;
    return elasticEnergyStationary(Ein, A);
}

double getDeltaE(const Material& m, double) {
    auto it = m.Qvals.find(4);
    if (it != m.Qvals.end()) return std::abs(it->second);
    return 0.5;
}

double sigmaTot(const std::vector<Material>& mats, double E, const std::vector<int>& mts_total){
    double S = 0.0;
    for (const auto& m : mats){
        double micro = 0.0;
        for (int mt : mts_total) micro += interpMT(m.mt, mt, E);
        S += std::max(0.0, m.proportion * m.rho * micro);
    }
    return S;
}

double neutronSpeed(double E) {  // Nonrelativistic, as max(E) around 2MeV
    const double E_J = E * MEV_TO_J;
    return std::sqrt(std::max(0.0, 2.0*E_J / M_N)); // Returning m/s
}

bool isThermal(double E, const Material& m){
    return E <= 3.0 * kB * m.T;   // 3kT band; tighten to 1kT if needed
}

std::string normSym(std::string s){ std::string t; t.reserve(s.size());
  for (unsigned char c: s) if (std::isalnum(c)) t.push_back((char)std::toupper(c));
  return t; }
bool isU235(const Material& m){ return normSym(m.sym)=="U235"; }
bool isU238(const Material& m){ return normSym(m.sym)=="U238"; }
bool isFuel(const Material& m){
  auto s=normSym(m.sym); return s=="U235"||s=="U238"||s=="PU239"||s=="PU241";
}

std::vector<double> computeColEnergy(const std::vector<Collisions>& cols) {
    size_t K = 0;
    for (const auto& c : cols) K = std::max(K, std::min(c.num.size(), c.sumEnergy.size()));

    std::vector<long long> N(K, 0);
    std::vector<double> S(K, 0.0);

    for (const auto& c : cols) {
        const size_t kmax = std::min(c.num.size(), c.sumEnergy.size());
        for (size_t k = 0; k < kmax; ++k) {
            N[k] += c.num[k];
            S[k] += c.sumEnergy[k];
        }
    }

    std::vector<double> avg(K, 0.0);
    for (size_t k = 0; k < K; ++k) avg[k] = (N[k] > 0) ? (S[k] / double(N[k])) : 0.0;
    return avg;
}

SimRes simulation(int nNeutrons, double energy, int iterations, int maxSteps, int inelastic, std::vector<Material>& mats) {
    // Data collection
    std::vector<Collisions> collisions;
    for (int i = 0; i < iterations; ++i)
        collisions.emplace_back(Collisions());
    std::vector<int> fNeutrons(iterations, 0);

    // Define modes of possible interaction
    const std::vector<int> MTs = inelastic ? std::vector<int>{2,4,18,102}
                                       : std::vector<int>{2,18,102};
    std::vector<std::vector<std::vector<int>>> statM(
            iterations, std::vector<std::vector<int>>(mats.size(), std::vector<int>(MTs.size(), 0)));

    std::vector<FourTally> fourFV(iterations);
    //TimeHist timeHist(1e-2, 1000);
    TimeHist timeHist(1e-11, 1000);
    std::deque<Neutron> bank;
    std::vector<double> matCum;
    std::vector<int> rxLabels;
    std::vector<double> rxCum;
    double matTotal = 0.0, rxTotal = 0.0;

    for (int iter = 0; iter < iterations; ++iter) {
        bank.clear();
        // Init neutrons
        for (int i = 0; i < nNeutrons; ++i) {
            Neutron n = {energy, 0, 0, false, true};
            bank.emplace_back(std::move(n));
        }
        auto& col = collisions[iter];
        fourFV[iter].started += nNeutrons;
        long long i = 0;
        while (i < 1LL * nNeutrons * maxSteps && !bank.empty()) {
            Neutron n = bank.front(); bank.pop_front();
            double E = n.energy; n.collisions++;
            // Sample time
            const double sTot = sigmaTot(mats, E, MTs);
            if (sTot <= 0.0) { ++i; continue; }
            const double xi = randomVal();

            const double path = -std::log(xi) / (sTot * 1e-2); // Correction to m instead of cm
            const double v = neutronSpeed(E);
            n.time += (v > 0.0 ? path / v : 0.0);
            buildMaterialCum(mats, E, MTs, matCum, matTotal);
            if (matTotal <= 0.0 || matCum.empty()) { ++i; continue; }
            size_t imat = pickIndex(matCum, randomVal() * matTotal);
            if (imat >= mats.size()) { ++i; continue; }

            buildReactionCum(mats[imat], E, MTs, rxLabels, rxCum, rxTotal);
            if (rxTotal <= 0.0 || rxCum.empty() || rxLabels.empty()) { ++i; continue; }
            size_t irx = pickIndex(rxCum, randomVal() * rxTotal);
            if (irx >= rxLabels.size()) { ++i; continue; }
            int mtChosen = rxLabels[irx];
            // Update Reaction stats
            auto it = std::find(MTs.begin(), MTs.end(), mtChosen);
            if (it != MTs.end()) {
                size_t mtpos = size_t(it - MTs.begin());
                ++statM[iter][imat][mtpos];
            }
            if (mtChosen == 2) {
                double Eout;
                elasticScatter(E, mats[imat].aw, mats[imat].T, 2e-4, Eout);
                if (n.isSource && !n.reachedTh && isThermal(Eout, mats[imat])) { n.reachedTh = true; ++fourFV[iter].reachedThermal; }

                n.energy = Eout;
                recordCollision(col, n.collisions, E);
                bank.emplace_back(n);
                i++;
                continue;
            }
            if (mtChosen == 4) {
                const double A = massRatioA(mats[imat]);
                const double Delta = getDeltaE(mats[imat], E);
                const double Eout = inelasticEnergyStationary(E, A, Delta);

                recordCollision(col, n.collisions, E);
                if (n.isSource && !n.reachedTh && isThermal(Eout, mats[imat])) { n.reachedTh = true; ++fourFV[iter].reachedThermal; }
                n.energy = Eout;
                bank.emplace_back(n);
                i++;
                continue;
            }
            if (mtChosen == 18) { // Fission                
                const int nEmit = sampleMultiplicity(valueInterp(mats[imat].neutrons, E));
                fNeutrons[iter] += nEmit;
                if (1) {
                    fourFV[iter].fissionBirthsTotal += nEmit;
                    if (isThermal(E, mats[imat])) {
                        ++fourFV[iter].absThTotal;
                        if (isFuel(mats[imat])) ++fourFV[iter].absThFuel;
                        fourFV[iter].fissionBirthsThermal += nEmit;
                    }
                }

                for (int k=0;k<nEmit;++k) {
                    timeHist.add(n.time);
                    const double T = 1.2895; // Standard value for U235
                    double Eout = sampleMaxwellE(T);

                    Neutron neu = {Eout, 0, 0, false, false};
                    bank.emplace_back(std::move(neu));
                }
                i++;
                continue;
            }
            if (mtChosen == 102) {
                if (n.isSource && !n.reachedTh && isU238(mats[imat]) && E>=1e-6 && E<=1e-2)
                    ++fourFV[iter].resAbsBeforeThermal;
                if (isThermal(E, mats[imat])) {
                    ++fourFV[iter].absThTotal;
                    if (isFuel(mats[imat])) ++fourFV[iter].absThFuel;
                }
                i++;
                continue;
            }
            i++;
        }
    }
    SimRes simRes = {statM, collisions, fNeutrons, fourFV, timeHist};
    return simRes;
}