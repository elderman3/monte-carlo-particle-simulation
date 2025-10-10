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

// cl /std:c++17 /Od /Zi /FS /RTC1 /EHsc /MDd mc2.cpp /Fe:mc2.exe

constexpr double kB = 8.617333262e-11;
constexpr double th = 2.5e-8;
constexpr double M_PI = 3.14159265358979323846;

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

bool readMaterialBlock(std::istream& in, Material& mat) {
    if (!(in >> mat.sym >> mat.z >> mat.a >> mat.aw >> mat.T))
        return false;

    int neu_num;
    if (!(in >> neu_num)) return false;
    mat.neutrons.reserve(neu_num);
    for (int i=0; i<neu_num; ++i) {
        double E, nu;
        if (!(in >> E >> nu)) return false;
        mat.neutrons.emplace_back(E, nu);
    }

    int mt;
    double Q;
    int nc;
    while (in >> mt >> Q >> nc) {
        mat.Qvals[mt] = Q;
        auto& xs = mat.mt[mt];
        xs.reserve(nc);
        for (int j=0; j<nc; ++j) {
            double E, sigma;
            if (!(in >> E >> sigma)) return false;
            xs.emplace_back(E, sigma);
        }
    }
    return true;
}

void fillData(std::vector<Material>& mats, std::vector<double>& x, int inelastic) {
    for (auto& m : mats) {
        auto& out1 = m.mt[1];
        auto& out4 = m.mt[4];
        out1.reserve(x.size());
        out4.reserve(x.size());
        // Get all availabe MT values
        std::vector<int> MTs;
        for (const auto& [k,v] : m.mt) if (k != 1 && k != 4) MTs.push_back(k);
        for (double d : x) {
            double sum1 = 0, sum4 = 0;
            for (auto mt : MTs) {
                double val = valueInterp(m.mt[mt], d);
                if (mt > 50 && mt < 92) {
                    sum4 += val;
                } else {
                    sum1 += val;
                }
            }
            if (!inelastic) {
                sum1 += sum4;
            }
        out1.emplace_back(d, sum1);
        out4.emplace_back(d, sum4);
        }
    }
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

double neutronSpeed(double E){ return std::sqrt(std::max(0.0, 2.0*E)); }

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

StatsOut computeStats(const std::vector<std::vector<std::vector<int>>>& statM) {
    const size_t I = statM.size();
    const size_t M = I? statM[0].size() : 0;
    const size_t R = (M? statM[0][0].size() : 0);
    StatsOut out;
    out.mean.assign(M, std::vector<double>(R, 0.0));
    out.relErr.assign(M, std::vector<double>(R, 0.0));
    out.sum.assign(M, std::vector<int>(R, 0));

    for (size_t m=0; m<M; ++m) {
        for (size_t r=0; r<R; ++r) {
            long long S = 0; long double Q = 0.0L;
            for (size_t i=0; i<I; ++i) {
                int c = statM[i][m][r];
                S += c;
                Q += 1.0L * c * c;
            }
            const double N = double(I);
            const double mu = (I? double(S)/N : 0.0);
            const double var = (I>1)? double((Q - (1.0L*S*S)/N) / (N - 1.0)) : 0.0;
            const double se = (I>0)? std::sqrt(std::max(0.0, var) / N) : 0.0;
            out.sum[m][r] = int(S);
            out.mean[m][r] = mu;
            out.relErr[m][r] = (mu>0.0)? se/mu : 0.0;
        }
    }
    return out;
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

std::string timePathCol() {
    using clock = std::chrono::system_clock;
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(
                        clock::now().time_since_epoch()).count();
    return "output/col_" + std::to_string(secs) + ".csv";
}

std::string timePathTime() {
    using clock = std::chrono::system_clock;
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(
                        clock::now().time_since_epoch()).count();
    return "output/time_" + std::to_string(secs) + ".csv";
}

std::string timePathkeff() {
    using clock = std::chrono::system_clock;
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(
                        clock::now().time_since_epoch()).count();
    return "output/keff_" + std::to_string(secs) + ".csv";
}

bool storeDatakeff(const std::vector<double>& data)
{
    std::ofstream os(timePathkeff());
    if (!os) return false;
    os << std::scientific << std::setprecision(6);
    for (size_t k = 0; k < data.size(); ++k)
        os << k << "," << data[k] << "\n";
    return true;
}

bool storeDataCol(const std::vector<double>& data)
{
    std::ofstream os(timePathCol());
    if (!os) return false;
    os << std::scientific << std::setprecision(6);
    for (size_t k = 0; k < data.size(); ++k)
        os << k << "," << data[k] << "\n";
    return true;
}

static void storeTimeHist(const TimeHist& H){
    std::ofstream os(timePathTime());
    os << std::scientific << std::setprecision(6);
    for (int k=0;k<H.nbins;++k){
        const double tmid = (k + 0.5) * H.dt;
        os << tmid << "," << H.counts[k] << "\n";
    }
}

std::string mtLabel(int mt) {
    switch (mt) {
        case 2:   return "MT2 Elastic";
        case 4:   return "MT4 Inelastic";
        case 18:  return "MT18 Fission";
        case 102: return "MT102 Capture";
        default:  return "MT" + std::to_string(mt);
    }
}

double calMeanF(const std::vector<int>& fNeutrons, int nNeutrons) {
    long long S = 0; for (int x : fNeutrons) S += x;
    const int I = int(fNeutrons.size());
    return (I > 0 && nNeutrons > 0) ? double(S) / (double(I) * nNeutrons) : 0.0;
}

FourTally sumFour(const std::vector<FourTally>& v) {
    FourTally s;
    for (const auto& x : v) {
        s.fissionBirthsTotal += x.fissionBirthsTotal;
        s.fissionBirthsThermal += x.fissionBirthsThermal;
        s.absThTotal += x.absThTotal;
        s.absThFuel += x.absThFuel;
        s.started += x.started;
        s.reachedThermal += x.reachedThermal;
        s.resAbsBeforeThermal+= x.resAbsBeforeThermal;
    }
    return s;
}

FourFactors factorsFrom(const FourTally& T) {
    FourFactors F;
    F.eta = (T.absThFuel > 0) ? double(T.fissionBirthsThermal)/double(T.absThFuel) : 0.0;
    F.eps = (T.fissionBirthsThermal > 0) ? double(T.fissionBirthsTotal)/double(T.fissionBirthsThermal) : 0.0;
    F.p = (T.started > 0) ? double(T.reachedThermal) / double(T.reachedThermal + T.resAbsBeforeThermal) : 0.0;
    F.f = (T.absThTotal > 0) ? double(T.absThFuel)/double(T.absThTotal) : 0.0;
    F.keff = F.eta * F.eps * F.p * F.f;
    return F;
}

FourFactors averageFour(const std::vector<FourTally>& perIter) {
    return factorsFrom(sumFour(perIter));
}

void printStatsOut(const StatsOut& S, const std::vector<std::string>& matNames, const std::vector<int>& MTs, std::ostream& os = std::cout) {
    const size_t M = S.sum.size();
    if (!M) { os << "(no data)\n"; return; }
    const size_t R = S.sum[0].size();

    std::vector<long long> rowTot(M,0), colTot(R,0);
    long long grand = 0;
    for (size_t i=0;i<M;++i)
        for (size_t j=0;j<R;++j){
            long long v = S.sum[i][j];
            rowTot[i]+=v; colTot[j]+=v; grand+=v;
        }
    auto pct=[&](long long x){ return grand? 100.0*double(x)/double(grand) : 0.0; };

    size_t wName = 4;
    for (size_t i=0;i<std::min(M,matNames.size());++i) wName = std::max(wName, matNames[i].size());
    std::vector<size_t> wCol(R, 14);
    for (size_t j=0;j<R;++j){
        wCol[j] = std::max(wCol[j], mtLabel(MTs[j]).size());
        for (size_t i=0;i<M;++i){
            if (S.sum[i][j]==0) continue;
            std::ostringstream ss;
            double rPct = std::isfinite(S.relErr[i][j]) ? 100.0*S.relErr[i][j] : 0.0;
            ss << std::setprecision(3) << std::scientific << S.mean[i][j]
               << " ± " << std::fixed << std::setprecision(1) << rPct << "% "
               << "(" << S.sum[i][j] << ")";
            wCol[j] = std::max<size_t>(wCol[j], ss.str().size());
        }
    }

    os << "\n=== Statistics ===\n";
    os << "Total events: " << grand << "\n";

    os << "\n-- By reaction (counts) --\n";
    for (size_t j=0;j<R;++j)
        os << std::left << std::setw(int(wCol[j])) << mtLabel(MTs[j]) << "  "
           << std::right << std::setw(10) << colTot[j] << "  "
           << std::fixed << std::setprecision(2) << std::setw(6) << pct(colTot[j]) << "%\n";

    os << "\n-- By material (counts) --\n";
    for (size_t i=0;i<M;++i){
        const std::string& name = (i<matNames.size()? matNames[i] : ("mat"+std::to_string(i)));
        os << std::left << std::setw(int(wName)) << name << "  "
           << std::right << std::setw(10) << rowTot[i] << "  "
           << std::fixed << std::setprecision(2) << std::setw(6) << pct(rowTot[i]) << "%\n";
    }

    os << "\n-- Matrix: mean +- relErr% (count) --\n";
    os << std::left << std::setw(int(wName)) << "" << "  ";
    for (size_t j=0;j<R;++j)
        os << std::left << std::setw(int(wCol[j])) << mtLabel(MTs[j]) << "  ";
    os << "\n";

    for (size_t i=0;i<M;++i){
        const std::string& name = (i<matNames.size()? matNames[i] : ("mat"+std::to_string(i)));
        os << std::left << std::setw(int(wName)) << name << "  ";
        for (size_t j=0;j<R;++j){
            if (S.sum[i][j]==0) {
                os << std::left << std::setw(int(wCol[j])) << "-" << "  ";
            } else {
                double rPct = std::isfinite(S.relErr[i][j]) ? 100.0*S.relErr[i][j] : 0.0;
                std::ostringstream cell;
                cell << std::setprecision(3) << std::scientific << S.mean[i][j]
                     << " +- " << std::fixed << std::setprecision(1) << rPct << "% "
                     << "(" << S.sum[i][j] << ")";
                os << std::left << std::setw(int(wCol[j])) << cell.str() << "  ";
            }
        }
        os << "\n";
    }
}

void printFF(const FourFactors& F, std::ostream& os = std::cout) {
    os << std::fixed << std::setprecision(6)
       << "FourFactors{eta=" << F.eta
       << ", eps=" << F.eps
       << ", p="   << F.p
       << ", f="   << F.f
       << ", keff="<< F.keff
       << "}\n";
}

void printFourTally(const FourTally& T, std::ostream& os = std::cout) {
    os << "FourTally{"
       << "fissionBirthsTotal=" << T.fissionBirthsTotal
       << ", fissionBirthsThermal=" << T.fissionBirthsThermal
       << ", absThTotal=" << T.absThTotal
       << ", absThFuel=" << T.absThFuel
       << ", started=" << T.started
       << ", reachedThermal=" << T.reachedThermal
       << ", resAbsBeforeThermal=" << T.resAbsBeforeThermal
       << "}\n";
}

void printMaterial(const Material& m, std::ostream& os = std::cout) {
    auto print_samples = [&](const std::vector<std::pair<double,double>>& v) {
        os << "size=" << v.size();
        if (!v.empty()) {
            os << ", first(E,σ)=(" << std::scientific << v.front().first
               << "," << v.front().second << ")"
               << ", last(E,σ)=(" << v.back().first
               << "," << v.back().second << ")" << std::defaultfloat;
        }
    };

    os << "Material {\n";
    os << "  sym=\"" << m.sym << "\""
       << ", Z=" << m.z
       << ", A=" << m.a
       << ", aw=" << m.aw
       << ", T=" << m.T
       << ", rho=" << m.rho
       << ", proportion=" << m.proportion << "\n";

    os << "  neutrons ν̄(E): ";
    if (m.neutrons.empty()) {
        os << "none\n";
    } else {
        os << "N=" << m.neutrons.size() << ", samples: ";
        const size_t kmax = std::min<size_t>(5, m.neutrons.size());
        for (size_t i = 0; i < kmax; ++i)
            os << "(" << m.neutrons[i].first << "," << m.neutrons[i].second << ")"
               << (i+1<kmax? ", ":"");
        if (m.neutrons.size() > kmax) os << ", ...";
        os << "\n";
    }

    os << "  Qvals (by MT): ";
    if (m.Qvals.empty()) {
        os << "none\n";
    } else {
        bool first = true;
        for (const auto& [mt,q] : m.Qvals) {
            os << (first? "" : ", ") << "MT" << mt << "=" << q;
            first = false;
        }
        os << "\n";
    }

    os << "  MT tables:\n";
    if (m.mt.empty()) {
        os << "    none\n";
    } else {
        for (const auto& [mt,vec] : m.mt) {
            os << "    MT" << mt << ": ";
            print_samples(vec);
            os << "\n";
        }
    }
    os << "}\n";
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
    TimeHist timeHist(1e-2, 1000);
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
            const double path = -std::log(xi) / sTot;
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

Geometry readGeometry() {
    // 
}

Universe readGeometryFile() {
    //
}



int main() {
    std::ifstream file("input.txt");
    if (!file) {
        std::cerr << "File opening failed\n";
        return 1;
    }

    auto x = logspace(-11.0, std::log10(20.0), 500); // same x values for all simulations

    std::string line;
    while (std::getline(file, line)) {
        char command[32];
        int nNeutrons, iterations, nMaterials, maxSteps, inelastic;
        float energy;
        if (std::sscanf(line.c_str(), "%31s %d %f %d %d %d %d", 
        command, &nNeutrons, &energy, &iterations, &nMaterials, &maxSteps, &inelastic) == 7) {
            // Standard simulation for standard parameters
            if (std::strcmp(command, "simulation") == 0) {
                // Read and store all material information for given material
                const std::vector<int> MTs = inelastic ? std::vector<int>{2,4,18,102}
                                       : std::vector<int>{2,18,102};
                std::vector<Material> mats;
                mats.reserve(nMaterials);
                for (int i = 0; i < nMaterials; ++i) {
                    if (!std::getline(file, line)) {
                        std::cerr << "Missing material spec line " << i << "\n";
                        return 1;
                    }
                    char fname[32];
                    double rho, relativeMoles;
                    if (std::sscanf(line.c_str(), "%31s %lf %lf", fname, &rho, &relativeMoles) == 3) {
                        Material mat;
                        mat.rho = rho;
                        mat.proportion = relativeMoles;
                        std::string path = std::string("data/") + fname + ".dat";
                        std::ifstream materialdata(path);
                        if (!readMaterialBlock(materialdata, mat)) {
                            std::cerr << "Block read fail " << i << "\n";
                        }
                        mats.push_back(std::move(mat));
                    }
                }
                std::vector<std::string> matNames; matNames.reserve(mats.size());
                for (auto& m : mats) matNames.push_back(m.sym);
                fillData(mats, x, inelastic);

                SimRes simRes = simulation(nNeutrons, energy, iterations, maxSteps, inelastic, mats);
                
                StatsOut matrixStats = computeStats(simRes.statM);
                std::vector<double> collisionEnergy = computeColEnergy(simRes.collisions);
                double meanF = calMeanF(simRes.fNeutrons, nNeutrons);
                
                FourFactors ff = averageFour(simRes.fVec);

                printStatsOut(matrixStats, matNames, MTs);
                std::cout << "Average Neutrons Produced " << meanF << "\n"; 
                printFF(ff);
                storeDataCol(collisionEnergy);
                storeTimeHist(simRes.timeHist);

            } else if (std::strcmp(command, "criticality") == 0) {
                const double N_H2O = 0.0334;
                const double N_UO2 = 0.0244672;

                const double enr = 0.05;
                const int nNeutrons = 100;
                const int iterations = 10;
                const int maxSteps = 1000;
                const int inelastic = 1;
                const double E0 = 0.5;

                std::vector<double> meankeff(200);

                for (int pctFuel = 2; pctFuel <= 20; ++pctFuel) {
                    const double phi_UO2 = pctFuel / 100.0;
                    const double phi_H2O = 1.0 - phi_UO2;

                    Material matH1, matO16_H2O, matU235, matU238, matO16_UO2;
                    {
                        std::ifstream dH1("data/H1.dat");
                        std::ifstream dO16a("data/O16.dat");
                        std::ifstream dU235("data/U235.dat");
                        std::ifstream dU238("data/U238.dat");
                        std::ifstream dO16b("data/O16.dat");
                        readMaterialBlock(dH1, matH1);
                        readMaterialBlock(dO16a, matO16_H2O);
                        readMaterialBlock(dU235, matU235);
                        readMaterialBlock(dU238, matU238);
                        readMaterialBlock(dO16b, matO16_UO2);
                    }
                    matO16_H2O.sym = "O16_H2O";
                    matO16_UO2.sym = "O16_UO2";

                    matH1.rho = N_H2O * phi_H2O; matH1.proportion = 2;
                    matO16_H2O.rho = N_H2O * phi_H2O; matO16_H2O.proportion = 1;

                    matU235.rho = N_UO2 * phi_UO2 * enr; matU235.proportion = 1;
                    matU238.rho = N_UO2 * phi_UO2 * (1-enr); matU238.proportion = 1;
                    matO16_UO2.rho = N_UO2 * phi_UO2; matO16_UO2.proportion = 2;

                    std::vector<Material> mats;
                    mats.reserve(5);
                    mats.push_back(std::move(matH1));
                    mats.push_back(std::move(matO16_H2O));
                    mats.push_back(std::move(matU235));
                    mats.push_back(std::move(matU238));
                    mats.push_back(std::move(matO16_UO2));

                    fillData(mats, x, inelastic);
                    SimRes simRes = simulation(nNeutrons, E0, iterations, maxSteps, inelastic, mats);

                    FourFactors ff = averageFour(simRes.fVec);
                    double meanF = calMeanF(simRes.fNeutrons, nNeutrons);

                    std::cout << "phi_UO2=" << phi_UO2
                            << "  phi_H2O=" << phi_H2O
                            << "  meanF=" << meanF << "\n";
                    printFF(ff);
                    printFourTally(simRes.fVec[0]);

                    meankeff[pctFuel] = ff.keff;
                }
                storeDatakeff(meankeff);
            }
        } else {
            std::cout << "Command fail: " << line << '\n';
        }
    }
    return 1;
}