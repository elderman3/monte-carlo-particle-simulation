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


constexpr double PI_2 = 1.57079632679489661923;
constexpr double kB = 8.617333262e-11;

inline size_t pickIndex(const std::vector<double>& cum, double u) {
    auto it = std::upper_bound(cum.begin(), cum.end(), u);
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

inline double interpMT(const std::map<int, std::vector<std::pair<double,double>>>& mt, int code, double E) {
    auto it = mt.find(code);
    return it == mt.end() ? 0.0 : valueInterp(it->second, E);
}

double randomVal(float min, float max) {
    static thread_local std::mt19937 gen(std::random_device{}()); // engine
    std::uniform_real_distribution<double> dist(min, max);        // [0,1)

    return dist(gen);
}

Result timing(float (*f)(int, float), int iter, float length) {
    auto start = std::chrono::steady_clock::now();
    float pi = f(iter, length);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    Result result{elapsed.count(), pi};
    return result;
}

Stats statistics(const std::vector<Result> inputV) {
    Stats result;
    int len = inputV.size();
    float sumTime = 0, sumPi = 0;
    for (int i = 0; i < len; ++i) {
        sumTime += inputV[i].runningTime;
        sumPi += inputV[i].pi;
    }
    result.time = sumTime; 
    result.piMean = sumPi / (float) len;
    result.timeMean = sumTime / (float) len;

    float stdSumTime = 0, stdSumPi = 0;
    float moment3 = 0, moment4 = 0;
    for (int i = 0; i < len; ++i) {
        stdSumTime += pow(result.timeMean - inputV[i].runningTime, 2);
        stdSumPi += pow(result.piMean - inputV[i].pi, 2);
        moment3 += pow(result.piMean - inputV[i].pi, 3);
        moment4 += pow(result.piMean - inputV[i].pi, 4);
    }
    result.piStd = sqrt(stdSumPi / len);
    result.timeStd = sqrt(stdSumTime / len);
    result.fom = 1 / (pow(result.piStd, 2) * result.time);

    float s = moment3 / (pow(result.piStd, 3) * len);
    float k = moment4 / (pow(result.piStd, 4) * len);
    result.jb = ((float) len / 6) * pow(s, 2) + ((float) len / 24) * pow((k - 3),2);
    if (result.jb > 5.991) result.normalized = 0;
    if (result.jb <= 5.991) result.normalized = 1;
    return result;
}

void printStats(const Stats& s) {
    std::cout << "Stats:\n"
              << "  piMean   = " << s.piMean     << "\n"
              << "  piStd    = " << s.piStd      << "\n"
              << "  time     = " << s.time       << "\n"
              << "  timeMean = " << s.timeMean   << "\n"
              << "  timeStd  = " << s.timeStd    << "\n"
              << "  fom      = " << s.fom        << "\n"
              << "  jb       = " << s.jb         << "\n"
              << "  normal   = " << s.normalized << "\n"
              << "\n" 
              << "  iterations = " << s.iterations << "\n"
              << "  trials     = " << s.trials     << "\n"
              << "  length     = " << s.length     << "\n";
}

bool readMaterialBlock(std::istream& in, Material& mat) {
    mat = Material();
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
    const size_t NE = x.size();
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
            if (inelastic) {
                sum1 += sum4;
            }
        out1.emplace_back(d, sum1);
        out4.emplace_back(d, sum4);
        }
    }
}

inline double materialWeight(const Material& m, double E, const std::vector<int>& mts_total) {
    double micro = 0.0;
    for (int mt : mts_total) micro += interpMT(m.mt, mt, E);
    return std::max(0.0, m.proportion * m.rho * micro);
}

inline void buildMaterialCum(const std::vector<Material>& mats, double E,
                               const std::vector<int>& mts_total,
                               std::vector<double>& cum, double& total) {
    cum.clear(); cum.reserve(mats.size());
    double run = 0.0;
    for (const auto& m : mats) { run += materialWeight(m, E, mts_total); cum.push_back(run); }
    total = run;
}

inline void buildReactionCum(const Material& m, double E,
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

inline int sampleMultiplicity(double nuBar) {
    const int n = (int)std::floor(nuBar);
    const double frac = nuBar - n;
    return n + (randomVal() < frac ? 1 : 0);
}

inline std::array<double,3> iso_dir() {
    const double u1 = randomVal(), u2 = randomVal();
    const double mu = 2.0*u1 - 1.0;
    const double phi = 2.0*M_PI*u2;
    const double s = std::sqrt(std::max(0.0, 1.0 - mu*mu));
    return { s*std::cos(phi), s*std::sin(phi), mu };
}

inline std::array<double,3> scale(const std::array<double,3>& a, double s){
    return {s*a[0], s*a[1], s*a[2]};
}

inline std::array<double,3> add(const std::array<double,3>& a, const std::array<double,3>& b){
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

inline double norm2(const std::array<double,3>& a){ return a[0]*a[0]+a[1]*a[1]+a[2]*a[2]; }

inline std::array<double,3> normalize(const std::array<double,3>& a){
    double n = std::sqrt(std::max(0.0, norm2(a)));
    return (n>0)? std::array<double,3>{a[0]/n, a[1]/n, a[2]/n} : std::array<double,3>{1,0,0};
}

inline double sampleMaxwellE(double T) {
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

inline void elasticScatter(double En, double A, double TK, double Efg, double& Eout)
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

inline void recordCollision(Collisions& c, int collIdx, double E) {
    if (collIdx >= (int)c.num.size()) {
        c.num.resize(collIdx + 1, 0);
        c.sumEnergy.resize(collIdx + 1, 0.0);
    }
    c.num[collIdx] += 1;
    c.sumEnergy[collIdx] += E;
}

SimRes simulation(int nNeutrons, double energy, int iterations, int maxSteps, int inelastic, std::vector<Material> mats, std::vector<double> x) {
    // Data collection
    std::vector<Collisions> collisions(iterations);
    for (int i = 0; i < iterations; ++i)
        collisions.emplace_back(Collisions());
    std::vector<int> fNeutrons(iterations, 0);

    // Define modes of possible interaction
    static const std::vector<int> MTs = {2, 18, 102};
    if (inelastic) {
        static const std::vector<int> MTs = {2, 4, 18, 102};
    }
    std::vector<std::vector<std::vector<int>>> statM(
            iterations, std::vector<std::vector<int>>(mats.size(), std::vector<int>(MTs.size(), 0)));

    std::deque<Neutron> bank;
    std::vector<double> matCum;
    std::vector<int> rxLabels;
    std::vector<double> rxCum;
    double matTotal = 0.0, rxTotal = 0.0;


    for (int iter = 0; iter < iterations; ++iter) {
        bank.clear();
        // Init neutrons
        for (int i = 0; i < nNeutrons; ++i) {
            Neutron n = {energy, 0};
            bank.emplace_back(std::move(n));
        }
        auto& col = collisions[iter];
        
        int i = 0;
        while (i < nNeutrons * maxSteps && !bank.empty()) {
            Neutron n = bank.front(); bank.pop_front();
            double E = n.energy; n.collisions++;
            
            buildMaterialCum(mats, E, MTs, matCum, matTotal);
            double u1 = randomVal() * matTotal;
            size_t imat = pickIndex(matCum, u1);

            buildReactionCum(mats[imat], E, MTs, rxLabels, rxCum, rxTotal);
            double u2 = randomVal() * rxTotal;
            size_t irx = pickIndex(rxCum, u2);
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
                n.energy = Eout;
                recordCollision(col, n.collisions, E);
            }
            if (mtChosen == 4) {

                
                recordCollision(col, n.collisions, E);
            }
            if (mtChosen == 18) { // Fission                
                const int nEmit = sampleMultiplicity(valueInterp(mats[imat].neutrons, E));
                fNeutrons[iter] += nEmit;
                for (int k=0;k<nEmit;++k) {
                    double Eout;
                    
                    const double T = 1.2895; // Standard value for U235
                    Eout = sampleMaxwellE(T);

                    Neutron neu = {Eout, 0};
                    bank.emplace_back(std::move(neu));
                }
            }
            if (mtChosen == 102) continue;
            i++;
        }
    }
    SimRes simRes = {statM, collisions, fNeutrons};
    return simRes;
}

inline StatsOut computeStats(const std::vector<std::vector<std::vector<int>>>& statM) {
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
            const double se  = (I>0)? std::sqrt(std::max(0.0, var) / N) : 0.0;
            out.sum[m][r]    = int(S);
            out.mean[m][r]   = mu;
            out.relErr[m][r]= (mu>0.0)? se/mu : 0.0;
        }
    }
    return out;
}

inline std::vector<double> computeColEnergy(const std::vector<Collisions>& cols) {
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


// TODO: Reimplement Statistical methods for simulation data
// TODO: Activate Inelastic Scattering function
// TODO: Criticality detection tool
// TODO: Write input.txt to perform analysis
// TODO: Create plotting functions for relevant cases M3,O3
// TODO: Collect necessary variables for four-factor-formula
// TODO: Time dependence
// TODO: Report


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
            fillData(mats, x, inelastic);
            SimRes simRes = simulation(nNeutrons, energy, iterations, maxSteps, inelastic, mats, x);
            StatsOut matrixStats = computeStats(simRes.statM);
            std::vector<double> collisionEnergy = computeColEnergy(simRes.collisions);



            // Simulation is used for 
            // M3 - H1 and H2 slowing down neutrons - plot neutron energy/collision number
            // Multiplication of neutrons in Uranium 
            

            // Simulation must record everything, bonus 1 and 2 are easyish?


            // Calculate macroscopic corss section vector for all x


            // call simulation and other end things
                

            }


        } else {
            std::cout << "Command fail: " << line << '\n';
        }
    }

}