#include "mc.h"
#include <algorithm>
#include <numeric>
#include <random>
#include <vector>
#include <array>
#include <deque>
#include <cmath>
#include <iostream>
#include <chrono>
#include <fstream>
#include <cuda_runtime.h>


using std::array; using std::vector;
#define M_PI 3.14159265358979323846264338327950288
#define BLOCK 32

constexpr double kB = 8.617333262e-11;
constexpr double MEV_TO_J = 1.602176634e-13;
constexpr double M_N = 1.674927498e-27;
constexpr double EFG = 2e-4;

// --- CUDA ERROR HANDLING ---
#define CHECK(x) check(x, #x)
static inline void check(cudaError_t err, const char* context) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: " << context << ": "
            << cudaGetErrorString(err) << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

// --- MATH UTILITIES ---

__host__ __device__ array<double,3> scaleByC(const array<double,3>& a, double s) {return {s*a[0], s*a[1], s*a[2]};}

__host__ __device__ vector<double> linspace(double a,double b,int n) {
    vector<double> x; x.reserve(n);
    if (n<=1) { if(n==1) x.push_back(a); return x; }
    double h=(b-a)/(n-1); for(int i=0;i<n;++i) x.push_back(a+i*h); return x;
}

__host__ __device__ vector<double> logspace(double ea,double eb,int n) {
    auto e=linspace(ea,eb,n); vector<double> x; x.reserve(n);
    for(double t:e) x.push_back(std::pow(10.0,t)); return x;
}

__host__ __device__ double randomVal(float min=0.f, float max=1.f) {
    thread_local std::mt19937_64 gen{0x9817981276389};
    std::uniform_real_distribution<double> d(min,max);
    return d(gen);
}

__host__ __device__ int sampleMultiplicity(double nuBar) {
    if (nuBar<=0) return 0;
    const int n = (int)std::floor(nuBar);
    const double frac = nuBar - n;
    return n + (randomVal() < frac ? 1 : 0);
}

__host__ __device__ array<double,3> iso_dir() {
    const double u1 = randomVal(), u2 = randomVal();
    const double mu = 2.0*u1 - 1.0;
    const double phi = 2.0*M_PI*u2;
    const double s = std::sqrt(std::max(0.0, 1.0 - mu*mu));
    return { s*std::cos(phi), s*std::sin(phi), mu };
}


// --- PHYSICS HELPERS ---

__host__ __device__ double neutronSpeed(double E) {
    const double E_J = E * MEV_TO_J;
    return std::sqrt(std::max(0.0, 2.0*E_J / M_N));
}

__host__ __device__ void elasticScatter(double En, double A, double TK, double Efg) {
    array<double,3> VL{0,0,0};
    const bool use_free_gas = (En < Efg && A <= 10.0);
    if (use_free_gas) {
        const double T_mev = kB * TK;
        const double u = randomVal(), v = randomVal();
        const double Et = -T_mev * std::log(std::max(1e-32, u*v));
        const double vmag = std::sqrt(std::max(0.0, 2.0*Et / A));
        VL = scaleByC(iso_dir(), vmag);
    }
    const array<double,3> n_hat_in = array<double,3>{0,0,1};
    const double vL_mag = std::sqrt(std::max(0.0, 2.0*En));
    const array<double,3> vL = scaleByC(n_hat_in, vL_mag);
    const array<double,3> VCM = scaleByC(add3(vL, scaleByC(VL, A)), 1.0/(1.0 + A));
    const array<double,3> vC  = add3(vL, scaleByC(VCM, -1.0));
    const double vC_mag = std::sqrt(std::max(0.0, normed(vC)));
    const array<double,3> vCprime = scaleByC(iso_dir(), vC_mag);
    const array<double,3> vLprime = add3(vCprime, VCM);
    return 0.5 * normed(vLprime);
}

__host__ __device__ double elasticEnergyStationary(double En, double A) {
    const double a = (A-1.0)/(A+1.0);
    const double alpha = a*a;
    const double u = randomVal();
    return (alpha + (1.0 - alpha)*u) * En;
}

__host__ __device__ double inelasticEnergyStationary(double En, double A, double Delta) {
    if (Delta <= 0.0) return elasticEnergyStationary(En, A);
    if (En <= Delta)  return 0.0;
    const double Ein = En - Delta;
    return elasticEnergyStationary(Ein, A);
}

__host__ __device__ double getDeltaE(const Material& m) {
    auto it = m.Qvals.find(4);
    if (it != m.Qvals.end()) return std::abs(it->second);
    return 0.5;
}

__host__ __device__ double macroscopicSigmaAt(const Universe& u, const array<double,3>& pWorld, double E, const vector<int>& mts_total, const Geometry** gOut) {
    if (gOut) *gOut = nullptr;
    const Geometry* g = findGeometryAt(u, pWorld);
    if (!g) return 0.0;
    if (gOut) *gOut = g;
    return sigmaTot(g->mats, E, mts_total);
}

__host__ __device__ double majorSigma(const Universe& u, double E, const vector<int>& mts_total) {
    double m = 0.0;
    for (const Geometry& g : u.geometries)
        m = std::max(m, sigmaTot(g.mats, E, mts_total));
    for (const Universe& su : u.subUniverse)
        m = std::max(m, majorSigma(su, E, mts_total));
    return m;
}

__host__ __device__ double sigmaTot(const vector<Material>& mats, double E, const vector<int>& mts_total) {
    double S = 0.0;
    for (const auto& m : mats) {
        double micro = 0.0;
        for (int mt : mts_total) micro += interpMT(m.mt, mt, E);
        S += std::max(0.0, m.proportion * m.rho * micro);
    }
    return S;
}

__host__ __device__ bool sampleReactionAtE(const Geometry& g, double E, const vector<int>& MTs_total, const vector<int>& MTs_sample, RxSample& out) {
    vector<double> matCum; double matTot=0.0;
    buildMaterialCum(g.mats, E, MTs_total, matCum, matTot);
    if (matTot<=0.0) return false;
    size_t imat = pickIndex(matCum, randomVal()*matTot);
    if (imat>=g.mats.size()) return false;

    vector<int> rxLabels; vector<double> rxCum; double rxTot=0.0;
    buildReactionCum(g.mats[imat], E, MTs_sample, rxLabels, rxCum, rxTot);
    if (rxTot<=0.0 || rxLabels.empty()) return false;
    size_t irx = pickIndex(rxCum, randomVal()*rxTot);
    if (irx>=rxLabels.size()) return false;

    out.mt = rxLabels[irx];
    out.matIndex = (int)imat;
    out.mat = &g.mats[imat];
    return true;
}

__host__ __device__ bool isFuelSym(const std::string& s) { std::string t; for(unsigned char c: s) { if(std::isalnum(c)) t.push_back((char)std::toupper(c)); } return t=="U235"||t=="U238"||t=="PU239"||t=="PU241"; }
__host__ __device__ bool isThermalE(double E, const Material& m) { return E <= 3.0 * kB * m.T; }


// --- SIMULATION HELPERS ---

__host__ __device__ size_t pickIndex(const vector<double>& cum, double u) {
    auto it = std::upper_bound(cum.begin(), cum.end(), u);
    if (it == cum.end()) return cum.empty() ? 0 : cum.size() - 1;
    return size_t(it - cum.begin());
}

__host__ __device__ double valueInterp(const vector<std::pair<double,double>>& data, double target) {
    const size_t n = data.size();
    if (n < 2) return 0.0;
    if (target < data.front().first || target > data.back().first) return 0.0;
    auto it = std::lower_bound(data.begin(), data.end(), target,
        [](const auto& p, double t) { return p.first < t; });
    if (it == data.begin()) return it->second;
    if (it == data.end())   return 0.0;
    const auto& [E2, s2] = *it; const auto& [E1, s1] = *(it - 1);
    const double denom = (E2 - E1); if (denom == 0.0) return s1;
    return s1 + (s2 - s1) * (target - E1) / denom;
}

__host__ __device__ double interpMT(const std::map<int, vector<std::pair<double,double>>>& mt, int code, double E) {
    auto it = mt.find(code); return it == mt.end() ? 0.0 : valueInterp(it->second, E);
}

__host__ __device__ double materialWeight(const Material& m, double E, const vector<int>& mts_total) {
    double micro = 0.0; for (int mt : mts_total) micro += interpMT(m.mt, mt, E);
    return std::max(0.0, m.proportion * m.rho * micro);
}

__host__ __device__ void buildMaterialCum(const vector<Material>& mats, double E, const vector<int>& mts_total, vector<double>& cum, double& total) {
    cum.clear(); cum.reserve(mats.size()); total = 0.0;
    for (const auto& m : mats) { total += materialWeight(m, E, mts_total); cum.push_back(total); }
}

__host__ __device__ void buildReactionCum(const Material& m, double E, const vector<int>& mts_sample, vector<int>& labels, vector<double>& cum, double& total) {
    labels.clear(); cum.clear(); total = 0.0;
    for (int mt : mts_sample) {
        double w = interpMT(m.mt, mt, E);
        if (w <= 0.0) continue;
        total += w; labels.push_back(mt); cum.push_back(total);
    }
}

__host__ __device__ int rxCol(const vector<int>& MTs, int mt) {
    auto it = std::find(MTs.begin(), MTs.end(), mt);
    return (it==MTs.end())? -1 : int(it - MTs.begin());
}

__host__ __device__ void recordCollision(Collisions& c, int collIdx, double E) {
    if (collIdx >= (int)c.num.size()) { c.num.resize(collIdx + 1, 0); c.sumEnergy.resize(collIdx + 1, 0.0); }
    c.num[collIdx] += 1; c.sumEnergy[collIdx] += E;
}

__host__ __device__ int pickMaterialIndexAtE(const Geometry& g, double E, const std::vector<int>& MTs_total) {
    std::vector<double> matCum; double matTot=0.0;
    buildMaterialCum(g.mats, E, MTs_total, matCum, matTot);
    if (matTot<=0.0) return -1;
    size_t imat = pickIndex(matCum, randomVal()*matTot);
    if (imat>=g.mats.size()) return -1;
    return (int)imat;
}

__host__ __device__ double defaultNuBar(const Material& m) {
    std::string s; for(unsigned char c: m.sym) { if(std::isalnum(c)) s.push_back((char)std::toupper(c)); }
    if (s=="U235") return 2.43;
    if (s=="U238") return 2.50;
    if (s=="PU239") return 2.90;
    return 2.40;
}


// --- FLIGHT FUNCTIONS ---

__host__ __device__ FlightResult surfaceFlight(const Universe& u, Neutron& n, double E, const vector<int>& MTs_total, Mesh3D* meshTLE) {
    const double EPS = 1e-8;
    FlightResult fr;
    const Geometry* g0=nullptr;
    double Sigma = macroscopicSigmaAt(u, n.pos, E, MTs_total, &g0);
    if (Sigma<0) Sigma=0;
    double l = (Sigma>0)? (-std::log(std::max(1e-32, randomVal()))/Sigma) : 1e300;

    double dSurf = nearestCollision(n, u);
    if (dSurf < 0) {
        fr.leaked = true; fr.traveled = 0.0; return fr;
    }
    if (l < dSurf) {
        fr.collided = true;
        fr.SigmaLocal = Sigma;
        fr.geom = g0;
        fr.traveled = l;
        n.pos = add3(n.pos, scaleByC(n.dir, l));
    } else {
        fr.collided = false;
        fr.traveled = dSurf;
        fr.geom = g0;
        if (meshTLE) {
            const double v = neutronSpeed(E);
            if (v>0) tallyTrackToMesh(*meshTLE, meshTLE->tle_density, n.pos, n.dir, dSurf, 1.0/v);
        }
        n.pos = add3(n.pos, scaleByC(n.dir, dSurf + EPS));
    }
    return fr;
}

__host__ __device__ FlightResult deltaFlight(const Universe& u, Neutron& n, double E, const vector<int>& MTs_total, double SigmaM, bool virtualAllowed) {
    FlightResult fr;
    if (SigmaM<=0.0) { fr.leaked=true; return fr; }
    const double l = -std::log(std::max(1e-32, randomVal()))/SigmaM;
    const array<double,3> newp = add3(n.pos, scaleByC(n.dir, l));

    const Geometry* g=nullptr;
    double SigmaLoc = macroscopicSigmaAt(u, newp, E, MTs_total, &g);
    if (!g) { fr.leaked = true; return fr; }

    const double Pacc = (SigmaLoc>0)? std::min(1.0, SigmaLoc/SigmaM) : 0.0;
    const bool accept = (randomVal() < Pacc);
    n.pos = newp; fr.traveled = l; fr.pos = n.pos; fr.geom = g; fr.SigmaLocal = SigmaLoc;
    if (accept) { fr.collided = true; }
    else { fr.virtualCollision = true; }
    return fr;
}


// --- NEUTRON ---

void walkHistory(const Universe& U, const RunParams& P, int batch, const vector<int>& MTs_total, const vector<int>& MTs_sample, TallyBook& T, 
                 Collisions& col, FourTally& four, std::deque<Neutron>& bankCur, std::deque<Neutron>& bankNext, int& fissionBirthsOut) {
    while (!bankCur.empty()) {
        Neutron n = bankCur.front(); bankCur.pop_front();
        n.dir = iso_dir();
        n.vel = neutronSpeed(n.energy);
        n.collisions = 0; n.time = 0.0; n.reachedTh = false;
        int steps=0;
        while (steps++ < P.maxSteps) {
            FlightResult fr;
            if (P.track==Tracking::Surface) {
                fr = surfaceFlight(U, n, n.energy, MTs_total, P.mesh);
                const double seg_m = fr.traveled * 0.01;
                n.time += seg_m / n.vel;
                if (fr.geom && fr.traveled > 0.0) {
                    scoreTLESegmentPerGeom(T, batch, *fr.geom, n.energy, fr.traveled, MTs_total);
                }
                if (!fr.collided) {
                    if (fr.leaked) { T.ensureBatchAll(batch, (int)T.matNames.size()); ++T.leaks[batch]; if (T.timeHist) T.timeHist->add3(n.time); break; }
                    continue;
                }
                scoreCFEDensityGlobal(T, batch, n.energy, fr.SigmaLocal);
            } else {
                const double SigmaM = T.SigmaM;
                fr = deltaFlight(U, n, n.energy, MTs_total, SigmaM, true);
                const double seg_m = fr.traveled * 0.01;
                n.time += seg_m / n.vel;
                if (fr.leaked) { T.ensureBatchAll(batch, (int)T.matNames.size()); ++T.leaks[batch]; if (T.timeHist) T.timeHist->add3(n.time); break; }
                if (fr.virtualCollision) {
                    scoreCFEDensityGlobal(T, batch, n.energy, SigmaM);
                    if (fr.geom) {
                        int mi = pickMaterialIndexAtE(*fr.geom, n.energy, MTs_total);
                        auto kjashd = fr.geom->mats[mi];
                        if (mi>=0) scoreCFERxPerMaterial(T, batch, mi, fr.geom->mats[mi], n.energy, MTs_total, SigmaM);
                    }
                    continue;
                }
                scoreCFEDensityGlobal(T, batch, n.energy, fr.SigmaLocal);
            }
            RxSample rx;
            if (!fr.geom || !sampleReactionAtE(*fr.geom, n.energy, MTs_total, MTs_sample, rx)) break;

            int mi = T.matIndex(*rx.mat);
            int rj = rxCol(MTs_sample, rx.mt);
            if (mi>=0 && rj>=0) {
                T.ensureBatchAll(batch, (int)T.matNames.size());
                T.statM[batch][mi][rj] += 1;
                scoreCFERxPerMaterial(T, batch, mi, *rx.mat, n.energy, MTs_total, (P.track==Tracking::Surface)? fr.SigmaLocal : T.SigmaM);
            }
            n.collisions++;
            recordCollision(col, n.collisions, n.energy);
            if (P.mesh) {
                const int vidx = vindex(*P.mesh, n.pos[0], n.pos[1], n.pos[2]);
                if (vidx >= 0) P.mesh->meshAnalogColl[vidx] += 1.0;
            }
            if (rx.mt==18) {
                double nu = valueInterp(rx.mat->neutrons, n.energy);
                if (!(nu>0.0)) nu = defaultNuBar(*rx.mat); 
                const int nEmit = sampleMultiplicity(nu);
                fissionBirthsOut += nEmit;
                four.fissionBirthsTotal += nEmit;
                if (isThermalE(n.energy, *rx.mat)) {
                    ++four.absThTotal;
                    if (isFuelSym(rx.mat->sym)) ++four.absThFuel;
                    four.fissionBirthsThermal += nEmit;
                }
                for (int k=0;k<nEmit;++k) {
                    Neutron child;
                    child.pos = n.pos; child.dir = iso_dir();
                    child.energy = std::max(0.0, -1.2895*std::log(std::max(1e-32, randomVal())));
                    child.isSource = false;
                    if (P.src==SourceMode::External) {
                        bankCur.push_back(child);
                    } else {
                        bankNext.push_back(child);
                    }
                }
                if (T.timeHist) T.timeHist->add3(n.time);
                break;
            } else if (rx.mt==2) {
                double A = rx.mat->a>0? (double)rx.mat->a : std::max(1, rx.mat->z);
                double Eout = elasticScatter(n.energy, A, rx.mat->T, EFG);
                if (n.isSource && !n.reachedTh && isThermalE(Eout, *rx.mat)) { n.reachedTh=true; ++four.reachedThermal; }
                n.energy = Eout; n.vel = neutronSpeed(n.energy);
                continue;
            } else if (rx.mt==4) {
                double A = rx.mat->a>0? (double)rx.mat->a : std::max(1, rx.mat->z);
                const double Delta = getDeltaE(*rx.mat);
                double Eout = inelasticEnergyStationary(n.energy, A, Delta);
                if (n.isSource && !n.reachedTh && isThermalE(Eout, *rx.mat)) { n.reachedTh=true; ++four.reachedThermal; }
                n.energy = Eout; n.vel = neutronSpeed(n.energy);
                continue;
            } else if (rx.mt==102) {
                if (n.isSource && !n.reachedTh && rx.mat->sym=="U238" && n.energy>=1e-6 && n.energy<=1e-2)
                    ++four.resAbsBeforeThermal;
                if (isThermalE(n.energy, *rx.mat)) {
                    ++four.absThTotal;
                    if (isFuelSym(rx.mat->sym)) ++four.absThFuel;
                }
                if (T.timeHist) T.timeHist->add3(n.time);
                break;
            } else {
                if (T.timeHist) T.timeHist->add3(n.time);
                break;
            }
        }
    }
}


// --- SIMULATION CUDA ---

typedef struct GPUInputParams {
    int neutrons;
    array<double, 3> sourcePos;
    double sourceE;
    Tracking tracking;
    
} GPUInputParams;


typedef struct GPUOutput {
    // Full simulation details
    // So 

    
    // TallyBook T;
    Collisions col;
    FourTally four;
    int fissionBirthsOut;

} GPUOutputParams;
    // Output has vectors for all variables, such that each 
__global__ void simulateOnGPUExternal(const GPUInputParams& inGPU, GPUOutput& outGPU) {
    // Singledimensional data, given thread ID 
    int th = threadIdx.x * blockDim.x + threadIdx.x;
    // Space to declare helper variables here
    

    // Each one of these does a single batch - Therefore I need data a
    // Prepare bank
    std::deque<Neutron> bank;
    double vel = neutronSpeed(inGPU.energy);
    for (int i=0;i<inGPU.neutrons;++i) {
            Neutron n; n.pos = P.sourcePos; n.dir = iso_dir();
            n.energy = P.sourceE; n.vel = vel; n.isSource = true; 
            bank.push_back(n);
        }
    // Stash fourfactor parameters etc
    simulateNeutrons();
    // Mesh, TLE, 
    



// const Universe& U, const RunParams& P, int batch, const vector<int>& MTs_total, const vector<int>& MTs_sample, TallyBook& T, 
// Collisions& col, FourTally& four, std::deque<Neutron>& bankCur, std::deque<Neutron>& bankNext, int& fissionBirthsOut

}

__global__ void simulateOnGPUCritical(const GPUInputParams& gpuP, GPUOutput& out) {
    // Singledimensional data, given thread ID 
    int th = threadIdx.x * blockDim.x + threadIdx.x;
    // Space to declare helper variables here
    

    // Each one of these does a single batch - Therefore I need data a
    int neutrons = GPUInputParams.neutrons; 
    std::deque<Neutron> bank;
    
    



// const Universe& U, const RunParams& P, int batch, const vector<int>& MTs_total, const vector<int>& MTs_sample, TallyBook& T, 
// Collisions& col, FourTally& four, std::deque<Neutron>& bankCur, std::deque<Neutron>& bankNext, int& fissionBirthsOut

}

RunOutputs runExternalCuda(const Universe& U, const RunParams& P) {
    auto ceil = [](int batches) {return int((batches - 1 + BLOCK) / BLOCK);};
    // Set up parameters, timing
    
    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();

    GPUInputParams inGPUdata;
    GPUOutput outGPUdata;
    // Somehow modify the data in these and alloc correct vector sizes
    {




    }

    float* inGPUptr = NULL;
    float* outGPUptr = NULL;
    CHECK(cudaMalloc((void**)&inGPUptr, sizeof(in)));
    CHECK(cudaMalloc((void**)&outGPUptr, sizeof(outGPUdata)));
    CHECK(cudaMemcpy(inGPUptr, inGPUdata, sizeof(inGPUdata), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(outGPUptr, outGPUdata, sizeof(outGPUdata), cudaMemcpyHostToDevice));

    { // Run kernels
        dim3 dimGrid(ceil(P.batches));
        dim3 dimBlock(BLOCK);
        simulateOnGPUExternal<<<dimGrid, dimBlock>>>()
        CHECK(cudaGetLastError());
    }

    cudaDeviceSynchronize();
    CHECK(cudaMemcpy(outGPUdata, outGPUptr, sizeof(outGPUdata), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));

    auto t1 = clk::now();
    // ANALYZE outGPUdata



    return outGPUdata

}

RunOutputs runCriticalityCuda(const Universe& U, const RunParams& P, int inactive=5) {
    // Parameter Setup
    GPUInputParams gpuP;
    gpuP.




    // Data to GPU
    // Launch kernels

    // Analyze
    // Return
}

// --- SIMULATION ---

RunOutputs runExternal(const Universe& U, const RunParams& P) {
    // Purpose: Observing system, eg shielding effectiveness
    // Timing start
    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();
    
    // Prepare run parameters and variables
    RunOutputs R; R.T.mesh = P.mesh; R.T.useTLE = (P.track==Tracking::Surface);
    R.T.useCFE = true; R.T.deltaMode = (P.track==Tracking::Delta);
    TimeHist lifeHist(1e-6, 200);
    R.T.timeHist = &lifeHist;
    const auto MTs = P.inelastic ? vector<int>{2,4,18,102} : vector<int>{2,18,102};
    R.T.Rcols = (int)MTs.size();
    const vector<int> MTs_total = MTs;
    const vector<int> MTs_sample = MTs;
    if (P.track==Tracking::Delta) R.T.SigmaM = majorSigma(U, P.sourceE, MTs_total);
    R.collisions.resize(P.batches);
    R.fissionChildren.assign(P.batches, 0);
    long long totalHist = 0;

    // Batch simulation
    // Have each batch behave as an independent simulation -> move to GPU with existing walkHistory function
    // Effectively here, move simParams to GPU and have a separate GPU running instance. Additionally, external/criticality can be done inside GPU so single function here
    // Have hard cap on how long a single iteration of walkHistory can be
    for (int b=0;b<P.batches;++b) {
        std::deque<Neutron> bank, next;

        // Initialize neutron bank of set size
        for (int i=0;i<P.historiesPerBatch;++i) {
            Neutron n; n.pos = P.sourcePos; n.dir = iso_dir();
            n.energy = P.sourceE; n.vel = neutronSpeed(n.energy);
            n.isSource = true; R.T.ensureBatch(b, (int)R.T.matNames.size(), (int)MTs.size());
            bank.push_back(n);
        }
        // Store stats
        totalHist += P.historiesPerBatch;
        FourTally four; four.started += P.historiesPerBatch;
        int births=0;
        // Simulate neutrons
        walkHistory(U, P, b, MTs_total, MTs_sample, R.T, R.collisions[b], four, bank, next, births);
        // More stats
        R.fissionChildren[b] = births;
    }
    // Timing end
    auto t1 = clk::now();
    R.perf.elapsed_s = std::chrono::duration<double>(t1 - t0).count();

    // Finalize stats
    R.perf.histories = totalHist;
    writeTimeHistCsv(lifeHist, (P.src==SourceMode::External? "output/timehist_external.csv" : "output/timehist_critical.csv"));
    R.T.timeHist = nullptr;
    return R;
}

RunOutputs runCriticality(const Universe& U, const RunParams& P, int inactive=5) {
    // Purpose: Find keff
    // Timing start
    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();

    // Prepare run parameters and variables
    RunOutputs R; R.T.mesh = P.mesh; R.T.useTLE = (P.track==Tracking::Surface);
    R.T.useCFE = true; R.T.deltaMode = (P.track==Tracking::Delta);
    TimeHist lifeHist(1e-6, 200);
    R.T.timeHist = &lifeHist;
    const auto MTs = P.inelastic ? vector<int>{2,4,18,102} : vector<int>{2,18,102};
    R.T.Rcols = (int)MTs.size();
    const vector<int> MTs_total = MTs;
    const vector<int> MTs_sample = MTs;
    std::deque<Neutron> bank, next;
    // Initialize neutron bank for simulation 
    for (int i=0;i<P.historiesPerBatch;++i) {
        Neutron n; n.pos = P.sourcePos; n.dir = iso_dir();
        n.energy = P.sourceE; n.vel = neutronSpeed(n.energy);
        n.isSource = true; bank.push_back(n);
    }
    double kLast = 1.0;
    long long totalHist = 0;
    for (int b=0;b<P.batches;++b) {
        if (P.track==Tracking::Delta) R.T.SigmaM = majorSigma(U, P.sourceE, MTs_total);
        FourTally four; four.started += (int)bank.size();
        R.collisions.emplace_back(Collisions{});
        int births=0;
        RunParams Pcycle = P; Pcycle.src = SourceMode::Criticality;
        const int Nk = (int)bank.size();
        walkHistory(U, Pcycle, b, MTs_total, MTs_sample, R.T, R.collisions.back(), four, bank, next, births);
        const int Nk1 = (int)next.size();
        totalHist += P.historiesPerBatch;
        double kcur = (Nk>0) ? double(Nk1)/double(Nk) : 0.0;
        R.keff_history.push_back(kcur);
        std::deque<Neutron> newBank;
        if (!next.empty()) {
            double scaleByC = (kLast>0)? 1.0/kLast : 1.0;
            int need = P.historiesPerBatch;
            while (need-- > 0) {
                size_t j = (size_t)std::floor(randomVal()*next.size());
                if (j>=next.size()) j = next.size()-1;
                newBank.push_back(next[j]);
            }
        }
        bank.swap(newBank);
        next.clear();
        kLast = (kcur>0? kcur : kLast);
    }
    auto t1 = clk::now();
    R.perf.elapsed_s = std::chrono::duration<double>(t1 - t0).count();
    R.perf.histories = totalHist;
    writeTimeHistCsv(lifeHist, (P.src==SourceMode::External? "output/timehist_external.csv" : "output/timehist_critical.csv"));
    R.T.timeHist = nullptr;

    return R;
}
