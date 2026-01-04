#include "mc.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef __CUDACC__
  #include <cuda_runtime.h>
  #define CUDA_CHECK(x) cuda_check((x), #x, __FILE__, __LINE__)
static inline void cuda_check(cudaError_t err, const char* expr, const char* file, int line) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: " << cudaGetErrorString(err)
                  << " | " << expr << " | " << file << ":" << line << "\n";
        std::abort();
    }
}
#endif

using std::array;
using std::vector;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

// --- CUDA ERROR HANDLING ---




// --- SHARED HELPERS ---

array<float,3> scaleByC(const array<float,3>& a, float s) {
    return {a[0]*s, a[1]*s, a[2]*s};
}

vector<float> logspace(float ea,float eb,int n) {
    vector<float> out;
    if (n<=0) return out;
    out.reserve((size_t)n);
    const float de = (n==1) ? 0.0f : (eb-ea)/(n-1);
    for (int i=0;i<n;++i) out.push_back(std::pow(10.0f, ea + de*i));
    return out;
}

// E in MeV -> v in m/s
__host__ __device__ float neutronSpeed(float E) {
    const float eJ = E * 1.602176634e-13;
    const float mn = 1.67492749804e-27;
    if (eJ <= 0.0f) return 0.0f;
    return std::sqrt(2.0 * eJ / mn);
}

float valueInterp(const vector<std::pair<float,float>>& data, float target) {
    if (data.empty()) return 0.0f;
    if (target <= data.front().first) return data.front().second;
    if (target >= data.back().first)  return data.back().second;

    // Binary search 
    size_t lo = 0, hi = data.size()-1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) >> 1;
        if (data[mid].first <= target) lo = mid;
        else hi = mid;
    }
    const float x0 = data[lo].first;
    const float x1 = data[hi].first;
    const float y0 = data[lo].second;
    const float y1 = data[hi].second;
    if (x1 == x0) return y0;
    const float t = (target - x0) / (x1 - x0);
    return y0 + t*(y1 - y0);
}


// --- CPU MISC HELPERS ---

// Faster random number generation than previous mt19937_64
uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

float randomVal() {
    thread_local uint64_t s = 0x1234abcd5678ef01ULL;
    uint64_t r = splitmix64(s);
    return (r >> 11) * (1.0f / 9007199254740992.0);
}

array<float,3> iso_dir() {
    const float u = 2.0*randomVal() - 1.0f;
    const float phi = 2.0*M_PI*randomVal();
    const float s = std::sqrt(std::max(0.0f, 1.0f - u*u));
    return {s*std::cos(phi), s*std::sin(phi), u};
}

int sampleMultiplicity(float nuBar) {
    if (nuBar <= 0.0f) return 0;
    int n = (int)std::floor(nuBar);
    const float frac = nuBar - n;
    if (randomVal() < frac) ++n;
    return std::max(0, n);
}

float sampleFissionEnergyMeV() {
    const float r = std::max(1e-14, (double)randomVal());
    return -std::log(r) * 1.0f; // MeV
}

bool isFuelSym(const std::string& sym) {
    std::string s;
    s.reserve(sym.size());
    for (char c : sym) {
        if (std::isalnum((unsigned char)c)) s.push_back((char)std::toupper((unsigned char)c));
    }
    return (s=="U235" || s=="U238" || s=="PU239" || s=="PU240" || s=="PU241");
}

bool isThermalE(float E, const Material& mat) {
    constexpr float kB = 1.380649e-23;
    const float T_K = (float)(mat.nuc ? mat.nuc->T : 0) + 273.15;
    const float eth_J = 2.0*kB*T_K;
    const float eth_MeV = eth_J / 1.602176634e-13;
    return E < eth_MeV;
}

float elasticEnergyStationary(float Ein, float A) {
    const float alpha = std::pow((A-1.0f)/(A+1.0f), 2.0);
    const float r = randomVal();
    return Ein * (alpha + (1.0f-alpha)*r);
}

float inelasticEnergyStationary(float Ein, float dE) {
    return std::max(0.0f, Ein - std::max(0.0f, dE));
}

float getDeltaE(const Material& m) {
    if (!m.nuc) return 1.0e-6;
    auto it = m.nuc->Qvals.find(4);
    if (it != m.nuc->Qvals.end()) return std::abs(it->second);
    return 1.0e-6;
}


// --- CPU CROSS SECTION HELPERS ---

float interpMT(const std::map<int, vector<std::pair<float,float>>>& mt, int MT, float E) {
    auto it = mt.find(MT);
    if (it == mt.end()) return 0.0f;
    return valueInterp(it->second, E);
}

float materialWeight(const Material& m, float E, const vector<int>& mts_total) {
    if (!m.nuc) return 0.0f;
    float sum = 0.0f;
    for (int mt : mts_total) sum += interpMT(m.nuc->mt, mt, E);
    return m.proportion * m.rho * sum;
}

float sigmaTot(const MatView& mats, float E, const vector<int>& mts_total) {
    float s = 0.0f;
    for (const auto& m : mats) s += materialWeight(m, E, mts_total);
    return s;
}

int pickMaterialIndexAtE(const Geometry& g, float E, const vector<int>& mts_total) {
    if (g.mats.empty()) return -1;
    float tot = 0.0f;
    for (const auto& m : g.mats) tot += materialWeight(m, E, mts_total);
    if (tot <= 0.0f) return 0;
    const float r = randomVal() * tot;
    float c = 0.0f;
    for (int i=0; i<(int)g.mats.size(); ++i) {
        c += materialWeight(g.mats[i], E, mts_total);
        if (r <= c) return i;
    }
    return (int)g.mats.size()-1;
}

// Reaction sampling among the enabled MTs.
struct RxEvent {
    int mt = 0;
    int mi = -1;
    const Material* mat = nullptr;
};

RxEvent sampleReactionAtE(const Geometry& g, float E, const vector<int>& mts_total, const vector<int>& mts_sample) {
    RxEvent rx{};
    if (g.mats.empty()) return rx;

    int mi = pickMaterialIndexAtE(g, E, mts_total);
    mi = std::max(0, std::min(mi, (int)g.mats.size()-1));
    const Material& m = g.mats[mi];

    float tot = 0.0f;
    float w[8];
    int   mtv[8];
    const int n = (int)std::min<size_t>(mts_sample.size(), 8);
    for (int i=0;i<n;++i){
        mtv[i] = mts_sample[i];
        w[i] = std::max(0.0f, interpMT(m.nuc->mt, mtv[i], E));
        tot += w[i];
    }
    if (tot<=0.0f) {
        rx.mt = mtv[0]; rx.mi = mi; rx.mat = &m;
        return rx;
    }
    float r = randomVal()*tot;
    float c = 0.0f;
    for (int i=0;i<n;++i){
        c += w[i];
        if (r <= c) { rx.mt = mtv[i]; break; }
    }
    rx.mi = mi;
    rx.mat = &m;
    return rx;
}

float macroscopicSigmaAt(const Universe& u, const array<float,3>& pos, float E, const vector<int>& mts_total, const Geometry** gOut) {
    const Geometry* g = findGeometryAt(u, pos);
    if (!g) {
        if (gOut) *gOut = nullptr;
        return -1.0f;
    }
    if (gOut) *gOut = g;
    return sigmaTot(g->mats, E, mts_total);
}

float majorSigma(const Universe& u, float E, const vector<int>& mts_total) {
    float m = 0.0f;
    for (const auto& g : u.geometries) m = std::max(m, sigmaTot(g.mats, E, mts_total));
    for (const auto& su : u.subUniverse) m = std::max(m, majorSigma(su, E, mts_total));
    return m;
}


// --- CPU SIMULATION ---

struct FlightResult {
    bool collided = false;
    bool leaked   = false;
    float traveled = 0.0f;
    float SigmaLocal = 0.0f;
    const Geometry* geom = nullptr;
};

int meshIndexHost(const Mesh3D& M, const array<float,3>& p) {
    return vindex(M, p[0], p[1], p[2]);
}

FlightResult surfaceFlight(const Universe& u, const Neutron& n, float E, const vector<int>& mts_total) {
    FlightResult fr{};

    const Geometry* g0 = nullptr;
    float Sigma = macroscopicSigmaAt(u, n.pos, E, mts_total, &g0);
    if (Sigma < 0.0f) { fr.leaked=true; return fr; }
    if (Sigma <= 0.0f) Sigma = 0.0f;

    const float r = std::max(1e-14, (double)randomVal());
    const float l = (Sigma>0.0f) ? (-std::log(r)/Sigma) : 1e300;

    const float dSurf = nearestCollision(n, u);
    if (dSurf < 0.0f) { fr.leaked = true; return fr; }

    fr.geom = g0;

    if (l < dSurf) {
        fr.collided = true;
        fr.traveled = l;
        fr.SigmaLocal = Sigma;
    } else {
        fr.collided = false;
        fr.traveled = dSurf;
        fr.SigmaLocal = Sigma;
    }
    return fr;
}

FlightResult deltaFlight(const Universe& u, const Neutron& n, float E, const vector<int>& mts_total, float SigmaM) {
    FlightResult fr{};

    const float r = std::max(1e-14, (double)randomVal());
    const float l = -std::log(r) / std::max(1e-30, (double)SigmaM);

    array<float,3> dHat = normed(n.dir);
    array<float,3> newp = madd(n.pos, dHat, l);

    const Geometry* g0 = nullptr;
    float SigmaLoc = macroscopicSigmaAt(u, newp, E, mts_total, &g0);
    if (SigmaLoc < 0.0f || !g0) {
        fr.leaked = true;
        fr.traveled = l;
        return fr;
    }

    fr.geom = g0;
    fr.traveled = l;
    fr.SigmaLocal = SigmaLoc;

    const float accept = std::min(1.0, (double)SigmaLoc / std::max(1e-30, (double)SigmaM));
    fr.collided = (randomVal() < accept);
    return fr;
}

void ensureMaterialIndex(TallyBook& T, const Material& m) {
    (void)T.matIndex(m);
}

void walkHistoryCPU(const Universe& u, int batch, int maxSteps,
                           std::deque<Neutron>& bankCur, std::vector<Neutron>& bankNext,
                           RunOutputs& R, const RunParams& P,
                           const vector<int>& MTs_total, const vector<int>& MTs_sample,
                           int* fissionChildrenOut) {

    while (!bankCur.empty()) {
        Neutron n = bankCur.front();
        bankCur.pop_front();

        n.dir = iso_dir();
        n.collisions = 0;
        n.reachedTh = false;
        n.time = 0.0f;
        n.vel = neutronSpeed(n.energy);

        for (int step=0; step<maxSteps; ++step) {
            FlightResult fr;

            if (P.track == Tracking::Surface) {
                fr = surfaceFlight(u, n, n.energy, MTs_total);
                if (fr.leaked) { R.T.leaks[batch]++; break; }

                if (P.mesh && fr.traveled > 0.0f) {
                    const array<float,3> dHat = normed(n.dir);
                    tallyTrackToMesh(*P.mesh, P.mesh->tle_density, n.pos, dHat, fr.traveled, (n.vel>0.0f)?(1.0f/n.vel):0.0f);
                }
                if (fr.geom && fr.traveled > 0.0f) {
                    scoreTLESegmentPerGeom(R.T, batch, *fr.geom, n.energy, fr.traveled, MTs_total);
                }

                n.pos = madd(n.pos, normed(n.dir), fr.traveled);

                if (!fr.collided) {
                    n.pos = madd(n.pos, normed(n.dir), 1e-6);
                    continue;
                }

                if (P.mesh) {
                    const int vid = meshIndexHost(*P.mesh, n.pos);
                    if (vid >= 0) P.mesh->meshAnalogColl[vid] += 1.0f;
                }

                const float SigmaRef = fr.SigmaLocal;
                scoreCFEDensityGlobal(R.T, batch, n.energy, SigmaRef);
                if (P.mesh) scoreCFE(R.T, batch, n.pos, n.energy, SigmaRef);

                RxEvent rx = sampleReactionAtE(*fr.geom, n.energy, MTs_total, MTs_sample);
                if (!rx.mat) break;

                ensureMaterialIndex(R.T, *rx.mat);
                const int mid = R.T.matIndex(*rx.mat);
                R.T.ensureBatchAll(batch, (int)R.T.matNames.size());

                int col = -1;
                for (int j=0;j<(int)MTs_total.size();++j) if (MTs_total[j]==rx.mt) { col=j; break; }
                if (col >= 0 && mid >= 0) R.T.statM[batch][mid][col] += 1;

                scoreCFERxPerMaterial(R.T, batch, mid, *rx.mat, n.energy, MTs_total, SigmaRef);

                if (isThermalE(n.energy, *rx.mat)) n.reachedTh = true;
                n.collisions += 1;

                if (rx.mt == 102) {
                    // capture
                    break;
                } else if (rx.mt == 18) {
                    // fission
                    const float nuBar = valueInterp(rx.mat->nuc->neutrons, n.energy);
                    int nEmit = sampleMultiplicity(nuBar>0?nuBar:2.43);
                    if (fissionChildrenOut) (*fissionChildrenOut) += nEmit;
                    for (int k=0;k<nEmit;++k) {
                        Neutron child;
                        child.pos = n.pos;
                        child.energy = sampleFissionEnergyMeV();
                        child.isSource = false;
                        if (P.src == SourceMode::External) bankCur.push_back(child);
                        else bankNext.push_back(child);
                    }
                    break;
                } else if (rx.mt == 2) {
                    // elastic
                    const float A = (rx.mat->nuc->aw > 0.0f) ? rx.mat->nuc->aw : (float)rx.mat->nuc->a;
                    n.energy = elasticEnergyStationary(n.energy, std::max(1.0f, A));
                } else if (rx.mt == 4) {
                    // inelastic
                    n.energy = inelasticEnergyStationary(n.energy, getDeltaE(*rx.mat));
                }
            } else {
                // Delta tracking
                const float SigmaM = R.T.SigmaM;
                fr = deltaFlight(u, n, n.energy, MTs_total, SigmaM);

                if (P.mesh && fr.traveled > 0.0f) {
                    const array<float,3> dHat = normed(n.dir);
                    tallyTrackToMesh(*P.mesh, P.mesh->tle_density, n.pos, dHat, fr.traveled, (n.vel>0.0f)?(1.0f/n.vel):0.0f);
                }

                n.pos = madd(n.pos, normed(n.dir), fr.traveled);

                if (fr.leaked) { R.T.leaks[batch]++; break; }

                scoreCFEDensityGlobal(R.T, batch, n.energy, SigmaM);
                if (P.mesh) scoreCFE(R.T, batch, n.pos, n.energy, SigmaM);

                if (fr.geom) {
                    int mi = pickMaterialIndexAtE(*fr.geom, n.energy, MTs_total);
                    mi = std::max(0, std::min(mi, (int)fr.geom->mats.size()-1));
                    const int mid = R.T.matIndex(fr.geom->mats[mi]);
                    R.T.ensureBatchAll(batch, (int)R.T.matNames.size());
                    scoreCFERxPerMaterial(R.T, batch, mid, fr.geom->mats[mi], n.energy, MTs_total, SigmaM);
                }

                if (!fr.collided) continue;

                if (P.mesh) {
                    const int vid = meshIndexHost(*P.mesh, n.pos);
                    if (vid >= 0) P.mesh->meshAnalogColl[vid] += 1.0f;
                }

                RxEvent rx = sampleReactionAtE(*fr.geom, n.energy, MTs_total, MTs_sample);
                if (!rx.mat) break;

                if (isThermalE(n.energy, *rx.mat)) n.reachedTh = true;
                n.collisions += 1;

                ensureMaterialIndex(R.T, *rx.mat);
                const int mid = R.T.matIndex(*rx.mat);
                R.T.ensureBatchAll(batch, (int)R.T.matNames.size());
                int col = -1;
                for (int j=0;j<(int)MTs_total.size();++j) if (MTs_total[j]==rx.mt) { col=j; break; }
                if (col >= 0 && mid >= 0) R.T.statM[batch][mid][col] += 1;

                if (rx.mt == 102) {
                    // capture
                    break;
                } else if (rx.mt == 18) {
                    // fission
                    const float nuBar = valueInterp(rx.mat->nuc->neutrons, n.energy);
                    int nEmit = sampleMultiplicity(nuBar>0?nuBar:2.43);
                    if (fissionChildrenOut) (*fissionChildrenOut) += nEmit;
                    for (int k=0;k<nEmit;++k) {
                        Neutron child;
                        child.pos = n.pos;
                        child.energy = sampleFissionEnergyMeV();
                        child.isSource = false;
                        if (P.src == SourceMode::External) bankCur.push_back(child);
                        else bankNext.push_back(child);
                    }
                    break;
                } else if (rx.mt == 2) {
                    // elastic
                    const float A = (rx.mat->nuc->aw > 0.0f) ? rx.mat->nuc->aw : (float)rx.mat->nuc->a;
                    n.energy = elasticEnergyStationary(n.energy, std::max(1.0f, A));
                } else if (rx.mt == 4) {
                    // inelastic
                    n.energy = inelasticEnergyStationary(n.energy, getDeltaE(*rx.mat));
                }
            }
        }
    }
}

RunOutputs runExternal(const Universe& U, const RunParams& P) {
    RunOutputs R;
    const vector<int> MTs_total  = P.inelastic ? vector<int>{2,4,18,102} : vector<int>{2,18,102};
    const vector<int> MTs_sample = MTs_total;

    const int Rcols = P.inelastic ? 4 : 3;

    if (P.mesh) P.mesh->zero();

    R.T.mesh = P.mesh;
    R.T.Rcols = Rcols;
    // R.T.nBatches = P.batches;
    // R.T.inelastic = P.inelastic;
    R.T.deltaMode = (P.track==Tracking::Delta);
    R.T.SigmaM = (P.track==Tracking::Delta) ? majorSigma(U, P.sourceE, MTs_total) : 0.0f;
    R.T.leaks.assign(P.batches, 0);
    R.T.cfe_global_time.assign(P.batches, 0.0f);

    auto t0 = std::chrono::steady_clock::now();

    int totalHist = 0;
    R.fissionChildren.assign(P.batches, 0);

    for (int b=0;b<P.batches;++b) {
        std::deque<Neutron> bankCur;
        std::vector<Neutron> bankNext;
        bankNext.reserve((size_t)P.historiesPerBatch*4);

        int fissCount = 0;

        for (int i=0;i<P.historiesPerBatch;++i) {
            Neutron n;
            n.pos = P.sourcePos;
            n.energy = P.sourceE;
            n.isSource = true;
            bankCur.push_back(n);
        }

        walkHistoryCPU(U, b, P.maxSteps, bankCur, bankNext, R, P, MTs_total, MTs_sample, &fissCount);
        totalHist += P.historiesPerBatch;

        R.fissionChildren[b] = fissCount;
    }

    auto t1 = std::chrono::steady_clock::now();
    R.perf.histories = totalHist;
    R.perf.elapsed_s = std::chrono::duration<float>(t1 - t0).count();

    return R;
}

RunOutputs runCriticality(const Universe& U, const RunParams& P, int inactive) {
    RunOutputs R;
    const vector<int> MTs_total  = P.inelastic ? vector<int>{2,4,18,102} : vector<int>{2,18,102};
    const vector<int> MTs_sample = MTs_total;
    const int Rcols = P.inelastic ? 4 : 3;

    if (P.mesh) P.mesh->zero();

    R.T.mesh = P.mesh;
    R.T.Rcols = Rcols;
    // R.T.nBatches = P.batches;
    // R.T.inelastic = P.inelastic;
    R.T.deltaMode = (P.track==Tracking::Delta);
    R.T.SigmaM = (P.track==Tracking::Delta) ? majorSigma(U, P.sourceE, MTs_total) : 0.0f;
    R.T.leaks.assign(P.batches, 0);
    R.T.cfe_global_time.assign(P.batches, 0.0f);

    R.fissionChildren.assign(P.batches, 0);

    std::vector<Neutron> bankCur;
    bankCur.reserve(P.historiesPerBatch);
    for (int i=0;i<P.historiesPerBatch;++i) {
        Neutron n;
        n.pos = P.sourcePos;
        n.energy = P.sourceE;
        n.isSource = true;
        bankCur.push_back(n);
    }

    auto t0 = std::chrono::steady_clock::now();

    for (int b=0;b<P.batches;++b) {
        std::deque<Neutron> q;
        for (auto& n : bankCur) q.push_back(n);

        std::vector<Neutron> bankNext;
        bankNext.reserve((size_t)P.historiesPerBatch*4);

        walkHistoryCPU(U, b, P.maxSteps, q, bankNext, R, P, MTs_total, MTs_sample, nullptr);

        const int Nk  = (int)bankCur.size();
        const int Nk1 = (int)bankNext.size();
        R.fissionChildren[b] = Nk1;
        const float kcur = (Nk>0) ? (float)Nk1/(float)Nk : 0.0f;
        R.keff_history.push_back(kcur);

        // Population control
        bankCur.clear();
        bankCur.reserve(P.historiesPerBatch);
        if (Nk1 <= 0) {
            for (int i=0;i<P.historiesPerBatch;++i) {
                Neutron n; n.pos = P.sourcePos; n.energy = P.sourceE; n.isSource = true;
                bankCur.push_back(n);
            }
        } else {
            for (int i=0;i<P.historiesPerBatch;++i) {
                const int j = (int)std::floor(randomVal() * Nk1);
                bankCur.push_back(bankNext[std::max(0, std::min(j, Nk1-1))]);
            }
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    R.perf.histories = P.historiesPerBatch * P.batches;
    R.perf.elapsed_s = std::chrono::duration<float>(t1 - t0).count();

    (void)inactive;
    return R;
}


// --- GPU ---

#ifdef __CUDACC__


// Device compatible vector
struct Vec3d { float x,y,z; };
__device__ Vec3d v3_make(float x,float y,float z){ return {x,y,z}; }
__device__ Vec3d v3_add(Vec3d a, Vec3d b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
__device__ Vec3d v3_sub(Vec3d a, Vec3d b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
__device__ Vec3d v3_mul(Vec3d a, float s){ return {a.x*s,a.y*s,a.z*s}; }
__device__ float v3_dot(Vec3d a, Vec3d b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
__device__ float v3_len(Vec3d a){ return sqrt(v3_dot(a,a)); }
__device__ Vec3d v3_norm(Vec3d a){ float L=v3_len(a); return (L>0)?v3_mul(a,1.0f/L):v3_make(0,0,0); }

// Threading compatible RNG (Fast)
struct RNG {
    unsigned long long s;
    __device__ unsigned long long nextU64() {
        unsigned long long x = s;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        s = x;
        return x * 2685821657736338717ULL;
    }
    __device__ float uni() {
        return ((nextU64() >> 11) * (1.0f / 9007199254740992.0));
    }
};

__device__ Vec3d isoDirGPU(RNG& rng) {
    const float u = 2.0*rng.uni() - 1.0f;
    const float phi = 2.0*M_PI*rng.uni();
    const float s = sqrt(fmax(0.0f, 1.0f - u*u));
    return v3_make(s*cos(phi), s*sin(phi), u);
}

__device__ int sampleMultiplicityGPU(RNG& rng, float nuBar) {
    if (nuBar <= 0.0f) return 0;
    int n = (int)floor(nuBar);
    float frac = nuBar - (float)n;
    if (rng.uni() < frac) ++n;
    return n < 0 ? 0 : n;
}

__device__ float sampleFissionEnergyMeVGPU(RNG& rng) {
    const float r = fmax(1e-14, (double)rng.uni());
    return -log(r) * 1.0f; // MeV
}


// --- GPU DATA ---
// Notably, one cannot use C++ std libraries in device code so new structs are required

struct GPUShape {
    float A,B,C,D,E,F,G,H,I,J;
    int torus;
};

struct GPURPNNode { int op; int shape; };

struct GPUGeometry {
    int shapeOff, nShapes;
    int nodeOff, nNodes;
    int mixId; // -1 for bounding-only geometry
};

struct GPUUniverse {
    Vec3d pos;
    int boundGeom;
    int childOff, nChild;
    int geomOff,  nGeom;
};

struct GPUMix { int compOff, nComp; };
struct GPUMixComp { int nuc; float rho; float prop; };

struct GPUNuclide {
    float A;      // mass ratio
    float T_K;    // temperature in Kelvin
    int isFuel;
};

// Memory layout:
//   xs[(nuc*4 + rx)*NE + i]
//   nuBar[nuc*NE + i]
// rx indices: 0=elastic(2), 1=inel(4), 2=fission(18), 3=capture(102)
struct GPUModel {
    const GPUUniverse* universes;
    const int* uChildren;
    const int* uGeoms;
    const GPUGeometry* geoms;
    const GPUShape* shapes;
    const GPURPNNode* nodes;
    const GPUMix* mixes;
    const GPUMixComp* comps;
    const GPUNuclide* nuclides;

    const float* xs;
    const float* nuBar;

    int NU;
    int NG;
    int NMix;
    int NComp;
    int NNuc;

    const float* Egrid;
    int NE;
    float logEmin;
    float invDlog;
    int rootUid;
};

struct GPURunParams {
    int historiesPerBatch;
    int batches;
    int maxSteps;
    int inelastic;
    int track; // 0=Surface, 1=Delta
    int src; // 0=External, 1=Criticality

    Vec3d sourcePos;
    float sourceE;
    float SigmaM;

    // Mesh
    int nx, ny, nz;
    Vec3d pmin;
    Vec3d pmax;
    Vec3d invH; // 1/dx, 1/dy, 1/dz
};

// Outputs
struct GPUOutputs {
    // Per batch
    int* leaks; // [batches]
    int* fissionChildren; // [batches]

    int* statM; // [batches * NNuc * Rcols]
    float* cfe_Rtot; // [batches * NNuc]
    float* cfe_Rabs; // -- || --
    float* tle_Rtot; // -- || --
    float* tle_Rabs; // -- || --
    float* cfe_global_time; // [batches]

    // Mesh fields (global)
    float* meshAnalogColl; // [nx*ny*nz]
    float* cfe_density;    // -- || --
    float* tle_density;    // -- || --

    int Rcols;
};

// Bank (criticality)
struct BankNeutron {
    Vec3d pos;
    float E;
    int isSource;
};


// --- GPU CROSS SECTION ---

__device__ float clampd(float x,float a,float b){ return x<a?a:(x>b?b:x); }

__device__ int clampi(int x,int a,int b){ return x<a?a:(x>b?b:x); }

__device__ float xsInterp(const GPUModel& M, int nuc, int rx, float E) {
    // rx: 0..3
    const float Emin = M.Egrid[0];
    const float Emax = M.Egrid[M.NE-1];
    float Ec = clampd(E, Emin, Emax);
    float u = (log10(Ec) - M.logEmin) * M.invDlog;
    int i = (int)floor(u);
    i = clampi(i, 0, M.NE-2);
    float t = u - (float)i;
    if (t < 0) t = 0; if (t>1) t = 1;
    const int base = (nuc*4 + rx)*M.NE;
    float y0 = M.xs[base + i];
    float y1 = M.xs[base + i + 1];
    return y0 + t*(y1 - y0);
}

__device__ float nuBarInterp(const GPUModel& M, int nuc, float E) {
    const float Emin = M.Egrid[0];
    const float Emax = M.Egrid[M.NE-1];
    float Ec = clampd(E, Emin, Emax);
    float u = (log10(Ec) - M.logEmin) * M.invDlog;
    int i = (int)floor(u);
    i = clampi(i, 0, M.NE-2);
    float t = u - (float)i;
    if (t < 0) t = 0; if (t>1) t = 1;
    const int base = nuc * M.NE;
    float y0 = M.nuBar[base + i];
    float y1 = M.nuBar[base + i + 1];
    return y0 + t*(y1 - y0);
}

// --- GPU GEOMETRY ---

__device__ bool insideLeafQuadric(const GPUShape& s, Vec3d p) {
    const float x=p.x,y=p.y,z=p.z;
    float f = s.A*x*x + s.B*y*y + s.C*z*z
             + s.D*x*y + s.E*y*z + s.F*x*z
             + s.G*x + s.H*y + s.I*z + s.J;
    return f <= 0.0f;
}

// No torus support
__device__ bool insideLeafGPU(const GPUShape& s, Vec3d p) {
    if (s.torus) return false;
    return insideLeafQuadric(s, p);
}

__device__ bool pointInGeomGPU(const GPUModel& M, int geomId, Vec3d p) {
    const GPUGeometry& g = M.geoms[geomId];
    bool stack[64];
    int sp = 0;
    const int nN = g.nNodes;
    for (int i=0;i<nN;++i) {
        const GPURPNNode nd = M.nodes[g.nodeOff + i];
        if (nd.op == 0) {
            // leaf
            const GPUShape& s = M.shapes[g.shapeOff + nd.shape];
            stack[sp++] = insideLeafGPU(s, p);
        } else if (nd.op == 1) {
            // NOT
            if (sp>0) stack[sp-1] = !stack[sp-1];
        } else if (nd.op == 2) {
            // OR
            if (sp>=2) {
                bool b = stack[--sp];
                bool a = stack[--sp];
                stack[sp++] = (a||b);
            }
        } else {
            // AND
            if (sp>=2) {
                bool b = stack[--sp];
                bool a = stack[--sp];
                stack[sp++] = (a&&b);
            }
        }
    }
    return (sp>0) ? stack[sp-1] : false;
}

__device__ int solveQuadraticGPU(float A,float B,float C,float& r0,float& r1) {
    // A t^2 + B t + C = 0, Returns number of roots, sorted
    const float EPS=1e-14;
    if (fabs(A) < EPS) {
        if (fabs(B) < EPS) return 0;
        r0 = -C/B;
        return 1;
    }
    float disc = B*B - 4*A*C;
    if (disc < 0.0f) return 0;
    float sdisc = sqrt(fmax(0.0f, disc));
    float q = -0.5*(B + (B>=0 ? sdisc : -sdisc));
    float t0 = q/A;
    float t1 = C/q;
    if (t0>t1) { float tmp=t0; t0=t1; t1=tmp; }
    r0=t0; r1=t1;
    return (disc==0.0f)?1:2;
}

__device__ int solveCubicMonicGPU(float a, float b, float c, float r[3]) {
    // Solve x^3 + a x^2 + b x + c = 0 (real roots only). Returns 1 or 3 roots.
    const float EPS = 1e-14;
    const float a2 = a*a;
    const float q = (a2 - 3.0*b)/9.0;
    const float r3 = (2.0*a2*a - 9.0*a*b + 27.0*c)/54.0;
    const float q3 = q*q*q;
    const float disc = r3*r3 - q3;

    if (disc > EPS) {
        float s = cbrt(r3 + sqrt(disc));
        float t = cbrt(r3 - sqrt(disc));
        r[0] = -a/3.0 + (s+t);
        return 1;
    }

    if (q < EPS) {
        r[0] = -a/3.0;
        r[1] = r[0];
        r[2] = r[0];
        return 3;
    }

    float cosarg = r3 / sqrt(fmax(EPS, q3));
    cosarg = fmax(-1.0f, fmin(1.0f, cosarg));
    float th = acos(cosarg);
    float sq = 2.0*sqrt(q);
    r[0] = -a/3.0 + sq*cos(th/3.0);
    r[1] = -a/3.0 + sq*cos((th + 2.0*M_PI)/3.0);
    r[2] = -a/3.0 + sq*cos((th + 4.0*M_PI)/3.0);

    if (r[0] > r[1]) { float t=r[0]; r[0]=r[1]; r[1]=t; }
    if (r[1] > r[2]) { float t=r[1]; r[1]=r[2]; r[2]=t; }
    if (r[0] > r[1]) { float t=r[0]; r[0]=r[1]; r[1]=t; }
    return 3;
}

__device__ int solveQuarticGPU(float a, float b, float c, float d, float e, float r[4]) {
    // Solve a x^4 + b x^3 + c x^2 + d x + e = 0 (real roots only)
    const float EPS = 1e-14;
    if (fabs(a) < EPS) {
        // degrade: cubic not implemented here
        return 0;
    }
    b /= a; c /= a; d /= a; e /= a;

    float bb = b*b;
    float p = c - 3.0*bb/8.0;
    float q = b*bb/8.0 - b*c/2.0 + d;
    float r0 = -3.0*bb*bb/256.0 + bb*c/16.0 - b*d/4.0 + e;

    float zroots[3];
    int nz = solveCubicMonicGPU(2.0*p, (p*p - 4.0*r0), (-q*q), zroots);
    float z = (nz==1) ? zroots[0] : zroots[nz-1];
    if (z < 0.0f) z = 0.0f;
    float alpha = sqrt(z);

    float m = p + z;
    float beta, gamma;
    if (fabs(alpha) < EPS) {
        float t0,t1;
        int n = solveQuadraticGPU(1.0f, m, r0, t0, t1);
        beta = (n>=1)?t0:0.0f;
        gamma = (n==2)?t1:beta;
    } else {
        beta = (m - q/alpha)/2.0;
        gamma = (m + q/alpha)/2.0;
    }

    int nr=0;
    float y0,y1;
    int n1 = solveQuadraticGPU(1.0f, alpha, beta, y0, y1);
    if (n1 >= 1) r[nr++] = y0;
    if (n1 == 2) r[nr++] = y1;
    float y2,y3;
    int n2 = solveQuadraticGPU(1.0f, -alpha, gamma, y2, y3);
    if (n2 >= 1) r[nr++] = y2;
    if (n2 == 2) r[nr++] = y3;

    float shift = b/4.0;
    for (int i=0;i<nr;++i) r[i] -= shift;

    for (int i=0;i<nr;++i) {
        for (int j=i+1;j<nr;++j) {
            if (r[j] < r[i]) { float tmp=r[i]; r[i]=r[j]; r[j]=tmp; }
        }
    }
    return nr;
}

__device__ int surfaceDistQuadricGPU(const GPUShape& s, Vec3d P, Vec3d D, float out[2]) {
    const float x=P.x,y=P.y,z=P.z;
    const float u=D.x,v=D.y,w=D.z;

    const float K = s.A*x*x + s.B*y*y + s.C*z*z
                   + s.D*x*y + s.E*y*z + s.F*x*z
                   + s.G*x + s.H*y + s.I*z + s.J;

    const float L = 2.0*(s.A*u*x + s.B*v*y + s.C*w*z)
                   + s.D*(v*x + u*y) + s.E*(w*y + v*z) + s.F*(w*x + u*z)
                   + s.G*u + s.H*v + s.I*w;

    const float Mq = s.A*u*u + s.B*v*v + s.C*w*w
                    + s.D*u*v + s.E*v*w + s.F*u*w;

    float r0=0.0f,r1=0.0f;
    int n = solveQuadraticGPU(Mq,L,K,r0,r1);
    int k=0;
    const float EPS = 1e-10;
    if (n>=1 && r0>EPS) out[k++]=r0;
    if (n==2 && r1>EPS) out[k++]=r1;
    if (k==2 && out[0]>out[1]) { float tmp=out[0]; out[0]=out[1]; out[1]=tmp; }
    return k;
}

__device__ void sortSmall(float* a, int n) {
    for (int i=1;i<n;++i) {
        float key=a[i];
        int j=i-1;
        while (j>=0 && a[j]>key) { a[j+1]=a[j]; --j; }
        a[j+1]=key;
    }
}

__device__ float geometryCollisionGPU(const GPUModel& M, int geomId, Vec3d pos, Vec3d dir) {
    const GPUGeometry& g = M.geoms[geomId];
    Vec3d dHat = v3_norm(dir);
    if (dHat.x==0 && dHat.y==0 && dHat.z==0) return -1.0f;

    float cand[128];
    int nCand=0;
    for (int si=0; si<g.nShapes; ++si) {
        const GPUShape& s = M.shapes[g.shapeOff + si];
        if (s.torus) {
            return -1.0f;
        } else {
            float roots[2];
            int nr = surfaceDistQuadricGPU(s, pos, dHat, roots);
            for (int k=0;k<nr && nCand<128;++k) cand[nCand++] = roots[k];
        }
    }
    if (nCand==0) return -1.0f;
    sortSmall(cand, nCand);

    float uniq[128];
    int nUniq=0;
    for (int i=0;i<nCand;++i) {
        float v = cand[i];
        if (v <= 1e-10) continue;
        if (nUniq==0 || fabs(v-uniq[nUniq-1]) > 1e-4) uniq[nUniq++] = v;
    }
    if (nUniq==0) return -1.0f;

    bool inside0 = pointInGeomGPU(M, geomId, pos);
    for (int i=0;i<nUniq;++i) {
        float sDist = uniq[i];
        float s0 = fmax(0.0f, sDist - 1e-4f);
        float s1 = sDist + 1e-4;
        Vec3d pA = v3_add(pos, v3_mul(dHat, s0));
        Vec3d pB = v3_add(pos, v3_mul(dHat, s1));
        bool inA = pointInGeomGPU(M, geomId, pA);
        bool inB = pointInGeomGPU(M, geomId, pB);
        if (inA != inB) {
            if (sDist <= 1e-10 && inB==inside0) continue;
            return sDist;
        }
    }
    return -1.0f;
}

__device__ bool pointInUniverseLocalGPU(const GPUModel& M, int uidx, Vec3d pLocal) {
    const GPUUniverse& u = M.universes[uidx];
    return pointInGeomGPU(M, u.boundGeom, pLocal);
}

__device__ int findGeometryAtGPU(const GPUModel& M, int rootUid, Vec3d pWorld) {
    // DFS with explicit stack.
    int uStack[64];
    Vec3d pStack[64];
    int childIt[64];
    int sp=0;

    Vec3d p0 = v3_sub(pWorld, M.universes[rootUid].pos);
    if (!pointInUniverseLocalGPU(M, rootUid, p0)) return -1;

    uStack[sp] = rootUid;
    pStack[sp] = p0;
    childIt[sp] = 0;
    ++sp;

    while (sp>0) {
        int uidx = uStack[sp-1];
        Vec3d pLocal = pStack[sp-1];
        int it = childIt[sp-1];
        const GPUUniverse& u = M.universes[uidx];

        if (it < u.nChild) {
            int cidx = M.uChildren[u.childOff + it];
            childIt[sp-1] = it + 1;
            Vec3d pc = v3_sub(pLocal, M.universes[cidx].pos);
            if (!pointInUniverseLocalGPU(M, cidx, pc)) continue;
            // descend
            if (sp < 64) {
                uStack[sp] = cidx;
                pStack[sp] = pc;
                childIt[sp] = 0;
                ++sp;
            }
            continue;
        }

        // No more children: search geometries (reverse order)
        for (int gi=u.nGeom-1; gi>=0; --gi) {
            int gidx = M.uGeoms[u.geomOff + gi];
            if (pointInGeomGPU(M, gidx, pLocal)) return gidx;
        }

        // pop
        --sp;
    }
    return -1;
}

__device__ float sigmaTotAtGPU(const GPUModel& M, int mixId, int inelastic, float E) {
    if (mixId < 0) return 0.0f;
    const GPUMix& mx = M.mixes[mixId];
    float sum = 0.0f;
    for (int i=0;i<mx.nComp;++i) {
        const GPUMixComp& c = M.comps[mx.compOff + i];
        const float w = c.rho * c.prop;
        float s2 = xsInterp(M, c.nuc, 0, E);
        float s4 = inelastic ? xsInterp(M, c.nuc, 1, E) : 0.0f;
        float s18 = xsInterp(M, c.nuc, 2, E);
        float s102 = xsInterp(M, c.nuc, 3, E);
        sum += w * (s2 + s4 + s18 + s102);
    }
    return sum;
}

__device__ float macroscopicSigmaAtGPU(const GPUModel& M, int rootUid, Vec3d pos, int inelastic, float E, int& outGeom, int& outMix) {
    int gid = findGeometryAtGPU(M, rootUid, pos);
    if (gid < 0) { outGeom=-1; outMix=-1; return -1.0f; }
    outGeom = gid;
    outMix  = M.geoms[gid].mixId;
    if (outMix < 0) return 0.0f;
    return sigmaTotAtGPU(M, outMix, inelastic, E);
}

__device__ float nearestCollisionGPU(const GPUModel& M, int rootUid, Vec3d posWorld, Vec3d dirWorld) {
    const float INF = 1e300;
    float best = INF;

    int uStack[2048];
    Vec3d pStack[2048];
    int sp=0;

    uStack[sp] = rootUid;
    pStack[sp] = v3_sub(posWorld, M.universes[rootUid].pos);
    ++sp;

    while (sp>0) {
        --sp;
        int uidx = uStack[sp];
        Vec3d pLocal = pStack[sp];
        const GPUUniverse& u = M.universes[uidx];

        // This universe
        for (int gi=0; gi<u.nGeom; ++gi) {
            int gidx = M.uGeoms[u.geomOff + gi];
            float d = geometryCollisionGPU(M, gidx, pLocal, dirWorld);
            if (d > 0.0f && d < best) best = d;
        }

        // Child universes
        for (int ci=0; ci<u.nChild; ++ci) {
            int cidx = M.uChildren[u.childOff + ci];
            Vec3d pc = v3_sub(pLocal, M.universes[cidx].pos);
            float dBound = geometryCollisionGPU(M, M.universes[cidx].boundGeom, pc, dirWorld);
            if (dBound > 0.0f && dBound <= best && sp < 2047) {
                uStack[sp] = cidx;
                pStack[sp] = pc;
                ++sp;
            }
        }
    }

    return (best<INF)? best : -1.0f;
}


// --- GPU MESH ---

__device__ int meshIndexGPU(const GPURunParams& P, Vec3d p) {
    if (p.x < P.pmin.x || p.x > P.pmax.x ||
        p.y < P.pmin.y || p.y > P.pmax.y ||
        p.z < P.pmin.z || p.z > P.pmax.z) return -1;

    float fx = (p.x - P.pmin.x) * P.invH.x;
    float fy = (p.y - P.pmin.y) * P.invH.y;
    float fz = (p.z - P.pmin.z) * P.invH.z;

    int ix = (int)floor(fx);
    int iy = (int)floor(fy);
    int iz = (int)floor(fz);

    ix = clampi(ix, 0, P.nx-1);
    iy = clampi(iy, 0, P.ny-1);
    iz = clampi(iz, 0, P.nz-1);

    return (iz*P.ny + iy)*P.nx + ix;
}

__device__ void tallyTrackToMeshGPU(const GPURunParams& P, float* dst, Vec3d p0, Vec3d dHat, float segLen, float wPerLen) {
    if (segLen <= 0.0f || wPerLen==0.0f) return;
    Vec3d p = p0;
    Vec3d p1 = v3_add(p0, v3_mul(dHat, segLen));

    if ((p.x < P.pmin.x && p1.x < P.pmin.x) || (p.x > P.pmax.x && p1.x > P.pmax.x) ||
        (p.y < P.pmin.y && p1.y < P.pmin.y) || (p.y > P.pmax.y && p1.y > P.pmax.y) ||
        (p.z < P.pmin.z && p1.z < P.pmin.z) || (p.z > P.pmax.z && p1.z > P.pmax.z)) {
        return;
    }

    float x = (p.x - P.pmin.x) * P.invH.x;
    float y = (p.y - P.pmin.y) * P.invH.y;
    float z = (p.z - P.pmin.z) * P.invH.z;

    int ix = clampi((int)floor(x), 0, P.nx-1);
    int iy = clampi((int)floor(y), 0, P.ny-1);
    int iz = clampi((int)floor(z), 0, P.nz-1);

    int stepX = (dHat.x > 0) ? 1 : (dHat.x < 0 ? -1 : 0);
    int stepY = (dHat.y > 0) ? 1 : (dHat.y < 0 ? -1 : 0);
    int stepZ = (dHat.z > 0) ? 1 : (dHat.z < 0 ? -1 : 0);

    float tMaxX, tMaxY, tMaxZ;
    float tDeltaX, tDeltaY, tDeltaZ;

    if (stepX != 0) {
        float nextBoundary = (stepX > 0) ? (floor(x) + 1.0f) : floor(x);
        tMaxX = (nextBoundary - x) / (dHat.x * P.invH.x);
        tDeltaX = (1.0f) / (fabs(dHat.x) * P.invH.x);
    } else {
        tMaxX = 1e300; tDeltaX = 1e300;
    }

    if (stepY != 0) {
        float nextBoundary = (stepY > 0) ? (floor(y) + 1.0f) : floor(y);
        tMaxY = (nextBoundary - y) / (dHat.y * P.invH.y);
        tDeltaY = (1.0f) / (fabs(dHat.y) * P.invH.y);
    } else {
        tMaxY = 1e300; tDeltaY = 1e300;
    }

    if (stepZ != 0) {
        float nextBoundary = (stepZ > 0) ? (floor(z) + 1.0f) : floor(z);
        tMaxZ = (nextBoundary - z) / (dHat.z * P.invH.z);
        tDeltaZ = (1.0f) / (fabs(dHat.z) * P.invH.z);
    } else {
        tMaxZ = 1e300; tDeltaZ = 1e300;
    }

    float t = 0.0f;
    while (t < segLen) {
        int idx = (iz*P.ny + iy)*P.nx + ix;
        float tNext = fmin(segLen, fmin(tMaxX, fmin(tMaxY, tMaxZ)));
        float dl = tNext - t;
        if (dl > 0.0f && idx >= 0) {
            atomicAdd(&dst[idx], dl * wPerLen);
        }
        t = tNext;
        if (t >= segLen) break;
        if (tMaxX <= tMaxY && tMaxX <= tMaxZ) {
            ix += stepX;
            tMaxX += tDeltaX;
            if (ix < 0 || ix >= P.nx) break;
        } else if (tMaxY <= tMaxX && tMaxY <= tMaxZ) {
            iy += stepY;
            tMaxY += tDeltaY;
            if (iy < 0 || iy >= P.ny) break;
        } else {
            iz += stepZ;
            tMaxZ += tDeltaZ;
            if (iz < 0 || iz >= P.nz) break;
        }
    }
}


// --- TALLIES ---

__device__ float CalculateCFEAddGPU(float E, float SigmaRef) {
    const float v = neutronSpeed(E);
    return (v>0.0f && SigmaRef>0.0f) ? (1.0f / (v * SigmaRef)) : 0.0f;
}

__device__ void scoreCFEDensityGlobalGPU(const GPURunParams& P, GPUOutputs& O, int batch, float E, float SigmaRef) {
    atomicAdd(&O.cfe_global_time[batch], CalculateCFEAddGPU(E, SigmaRef));
}

__device__ void scoreCFEGPU(const GPURunParams& P, GPUOutputs& O, Vec3d pos, float E, float SigmaRef) {
    const int idx = meshIndexGPU(P, pos);
    if (idx < 0) return;
    atomicAdd(&O.cfe_density[idx], CalculateCFEAddGPU(E, SigmaRef));
}

__device__ void scoreTLESegmentPerMixGPU(const GPUModel& M, const GPURunParams& P, GPUOutputs& O,
                                                         int batch, int mixId, float E, float segLen) {
    if (mixId < 0 || segLen<=0.0f) return;
    const GPUMix& mx = M.mixes[mixId];
    for (int i=0;i<mx.nComp;++i) {
        const GPUMixComp& c = M.comps[mx.compOff + i];
        const float w = c.rho*c.prop;
        const float s2 = xsInterp(M, c.nuc, 0, E);
        const float s4 = (P.inelastic?xsInterp(M, c.nuc, 1, E):0.0f);
        const float s18 = xsInterp(M, c.nuc, 2, E);
        const float s102 = xsInterp(M, c.nuc, 3, E);
        const float Stot = w*(s2+s4+s18+s102);
        const float Sabs = w*(s18+s102);
        atomicAdd(&O.tle_Rtot[batch*M.NNuc + c.nuc], segLen*Stot);
        atomicAdd(&O.tle_Rabs[batch*M.NNuc + c.nuc], segLen*Sabs);
    }
}

__device__ void scoreCFERxPerMaterialGPU(const GPUModel& M, const GPURunParams& P, GPUOutputs& O,
                                                         int batch, int nuc, const GPUMixComp& c, float E, float SigmaRef) {
    const float w = c.rho*c.prop;
    const float s2 = xsInterp(M, nuc, 0, E);
    const float s4 = (P.inelastic?xsInterp(M, nuc, 1, E):0.0f);
    const float s18 = xsInterp(M, nuc, 2, E);
    const float s102 = xsInterp(M, nuc, 3, E);
    const float Stot = w*(s2+s4+s18+s102);
    const float Sabs = w*(s18+s102);
    const float inv = (SigmaRef>0.0f) ? (1.0f/SigmaRef) : 0.0f;
    atomicAdd(&O.cfe_Rtot[batch*M.NNuc + nuc], Stot*inv);
    atomicAdd(&O.cfe_Rabs[batch*M.NNuc + nuc], Sabs*inv);
}

__device__ int pickCompInMixGPU(const GPUModel& M, const GPURunParams& P, RNG& rng, int mixId, float E) {
    const GPUMix& mx = M.mixes[mixId];
    float tot = 0.0f;
    for (int i=0;i<mx.nComp;++i) {
        const GPUMixComp& c = M.comps[mx.compOff + i];
        const float w = c.rho*c.prop;
        float s2 = xsInterp(M, c.nuc, 0, E);
        float s4 = (P.inelastic?xsInterp(M, c.nuc, 1, E):0.0f);
        float s18 = xsInterp(M, c.nuc, 2, E);
        float s102 = xsInterp(M, c.nuc, 3, E);
        tot += w*(s2+s4+s18+s102);
    }
    if (tot <= 0.0f) return 0;
    float r = rng.uni()*tot;
    float csum = 0.0f;
    for (int i=0;i<mx.nComp;++i) {
        const GPUMixComp& c = M.comps[mx.compOff + i];
        const float w = c.rho*c.prop;
        float s2 = xsInterp(M, c.nuc, 0, E);
        float s4 = (P.inelastic?xsInterp(M, c.nuc, 1, E):0.0f);
        float s18 = xsInterp(M, c.nuc, 2, E);
        float s102 = xsInterp(M, c.nuc, 3, E);
        csum += w*(s2+s4+s18+s102);
        if (r <= csum) return i;
    }
    return mx.nComp-1;
}

__device__ int sampleRxInNuclideGPU(const GPUModel& M, const GPURunParams& P, RNG& rng, int nuc, float E) {
    float w0 = fmax(0.0f, xsInterp(M, nuc, 0, E));
    float w1 = P.inelastic ? fmax(0.0f, xsInterp(M, nuc, 1, E)) : 0.0f;
    float w2 = fmax(0.0f, xsInterp(M, nuc, 2, E));
    float w3 = fmax(0.0f, xsInterp(M, nuc, 3, E));
    float tot = w0+w1+w2+w3;
    if (tot <= 0.0f) return 0;
    float r = rng.uni()*tot;
    if ((r -= w0) <= 0) return 0;
    if ((r -= w1) <= 0) return 1;
    if ((r -= w2) <= 0) return 2;
    return 3;
}

__device__ int rxToColGPU(const GPUOutputs& O, int rx) {
    // rx codes: 0 elastic, 1 inelastic, 2 fission, 3 capture
    if (O.Rcols == 4) return rx;
    // Rcols==3 (no inelastic): columns are [elastic, fission, capture]
    if (rx == 0) return 0;
    if (rx == 2) return 1;
    return 2;
}

// --- KERNELS ---

__global__ void simulateOnGPUExternal(const GPUModel M, const GPURunParams P, GPUOutputs O) {
    int batch = blockIdx.x;
    if (batch >= P.batches) return;

    const int tid = threadIdx.x;
    const int stride = blockDim.x;

    // Seed per thread and batch.
    RNG rng;
    rng.s = 0x9e3779b97f4a7c15ULL ^ ((unsigned long long)batch<<32) ^ (unsigned long long)(tid+1);

    // Each thread processes multiple batches via striding
    for (int h = tid; h < P.historiesPerBatch; h += stride) {
        // Local stack for fission descendants within this source history.
        BankNeutron stack[64];
        int sp = 0;
        stack[sp++] = BankNeutron{P.sourcePos, P.sourceE, 1};

        while (sp > 0) {
            BankNeutron bn = stack[--sp];
            Vec3d pos = bn.pos;
            float E = bn.E;
            Vec3d dir = isoDirGPU(rng);

            float vel = neutronSpeed(E);

            for (int step=0; step<P.maxSteps; ++step) {
                Vec3d dHat = v3_norm(dir);
                if (P.track == 0) {
                    // Surface tracking 
                    int gid=-1, mix=-1;
                    float Sigma = macroscopicSigmaAtGPU(M, M.rootUid, pos, P.inelastic, E, gid, mix);
                    if (Sigma < 0.0f || mix < 0) { atomicAdd(&O.leaks[batch], 1); break; }

                    float r = fmax(1e-14, (double)rng.uni());
                    float l = (Sigma>0.0f) ? (-log(r)/Sigma) : 1e300;
                    float dSurf = nearestCollisionGPU(M, M.rootUid, pos, dir);
                    if (dSurf < 0.0f) { atomicAdd(&O.leaks[batch], 1); break; }

                    float travel = (l < dSurf) ? l : dSurf;

                    tallyTrackToMeshGPU(P, O.tle_density, pos, dHat, travel, (vel>0.0f)?(1.0f/vel):0.0f);
                    scoreTLESegmentPerMixGPU(M, P, O, batch, mix, E, travel);

                    pos = v3_add(pos, v3_mul(dHat, travel));

                    // Collision?
                    if (l >= dSurf) {
                        pos = v3_add(pos, v3_mul(dHat, 1e-4));
                        continue;
                    }

                    // Real collision
                    int vid = meshIndexGPU(P, pos);
                    if (vid >= 0) atomicAdd(&O.meshAnalogColl[vid], 1.0f);

                    const float SigmaRef = Sigma;
                    scoreCFEDensityGlobalGPU(P, O, batch, E, SigmaRef);
                    scoreCFEGPU(P, O, pos, E, SigmaRef);

                    // Sample component and rx
                    int ci = pickCompInMixGPU(M, P, rng, mix, E);
                    const GPUMixComp& c = M.comps[M.mixes[mix].compOff + ci];
                    int nuc = c.nuc;
                    scoreCFERxPerMaterialGPU(M, P, O, batch, nuc, c, E, SigmaRef);

                    int rx = sampleRxInNuclideGPU(M, P, rng, nuc, E);

                    // Analog reaction counts per nuclide
                    {
                        const int col = rxToColGPU(O, rx);
                        if (col >= 0) {
                            const int idx = ((batch * M.NNuc + nuc) * O.Rcols) + col;
                            atomicAdd(&O.statM[idx], 1);
                        }
                    }

                    if (rx == 3) {
                        // capture
                        break;
                    } else if (rx == 2) {
                        // fission
                        float nuBar = nuBarInterp(M, nuc, E);
                        if (nuBar <= 0.0f) nuBar = 2.43;
                        int nEmit = sampleMultiplicityGPU(rng, nuBar);
                        atomicAdd(&O.fissionChildren[batch], nEmit);
                        for (int k=0;k<nEmit;++k) {
                            if (sp < 64) {
                                BankNeutron child;
                                child.pos = pos;
                                child.E = sampleFissionEnergyMeVGPU(rng);
                                child.isSource = 0;
                                stack[sp++] = child;
                            }
                        }
                        break;
                    } else if (rx == 0) {
                        // elastic
                        float A = fmax(1.0f, M.nuclides[nuc].A);
                        float alpha = ((A-1.0f)/(A+1.0f)); alpha = alpha*alpha;
                        float rr = rng.uni();
                        E = E*(alpha + (1.0f-alpha)*rr);
                        vel = neutronSpeed(E);
                    } else {
                        // inelastic
                        E = fmax(0.0f, E - 1e-4f);
                        vel = neutronSpeed(E);
                    }

                } else {
                    // Delta tracking
                    float r = fmax(1e-14f, rng.uni());
                    float l = -log(r) / fmax(1e-30, (double)P.SigmaM);

                    tallyTrackToMeshGPU(P, O.tle_density, pos, dHat, l, (vel>0.0f)?(1.0f/vel):0.0f);

                    Vec3d newp = v3_add(pos, v3_mul(dHat, l));

                    int gid=-1, mix=-1;
                    float SigmaLoc = macroscopicSigmaAtGPU(M, M.rootUid, newp, P.inelastic, E, gid, mix);
                    if (SigmaLoc < 0.0f || mix < 0) { atomicAdd(&O.leaks[batch], 1); break; }

                    pos = newp;

                    // Collision estimator
                    scoreCFEDensityGlobalGPU(P, O, batch, E, P.SigmaM);
                    scoreCFEGPU(P, O, pos, E, P.SigmaM);

                    int ci = pickCompInMixGPU(M, P, rng, mix, E);
                    const GPUMixComp& c = M.comps[M.mixes[mix].compOff + ci];
                    scoreCFERxPerMaterialGPU(M, P, O, batch, c.nuc, c, E, P.SigmaM);

                    float accept = fmin(1.0, SigmaLoc / fmax(1e-30, (double)P.SigmaM));
                    if (rng.uni() >= accept) {
                        // virtual
                        continue;
                    }

                    // Real collision
                    int vid = meshIndexGPU(P, pos);
                    if (vid >= 0) atomicAdd(&O.meshAnalogColl[vid], 1.0f);

                    int nuc = c.nuc;
                    int rx = sampleRxInNuclideGPU(M, P, rng, nuc, E);

                    // Analog reaction counts per nuclide
                    {
                        const int col = rxToColGPU(O, rx);
                        if (col >= 0) {
                            const int idx = ((batch * M.NNuc + nuc) * O.Rcols) + col;
                            atomicAdd(&O.statM[idx], 1);
                        }
                    }

                    if (rx == 3) {
                        // capture
                        break; 
                    } else if (rx == 2) {
                        // fission
                        float nuBar = nuBarInterp(M, nuc, E);
                        if (nuBar <= 0.0f) nuBar = 2.43;
                        int nEmit = sampleMultiplicityGPU(rng, nuBar);
                        atomicAdd(&O.fissionChildren[batch], nEmit);
                        for (int k=0;k<nEmit;++k) {
                            if (sp < 64) {
                                BankNeutron child{pos, sampleFissionEnergyMeVGPU(rng), 0};
                                stack[sp++] = child;
                            }
                        }
                        break;
                    }
                    if (rx == 0) {
                        // elastic
                        float A = fmax(1.0f, M.nuclides[nuc].A);
                        float alpha = ((A-1.0f)/(A+1.0f)); alpha = alpha*alpha;
                        float rr = rng.uni();
                        E = E*(alpha + (1.0f-alpha)*rr);
                        vel = neutronSpeed(E);
                    } else {
                        // inelastic
                        E = fmax(0.0f, E - 1e-4f);
                        vel = neutronSpeed(E);
                    }
                }
            }
        }
    }
}

__global__ void simulateOnGPUCritical(const GPUModel M, const GPURunParams P, GPUOutputs O,
                                     const BankNeutron* bankCur, int nCur,
                                     BankNeutron* bankNext, int maxNext, int* nNextOut,
                                     int cycle) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nCur) return;
    
    // One history per thread
    RNG rng;
    rng.s = 0x9e3779b97f4a7c15ULL ^ ((unsigned long long)cycle<<32) ^ (unsigned long long)(i+1);

    BankNeutron bn = bankCur[i];
    Vec3d pos = bn.pos;
    float E  = bn.E;
    Vec3d dir = isoDirGPU(rng);

    float vel = neutronSpeed(E);

    for (int step=0; step<P.maxSteps; ++step) {
        Vec3d dHat = v3_norm(dir);
        if (P.track == 0) {
            // Surface tracking
            int gid=-1, mix=-1;
            float Sigma = macroscopicSigmaAtGPU(M, M.rootUid, pos, P.inelastic, E, gid, mix);
            if (Sigma < 0.0f || mix < 0) { atomicAdd(&O.leaks[cycle], 1); break; }

            float r = fmax(1e-14f, rng.uni());
            float l = (Sigma>0.0f) ? (-log(r)/Sigma) : 1e300;
            float dSurf = nearestCollisionGPU(M, M.rootUid, pos, dir);
            if (dSurf < 0.0f) { atomicAdd(&O.leaks[cycle], 1); break; }

            float travel = (l < dSurf) ? l : dSurf;

            tallyTrackToMeshGPU(P, O.tle_density, pos, dHat, travel, (vel>0.0f)?(1.0f/vel):0.0f);
            scoreTLESegmentPerMixGPU(M, P, O, cycle, mix, E, travel);

            pos = v3_add(pos, v3_mul(dHat, travel));

            if (l >= dSurf) {
                pos = v3_add(pos, v3_mul(dHat, 1e-4));
                continue;
            }

            // Real collision
            int vid = meshIndexGPU(P, pos);
            if (vid >= 0) atomicAdd(&O.meshAnalogColl[vid], 1.0f);

            const float SigmaRef = Sigma;
            scoreCFEDensityGlobalGPU(P, O, cycle, E, SigmaRef);
            scoreCFEGPU(P, O, pos, E, SigmaRef);

            // Component + reaction
            int ci = pickCompInMixGPU(M, P, rng, mix, E);
            const GPUMixComp& c = M.comps[M.mixes[mix].compOff + ci];
            int nuc = c.nuc;
            scoreCFERxPerMaterialGPU(M, P, O, cycle, nuc, c, E, SigmaRef);
            int rx = sampleRxInNuclideGPU(M, P, rng, nuc, E);

            // Analog reaction counts per nuclide
            {
                const int col = rxToColGPU(O, rx);
                if (col >= 0) {
                    const int idx = ((cycle * M.NNuc + nuc) * O.Rcols) + col;
                    atomicAdd(&O.statM[idx], 1);
                }
            }

            if (rx == 3) {
                // capture
                break;
            } else if (rx == 2) {
                // fission 
                float nuBar = nuBarInterp(M, nuc, E);
                if (nuBar <= 0.0f) nuBar = 2.43;
                int nEmit = sampleMultiplicityGPU(rng, nuBar);
                atomicAdd(&O.fissionChildren[cycle], nEmit);

                int base = atomicAdd(nNextOut, nEmit);
                for (int k=0;k<nEmit;++k) {
                    int idx = base + k;
                    if (idx < maxNext) {
                        bankNext[idx] = BankNeutron{pos, sampleFissionEnergyMeVGPU(rng), 0};
                    }
                }
                break;
            }
            if (rx == 0) {
                // elastic
                float A = fmax(1.0f, M.nuclides[nuc].A);
                float alpha = ((A-1.0f)/(A+1.0f)); alpha = alpha*alpha;
                float rr = rng.uni();
                E = E*(alpha + (1.0f-alpha)*rr);
                vel = neutronSpeed(E);
            } else {
                // inelastic
                E = fmax(0.0f, E - 1e-4f);
                vel = neutronSpeed(E);
            }

        } else {
            // Delta tracking
            float r = fmax(1e-14f, rng.uni());
            float l = -log(r) / fmax(1e-30, (double)P.SigmaM) * 0.1;

            tallyTrackToMeshGPU(P, O.tle_density, pos, dHat, l, (vel>0.0f)?(1.0f/vel):0.0f);

            Vec3d newp = v3_add(pos, v3_mul(dHat, l));
            int gid=-1, mix=-1;
            float SigmaLoc = macroscopicSigmaAtGPU(M, M.rootUid, newp, P.inelastic, E, gid, mix);
            if (SigmaLoc < 0.0f || mix < 0) { atomicAdd(&O.leaks[cycle], 1); break; }
            pos = newp;

            scoreCFEDensityGlobalGPU(P, O, cycle, E, P.SigmaM);
            scoreCFEGPU(P, O, pos, E, P.SigmaM);

            int ci = pickCompInMixGPU(M, P, rng, mix, E);
            const GPUMixComp& c = M.comps[M.mixes[mix].compOff + ci];
            scoreCFERxPerMaterialGPU(M, P, O, cycle, c.nuc, c, E, P.SigmaM);

            float accept = fmin(1.0, SigmaLoc / fmax(1e-30, (double)P.SigmaM));
            if (rng.uni() >= accept) continue;

            int vid = meshIndexGPU(P, pos);
            if (vid >= 0) atomicAdd(&O.meshAnalogColl[vid], 1.0f);

            int nuc = c.nuc;
            int rx = sampleRxInNuclideGPU(M, P, rng, nuc, E);

            // Analog reaction counts per nuclide
            {
                const int col = rxToColGPU(O, rx);
                if (col >= 0) {
                    const int idx = ((cycle * M.NNuc + nuc) * O.Rcols) + col;
                    atomicAdd(&O.statM[idx], 1);
                }
            }
            if (rx == 3) {
                // capture
                break;
            } else if (rx == 2) {
                // fission
                float nuBar = nuBarInterp(M, nuc, E);
                if (nuBar <= 0.0f) nuBar = 2.43;
                int nEmit = sampleMultiplicityGPU(rng, nuBar);
                atomicAdd(&O.fissionChildren[cycle], nEmit);

                int base = atomicAdd(nNextOut, nEmit);
                for (int k=0;k<nEmit;++k) {
                    int idx = base + k;
                    if (idx < maxNext) {
                        bankNext[idx] = BankNeutron{pos, sampleFissionEnergyMeVGPU(rng), 0};
                    }
                }
                break;
            }
            if (rx == 0) {
                // elastic
                float A = fmax(1.0f, M.nuclides[nuc].A);
                float alpha = ((A-1.0f)/(A+1.0f)); alpha = alpha*alpha;
                float rr = rng.uni();
                E = E*(alpha + (1.0f-alpha)*rr);
                vel = neutronSpeed(E);
            } else {
                // inelastic
                E = fmax(0.0f, E - 1e-4f);
                vel = neutronSpeed(E);
            }
        }
    }
}

__global__ void resampleBankKernel(const BankNeutron* bankNext, int nNext,
                                  BankNeutron* bankCur, int nCur,
                                  Vec3d srcPos, float srcE, unsigned long long seed) {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nCur) return;
    RNG rng; rng.s = seed ^ (unsigned long long)(i+1);

    if (nNext <= 0) {
        bankCur[i] = BankNeutron{srcPos, srcE, 1};
        return;
    }
    int j = (int)floor(rng.uni() * (float)nNext);
    if (j < 0) j = 0; if (j >= nNext) j = nNext-1;
    BankNeutron pick = bankNext[j];
    pick.isSource = 0;
    bankCur[i] = pick;
}

// --- GPU START FUNCTIONS ---

struct GPUModelHost {
    // Host arrays
    vector<GPUUniverse> hU;
    vector<int> hUChildren;
    vector<int> hUGeoms;
    vector<GPUGeometry> hG;
    vector<GPUShape> hS;
    vector<GPURPNNode> hN;
    vector<GPUMix> hMix;
    vector<GPUMixComp> hComp;
    vector<GPUNuclide> hNuc;
    vector<std::string> nucNames;
    vector<const Nuclide*> hNucData; // host-only: points into cached nuclide database

    vector<float> hEgrid;
    int NE = 0;
    float logEmin = 0.0f;
    float invDlog = 1.0f;

    vector<float> hXS;    // nuc*4*NE
    vector<float> hNuBar; // nuc*NE

    // Device pointers
    GPUUniverse* dU=nullptr;
    int* dUChildren=nullptr;
    int* dUGeoms=nullptr;
    GPUGeometry* dG=nullptr;
    GPUShape* dS=nullptr;
    GPURPNNode* dN=nullptr;
    GPUMix* dMix=nullptr;
    GPUMixComp* dComp=nullptr;
    GPUNuclide* dNuc=nullptr;
    float* dEgrid=nullptr;
    float* dXS=nullptr;
    float* dNuBar=nullptr;

    GPUModel GPUView{};

    void freeGPU() {
        if (dU) CUDA_CHECK(cudaFree(dU));
        if (dUChildren) CUDA_CHECK(cudaFree(dUChildren));
        if (dUGeoms) CUDA_CHECK(cudaFree(dUGeoms));
        if (dG) CUDA_CHECK(cudaFree(dG));
        if (dS) CUDA_CHECK(cudaFree(dS));
        if (dN) CUDA_CHECK(cudaFree(dN));
        if (dMix) CUDA_CHECK(cudaFree(dMix));
        if (dComp) CUDA_CHECK(cudaFree(dComp));
        if (dNuc) CUDA_CHECK(cudaFree(dNuc));
        if (dEgrid) CUDA_CHECK(cudaFree(dEgrid));
        if (dXS) CUDA_CHECK(cudaFree(dXS));
        if (dNuBar) CUDA_CHECK(cudaFree(dNuBar));
        dU=nullptr; dUChildren=nullptr; dUGeoms=nullptr; dG=nullptr; dS=nullptr; dN=nullptr; 
        dMix=nullptr; dComp=nullptr;dNuc=nullptr; dEgrid=nullptr; dXS=nullptr; dNuBar=nullptr;
    }
};

int opToInt(Op op) {
    return (int)op;
}

int addNuclide(GPUModelHost& H, const Material& m) {
    // Deduplicate by nuclide symbol.
    if (!m.nuc) return -1;
    const std::string& sym = m.nuc->sym;
    auto it = std::find(H.nucNames.begin(), H.nucNames.end(), sym);
    if (it != H.nucNames.end()) return (int)std::distance(H.nucNames.begin(), it);

    int idx = (int)H.nucNames.size();
    H.nucNames.push_back(sym);
    H.hNucData.push_back(m.nuc);

    GPUNuclide nu;
    nu.A = (m.nuc->aw > 0.0f) ? m.nuc->aw : (float)m.nuc->a;
    nu.T_K = (float)m.nuc->T + 273.15;
    nu.isFuel = isFuelSym(sym) ? 1 : 0;
    H.hNuc.push_back(nu);

    return idx;
}

void ensureXSStorage(GPUModelHost& H) {
    if (H.NE > 0) return;
    H.hEgrid = logspace(-11.0f, std::log10(20.0f), 500);
    H.NE = (int)H.hEgrid.size();
    H.logEmin = std::log10(H.hEgrid.front());
    const float logEmax = std::log10(H.hEgrid.back());
    const float dlog = (H.NE>1) ? ((logEmax - H.logEmin)/(H.NE-1)) : 1.0f;
    H.invDlog = (dlog>0.0f) ? (1.0f/dlog) : 1.0f;
}

void buildNuclideTables(GPUModelHost& H, bool inelastic) {
    ensureXSStorage(H);

    const int NNuc = (int)H.hNuc.size();
    H.hXS.assign((size_t)NNuc * 4 * (size_t)H.NE, 0.0f);
    H.hNuBar.assign((size_t)NNuc * (size_t)H.NE, 0.0f);

    for (int ni = 0; ni < NNuc; ++ni) {
        const Nuclide* n = (ni < (int)H.hNucData.size()) ? H.hNucData[ni] : nullptr;
        if (!n) continue;
        const std::string& sym = H.nucNames[ni];

        for (int ei = 0; ei < H.NE; ++ei) {
            float E = H.hEgrid[ei];
            float s2   = interpMT(n->mt, 2, E);
            float s4   = inelastic ? interpMT(n->mt, 4, E) : 0.0f;
            float s18  = interpMT(n->mt, 18, E);
            float s102 = interpMT(n->mt, 102, E);

            H.hXS[((ni*4 + 0)*H.NE) + ei] = s2;
            H.hXS[((ni*4 + 1)*H.NE) + ei] = s4;
            H.hXS[((ni*4 + 2)*H.NE) + ei] = s18;
            H.hXS[((ni*4 + 3)*H.NE) + ei] = s102;

            float nuBar = valueInterp(n->neutrons, E);
            if (nuBar <= 0.0f && isFuelSym(sym)) nuBar = 2.43;
            if (nuBar <= 0.0f) nuBar = 0.0f;
            H.hNuBar[(ni*H.NE) + ei] = nuBar;
        }
    }
}

int addMix(GPUModelHost& H, const Geometry& g, vector<Material>& outAllNuclides) {
    if (g.mats.empty()) return -1;

    int compOff = (int)H.hComp.size();
    int nComp = (int)g.mats.size();

    for (const auto& m : g.mats) {
        int nuc = addNuclide(H, m);
        outAllNuclides.push_back(m);

        GPUMixComp c;
        c.nuc = nuc;
        c.rho = m.rho;
        c.prop = m.proportion;
        H.hComp.push_back(c);
    }

    GPUMix mx;
    mx.compOff = compOff;
    mx.nComp = nComp;
    H.hMix.push_back(mx);
    return (int)H.hMix.size()-1;
}

int addGeometryFlat(GPUModelHost& H, const Geometry& g, int mixId) {
    GPUGeometry gg;
    gg.shapeOff = (int)H.hS.size();
    gg.nShapes  = (int)g.shapes.size();
    gg.nodeOff  = (int)H.hN.size();
    gg.nNodes   = (int)g.nodes.size();
    gg.mixId    = mixId;

    for (const auto& s : g.shapes) {
        GPUShape ds;
        ds.A=s.A; ds.B=s.B; ds.C=s.C; ds.D=s.D; ds.E=s.E; ds.F=s.F;
        ds.G=s.G; ds.H=s.H; ds.I=s.I; ds.J=s.J; ds.torus = (int)s.torus;
        H.hS.push_back(ds);
    }

    for (const auto& nd : g.nodes) {
        GPURPNNode dn;
        dn.op = opToInt(nd.op);
        dn.shape = nd.shape;
        H.hN.push_back(dn);
    }

    H.hG.push_back(gg);
    return (int)H.hG.size()-1;
}

int addUniverseRec(GPUModelHost& H, const Universe& U, vector<Material>& allNuclides)
{
    int uidx = (int)H.hU.size();
    H.hU.push_back(GPUUniverse{});

    GPUUniverse du{};
    du.pos = Vec3d{U.pos[0], U.pos[1], U.pos[2]};
    du.boundGeom = addGeometryFlat(H, U.boundingGeometry, -1);

    // ---- FIX: build children first, but don't write to parent's child list yet
    vector<int> directChildIdx;
    directChildIdx.reserve(U.subUniverse.size());
    for (auto& su : U.subUniverse) {
        directChildIdx.push_back(addUniverseRec(H, su, allNuclides));
    }

    du.childOff = (int)H.hUChildren.size();
    du.nChild = (int)directChildIdx.size();
    H.hUChildren.insert(H.hUChildren.end(), directChildIdx.begin(), directChildIdx.end());

    // geometries (already safe because you do it after recursion)
    du.geomOff = (int)H.hUGeoms.size();
    du.nGeom = (int)U.geometries.size();
    for (auto& g : U.geometries) {
        int mixId = addMix(H, g, allNuclides);
        int gidx = addGeometryFlat(H, g, mixId);
        H.hUGeoms.push_back(gidx);
    }

    H.hU[uidx] = du;
    return uidx;
}

// All necessary transformations to CPU data to start GPU processing
GPUModelHost buildGPUModelHost(const Universe& U, bool inelastic) {
    GPUModelHost H;
    vector<Material> allNuclides;
    allNuclides.reserve(1024);
    const int rootUid = addUniverseRec(H, U, allNuclides);

    std::unordered_map<std::string, Material> unique;
    unique.reserve(allNuclides.size()*2);
    for (const auto& m : allNuclides) {
        if (unique.find(m.nuc->sym) == unique.end()) unique.emplace(m.nuc->sym, m);
    }
    vector<Material> uniqList;
    uniqList.reserve(unique.size());
    for (auto& kv : unique) uniqList.push_back(kv.second);

    // Build per-nuclide XS/nuBar tables from cached nuclide data (no copies)
    buildNuclideTables(H, inelastic);

    // Upload to device
    CUDA_CHECK(cudaMalloc(&H.dU, H.hU.size()*sizeof(GPUUniverse)));
    CUDA_CHECK(cudaMalloc(&H.dUChildren, H.hUChildren.size()*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&H.dUGeoms, H.hUGeoms.size()*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&H.dG, H.hG.size()*sizeof(GPUGeometry)));
    CUDA_CHECK(cudaMalloc(&H.dS, H.hS.size()*sizeof(GPUShape)));
    CUDA_CHECK(cudaMalloc(&H.dN, H.hN.size()*sizeof(GPURPNNode)));
    CUDA_CHECK(cudaMalloc(&H.dMix, H.hMix.size()*sizeof(GPUMix)));
    CUDA_CHECK(cudaMalloc(&H.dComp, H.hComp.size()*sizeof(GPUMixComp)));
    CUDA_CHECK(cudaMalloc(&H.dNuc, H.hNuc.size()*sizeof(GPUNuclide)));
    CUDA_CHECK(cudaMalloc(&H.dEgrid, H.hEgrid.size()*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&H.dXS, H.hXS.size()*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&H.dNuBar, H.hNuBar.size()*sizeof(float)));

    CUDA_CHECK(cudaMemcpy(H.dU, H.hU.data(), H.hU.size()*sizeof(GPUUniverse), cudaMemcpyHostToDevice));
    if (!H.hUChildren.empty()) CUDA_CHECK(cudaMemcpy(H.dUChildren, H.hUChildren.data(), H.hUChildren.size()*sizeof(int), cudaMemcpyHostToDevice));
    if (!H.hUGeoms.empty())    CUDA_CHECK(cudaMemcpy(H.dUGeoms, H.hUGeoms.data(), H.hUGeoms.size()*sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(H.dG, H.hG.data(), H.hG.size()*sizeof(GPUGeometry), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(H.dS, H.hS.data(), H.hS.size()*sizeof(GPUShape), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(H.dN, H.hN.data(), H.hN.size()*sizeof(GPURPNNode), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(H.dMix, H.hMix.data(), H.hMix.size()*sizeof(GPUMix), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(H.dComp, H.hComp.data(), H.hComp.size()*sizeof(GPUMixComp), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(H.dNuc, H.hNuc.data(), H.hNuc.size()*sizeof(GPUNuclide), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(H.dEgrid, H.hEgrid.data(), H.hEgrid.size()*sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(H.dXS, H.hXS.data(), H.hXS.size()*sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(H.dNuBar, H.hNuBar.data(), H.hNuBar.size()*sizeof(float), cudaMemcpyHostToDevice));

    H.GPUView.universes = H.dU;
    H.GPUView.uChildren = H.dUChildren;
    H.GPUView.uGeoms = H.dUGeoms;
    H.GPUView.geoms = H.dG;
    H.GPUView.shapes = H.dS;
    H.GPUView.nodes = H.dN;
    H.GPUView.mixes = H.dMix;
    H.GPUView.comps = H.dComp;
    H.GPUView.nuclides = H.dNuc;
    H.GPUView.xs = H.dXS;
    H.GPUView.nuBar = H.dNuBar;
    H.GPUView.NU = (int)H.hU.size();
    H.GPUView.rootUid = rootUid;
    H.GPUView.NG = (int)H.hG.size();
    H.GPUView.NMix = (int)H.hMix.size();
    H.GPUView.NComp = (int)H.hComp.size();
    H.GPUView.NNuc = (int)H.hNuc.size();
    H.GPUView.Egrid = H.dEgrid;
    H.GPUView.NE = H.NE;
    H.GPUView.logEmin = H.logEmin;
    H.GPUView.invDlog = H.invDlog;

    return H;
}

// Convert run parameters to GPU compatible
GPURunParams buildRunParamsGPU(const RunParams& P, float SigmaM) {
    GPURunParams G{};
    G.historiesPerBatch = P.historiesPerBatch;
    G.batches = P.batches;
    G.maxSteps = P.maxSteps;
    G.inelastic = P.inelastic ? 1 : 0;
    G.track = (P.track==Tracking::Delta) ? 1 : 0;
    G.src = (P.src==SourceMode::Criticality) ? 1 : 0;
    G.sourcePos = {P.sourcePos[0], P.sourcePos[1], P.sourcePos[2]};
    G.sourceE = P.sourceE;
    G.SigmaM = SigmaM;

    if (P.mesh) {
        G.nx = P.mesh->nx; G.ny = P.mesh->ny; G.nz = P.mesh->nz;
        G.pmin = {P.mesh->pmin[0], P.mesh->pmin[1], P.mesh->pmin[2]};
        G.pmax = {P.mesh->pmax[0], P.mesh->pmax[1], P.mesh->pmax[2]};
        float dx = (P.mesh->pmax[0]-P.mesh->pmin[0]) / std::max(1, P.mesh->nx);
        float dy = (P.mesh->pmax[1]-P.mesh->pmin[1]) / std::max(1, P.mesh->ny);
        float dz = (P.mesh->pmax[2]-P.mesh->pmin[2]) / std::max(1, P.mesh->nz);
        G.invH = { (dx>0)?(1.0f/dx):0.0f, (dy>0)?(1.0f/dy):0.0f, (dz>0)?(1.0f/dz):0.0f };
    } else {
        G.nx=G.ny=G.nz=0;
        G.pmin = {0,0,0}; G.pmax={0,0,0};
        G.invH = {0,0,0};
    }
    return G;
}

// Allocate GPU output
void allocOutputsGPU(const GPUModelHost& H, const GPURunParams& GP, int Rcols, GPUOutputs& O,
                            int*& dLeaks, int*& dFiss, int*& dStat,
                            float*& dCfeRtot, float*& dCfeRabs, float*& dTleRtot, float*& dTleRabs, float*& dCfeTau,
                            float*& dMeshA, float*& dMeshC, float*& dMeshT) {

    const int B = GP.batches;
    const int M = (int)H.hNuc.size();
    const int Nmesh = std::max(1, GP.nx*GP.ny*GP.nz);

    CUDA_CHECK(cudaMalloc(&dLeaks, B*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dFiss,  B*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dStat,  (size_t)B*M*Rcols*sizeof(int)));
    CUDA_CHECK(cudaMalloc(&dCfeRtot, (size_t)B*M*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dCfeRabs, (size_t)B*M*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dTleRtot, (size_t)B*M*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dTleRabs, (size_t)B*M*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dCfeTau,  B*sizeof(float)));

    CUDA_CHECK(cudaMalloc(&dMeshA, (size_t)Nmesh*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dMeshC, (size_t)Nmesh*sizeof(float)));
    CUDA_CHECK(cudaMalloc(&dMeshT, (size_t)Nmesh*sizeof(float)));

    CUDA_CHECK(cudaMemset(dLeaks, 0, B*sizeof(int)));
    CUDA_CHECK(cudaMemset(dFiss,  0, B*sizeof(int)));
    CUDA_CHECK(cudaMemset(dStat,  0, (size_t)B*M*Rcols*sizeof(int)));
    CUDA_CHECK(cudaMemset(dCfeRtot,0,(size_t)B*M*sizeof(float)));
    CUDA_CHECK(cudaMemset(dCfeRabs,0,(size_t)B*M*sizeof(float)));
    CUDA_CHECK(cudaMemset(dTleRtot,0,(size_t)B*M*sizeof(float)));
    CUDA_CHECK(cudaMemset(dTleRabs,0,(size_t)B*M*sizeof(float)));
    CUDA_CHECK(cudaMemset(dCfeTau, 0, B*sizeof(float)));

    CUDA_CHECK(cudaMemset(dMeshA, 0, (size_t)Nmesh*sizeof(float)));
    CUDA_CHECK(cudaMemset(dMeshC, 0, (size_t)Nmesh*sizeof(float)));
    CUDA_CHECK(cudaMemset(dMeshT, 0, (size_t)Nmesh*sizeof(float)));

    O.leaks = dLeaks;
    O.fissionChildren = dFiss;
    O.statM = dStat;
    O.cfe_Rtot = dCfeRtot;
    O.cfe_Rabs = dCfeRabs;
    O.tle_Rtot = dTleRtot;
    O.tle_Rabs = dTleRabs;
    O.cfe_global_time = dCfeTau;
    O.meshAnalogColl = dMeshA;
    O.cfe_density = dMeshC;
    O.tle_density = dMeshT;
    O.Rcols = Rcols;
}

// Copy GPU tallies back into RunOutputs
void downloadOutputsGPU(const GPUModelHost& H, const GPURunParams& GP, int Rcols, const GPUOutputs& O,
                               int* dLeaks, int* dFiss, int* dStat,
                               float* dCfeRtot, float* dCfeRabs, float* dTleRtot, float* dTleRabs, float* dCfeTau,
                               float* dMeshA, float* dMeshC, float* dMeshT,
                               RunOutputs& R, Mesh3D* mesh) {

    const int B = GP.batches;
    const int M = (int)H.hNuc.size();
    const int Nmesh = std::max(1, GP.nx*GP.ny*GP.nz);

    vector<int> hLeaks(B,0), hFiss(B,0);
    vector<int> hStat((size_t)B*M*Rcols, 0);
    vector<float> hCfeRtot((size_t)B*M,0.0f), hCfeRabs((size_t)B*M,0.0f);
    vector<float> hTleRtot((size_t)B*M,0.0f), hTleRabs((size_t)B*M,0.0f);
    vector<float> hTau(B,0.0f);
    vector<float> hMeshA(Nmesh,0.0f), hMeshC(Nmesh,0.0f), hMeshT(Nmesh,0.0f);

    CUDA_CHECK(cudaMemcpy(hLeaks.data(), dLeaks, B*sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hFiss.data(),  dFiss,  B*sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hStat.data(),  dStat,  (size_t)B*M*Rcols*sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hCfeRtot.data(), dCfeRtot, (size_t)B*M*sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hCfeRabs.data(), dCfeRabs, (size_t)B*M*sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hTleRtot.data(), dTleRtot, (size_t)B*M*sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hTleRabs.data(), dTleRabs, (size_t)B*M*sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hTau.data(), dCfeTau, B*sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hMeshA.data(), dMeshA, (size_t)Nmesh*sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hMeshC.data(), dMeshC, (size_t)Nmesh*sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hMeshT.data(), dMeshT, (size_t)Nmesh*sizeof(float), cudaMemcpyDeviceToHost));

    // Fill RunOutputs
    R.T.mesh = mesh;
    R.T.Rcols = Rcols;
    // R.T.nBatches = B;
    // R.T.inelastic = (GP.inelastic != 0);
    R.T.deltaMode = (GP.track != 0);
    R.T.leaks = hLeaks;
    R.fissionChildren = hFiss;
    R.T.cfe_global_time = hTau;

    R.T.matNames = H.nucNames;
    R.T.mat2idx.clear();
    for (int i=0;i<(int)H.nucNames.size();++i) R.T.mat2idx[H.nucNames[i]] = i;

    // statM
    R.T.statM.assign(B, vector<vector<int>>(M, vector<int>(Rcols, 0)));
    for (int b=0;b<B;++b) {
        for (int m=0;m<M;++m) {
            for (int c=0;c<Rcols;++c) {
                R.T.statM[b][m][c] = hStat[((b*M + m)*Rcols) + c];
            }
        }
    }

    R.T.cfe_Rtot.assign(B, vector<float>(M,0.0f));
    R.T.cfe_Rabs.assign(B, vector<float>(M,0.0f));
    R.T.tle_Rtot.assign(B, vector<float>(M,0.0f));
    R.T.tle_Rabs.assign(B, vector<float>(M,0.0f));

    for (int b=0;b<B;++b) {
        for (int m=0;m<M;++m) {
            R.T.cfe_Rtot[b][m] = hCfeRtot[b*M + m];
            R.T.cfe_Rabs[b][m] = hCfeRabs[b*M + m];
            R.T.tle_Rtot[b][m] = hTleRtot[b*M + m];
            R.T.tle_Rabs[b][m] = hTleRabs[b*M + m];
        }
    }

    if (mesh) {
        mesh->meshAnalogColl = std::move(hMeshA);
        mesh->cfe_density = std::move(hMeshC);
        mesh->tle_density = std::move(hMeshT);
    }
}

RunOutputs runExternalCuda(const Universe& U, const RunParams& P) {
    RunOutputs R;

    const int Rcols = P.inelastic ? 4 : 3;

    // Prepare mesh
    if (P.mesh) P.mesh->zero();

    auto t0 = std::chrono::steady_clock::now();

    // Build and upload model
    GPUModelHost H = buildGPUModelHost(U, P.inelastic);

    // Delta majorant
    const vector<int> MTs_total = P.inelastic ? vector<int>{2,4,18,102} : vector<int>{2,18,102};
    const float SigmaM = (P.track==Tracking::Delta) ? majorSigma(U, P.sourceE, MTs_total) : 0.0f;

    GPURunParams GP = buildRunParamsGPU(P, SigmaM);

    // Outputs
    GPUOutputs O{};
    int *dLeaks=nullptr, *dFiss=nullptr, *dStat=nullptr;
    float *dCfeRtot=nullptr,*dCfeRabs=nullptr,*dTleRtot=nullptr,*dTleRabs=nullptr,*dCfeTau=nullptr;
    float *dMeshA=nullptr,*dMeshC=nullptr,*dMeshT=nullptr;
    allocOutputsGPU(H, GP, Rcols, O, dLeaks, dFiss, dStat, dCfeRtot, dCfeRabs, dTleRtot, dTleRabs, dCfeTau, dMeshA, dMeshC, dMeshT);

    // Launch
    const int block = 32;
    dim3 grid(GP.batches, 1, 1);
    simulateOnGPUExternal<<<grid, block>>>(H.GPUView, GP, O);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Download
    downloadOutputsGPU(H, GP, Rcols, O, dLeaks, dFiss, dStat, dCfeRtot, dCfeRabs, dTleRtot, dTleRabs, dCfeTau,
                       dMeshA, dMeshC, dMeshT, R, P.mesh);

    // Fill perf
    auto t1 = std::chrono::steady_clock::now();
    R.perf.histories = P.historiesPerBatch * P.batches;
    R.perf.elapsed_s = std::chrono::duration<float>(t1 - t0).count();

    // Cleanup
    CUDA_CHECK(cudaFree(dLeaks));
    CUDA_CHECK(cudaFree(dFiss));
    CUDA_CHECK(cudaFree(dStat));
    CUDA_CHECK(cudaFree(dCfeRtot));
    CUDA_CHECK(cudaFree(dCfeRabs));
    CUDA_CHECK(cudaFree(dTleRtot));
    CUDA_CHECK(cudaFree(dTleRabs));
    CUDA_CHECK(cudaFree(dCfeTau));
    CUDA_CHECK(cudaFree(dMeshA));
    CUDA_CHECK(cudaFree(dMeshC));
    CUDA_CHECK(cudaFree(dMeshT));

    H.freeGPU();
    return R;
}

RunOutputs runCriticalityCuda(const Universe& U, const RunParams& P, int inactive) {
    RunOutputs R;
    const int Rcols = P.inelastic ? 4 : 3;

    if (P.mesh) P.mesh->zero();

    auto t0 = std::chrono::steady_clock::now();

    GPUModelHost H = buildGPUModelHost(U, P.inelastic);

    const vector<int> MTs_total = P.inelastic ? vector<int>{2,4,18,102} : vector<int>{2,18,102};
    const float SigmaM = (P.track==Tracking::Delta) ? majorSigma(U, P.sourceE, MTs_total) : 0.0f;

    GPURunParams GP = buildRunParamsGPU(P, SigmaM);

    // Outputs 
    GPUOutputs O{};
    int *dLeaks=nullptr, *dFiss=nullptr, *dStat=nullptr;
    float *dCfeRtot=nullptr,*dCfeRabs=nullptr,*dTleRtot=nullptr,*dTleRabs=nullptr,*dCfeTau=nullptr;
    float *dMeshA=nullptr,*dMeshC=nullptr,*dMeshT=nullptr;
    allocOutputsGPU(H, GP, Rcols, O, dLeaks, dFiss, dStat, dCfeRtot, dCfeRabs, dTleRtot, dTleRabs, dCfeTau, dMeshA, dMeshC, dMeshT);

    // Banks
    const int nCur = P.historiesPerBatch;
    const int maxNext = std::max(1, nCur*8);
    BankNeutron* dBankCur=nullptr;
    BankNeutron* dBankNext=nullptr;
    int* dNnext=nullptr;
    CUDA_CHECK(cudaMalloc(&dBankCur, (size_t)nCur*sizeof(BankNeutron)));
    CUDA_CHECK(cudaMalloc(&dBankNext,(size_t)maxNext*sizeof(BankNeutron)));
    CUDA_CHECK(cudaMalloc(&dNnext, sizeof(int)));

    // Init bankCur on host
    vector<BankNeutron> hBank(nCur);
    for (int i=0;i<nCur;++i) hBank[i] = BankNeutron{{P.sourcePos[0],P.sourcePos[1],P.sourcePos[2]}, P.sourceE, 1};
    CUDA_CHECK(cudaMemcpy(dBankCur, hBank.data(), (size_t)nCur*sizeof(BankNeutron), cudaMemcpyHostToDevice));

    // Per-cycle loop
    R.keff_history.clear();
    R.keff_history.reserve(P.batches);

    for (int cyc=0; cyc<P.batches; ++cyc) {
        // Reset next counter
        CUDA_CHECK(cudaMemset(dNnext, 0, sizeof(int)));

        // Launch kernel
        int block = 32;
        int grid = (nCur + block - 1)/block;
        simulateOnGPUCritical<<<grid, block>>>(H.GPUView, GP, O, dBankCur, nCur, dBankNext, maxNext, dNnext, cyc);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        // Get Nnext
        int hNnext=0;
        CUDA_CHECK(cudaMemcpy(&hNnext, dNnext, sizeof(int), cudaMemcpyDeviceToHost));
        float kcur = (nCur>0) ? ((float)hNnext / (float)nCur) : 0.0f;
        R.keff_history.push_back(kcur);

        // Resample
        const int nStore = std::min(hNnext, maxNext);
        unsigned long long seed = 0x1234abcdULL ^ ((unsigned long long)cyc<<32);
        int grid2 = (nCur + block - 1)/block;
        resampleBankKernel<<<grid2, block>>>(dBankNext, nStore, dBankCur, nCur,
                                             Vec3d{P.sourcePos[0],P.sourcePos[1],P.sourcePos[2]}, P.sourceE, seed);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    // Download
    downloadOutputsGPU(H, GP, Rcols, O, dLeaks, dFiss, dStat, dCfeRtot, dCfeRabs, dTleRtot, dTleRabs, dCfeTau,
                       dMeshA, dMeshC, dMeshT, R, P.mesh);

    auto t1 = std::chrono::steady_clock::now();
    R.perf.histories = P.historiesPerBatch * P.batches;
    R.perf.elapsed_s = std::chrono::duration<float>(t1 - t0).count();

    (void)inactive;

    // Cleanup
    CUDA_CHECK(cudaFree(dLeaks));
    CUDA_CHECK(cudaFree(dFiss));
    CUDA_CHECK(cudaFree(dStat));
    CUDA_CHECK(cudaFree(dCfeRtot));
    CUDA_CHECK(cudaFree(dCfeRabs));
    CUDA_CHECK(cudaFree(dTleRtot));
    CUDA_CHECK(cudaFree(dTleRabs));
    CUDA_CHECK(cudaFree(dCfeTau));
    CUDA_CHECK(cudaFree(dMeshA));
    CUDA_CHECK(cudaFree(dMeshC));
    CUDA_CHECK(cudaFree(dMeshT));
    CUDA_CHECK(cudaFree(dBankCur));
    CUDA_CHECK(cudaFree(dBankNext));
    CUDA_CHECK(cudaFree(dNnext));

    H.freeGPU();
    return R;
}

RunOutputs runCuda(const Universe& U, const RunParams& P, int inactive) {
    if (P.src == SourceMode::Criticality) return runCriticalityCuda(U, P, inactive);
    return runExternalCuda(U, P);
}


#else

// If compiled without nvcc, use CPU code 
RunOutputs runExternalCuda(const Universe& U, const RunParams& P) { return runExternal(U,P); }
RunOutputs runCriticalityCuda(const Universe& U, const RunParams& P, int inactive) { return runCriticality(U,P,inactive); }
RunOutputs runCuda(const Universe& U, const RunParams& P, int inactive) {
    if (P.src == SourceMode::Criticality) return runCriticality(U,P,inactive);
    return runExternal(U,P);
}

#endif


bool compiledWithCuda() {
#ifdef __CUDACC__
    return true;
#else
    return false;
#endif
}
