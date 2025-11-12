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

using std::array; using std::vector;
#define M_PI 3.14159265358979323846264338327950288

constexpr double kB = 8.617333262e-11;
constexpr double MEV_TO_J = 1.602176634e-13;
constexpr double M_N = 1.674927498e-27;

size_t pickIndex(const vector<double>& cum, double u) {
    auto it = std::upper_bound(cum.begin(), cum.end(), u);
    if (it == cum.end()) return cum.empty() ? 0 : cum.size() - 1;
    return size_t(it - cum.begin());
}

vector<double> linspace(double a,double b,int n) {
    vector<double> x; x.reserve(n);
    if (n<=1) { if(n==1) x.push_back(a); return x; }
    double h=(b-a)/(n-1); for(int i=0;i<n;++i) x.push_back(a+i*h); return x;
}

vector<double> logspace(double ea,double eb,int n) {
    auto e=linspace(ea,eb,n); vector<double> x; x.reserve(n);
    for(double t:e) x.push_back(std::pow(10.0,t)); return x;
}

double valueInterp(const vector<std::pair<double,double>>& data, double target) {
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

double interpMT(const std::map<int, vector<std::pair<double,double>>>& mt, int code, double E) {
    auto it = mt.find(code); return it == mt.end() ? 0.0 : valueInterp(it->second, E);
}

double randomVal(float min=0.f, float max=1.f) {
    thread_local std::mt19937_64 gen{0x9817981276389};
    std::uniform_real_distribution<double> d(min,max);
    return d(gen);
}

double materialWeight(const Material& m, double E, const vector<int>& mts_total) {
    double micro = 0.0; for (int mt : mts_total) micro += interpMT(m.mt, mt, E);
    return std::max(0.0, m.proportion * m.rho * micro);
}

void buildMaterialCum(const vector<Material>& mats, double E, const vector<int>& mts_total, vector<double>& cum, double& total) {
    cum.clear(); cum.reserve(mats.size()); total = 0.0;
    for (const auto& m : mats) { total += materialWeight(m, E, mts_total); cum.push_back(total); }
}

void buildReactionCum(const Material& m, double E, const vector<int>& mts_sample, vector<int>& labels, vector<double>& cum, double& total) {
    labels.clear(); cum.clear(); total = 0.0;
    for (int mt : mts_sample) {
        double w = interpMT(m.mt, mt, E);
        if (w <= 0.0) continue;
        total += w; labels.push_back(mt); cum.push_back(total);
    }
}

int sampleMultiplicity(double nuBar) {
    if (nuBar<=0) return 0;
    const int n = (int)std::floor(nuBar);
    const double frac = nuBar - n;
    return n + (randomVal() < frac ? 1 : 0);
}

array<double,3> iso_dir() {
    const double u1 = randomVal(), u2 = randomVal();
    const double mu = 2.0*u1 - 1.0;
    const double phi = 2.0*M_PI*u2;
    const double s = std::sqrt(std::max(0.0, 1.0 - mu*mu));
    return { s*std::cos(phi), s*std::sin(phi), mu };
}

inline array<double,3> add(const array<double,3>& a, const array<double,3>& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

array<double,3> scaleByC(const array<double,3>& a, double s) {
    return {s*a[0], s*a[1], s*a[2]};
}

double norm2(const array<double,3>& a) { return a[0]*a[0]+a[1]*a[1]+a[2]*a[2]; }

double neutronSpeed(double E) {
    const double E_J = E * MEV_TO_J;
    return std::sqrt(std::max(0.0, 2.0*E_J / M_N));
}

constexpr double EFG = 2e-4;
void elasticScatter(double En, double A, double TK, double Efg, double& Eout) {
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
    const array<double,3> VCM = scaleByC(add(vL, scaleByC(VL, A)), 1.0/(1.0 + A));
    const array<double,3> vC  = add(vL, scaleByC(VCM, -1.0));
    const double vC_mag = std::sqrt(std::max(0.0, norm2(vC)));
    const array<double,3> vCprime = scaleByC(iso_dir(), vC_mag);
    const array<double,3> vLprime = add(vCprime, VCM);
    Eout = 0.5 * norm2(vLprime);
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

double getDeltaE(const Material& m) {
    auto it = m.Qvals.find(4);
    if (it != m.Qvals.end()) return std::abs(it->second);
    return 0.5;
}

double sigmaTot(const vector<Material>& mats, double E, const vector<int>& mts_total) {
    double S = 0.0;
    for (const auto& m : mats) {
        double micro = 0.0;
        for (int mt : mts_total) micro += interpMT(m.mt, mt, E);
        S += std::max(0.0, m.proportion * m.rho * micro);
    }
    return S;
}

void tallyTrackToMesh(const Mesh3D& M, vector<double>& dst, const array<double,3>& p0, const array<double,3>& dHat, double segLen, double weightPerLength) {
    if (dst.empty() || segLen<=0) return;
    const array<double,3> p1 = add(p0, scaleByC(dHat, segLen));
    array<double,3> t0={0,0,0}, t1={segLen,segLen,segLen};
    array<double,3> inv = { dHat[0]!=0?1.0/dHat[0]:1e300, dHat[1]!=0?1.0/dHat[1]:1e300, dHat[2]!=0?1.0/dHat[2]:1e300 };
    double tmin=0.0, tmax=segLen;
    auto clipAxis=[&](int a, double minv, double maxv) {
        double tA = (minv - p0[a]) * inv[a];
        double tB = (maxv - p0[a]) * inv[a];
        if (tA>tB) std::swap(tA,tB);
        tmin = std::max(tmin, tA);
        tmax = std::min(tmax, tB);
    };
    clipAxis(0,M.pmin[0],M.pmax[0]); clipAxis(1,M.pmin[1],M.pmax[1]); clipAxis(2,M.pmin[2],M.pmax[2]);
    if (tmax<=tmin) return;

    auto clampi=[&](int i,int n) { return std::min(std::max(i,0), n-1); };
    const double dx=(M.pmax[0]-M.pmin[0])/M.nx, dy=(M.pmax[1]-M.pmin[1])/M.ny, dz=(M.pmax[2]-M.pmin[2])/M.nz;
    array<double,3> p = add(p0, scaleByC(dHat, tmin));
    int i = clampi(int((p[0]-M.pmin[0])/dx), M.nx);
    int j = clampi(int((p[1]-M.pmin[1])/dy), M.ny);
    int k = clampi(int((p[2]-M.pmin[2])/dz), M.nz);
    int stepX = (dHat[0]>0)?+1: (dHat[0]<0?-1:0);
    int stepY = (dHat[1]>0)?+1: (dHat[1]<0?-1:0);
    int stepZ = (dHat[2]>0)?+1: (dHat[2]<0?-1:0);

    double nextX = (stepX>0)? (M.pmin[0] + (i+1)*dx) : (M.pmin[0] + i*dx);
    double nextY = (stepY>0)? (M.pmin[1] + (j+1)*dy) : (M.pmin[1] + j*dy);
    double nextZ = (stepZ>0)? (M.pmin[2] + (k+1)*dz) : (M.pmin[2] + k*dz);

    double t = tmin;
    double tMaxX = (stepX==0)? 1e300 : (nextX - p[0]) * inv[0];
    double tMaxY = (stepY==0)? 1e300 : (nextY - p[1]) * inv[1];
    double tMaxZ = (stepZ==0)? 1e300 : (nextZ - p[2]) * inv[2];
    const double tDeltaX = (stepX==0)? 1e300 : std::abs(dx * inv[0]);
    const double tDeltaY = (stepY==0)? 1e300 : std::abs(dy * inv[1]);
    const double tDeltaZ = (stepZ==0)? 1e300 : std::abs(dz * inv[2]);

    while (t < tmax && i>=0 && j>=0 && k>=0 && i<M.nx && j<M.ny && k<M.nz) {
        double tNext = std::min({tMaxX, tMaxY, tMaxZ, tmax});
        double seg = std::max(0.0, tNext - t);
        if (seg>0) {
            dst[M.idx(i,j,k)] += weightPerLength * seg;
        }
        t = tNext;
        if (t>=tmax) break;
        if (tNext==tMaxX) { i += stepX; tMaxX += tDeltaX; }
        else if (tNext==tMaxY) { j += stepY; tMaxY += tDeltaY; }
        else { k += stepZ; tMaxZ += tDeltaZ; }
    }
}

bool sampleReactionAtE(const Geometry& g, double E, const vector<int>& MTs_total, const vector<int>& MTs_sample, RxSample& out) {
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

FlightResult surfaceFlight(const Universe& u, Neutron& n, double E, const vector<int>& MTs_total, Mesh3D* meshTLE) {
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
        n.pos = add(n.pos, scaleByC(n.dir, l));
    } else {
        fr.collided = false;
        fr.traveled = dSurf;
        fr.geom = g0;
        if (meshTLE) {
            const double v = neutronSpeed(E);
            if (v>0) tallyTrackToMesh(*meshTLE, meshTLE->tle_density, n.pos, n.dir, dSurf, 1.0/v);
        }
        n.pos = add(n.pos, scaleByC(n.dir, dSurf + EPS));
    }
    return fr;
}

FlightResult deltaFlight(const Universe& u, Neutron& n, double E, const vector<int>& MTs_total, double SigmaM, bool virtualAllowed) {
    FlightResult fr;
    if (SigmaM<=0.0) { fr.leaked=true; return fr; }
    const double l = -std::log(std::max(1e-32, randomVal()))/SigmaM;
    const array<double,3> newp = add(n.pos, scaleByC(n.dir, l));

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

int rxCol(const vector<int>& MTs, int mt) {
    auto it = std::find(MTs.begin(), MTs.end(), mt);
    return (it==MTs.end())? -1 : int(it - MTs.begin());
}

void recordCollision(Collisions& c, int collIdx, double E) {
    if (collIdx >= (int)c.num.size()) { c.num.resize(collIdx + 1, 0); c.sumEnergy.resize(collIdx + 1, 0.0); }
    c.num[collIdx] += 1; c.sumEnergy[collIdx] += E;
}
bool isFuelSym(const std::string& s) { std::string t; for(unsigned char c: s) { if(std::isalnum(c)) t.push_back((char)std::toupper(c)); } return t=="U235"||t=="U238"||t=="PU239"||t=="PU241"; }
bool isThermalE(double E, const Material& m) { return E <= 3.0 * kB * m.T; }

void scoreCFE(TallyBook& T, int batch, const array<double,3>& pos, double E, double SigmaRef) {
    if (!T.mesh || !T.useCFE) return;
    if (!T.mesh->inside(pos)) return;
    const double v = neutronSpeed(E); if (v<=0 || SigmaRef<=0) return;
    const int nx=T.mesh->nx, ny=T.mesh->ny, nz=T.mesh->nz;
    const double dx=(T.mesh->pmax[0]-T.mesh->pmin[0])/nx;
    const double dy=(T.mesh->pmax[1]-T.mesh->pmin[1])/ny;
    const double dz=(T.mesh->pmax[2]-T.mesh->pmin[2])/nz;
    int i = std::min(std::max(int((pos[0]-T.mesh->pmin[0])/dx),0), nx-1);
    int j = std::min(std::max(int((pos[1]-T.mesh->pmin[1])/dy),0), ny-1);
    int k = std::min(std::max(int((pos[2]-T.mesh->pmin[2])/dz),0), nz-1);
    T.mesh->cfe_density[T.mesh->idx(i,j,k)] += 1.0/(v*SigmaRef);
}

double microXS(const Material& m, int mt, double E) {
    auto it = m.mt.find(mt);
    return (it==m.mt.end()) ? 0.0 : valueInterp(it->second, E);
}

double macroXS_comp(const Material& m, int mt, double E) {
    return m.rho * m.proportion * microXS(m, mt, E);
}

double macroXS_comp_total(const Material& m, double E, const std::vector<int>& mts_total) {
    double s=0.0; for(int mt: mts_total) s += macroXS_comp(m, mt, E); return s;
}

double macroXS_comp_abs(const Material& m, double E) {
    return macroXS_comp(m, 102, E) + macroXS_comp(m, 18, E);
}

void scoreCFE_rx_per_material(TallyBook& T, int batch,
                                            int mi, const Material& m, double E,
                                            const std::vector<int>& mts_total, double SigmaRef) {
    if (SigmaRef <= 0.0) return;
    const double Stot = macroXS_comp_total(m, E, mts_total);
    const double Sabs = macroXS_comp_abs(m, E);
    const int mit = T.matIndex(m);
    if (mit >= 0) {
        T.ensureBatchAll(batch, (int)T.matNames.size());
        T.cfe_Rtot[batch][mit] += Stot / SigmaRef;
        T.cfe_Rabs[batch][mit] += Sabs / SigmaRef;
    }
}

void scoreTLE_segment_per_geom(TallyBook& T, int batch, const Geometry& g0, double E, double segLen, const std::vector<int>& mts_total) {
    if (segLen <= 0.0) return;
    T.ensureBatchAll(batch, (int)T.matNames.size());

    const bool isMixture = (g0.mats.size() > 1);
    double Stot_mix = 0.0, Sabs_mix = 0.0;
    if (isMixture) {
        for (const Material& m : g0.mats) {
            Stot_mix += macroXS_comp_total(m, E, mts_total);
            Sabs_mix += macroXS_comp_abs(m, E);
        }
        if (Stot_mix <= 0.0 && Sabs_mix <= 0.0) return;
    }

    const double L = segLen * 1;

    if (!isMixture) {
        const Material& m = g0.mats[0];
        const int mi = T.matIndex(m);
        const double Stot = macroXS_comp_total(m, E, mts_total);
        const double Sabs = macroXS_comp_abs(m, E);
        T.tle_Rtot[batch][mi] += L * Stot;
        T.tle_Rabs[batch][mi] += L * Sabs;
        return;
    }

    for (const Material& m : g0.mats) {
        const int mi = T.matIndex(m);
        const double Stot_c = macroXS_comp_total(m, E, mts_total);
        const double Sabs_c = macroXS_comp_abs(m, E);
        if (Stot_mix > 0.0) T.tle_Rtot[batch][mi] += L * Stot_c;
        if (Sabs_mix > 0.0) T.tle_Rabs[batch][mi] += L * Sabs_c;
    }
}

void scoreCFE_density_global(TallyBook& T, int batch, double E, double SigmaRef) {
    if (SigmaRef <= 0.0) return;
    const double v = neutronSpeed(E);
    if (v>0.0) { T.ensureBatchAll(batch, (int)T.matNames.size()); T.cfe_global_time[batch] += 1.0/(v*SigmaRef); }
}

int pickMaterialIndexAtE(const Geometry& g, double E, const std::vector<int>& MTs_total) {
    std::vector<double> matCum; double matTot=0.0;
    buildMaterialCum(g.mats, E, MTs_total, matCum, matTot);
    if (matTot<=0.0) return -1;
    size_t imat = pickIndex(matCum, randomVal()*matTot);
    if (imat>=g.mats.size()) return -1;
    return (int)imat;
}

double defaultNuBar(const Material& m) {
    std::string s; for(unsigned char c: m.sym) { if(std::isalnum(c)) s.push_back((char)std::toupper(c)); }
    if (s=="U235") return 2.43;
    if (s=="U238") return 2.50;
    if (s=="PU239") return 2.90;
    return 2.40;
}

int vindex(const Mesh3D& M, double x, double y, double z) {
    auto toI = [&](double X, double a, double b, int n) {
        const double t = (X - a) / (b - a);
        int i = (int)std::floor(t * n);
        if (i < 0) i = 0;
        if (i >= n) i = n - 1;
        return i;
    };
    const int ix = toI(x, M.pmin[0], M.pmax[0], M.nx);
    const int iy = toI(y, M.pmin[1], M.pmax[1], M.ny);
    const int iz = toI(z, M.pmin[2], M.pmax[2], M.nz);
    return (iz * M.ny + iy) * M.nx + ix;
}

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
                    scoreTLE_segment_per_geom(T, batch, *fr.geom, n.energy, fr.traveled, MTs_total);
                }
                if (!fr.collided) {
                    if (fr.leaked) { T.ensureBatchAll(batch, (int)T.matNames.size()); ++T.leaks[batch]; if (T.timeHist) T.timeHist->add(n.time); break; }
                    continue;
                }
                scoreCFE_density_global(T, batch, n.energy, fr.SigmaLocal);
            } else {
                const double SigmaM = T.SigmaM;
                fr = deltaFlight(U, n, n.energy, MTs_total, SigmaM, true);
                const double seg_m = fr.traveled * 0.01;
                n.time += seg_m / n.vel;
                if (fr.leaked) { T.ensureBatchAll(batch, (int)T.matNames.size()); ++T.leaks[batch]; if (T.timeHist) T.timeHist->add(n.time); break; }
                if (fr.virtualCollision) {
                    scoreCFE_density_global(T, batch, n.energy, SigmaM);
                    if (fr.geom) {
                        int mi = pickMaterialIndexAtE(*fr.geom, n.energy, MTs_total);
                        auto kjashd = fr.geom->mats[mi];
                        if (mi>=0) scoreCFE_rx_per_material(T, batch, mi, fr.geom->mats[mi], n.energy, MTs_total, SigmaM);
                    }
                    continue;
                }
                scoreCFE_density_global(T, batch, n.energy, fr.SigmaLocal);
            }
            RxSample rx;
            if (!fr.geom || !sampleReactionAtE(*fr.geom, n.energy, MTs_total, MTs_sample, rx)) break;

            int mi = T.matIndex(*rx.mat);
            int rj = rxCol(MTs_sample, rx.mt);
            if (mi>=0 && rj>=0) {
                T.ensureBatchAll(batch, (int)T.matNames.size());
                T.statM[batch][mi][rj] += 1;
                scoreCFE_rx_per_material(T, batch, mi, *rx.mat, n.energy, MTs_total, (P.track==Tracking::Surface)? fr.SigmaLocal : T.SigmaM);
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
                if (T.timeHist) T.timeHist->add(n.time);
                break;
            } else if (rx.mt==2) {
                double A = rx.mat->a>0? (double)rx.mat->a : std::max(1, rx.mat->z);
                double Eout; elasticScatter(n.energy, A, rx.mat->T, EFG, Eout);
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
                if (T.timeHist) T.timeHist->add(n.time);
                break;
            } else {
                if (T.timeHist) T.timeHist->add(n.time);
                break;
            }
        }
    }
}

RunOutputs runExternal(const Universe& U, const RunParams& P) {
    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();
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
    for (int b=0;b<P.batches;++b) {
        std::deque<Neutron> bank, next;
        for (int i=0;i<P.historiesPerBatch;++i) {
            Neutron n; n.pos = P.sourcePos; n.dir = iso_dir();
            n.energy = P.sourceE; n.vel = neutronSpeed(n.energy);
            n.isSource = true; R.T.ensureBatch(b, (int)R.T.matNames.size(), (int)MTs.size());
            bank.push_back(n);
        }
        totalHist += P.historiesPerBatch;
        FourTally four; four.started += P.historiesPerBatch;
        int births=0;
        walkHistory(U, P, b, MTs_total, MTs_sample, R.T, R.collisions[b], four, bank, next, births);
        R.fissionChildren[b] = births;
    }
    auto t1 = clk::now();
    R.perf.elapsed_s = std::chrono::duration<double>(t1 - t0).count();
    R.perf.histories = totalHist;
    auto writeTimeHistCsv = [](const TimeHist& H, const char* path){
        std::ofstream os(path);
        os << "t_center_s,count\n";
        for (int k=0;k<H.nbins;++k){
            double t0 = k*H.dt, t1 = (k+1)*H.dt;
            os << std::scientific << 0.5*(t0+t1) << "," << H.counts[k] << "\n";
        }
    };
    writeTimeHistCsv(lifeHist, (P.src==SourceMode::External? "output/timehist_external.csv" : "output/timehist_critical.csv"));
    R.T.timeHist = nullptr;
    return R;
}

void printFour(const FourTally&F) {
    std::cout << F.fissionBirthsTotal << " BirthsTotal" << "\n"
    << F.fissionBirthsThermal << " BirthsThermal" << "\n"
    << F.absThTotal << " AbsThermal" << "\n"
    << F.absThFuel << " AbsFuel" << "\n"
    << F.started << " Started" << "\n"
    << F.reachedThermal << " ReachThermal" << "\n"
    << F.resAbsBeforeThermal << " AbsBefTherm" << "\n";
}
RunOutputs runCriticality(const Universe& U, const RunParams& P, int inactive=5) {
    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();
    RunOutputs R; R.T.mesh = P.mesh; R.T.useTLE = (P.track==Tracking::Surface);
    R.T.useCFE = true; R.T.deltaMode = (P.track==Tracking::Delta);
    TimeHist lifeHist(1e-6, 200);
    R.T.timeHist = &lifeHist;
    const auto MTs = P.inelastic ? vector<int>{2,4,18,102} : vector<int>{2,18,102};
    R.T.Rcols = (int)MTs.size();
    const vector<int> MTs_total = MTs;
    const vector<int> MTs_sample = MTs;
    std::deque<Neutron> bank, next;
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
        //printFour(four);
        const int Nk1 = (int)next.size();
        totalHist += P.historiesPerBatch;
        //double kcur = (double)births / std::max(1, P.historiesPerBatch);
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
    auto writeTimeHistCsv = [](const TimeHist& H, const char* path){
        std::ofstream os(path);
        os << "t_center_s,count\n";
        for (int k=0;k<H.nbins;++k){
            double t0 = k*H.dt, t1 = (k+1)*H.dt;
            os << std::scientific << 0.5*(t0+t1) << "," << H.counts[k] << "\n";
        }
    };
    writeTimeHistCsv(lifeHist, (P.src==SourceMode::External? "output/timehist_external.csv" : "output/timehist_critical.csv"));
    R.T.timeHist = nullptr;

    return R;
}
