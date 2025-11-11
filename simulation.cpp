#include "mc.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <map>
#include <array>
#include <algorithm>

using std::array; using std::vector;

constexpr double kB = 8.617333262e-11;
constexpr double MEV_TO_J = 1.602176634e-13;
constexpr double M_N = 1.674927498e-27;

size_t pickIndex(const std::vector<double>& cum, double u) {
    auto it = std::upper_bound(cum.begin(), cum.end(), u);
    if (it == cum.end()) return cum.empty() ? 0 : cum.size() - 1;
    return size_t(it - cum.begin());
}

std::vector<double> linspace(double a,double b,int n) {
    std::vector<double> x; x.reserve(n);
    if (n<=1) { if(n==1) x.push_back(a); return x; }
    double h=(b-a)/(n-1); for(int i=0;i<n;++i) x.push_back(a+i*h); return x;
}

std::vector<double> logspace(double ea,double eb,int n) {
    auto e=linspace(ea,eb,n); std::vector<double> x; x.reserve(n);
    for(double t:e) x.push_back(std::pow(10.0,t)); return x;
}

double valueInterp(const std::vector<std::pair<double,double>>& data, double target) {
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

double interpMT(const std::map<int, std::vector<std::pair<double,double>>>& mt, int code, double E) {
    auto it = mt.find(code); return it == mt.end() ? 0.0 : valueInterp(it->second, E);
}

double randomVal(float min=0.f, float max=1.f) {
    thread_local std::mt19937_64 gen{0x88123764};
    std::uniform_real_distribution<double> d(min,max);
    return d(gen);
}

double materialWeight(const Material& m, double E, const std::vector<int>& mts_total) {
    double micro = 0.0; for (int mt : mts_total) micro += interpMT(m.mt, mt, E);
    return std::max(0.0, m.proportion * m.rho * micro);
}

void buildMaterialCum(const std::vector<Material>& mats, double E, const std::vector<int>& mts_total, std::vector<double>& cum, double& total) {
    cum.clear(); cum.reserve(mats.size()); total = 0.0;
    for (const auto& m : mats) { total += materialWeight(m, E, mts_total); cum.push_back(total); }
}

void buildReactionCum(const Material& m, double E, const std::vector<int>& mts_sample, std::vector<int>& labels, std::vector<double>& cum, double& total) {
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

std::array<double,3> iso_dir() {
    const double u1 = randomVal(), u2 = randomVal();
    const double mu = 2.0*u1 - 1.0;
    const double phi = 2.0*M_PI*u2;
    const double s = std::sqrt(std::max(0.0, 1.0 - mu*mu));
    return { s*std::cos(phi), s*std::sin(phi), mu };
}

std::array<double,3> add(const std::array<double,3>& a, const std::array<double,3>& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

std::array<double,3> scale(const std::array<double,3>& a, double s) {
    return {s*a[0], s*a[1], s*a[2]};
}

double dot3(const array<double,3>& a, const array<double,3>& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
double norm2(const std::array<double,3>& a) { return a[0]*a[0]+a[1]*a[1]+a[2]*a[2]; }
std::array<double,3> normalize(const std::array<double,3>& a) {
    double n = std::sqrt(std::max(0.0, norm2(a)));
    return (n>0)? std::array<double,3>{a[0]/n, a[1]/n, a[2]/n} : std::array<double,3>{1,0,0};
}

double neutronSpeed(double E) { // non-relativistic
    const double E_J = E * MEV_TO_J;
    return std::sqrt(std::max(0.0, 2.0*E_J / M_N));
}

constexpr double EFG = 2e-4; // Free-gas upper bound (MeV)
void elasticScatter(double En, double A, double TK, double Efg, double& Eout) {
    std::array<double,3> VL{0,0,0};
    const bool use_free_gas = (En < Efg && A <= 10.0);
    if (use_free_gas) {
        const double T_mev = kB * TK;
        const double u = randomVal(), v = randomVal();
        const double Et = -T_mev * std::log(std::max(1e-32, u*v));
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

double sigmaTot(const std::vector<Material>& mats, double E, const std::vector<int>& mts_total) {
    double S = 0.0;
    for (const auto& m : mats) {
        double micro = 0.0;
        for (int mt : mts_total) micro += interpMT(m.mt, mt, E);
        S += std::max(0.0, m.proportion * m.rho * micro);
    }
    return S;
}

void tallyTrackToMesh(const Mesh3D& M, std::vector<double>& dst, const array<double,3>& p0, const array<double,3>& dHat, double segLen, double weightPerLength) {
    if (dst.empty() || segLen<=0) return;
    const array<double,3> p1 = add(p0, scale(dHat, segLen));
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
    array<double,3> p = add(p0, scale(dHat, tmin));
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

bool sampleReactionAtE(const Geometry& g, double E, const std::vector<int>& MTs_total, const std::vector<int>& MTs_sample, RxSample& out) {
    std::vector<double> matCum; double matTot=0.0;
    buildMaterialCum(g.mats, E, MTs_total, matCum, matTot);
    if (matTot<=0.0) return false;
    size_t imat = pickIndex(matCum, randomVal()*matTot);
    if (imat>=g.mats.size()) return false;

    std::vector<int> rxLabels; std::vector<double> rxCum; double rxTot=0.0;
    buildReactionCum(g.mats[imat], E, MTs_sample, rxLabels, rxCum, rxTot);
    if (rxTot<=0.0 || rxLabels.empty()) return false;
    size_t irx = pickIndex(rxCum, randomVal()*rxTot);
    if (irx>=rxLabels.size()) return false;

    out.mt = rxLabels[irx];
    out.matIndex = (int)imat;
    out.mat = &g.mats[imat];
    return true;
}

extern const Geometry* findGeometryAt(const Universe&, const std::array<double,3>& pWorld);
extern double macroscopicSigmaAt(const Universe&, const std::array<double,3>&, double, const std::vector<int>&, const Geometry** gOut);
extern double majorSigma(const Universe&, double, const std::vector<int>&);
extern double nearestCollision(const Neutron& n, const Universe& u);

FlightResult surfaceFlight(const Universe& u, Neutron& n, double E, const std::vector<int>& MTs_total, Mesh3D* meshTLE) {
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
        n.pos = add(n.pos, scale(n.dir, l));
    } else {
        fr.collided = false;
        fr.traveled = dSurf;
        if (meshTLE) {
            const double v = neutronSpeed(E);
            if (v>0) tallyTrackToMesh(*meshTLE, meshTLE->tle_density, n.pos, n.dir, dSurf, 1.0/v);
        }
        n.pos = add(n.pos, scale(n.dir, dSurf + EPS));
    }
    return fr;
}

FlightResult deltaFlight(const Universe& u, Neutron& n, double E, const std::vector<int>& MTs_total, double SigmaM, bool virtualAllowed) {
    FlightResult fr;
    if (SigmaM<=0.0) { fr.leaked=true; return fr; }
    const double l = -std::log(std::max(1e-32, randomVal()))/SigmaM;
    const array<double,3> newp = add(n.pos, scale(n.dir, l));

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

const std::vector<int>& rxList(bool inelastic) {
    const std::vector<int> A{2,18,102};
    const std::vector<int> B{2,4,18,102};
    return inelastic?B:A;
}

int rxCol(const std::vector<int>& MTs, int mt) {
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

void walkHistory(const Universe& U, const RunParams& P, int batch, const std::vector<int>& MTs_total, const std::vector<int>& MTs_sample, TallyBook& T, 
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
                if (!fr.collided) {
                    if (fr.leaked) break;
                    continue;
                }
                if (P.mesh && P.mesh->tle_density.size()) {
                    const double v = neutronSpeed(n.energy);
                    if (v>0) tallyTrackToMesh(*P.mesh, P.mesh->tle_density, n.pos, n.dir, 0.0, 0.0);
                }
                scoreCFE(T, batch, n.pos, n.energy, fr.SigmaLocal);
            } else {
                const double SigmaM = T.SigmaM;
                fr = deltaFlight(U, n, n.energy, MTs_total, SigmaM, true);
                if (fr.leaked) break;
                if (fr.virtualCollision) {
                    scoreCFE(T, batch, n.pos, n.energy, SigmaM);
                    continue;
                }
                scoreCFE(T, batch, n.pos, n.energy, fr.SigmaLocal);
            }
            RxSample rx;
            if (!fr.geom || !sampleReactionAtE(*fr.geom, n.energy, MTs_total, MTs_sample, rx)) {
                break;
            }
            int mi = T.matIndex(*rx.mat);
            int rj = rxCol(MTs_sample, rx.mt);
            if (mi>=0 && rj>=0) {
                T.ensure_batch(batch, (int)T.matNames.size(), (int)MTs_sample.size());
                T.statM[batch][mi][rj] += 1;
            }
            n.collisions++;
            recordCollision(col, n.collisions, n.energy);
            if (rx.mt==18) {
                const int nEmit = sampleMultiplicity(valueInterp(rx.mat->neutrons, n.energy));
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
                break;
            } else {
                break;
            }
        }
    }
}

RunOutputs runExternal(const Universe& U, const RunParams& P) {
    RunOutputs R; R.T.mesh = P.mesh; R.T.useTLE = (P.track==Tracking::Surface);
    R.T.useCFE = true; R.T.deltaMode = (P.track==Tracking::Delta);
    const auto& MTs = rxList(P.inelastic);
    const std::vector<int> MTs_total = MTs;
    const std::vector<int> MTs_sample = MTs;
    if (P.track==Tracking::Delta) R.T.SigmaM = majorSigma(U, P.sourceE, MTs_total);
    R.collisions.resize(P.batches);
    R.fissionChildren.assign(P.batches, 0);
    for (int b=0;b<P.batches;++b) {
        std::deque<Neutron> bank, next;
        for (int i=0;i<P.historiesPerBatch;++i) {
            Neutron n; n.pos = P.sourcePos; n.dir = iso_dir();
            n.energy = P.sourceE; n.vel = neutronSpeed(n.energy);
            n.isSource = true; R.T.ensure_batch(b, (int)R.T.matNames.size(), (int)MTs.size());
            bank.push_back(n);
        }
        FourTally four; four.started += P.historiesPerBatch;
        int births=0;
        walkHistory(U, P, b, MTs_total, MTs_sample, R.T, R.collisions[b], four, bank, next, births);
        R.fissionChildren[b] = births;
    }
    return R;
}

RunOutputs runCriticality(const Universe& U, const RunParams& P, int inactive=5) {
    RunOutputs R; R.T.mesh = P.mesh; R.T.useTLE = (P.track==Tracking::Surface);
    R.T.useCFE = true; R.T.deltaMode = (P.track==Tracking::Delta);
    const auto& MTs = rxList(P.inelastic);
    const std::vector<int> MTs_total = MTs;
    const std::vector<int> MTs_sample = MTs;
    std::deque<Neutron> bank, next;
    for (int i=0;i<P.historiesPerBatch;++i) {
        Neutron n; n.pos = P.sourcePos; n.dir = iso_dir();
        n.energy = P.sourceE; n.vel = neutronSpeed(n.energy);
        n.isSource = true; bank.push_back(n);
    }
    double kLast = 1.0;
    for (int b=0;b<P.batches;++b) {
        if (P.track==Tracking::Delta) R.T.SigmaM = majorSigma(U, P.sourceE, MTs_total);
        FourTally four; four.started += (int)bank.size();
        R.collisions.emplace_back(Collisions{});
        int births=0;
        RunParams Pcycle = P; Pcycle.src = SourceMode::Criticality;
        walkHistory(U, Pcycle, b, MTs_total, MTs_sample, R.T, R.collisions.back(), four, bank, next, births);
        double kcur = (double)births / std::max(1, P.historiesPerBatch);
        R.keff_history.push_back(kcur);
        std::deque<Neutron> newBank;
        if (!next.empty()) {
            double scale = (kLast>0)? 1.0/kLast : 1.0;
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
    return R;
}
