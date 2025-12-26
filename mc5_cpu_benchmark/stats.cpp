#include "mc.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <array>
#include <string>
#include <deque>


// --- COMMON FUNCTIONS ---

std::vector<RegionMap> defaultRegions() {
    return {
        {"Water",     {"H-1","O-16"}},
        {"Air",       {"N-14"}},
        {"Container", {"Fe-56"}}
    };
}


// --- STATS FUNCTIONS ---

StatsOut computeStats(const vector<vector<vector<int>>>& statM) {
    const size_t I = statM.size();
    const size_t M = I? statM[0].size() : 0;
    const size_t R = (M? statM[0][0].size() : 0);
    StatsOut out;
    out.mean.assign(M, vector<double>(R, 0.0));
    out.relErr.assign(M, vector<double>(R, 0.0));
    out.sum.assign(M, vector<int>(R, 0));

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

std::string mtLabel(int mt) {
    switch (mt) {
        case 2:   return "MT2 Elastic";
        case 4:   return "MT4 Inelastic";
        case 18:  return "MT18 Fission";
        case 102: return "MT102 Capture";
        default:  return "MT" + std::to_string(mt);
    }
}

void printStatsOut(const StatsOut& S, const vector<std::string>& matNames, const vector<int>& MTs, std::ostream& os = std::cout) {
    const size_t M = S.sum.size();
    if (!M) { os << "(no data)\n"; return; }
    const size_t R = S.sum[0].size();

    vector<long long> rowTot(M,0), colTot(R,0);
    long long grand = 0;
    for (size_t i=0;i<M;++i)
        for (size_t j=0;j<R;++j) {
            long long v = S.sum[i][j];
            rowTot[i]+=v; colTot[j]+=v; grand+=v;
        }
    auto pct=[&](long long x) { return grand? 100.0*double(x)/double(grand) : 0.0; };

    size_t wName = 4;
    for (size_t i=0;i<std::min(M,matNames.size());++i) wName = std::max(wName, matNames[i].size());
    vector<size_t> wCol(R, 14);
    for (size_t j=0;j<R;++j) {
        wCol[j] = std::max(wCol[j], mtLabel(MTs[j]).size());
        for (size_t i=0;i<M;++i) {
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
    for (size_t i=0;i<M;++i) {
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

    for (size_t i=0;i<M;++i) {
        const std::string& name = (i<matNames.size()? matNames[i] : ("mat"+std::to_string(i)));
        os << std::left << std::setw(int(wName)) << name << "  ";
        for (size_t j=0;j<R;++j) {
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

ScalarStat statsFromBatches(const std::vector<double>& x, double elapsed_s) {
    ScalarStat s; s.batches = (int)x.size();
    if (s.batches==0) return s;
    double mu=0; for(double v:x) mu+=v; mu/=x.size();
    double var=0; for(double v:x) { double d=v-mu; var+=d*d; }
    var = (x.size()>1)? var/(x.size()-1) : 0.0;
    double se = (x.size()>0)? std::sqrt(var / x.size()) : 0.0;
    s.mean = mu;
    s.stddev = std::sqrt(var);
    s.relErr = (mu!=0.0)? (se/std::abs(mu)) : 0.0;
    s.fom = (elapsed_s>0.0 && s.relErr>0.0)? 1.0/(s.relErr*s.relErr*elapsed_s) : 0.0;
    return s;
}

void estimateNeutronDensity1W(const RunParams& P, const RunOutputs& R, double mixtureVolume_m3, double meanNuBar = 2.43) {
    const double Ef = 200.0e6 * 1.602176634e-19;
    const double P_w = 1.0;
    const double Fdot = P_w / Ef;
    const double Qn   = Fdot * meanNuBar;
    double sumTime = 0.0;
    for (double x : R.T.cfe_global_time) sumTime += x;
    const double Nhist = double(P.historiesPerBatch) * double(P.batches);
    const double tau_per_hist = (Nhist>0)? (sumTime/Nhist) : 0.0;
    const double N_neutrons = Qn * tau_per_hist;
    const double n_density  = (mixtureVolume_m3>0.0)? (N_neutrons / mixtureVolume_m3) : 0.0;

    std::cout << "\n[1 W] fission rate: " << std::scientific << Fdot
              << "  neutrons emitted/s: " << Qn
              << "  <lifetime>: " << tau_per_hist << " s"
              << "  <n> in system: " << N_neutrons
              << "  avg neutron density: " << n_density << " 1/m^3\n";
    std::ofstream os("output/neutron_density_1W.csv");
    os << "Fdot_n_per_s,nu_bar,tau_s,N_neutrons,volume_m3,n_density_per_m3\n"
       << Fdot << "," << meanNuBar << "," << tau_per_hist << "," << N_neutrons << "," << mixtureVolume_m3 << "," << n_density << "\n";
}


// --- BATCH FUNCTIONS ---


double normPerBatch(const RunParams& P) {
    return (P.historiesPerBatch>0)? (P.sourceRate / double(P.historiesPerBatch)) : 0.0;
}

std::map<std::string,int> matIndexMap(const TallyBook& T) {
    std::map<std::string,int> m2i; for (int i=0;i<(int)T.matNames.size();++i) m2i[T.matNames[i]] = i; return m2i;
}

int colIdx(const std::vector<int>& mts, int mt) {
    for (int j=0;j<(int)mts.size();++j) if (mts[j]==mt) return j; return -1;
}

BatchMetrics buildBatchMetrics(const RunParams& P, const RunOutputs& R, const std::vector<int>& MTs) {
    const auto regs = defaultRegions();
    const auto m2i = matIndexMap(R.T);
    const double norm = normPerBatch(P);
    const int bN = (int)R.T.statM.size();
    BatchMetrics B;
    B.leak_nps.resize(bN,0);
    B.abs_water_nps.resize(bN,0);
    B.abs_steel_nps.resize(bN,0);
    B.tle_abs_water_nps.resize(bN,0);
    B.cfe_abs_water_nps.resize(bN,0);
    B.cfe_tau_s.resize(bN,0);
    const int jCap = colIdx(MTs,102);
    const int jFis = colIdx(MTs,18);
    auto sum_region_analog = [&](int b, const std::vector<std::string>& mats)->int {
        long long cnt = 0;
        for (const auto& name : mats) {
            auto it = m2i.find(name); if (it==m2i.end()) continue;
            int i = it->second;
            if (b >= (int)R.T.statM.size()) continue;
            if (i >= (int)R.T.statM[b].size()) continue;
            if (jCap>=0 && jCap<(int)R.T.statM[b][i].size()) cnt += R.T.statM[b][i][jCap];
            if (jFis>=0 && jFis<(int)R.T.statM[b][i].size()) cnt += R.T.statM[b][i][jFis];
        }
        return (int)cnt;
    };

    for (int b=0;b<bN;++b) {
        long long l = (b<(int)R.T.leaks.size() ? R.T.leaks[b] : 0);
        B.leak_nps[b] = norm * double(l);
        int absW = sum_region_analog(b, regs[0].mats);
        int absFe= sum_region_analog(b, regs[2].mats);
        B.abs_water_nps[b] = norm * double(absW);
        B.abs_steel_nps[b] = norm * double(absFe);
        double tleW=0, cfeW=0;
        for (const auto& name : regs[0].mats) {
            auto it = m2i.find(name); if (it==m2i.end()) continue;
            int mi = it->second;
            if (b<(int)R.T.tle_Rabs.size() && mi<(int)R.T.tle_Rabs[b].size()) tleW += R.T.tle_Rabs[b][mi];
            if (b<(int)R.T.cfe_Rabs.size() && mi<(int)R.T.cfe_Rabs[b].size()) cfeW += R.T.cfe_Rabs[b][mi];
        }
        B.tle_abs_water_nps[b] = norm * tleW;
        B.cfe_abs_water_nps[b] = norm * cfeW;
        if (b<(int)R.T.cfe_global_time.size()) B.cfe_tau_s[b] = R.T.cfe_global_time[b];
    }
    return B;
}

void printFoms(const RunParams& P, const RunOutputs& R, const std::vector<int>& MTs, const std::string& tag) {
    const auto B = buildBatchMetrics(P, R, MTs);
    auto pr = [&](const char* name, const std::vector<double>& x) {
        auto st = statsFromBatches(x, R.perf.elapsed_s);
        std::cout << std::left << std::setw(20) << name
                  << " mean=" << std::scientific << st.mean
                  << " relErr=" << st.relErr
                  << " FOM=" << st.fom
                  << " time[s]=" << R.perf.elapsed_s
                  << " batches=" << st.batches << "\n";
        return st;
    };
    std::cout << "\n=== FOM (normalized to " << std::scientific << P.sourceRate << " n/s) ["<<tag<<"] ===\n";
    auto s1 = pr("leak_rate_nps",        B.leak_nps);
    auto s2 = pr("abs_water_nps(analog)",B.abs_water_nps);
    auto s3 = pr("abs_steel_nps(analog)",B.abs_steel_nps);
    auto s4 = pr("abs_water_nps(TLE)",   B.tle_abs_water_nps);
    auto s5 = pr("abs_water_nps(CFE)",   B.cfe_abs_water_nps);
    auto s6 = pr("cfe_tau_per_hist[s]",  B.cfe_tau_s);

    std::ofstream os("output/fom_"+tag+".csv");
    os << "metric,mean,relerr,fom,time_s,batches\n";
    os << "leak,"           << s1.mean << "," << s1.relErr << "," << s1.fom << "," << R.perf.elapsed_s << "," << s1.batches << "\n";
    os << "abs_water_analog,"<< s2.mean << "," << s2.relErr << "," << s2.fom << "," << R.perf.elapsed_s << "," << s2.batches << "\n";
    os << "abs_steel_analog,"<< s3.mean << "," << s3.relErr << "," << s3.fom << "," << R.perf.elapsed_s << "," << s3.batches << "\n";
    os << "abs_water_TLE,"   << s4.mean << "," << s4.relErr << "," << s4.fom << "," << R.perf.elapsed_s << "," << s4.batches << "\n";
    os << "abs_water_CFE,"   << s5.mean << "," << s5.relErr << "," << s5.fom << "," << R.perf.elapsed_s << "," << s5.batches << "\n";
    os << "cfe_tau_per_hist,"<< s6.mean << "," << s6.relErr << "," << s6.fom << "," << R.perf.elapsed_s << "," << s6.batches << "\n";
}


// --- MESH FUNCTIONS ---

Mesh3D makeMeshFromUniverse(const Universe& U, int nx,int ny,int nz) {
    Universe tmp = U;
    auto ext = boundingBox(tmp);
    Mesh3D M; M.nx=nx; M.ny=ny; M.nz=nz;
    M.pmin = {-0.5*ext[0], -0.5*ext[1], -0.5*ext[2]};
    M.pmax = {+0.5*ext[0], +0.5*ext[1], +0.5*ext[2]};
    M.zero();
    return M;
}

void accumulateByRegion(const TallyBook& T, const StatsOut& S, std::map<std::string,int>& name2idx, std::vector<double>& analogColl, std::vector<double>& analogAbs) {
    const auto regs = defaultRegions();
    analogColl.assign(regs.size(), 0.0);
    analogAbs.assign(regs.size(), 0.0);
    auto findCol = [&](int mt)->int {
        std::map<int,int> col;
        if (col.empty()) {
            col[2]=0; col[4]=1; col[18]=2; col[102]=3;
        }
        return col.count(mt)? col[mt] : -1;
    };
    const int M = (int)T.matNames.size();
    for (size_t r=0;r<regs.size();++r) {
        for (const auto& mname : regs[r].mats) {
            auto it = name2idx.find(mname);
            if (it==name2idx.end()) continue;
            int i = it->second; if (i<0 || i>=M) continue;
            double tot = 0.0;
            for (int mt : {2,4,18,102}) { int j=findCol(mt); if (j>=0 && j<(int)S.sum[i].size()) tot += S.sum[i][j]; }
            analogColl[r] += tot;
            double abs = 0.0;
            int jc = findCol(102); if (jc>=0) abs += S.sum[i][jc];
            int jf = findCol(18);  if (jf>=0) abs += S.sum[i][jf];
            analogAbs[r] += abs;
        }
    }
}

void printAndStoreRates(const RunParams& P, const RunOutputs& R, const StatsOut& S, const std::string& tag) {
    const double srcRate = P.sourceRate;
    const double Nhist   = double(P.historiesPerBatch) * double(P.batches);
    const double norm    = (Nhist>0)? (srcRate/Nhist) : 0.0;
    std::map<std::string,int> m2i;
    for (int i=0;i<(int)R.T.matNames.size();++i) m2i[R.T.matNames[i]] = i;
    const auto regs = defaultRegions();
    std::vector<double> analogColl, analogAbs;
    accumulate_by_region(R.T, S, m2i, analogColl, analogAbs);
    std::vector<double> cfeColl(regs.size(),0.0), cfeAbs(regs.size(),0.0),
                        tleColl(regs.size(),0.0), tleAbs(regs.size(),0.0);

    for (size_t b=0;b<R.T.cfe_Rtot.size();++b) {
        for (size_t r=0;r<regs.size();++r) {
            for (const auto& mname : regs[r].mats) {
                auto it = m2i.find(mname); if (it==m2i.end()) continue;
                int mi = it->second;
                if (mi<(int)R.T.cfe_Rtot[b].size()) {
                    cfeColl[r] += R.T.cfe_Rtot[b][mi];
                    cfeAbs [r] += R.T.cfe_Rabs [b][mi];
                    tleColl[r] += R.T.tle_Rtot[b][mi];
                    tleAbs [r] += R.T.tle_Rabs [b][mi];
                }
            }
        }
    }
    long long leakCounts = 0;
    for (int x : R.T.leaks) leakCounts += x;
    const double leakRate = norm * double(leakCounts);
    std::cout << "\n=== Rates normalized to " << std::scientific << srcRate << " n/s (" << tag << ") ===\n";
    std::cout << "Total leak rate [n/s]: " << leakRate << "\n";
    auto prRow = [&](const char* label, double aC, double aA, double cC, double cA, double tC, double tA) {
        std::cout << std::left << std::setw(12) << label
                  << " analog Coll: " << std::scientific << aC*norm
                  << "  analog Abs: " << aA*norm
                  << "  CFE Coll: "  << cC*norm
                  << "  CFE Abs: "   << cA*norm
                  << "  TLE Coll: "  << tC*norm
                  << "  TLE Abs: "   << tA*norm << "\n";
    };
    for (size_t r=0;r<regs.size();++r) {
        prRow(regs[r].label, analogColl[r], analogAbs[r], cfeColl[r], cfeAbs[r], tleColl[r], tleAbs[r]);
    }
    std::ofstream os("output/rates_"+tag+".csv");
    os << "region,analog_coll,analog_abs,cfe_coll,cfe_abs,tle_coll,tle_abs,units\n";
    for (size_t r=0;r<regs.size();++r) {
        os << regs[r].label << ","
           << (analogColl[r]*norm) << "," << (analogAbs[r]*norm) << ","
           << (cfeColl[r]*norm)    << "," << (cfeAbs[r]*norm)    << ","
           << (tleColl[r]*norm)    << "," << (tleAbs[r]*norm)    << ",1/s\n";
    }
    std::ofstream os2("output/leak_"+tag+".csv");
    os2 << "leak_rate_n_per_s\n" << leakRate << "\n";
}

void tallyTrackToMesh(const Mesh3D& M, vector<double>& dst, const array<double,3>& p0, const array<double,3>& dHat, double segLen, double weightPerLength) {
    if (dst.empty() || segLen<=0) return;
    const array<double,3> p1 = add3(p0, scaleByC(dHat, segLen));
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
    array<double,3> p = add3(p0, scaleByC(dHat, tmin));
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

double macroXSComp(const Material& m, int mt, double E) {
    return m.rho * m.proportion * microXS(m, mt, E);
}

double macroXSCompTotal(const Material& m, double E, const std::vector<int>& mts_total) {
    double s=0.0; for(int mt: mts_total) s += macroXSComp(m, mt, E); return s;
}

double macroXSCompAbs(const Material& m, double E) {
    return macroXSComp(m, 102, E) + macroXSComp(m, 18, E);
}

void scoreCFERxPerMaterial(TallyBook& T, int batch, int mi, const Material& m, double E, const std::vector<int>& mts_total, double SigmaRef) {
    if (SigmaRef <= 0.0) return;
    const double Stot = macroXSCompTotal(m, E, mts_total);
    const double Sabs = macroXSCompAbs(m, E);
    const int mit = T.matIndex(m);
    if (mit >= 0) {
        T.ensureBatchAll(batch, (int)T.matNames.size());
        T.cfe_Rtot[batch][mit] += Stot / SigmaRef;
        T.cfe_Rabs[batch][mit] += Sabs / SigmaRef;
    }
}

void scoreTLESegmentPerGeom(TallyBook& T, int batch, const Geometry& g0, double E, double segLen, const std::vector<int>& mts_total) {
    if (segLen <= 0.0) return;
    T.ensureBatchAll(batch, (int)T.matNames.size());

    const bool isMixture = (g0.mats.size() > 1);
    double Stot_mix = 0.0, Sabs_mix = 0.0;
    if (isMixture) {
        for (const Material& m : g0.mats) {
            Stot_mix += macroXSCompTotal(m, E, mts_total);
            Sabs_mix += macroXSCompAbs(m, E);
        }
        if (Stot_mix <= 0.0 && Sabs_mix <= 0.0) return;
    }

    const double L = segLen * 1;

    if (!isMixture) {
        const Material& m = g0.mats[0];
        const int mi = T.matIndex(m);
        const double Stot = macroXSCompTotal(m, E, mts_total);
        const double Sabs = macroXSCompAbs(m, E);
        T.tle_Rtot[batch][mi] += L * Stot;
        T.tle_Rabs[batch][mi] += L * Sabs;
        return;
    }

    for (const Material& m : g0.mats) {
        const int mi = T.matIndex(m);
        const double Stot_c = macroXSCompTotal(m, E, mts_total);
        const double Sabs_c = macroXSCompAbs(m, E);
        if (Stot_mix > 0.0) T.tle_Rtot[batch][mi] += L * Stot_c;
        if (Sabs_mix > 0.0) T.tle_Rabs[batch][mi] += L * Sabs_c;
    }
}

void scoreCFEDensityGlobal(TallyBook& T, int batch, double E, double SigmaRef) {
    if (SigmaRef <= 0.0) return;
    const double v = neutronSpeed(E);
    if (v>0.0) { T.ensureBatchAll(batch, (int)T.matNames.size()); T.cfe_global_time[batch] += 1.0/(v*SigmaRef); }
}
