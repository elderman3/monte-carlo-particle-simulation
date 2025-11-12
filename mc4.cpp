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

using std::array; using std::vector;
#define M_PI 3.14159265358979323846264338327950288

std::string mtLabel(int mt) {
    switch (mt) {
        case 2:   return "MT2 Elastic";
        case 4:   return "MT4 Inelastic";
        case 18:  return "MT18 Fission";
        case 102: return "MT102 Capture";
        default:  return "MT" + std::to_string(mt);
    }
}

Mesh3D makeMeshFromUniverse(const Universe& U, int nx,int ny,int nz) {
    Universe tmp = U;
    auto ext = boundingBox(tmp);
    Mesh3D M; M.nx=nx; M.ny=ny; M.nz=nz;
    M.pmin = {-0.5*ext[0], -0.5*ext[1], -0.5*ext[2]};
    M.pmax = {+0.5*ext[0], +0.5*ext[1], +0.5*ext[2]};
    M.zero();
    return M;
}

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

double calMeanF(const vector<int>& fNeutrons, int nNeutrons) {
    long long S = 0; for (int x : fNeutrons) S += x;
    const int I = int(fNeutrons.size());
    return (I > 0 && nNeutrons > 0) ? double(S) / (double(I) * nNeutrons) : 0.0;
}

FourTally sumFour(const vector<FourTally>& v) {
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

FourFactors averageFour(const vector<FourTally>& perIter) {
    return factorsFrom(sumFour(perIter));
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

struct ScalarStat {
    double mean=0, stddev=0, relErr=0, fom=0;
    int batches=0;
};

ScalarStat stats_from_batches(const std::vector<double>& x, double elapsed_s) {
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

struct RegionMap { const char* label; std::vector<std::string> mats; };

std::vector<RegionMap> defaultRegions() {
    return {
        {"Water",     {"H-1","O-16"}},
        {"Air",       {"N-14"}},
        {"Container", {"Fe-56"}}
    };
}

double norm_per_batch(const RunParams& P) {
    return (P.historiesPerBatch>0)? (P.sourceRate / double(P.historiesPerBatch)) : 0.0;
}

std::map<std::string,int> matIndexMap(const TallyBook& T) {
    std::map<std::string,int> m2i; for (int i=0;i<(int)T.matNames.size();++i) m2i[T.matNames[i]] = i; return m2i;
}

int colIdx(const std::vector<int>& mts, int mt) {
    for (int j=0;j<(int)mts.size();++j) if (mts[j]==mt) return j; return -1;
}
struct BatchMetrics {
    std::vector<double> leak_nps;
    std::vector<double> abs_water_nps;
    std::vector<double> abs_steel_nps;
    std::vector<double> tle_abs_water_nps;
    std::vector<double> cfe_abs_water_nps;
    std::vector<double> cfe_tau_s;
};

BatchMetrics build_batch_metrics(const RunParams& P, const RunOutputs& R, const std::vector<int>& MTs) {
    const auto regs = defaultRegions();
    const auto m2i = matIndexMap(R.T);
    const double norm = norm_per_batch(P);
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

void print_foms(const RunParams& P, const RunOutputs& R, const std::vector<int>& MTs, const std::string& tag) {
    const auto B = build_batch_metrics(P, R, MTs);
    auto pr = [&](const char* name, const std::vector<double>& x) {
        auto st = stats_from_batches(x, R.perf.elapsed_s);
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

void accumulate_by_region(const TallyBook& T, const StatsOut& S, std::map<std::string,int>& name2idx, std::vector<double>& analogColl, std::vector<double>& analogAbs) {
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

void print_and_store_rates(const RunParams& P, const RunOutputs& R, const StatsOut& S, const std::string& tag) {
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

void estimate_neutron_density_1W(const RunParams& P, const RunOutputs& R, double mixtureVolume_m3, double meanNuBar = 2.43) {
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

int main() {
    std::ifstream gfile("geometry/geometry.txt");
    if (!gfile) { std::cout << "File opening failed\n"; return 1; }
    Universe u; // This universe is used in all further calculations, where the possible lattice will be constructed 
    Universe singleUniverse;
    std::string line;
    if (!readUniverseFile(gfile, singleUniverse)) { std::cout << "Failed reading Universe"; return 1; }
    if (singleUniverse.latticeType) {
        if (singleUniverse.latticeType == 1) {
            const int cols = singleUniverse.lattice[0];
            const int rows = singleUniverse.lattice[1];
            if (cols <= 0 || rows <= 0) return 1;
            u.latticeType = singleUniverse.latticeType;
            u.lattice = {cols, rows};
            u.boundDim = singleUniverse.boundDim;
            u.universeShape = singleUniverse.universeShape;
            for (int i = 0; i < cols; ++i) {
                for (int j = 0; j < rows; ++j) {
                    Universe temp = singleUniverse;
                    temp.latticeType = 1;
                    auto C = squareCellCenter(u, i, j);
                    temp.pos = { C[0], C[1], 0.0 };
                    u.subUniverse.push_back(std::move(temp));
                }
            }
        } else if (singleUniverse.latticeType == 2) {
            const int cols = singleUniverse.lattice[0];
            const int rows = singleUniverse.lattice[1];
            if (cols <= 0 || rows <= 0) return 1;
            double t = singleUniverse.boundDim[0];
            double d = t * std::sqrt(3);
            double lat = 1.5 * t; 
            u.latticeType = singleUniverse.latticeType;
            u.lattice = {cols, rows};
            u.boundDim = singleUniverse.boundDim;
            u.universeShape = singleUniverse.universeShape;
            for (int i = 0; i < cols; ++i) {
                for (int j = 0; j < rows; ++j) {
                    Universe temp = singleUniverse;
                    temp.latticeType = 2;
                    auto C = hexCellCenter(u, i, j);
                    temp.pos = { C[0], C[1], 0.0 };
                    u.subUniverse.push_back(std::move(temp));
                }
            }
        }
    } else {
        u = singleUniverse;
    }
    /*
    // Size for debugging
    std::cout << "\n Point Method\n";
    volumePointMethod(u, 1000);
    std::cout << "\n Line Method\n";
    volumeLineMethodTorus(u, 1000);
    return 1;
    */
    // Universe Fully built
    // Start reading commands to execute simulation
    std::ifstream cfile("inputs/input.txt");
    if (!cfile) { std::cerr << "File opening failed\n"; return 1; }
    int nCommands;
    if (!(cfile >> nCommands)) {std::cout << "Failed reading commCount"; return 1; }
    std::string command; 
    int nCount, bCount;
    bool inelastic, viz;
    double energy; double sx=0,sy=0,sz=0; // Source positions if applicable
    for (int i = 0; i < nCommands; ++i) {
        if (!(cfile >> command >> nCount >> bCount >> inelastic >> viz >> energy)) { std::cout << "Failed reading parameters"; return 1; }
        if (cfile.peek()==' ') (void)(cfile >> sx >> sy >> sz);
        RunParams P;
        const auto MTs = inelastic ? std::vector<int>{2,4,18,102} : std::vector<int>{2,18,102};
        P.historiesPerBatch = nCount; P.batches = bCount; P.inelastic = inelastic; P.sourceE = energy;
        P.sourcePos = {sx,sy,sz}; P.maxSteps = 1000000;
        Mesh3D mesh = makeMeshFromUniverse(u, 60, 60, 60);
        P.mesh = &mesh;
        mesh.meshAnalogColl.assign(mesh.nx*mesh.ny*mesh.nz, 0.0);
        if (command == "delta") {
            P.track = Tracking::Delta; P.src = SourceMode::External;
            auto R = runExternal(u, P);
            auto S = computeStats(R.T.statM);
            print_and_store_rates(P, R, S, (command=="delta" ? "delta" : "surface"));
            std::cout << "\n[External | Delta-tracking] results\n"; printStatsOut(S, R.T.matNames, MTs, std::cout);
            print_foms(P, R, MTs, (command=="delta" ? "delta" : "surface"));
            if (viz) {
                writeVTKStructuredPoints(mesh, mesh.cfe_density, "density_cfe_d", "density_cfe_d");
                writeVTKStructuredPoints(mesh, mesh.meshAnalogColl, "analog_collisions", "analog_collisions");
            }
        } else if (command == "surface") {
            P.track = Tracking::Surface; P.src = SourceMode::External;
            auto R = runExternal(u, P);
            auto S = computeStats(R.T.statM);
            print_and_store_rates(P, R, S, (command=="delta" ? "delta" : "surface"));
            std::cout << "\n[External | Surface-tracking] results\n"; printStatsOut(S, R.T.matNames, MTs, std::cout);
            print_foms(P, R, MTs, (command=="delta" ? "delta" : "surface"));
            if (viz) {
                writeVTKStructuredPoints(mesh, mesh.cfe_density, "density_cfe", "density_cfe");
                writeVTKStructuredPoints(mesh, mesh.tle_density, "density_tle", "density_tle");
                writeVTKStructuredPoints(mesh, mesh.meshAnalogColl, "analog_collisions", "analog_collisions");
            }
        } else if (command == "criticality") {
            P.track = Tracking::Surface; P.src = SourceMode::Criticality;
            auto R = runCriticality(u, P, 5);
            auto S = computeStats(R.T.statM);
            std::cout << "\n[Criticality | Surface-Tracking] results\n"; printStatsOut(S, R.T.matNames, MTs, std::cout);
            storeDatakeff(R.keff_history);
            std::cout << "\n[Criticality] keff history written. Mean ~ "
                      << std::accumulate(R.keff_history.begin(), R.keff_history.end(), 0.0)/std::max<size_t>(1, R.keff_history.size()) << "\n";
            double h = 0.50;
            double Vmix = M_PI * 0.20 * 0.20 * h;
            estimate_neutron_density_1W(P, R, Vmix, 2.43);
        } else {
            std::cout << "Command not valid";
        }
    }
}
