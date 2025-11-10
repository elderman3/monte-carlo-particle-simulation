#include "mc.h"
#include "geomops.cpp"
#include "simulation.cpp"
#include "io.cpp"
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

// surface and delta-tracking algorithm for transport
// Statistical estimate of neutron flux: track-length and collision
// External source simulation algorithm

/*
Demonstration:

    Cylinder
    - Hollow cylinder, air + water
    - Point neutron-source in center
    - Total leak rate + analog + implicit estimates of collisions + absorptions

    - Collision visualization
    - Criticality source simulation
    - Average neutron density



*/


// Surface tracking = effectively use neutron direction to find next collision surface and collide there.
// Delta tracking = take max S and use that to sample 
// These both need effectively 1 function inside simulation to do the interaction sampling. 

// Simulation is done either in batches (External source) or in cycles (Criticality source) [Removes neutrons if too many, dupes if few]

// Fix 2 old problems


// For stats:
/*
It is assumed that the simulated population size is divided into a number of equal size
batches. The statistical estimators (mean + standard deviation) are collected by averaging over the batch-wise results.
*/

std::string mtLabel(int mt) {
    switch (mt) {
        case 2:   return "MT2 Elastic";
        case 4:   return "MT4 Inelastic";
        case 18:  return "MT18 Fission";
        case 102: return "MT102 Capture";
        default:  return "MT" + std::to_string(mt);
    }
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






int main() {
    std::ifstream file("geometry/geometry.txt");
    if (!file) { std::cout << "File opening failed\n"; return 1; }
    Universe u; // This universe is used in all further calculations, where the possible lattice will be constructed 
    Universe singleUniverse;
    std::string line;
    if (!readUniverseFile(file, singleUniverse)) { std::cout << "Failed reading Universe"; return 1; }
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
    // Universe Fully built
    // Start reading commands to execute simulation
    std::ifstream file("inputs/input.txt");
    if (!file) {
        std::cerr << "File opening failed\n";
        return 1;
    }
    int nCommands;
    if (!(file >> nCommands)) {std::cout << "Failed reading commCount"; }
    std::string command; 
    int nCount, bCount;
    bool elastic, viz;
    double energy;
    for (int i = 0; i < nCommands; ++i) {
        if (!(file >> command >> nCount >> bCount >> elastic >> viz >> energy)) { std::cout << "Failed reading parameters"; }
        if (command == "delta") {

        } else if (command == "surface") {

        } else if (command == "criticality") {

        } else if (command == "density") {
            // Run basic simulation but now calculate density statistics? -> How granular should the grid be and could I use parameters to control this? 
            
        } else {
            std::cout << "Command not valid";
        }
        // Printout Statistics, for both estimates etc
        if (viz) {
            // Run visualization of density??
        }

    }


    // Simulation must return 
        // 3d Density 
        // k-eff values
        // Estimate models 1 and 2
            // These effectively require what

    
/*
Later things
    // After IO-reading

    //printUniverse(u);
    int iter = 100000; 
    volumePointMethod(u, iter);
    std::cout << "\n Line Method\n";
    volumeLineMethod(u, iter);
    renderSliceASCII(u, SliceAxis::Z, 0.0, 100, 50);

*/

}
