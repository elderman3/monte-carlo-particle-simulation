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
            printAndStoreRates(P, R, S, (command=="delta" ? "delta" : "surface"));
            std::cout << "\n[External | Delta-tracking] results\n"; printStatsOut(S, R.T.matNames, MTs, std::cout);
            printFoms(P, R, MTs, (command=="delta" ? "delta" : "surface"));
            if (viz) {
                writeVTKStructuredPoints(mesh, mesh.cfe_density, "density_cfe_d", "density_cfe_d");
                writeVTKStructuredPoints(mesh, mesh.meshAnalogColl, "analog_collisions", "analog_collisions");
            }

        } else if (command == "surface") {
            P.track = Tracking::Surface; P.src = SourceMode::External;
            auto R = runExternal(u, P);
            auto S = computeStats(R.T.statM);
            printAndStoreRates(P, R, S, (command=="delta" ? "delta" : "surface"));
            std::cout << "\n[External | Surface-tracking] results\n"; printStatsOut(S, R.T.matNames, MTs, std::cout);
            printFoms(P, R, MTs, (command=="delta" ? "delta" : "surface"));
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
            estimateNeutronDensity1W(P, R, Vmix, 2.43);

        } else {
            std::cout << "Command not valid";
        }
    }
}