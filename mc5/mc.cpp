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
#include <chrono>

using std::array; using std::vector;
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif


bool endsWithCaseInsensitive(const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) return false;
    for (size_t i=0;i<suffix.size();++i) {
        const char a = (char)std::tolower((unsigned char)s[s.size()-suffix.size()+i]);
        const char b = (char)std::tolower((unsigned char)suffix[i]);
        if (a != b) return false;
    }
    return true;
}

std::string stripCudaSuffix(const std::string& cmd, bool* useCuda) {
    *useCuda = false;
    if (endsWithCaseInsensitive(cmd, "cuda")) {
        *useCuda = true;
        return cmd.substr(0, cmd.size() - 4);
    }
    return cmd;
}

int main(int argc, char** argv) {
    auto t0 = std::chrono::steady_clock::now();

    std::string geomFile = "geometry/geometry.txt";
    std::string inputFile = "inputs/input.txt";

    for (int ai = 1; ai < argc; ++ai) {
        const std::string a = argv[ai];
        auto needVal = [&](const char* opt) -> const char* {
            if (ai + 1 >= argc) {
                std::cerr << "Missing value after " << opt << "\n";
                std::exit(2);
            }
            return argv[++ai];
        };

        if (a == "--help" || a == "-h") {
            return 0;
        } else if (a == "--geom" || a == "-g") {
            geomFile = needVal("--geom");
        } else if (a == "--input" || a == "-i") {
            inputFile = needVal("--input");
        } else {
            std::cerr << "Unknown argument: " << a << "\n";
            return 2;
        }
    }


    std::cout << "[Build mode] " << (compiledWithCuda() ? "CUDA" : "CPU") << "\n";
    std::ifstream gfile(geomFile);
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
                    temp.pos = { C[0], C[1], 0.0f };
                    u.subUniverse.push_back(std::move(temp));
                }
            }
        } else if (singleUniverse.latticeType == 2) {
            const int cols = singleUniverse.lattice[0];
            const int rows = singleUniverse.lattice[1];
            if (cols <= 0 || rows <= 0) return 1;
            float t = singleUniverse.boundDim[0];
            float d = t * std::sqrt(3);
            float lat = 1.5 * t; 
            u.latticeType = singleUniverse.latticeType;
            u.lattice = {cols, rows};
            u.boundDim = singleUniverse.boundDim;
            u.universeShape = singleUniverse.universeShape;
            for (int i = 0; i < cols; ++i) {
                for (int j = 0; j < rows; ++j) {
                    Universe temp = singleUniverse;
                    temp.latticeType = 2;
                    auto C = hexCellCenter(u, i, j);
                    temp.pos = { C[0], C[1], 0.0f };
                    u.subUniverse.push_back(std::move(temp));
                }
            }
        }
    } else {
        u = singleUniverse;
    }
    auto t1 = std::chrono::steady_clock::now();
    auto t = std::chrono::duration<float>(t1 - t0).count();
    std::cout << "Begin simulation | Input time: " << t << "\n";


    std::ifstream cfile(inputFile);
    if (!cfile) { std::cerr << "File opening failed\n"; return 1; }
    int nCommands;
    if (!(cfile >> nCommands)) {std::cout << "Failed reading commCount"; return 1; }
    std::string command; 
    int nCount, bCount;
    bool inelastic, viz;
    float energy; float sx=0,sy=0,sz=0; // Source positions if applicable


    for (int i = 0; i < nCommands; ++i) {
        if (!(cfile >> command >> nCount >> bCount >> inelastic >> viz >> energy)) { std::cout << "Failed reading parameters"; return 1; }
        if (cfile.peek()==' ') (void)(cfile >> sx >> sy >> sz);
        RunParams P;
        const auto MTs = inelastic ? std::vector<int>{2,4,18,102} : std::vector<int>{2,18,102};
        P.historiesPerBatch = nCount; P.batches = bCount; P.inelastic = inelastic; P.sourceE = energy;
        P.sourcePos = {sx,sy,sz}; P.maxSteps = 100000;
        Mesh3D mesh = makeMeshFromUniverse(u, 300, 300, 300);
        P.mesh = &mesh;
        mesh.meshAnalogColl.assign(mesh.nx*mesh.ny*mesh.nz, 0.0f);

        bool useCuda = false;
        std::string cmdBase = command;
        if (endsWithCaseInsensitive(cmdBase, "cuda")) {
            useCuda = true;
            cmdBase = cmdBase.substr(0, cmdBase.size() - 4);
        }
        std::string outTag = cmdBase + (useCuda ? "_cuda" : "");
        std::cout << "Using " << ((useCuda) ? "CUDA" : "CPU") << "\n";

        if (cmdBase == "delta") {
            P.track = Tracking::Delta; P.src = SourceMode::External;
            auto R = useCuda ? runExternalCuda(u, P) : runExternal(u, P);
            auto S = computeStats(R.T.statM);
            printAndStoreRates(P, R, S, outTag);
            std::cout << "\n[External | Delta-tracking] results\n"; printStatsOut(S, R.T.matNames, MTs, std::cout);
            printFoms(P, R, MTs, outTag);
            if (viz) {
                writeVTKStructuredPoints(mesh, mesh.cfe_density, "density_cfe", "density_cfe_d");
                writeVTKStructuredPoints(mesh, mesh.meshAnalogColl, "analog_collisions", "analog_collisions");
                auto gtag = buildGeometryTagField(mesh, u);
                writeVTKStructuredPoints(mesh, gtag, "geometry_tag", "geometry_tag");
            }

        } else if (cmdBase == "surface") {
            P.track = Tracking::Surface; P.src = SourceMode::External;
            auto R = useCuda ? runExternalCuda(u, P) : runExternal(u, P);
            auto S = computeStats(R.T.statM);
            printAndStoreRates(P, R, S, outTag);
            std::cout << "\n[External | Surface-tracking] results\n"; printStatsOut(S, R.T.matNames, MTs, std::cout);
            printFoms(P, R, MTs, outTag);
            if (viz) {
                writeVTKStructuredPoints(mesh, mesh.cfe_density, "density_cfe", "density_cfe");
                writeVTKStructuredPoints(mesh, mesh.tle_density, "density_tle", "density_tle");
                writeVTKStructuredPoints(mesh, mesh.meshAnalogColl, "analog_collisions", "analog_collisions");
                auto gtag = buildGeometryTagField(mesh, u);
                writeVTKStructuredPoints(mesh, gtag, "geometry_tag", "geometry_tag");
            }

        } else if (cmdBase == "criticality" || cmdBase == "criticalitydelta" || cmdBase == "criticality_delta") {
            P.src = SourceMode::Criticality;
            P.track = (cmdBase == "criticality") ? Tracking::Surface : Tracking::Delta;
            auto R = useCuda ? runCriticalityCuda(u, P, 5) : runCriticality(u, P, 5);
            auto S = computeStats(R.T.statM);
            std::cout << "\n[Criticality | " << (P.track == Tracking::Delta ? "Delta" : "Surface") << "-Tracking] results\n";
            printStatsOut(S, R.T.matNames, MTs, std::cout);
            storeDatakeff(R.keff_history);
            std::cout << "\n[Criticality] keff history written. Mean ~ "
                      << std::accumulate(R.keff_history.begin(), R.keff_history.end(), 0.0f)/std::max<size_t>(1, R.keff_history.size()) << "\n";
            float h = 0.50;
            float Vmix = M_PI * 0.20 * 0.20 * h;
            estimateNeutronDensity1W(P, R, Vmix, 2.43);

            if (viz) {
                writeVTKStructuredPoints(mesh, mesh.cfe_density, "density_cfe", "criticality_density_cfe");
                writeVTKStructuredPoints(mesh, mesh.tle_density, "density_tle", "criticality_density_tle");
                writeVTKStructuredPoints(mesh, mesh.meshAnalogColl, "analog_collisions", "criticality_analog_collisions");
                auto gtag = buildGeometryTagField(mesh, u);
                writeVTKStructuredPoints(mesh, gtag, "geometry_tag", "criticality_geometry_tag");
            }

        } else {
            std::cout << "Command not valid";
        }
    }
}