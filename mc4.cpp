#include "mc.h"
#include "mc3m.cpp"
#include "mc2m.cpp"
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

int main() {
    std::ifstream file("inputs/input7.txt");
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
                const std::vector<int> MTs = inelastic ? std::vector<int>{2,4,18,102}
                                       : std::vector<int>{2,18,102};
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
                std::vector<std::string> matNames; matNames.reserve(mats.size());
                for (auto& m : mats) matNames.push_back(m.sym);
                fillData(mats, x, inelastic);

                SimRes simRes = simulation(nNeutrons, energy, iterations, maxSteps, inelastic, mats);
                
                StatsOut matrixStats = computeStats(simRes.statM);
                std::vector<double> collisionEnergy = computeColEnergy(simRes.collisions);
                double meanF = calMeanF(simRes.fNeutrons, nNeutrons);
                
                FourFactors ff = averageFour(simRes.fVec);

                printStatsOut(matrixStats, matNames, MTs);
                std::cout << "Average Neutrons Produced " << meanF << "\n"; 
                printFF(ff);
                storeDataCol(collisionEnergy);
                storeTimeHist(simRes.timeHist);

            } else {
                std::cout << "Command fail: " << line << '\n';
            }
    }}

    
    return 1;
}

// Using this as the template 

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
    //printUniverse(u);
    int iter = 100000; 
    volumePointMethod(u, iter);
    std::cout << "\n Line Method\n";
    volumeLineMethod(u, iter);
    renderSliceASCII(u, SliceAxis::Z, 0.0, 100, 50);
    
    // number indicates geometry file
    // volume calculation for:
    // 1     - Fuel pin slide 13 -- Success/Documented : All volumes slightly too large
    // 2     - Hollow cylinder slide 25 -- Success/Documented
    // 3     - Square lattice -- Success
    // 4     - Hex prism -- Success
    // 5     - Elliptical Torus -- Success : Line gives bad Size

    // Visualization for:
    // 6     - Translations and rotations -- Success
    // 7     - Hex lattices -- Success

}
