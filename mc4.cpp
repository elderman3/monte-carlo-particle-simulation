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

            } else if (std::strcmp(command, "criticality") == 0) {
                const double N_H2O = 0.0334;
                const double N_UO2 = 0.0244672;

                const double enr = 0.05;
                const int nNeutrons = 100;
                const int iterations = 10;
                const int maxSteps = 1000;
                const int inelastic = 1;
                const double E0 = 0.5;

                std::vector<double> meankeff(200);

                for (int pctFuel = 2; pctFuel <= 20; ++pctFuel) {
                    const double phi_UO2 = pctFuel / 100.0;
                    const double phi_H2O = 1.0 - phi_UO2;

                    Material matH1, matO16_H2O, matU235, matU238, matO16_UO2;
                    {
                        std::ifstream dH1("data/H1.dat");
                        std::ifstream dO16a("data/O16.dat");
                        std::ifstream dU235("data/U235.dat");
                        std::ifstream dU238("data/U238.dat");
                        std::ifstream dO16b("data/O16.dat");
                        readMaterialBlock(dH1, matH1);
                        readMaterialBlock(dO16a, matO16_H2O);
                        readMaterialBlock(dU235, matU235);
                        readMaterialBlock(dU238, matU238);
                        readMaterialBlock(dO16b, matO16_UO2);
                    }
                    matO16_H2O.sym = "O16_H2O";
                    matO16_UO2.sym = "O16_UO2";

                    matH1.rho = N_H2O * phi_H2O; matH1.proportion = 2;
                    matO16_H2O.rho = N_H2O * phi_H2O; matO16_H2O.proportion = 1;

                    matU235.rho = N_UO2 * phi_UO2 * enr; matU235.proportion = 1;
                    matU238.rho = N_UO2 * phi_UO2 * (1-enr); matU238.proportion = 1;
                    matO16_UO2.rho = N_UO2 * phi_UO2; matO16_UO2.proportion = 2;

                    std::vector<Material> mats;
                    mats.reserve(5);
                    mats.push_back(std::move(matH1));
                    mats.push_back(std::move(matO16_H2O));
                    mats.push_back(std::move(matU235));
                    mats.push_back(std::move(matU238));
                    mats.push_back(std::move(matO16_UO2));

                    fillData(mats, x, inelastic);
                    SimRes simRes = simulation(nNeutrons, E0, iterations, maxSteps, inelastic, mats);

                    FourFactors ff = averageFour(simRes.fVec);
                    double meanF = calMeanF(simRes.fNeutrons, nNeutrons);

                    std::cout << "phi_UO2=" << phi_UO2
                            << "  phi_H2O=" << phi_H2O
                            << "  meanF=" << meanF << "\n";
                    printFF(ff);
                    printFourTally(simRes.fVec[0]);

                    meankeff[pctFuel] = ff.keff;
                }
                storeDatakeff(meankeff);
            }
        } else {
            std::cout << "Command fail: " << line << '\n';
        }
    }

    
    return 1;
}


