#include "mc.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <chrono>
#include <random>
#include <cmath>
#include "mathplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

constexpr double PI_2 = 1.57079632679489661923;

double randomVal(float min, float max) {
    static thread_local std::mt19937 gen(std::random_device{}()); // engine
    std::uniform_real_distribution<double> dist(min, max);        // [0,1)

    return dist(gen);
}

Result timing(float (*f)(int, float), int iter, float length) {
    auto start = std::chrono::steady_clock::now();
    float pi = f(iter, length);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    Result result{elapsed.count(), pi};
    return result;
}

Stats statistics(const std::vector<Result> inputV) {
    Stats result;
    int len = inputV.size();
    float sumTime = 0, sumPi = 0;
    for (int i = 0; i < len; ++i) {
        sumTime += inputV[i].runningTime;
        sumPi += inputV[i].pi;
    }
    result.time = sumTime; 
    result.piMean = sumPi / (float) len;
    result.timeMean = sumTime / (float) len;

    float stdSumTime = 0, stdSumPi = 0;
    float moment3 = 0, moment4 = 0;
    for (int i = 0; i < len; ++i) {
        stdSumTime += pow(result.timeMean - inputV[i].runningTime, 2);
        stdSumPi += pow(result.piMean - inputV[i].pi, 2);
        moment3 += pow(result.piMean - inputV[i].pi, 3);
        moment4 += pow(result.piMean - inputV[i].pi, 4);
    }
    result.piStd = sqrt(stdSumPi / len);
    result.timeStd = sqrt(stdSumTime / len);
    result.fom = 1 / (pow(result.piStd, 2) * result.time);

    float s = moment3 / (pow(result.piStd, 3) * len);
    float k = moment4 / (pow(result.piStd, 4) * len);
    result.jb = ((float) len / 6) * pow(s, 2) + ((float) len / 24) * pow((k - 3),2);
    if (result.jb > 5.991) result.normalized = 0;
    if (result.jb <= 5.991) result.normalized = 1;
    return result;
}

void printStats(const Stats& s) {
    std::cout << "Stats:\n"
              << "  piMean   = " << s.piMean     << "\n"
              << "  piStd    = " << s.piStd      << "\n"
              << "  time     = " << s.time       << "\n"
              << "  timeMean = " << s.timeMean   << "\n"
              << "  timeStd  = " << s.timeStd    << "\n"
              << "  fom      = " << s.fom        << "\n"
              << "  jb       = " << s.jb         << "\n"
              << "  normal   = " << s.normalized << "\n"
              << "\n" 
              << "  iterations = " << s.iterations << "\n"
              << "  trials     = " << s.trials     << "\n"
              << "  length     = " << s.length     << "\n";
}

int main() {
    std::ifstream file("input.txt");
    if (!file) {
        std::cerr << "File opening failed\n";
        return 1;
    }
    std::string line;
    while (std::getline(file, line)) {
        char command[32];
        int nNeutrons, iterations, nNmaterials;
        float energy;
        if (std::sscanf(line.c_str(), "%31s %d %f %d %d", command, &nNeutrons, &energy, &iterations, &nMaterials) == 5) {
            if (std::strcmp(command, "simulation") == 0) {
                std::vector<Material> mats;
                for (int i = 0; i < nMaterials; ++i) {
                    if (std::sscanf(line.c_str(), "%31s %d %f %d %d", command, &nNeutrons, &energy, &iterations, &nMaterials) == 5) {
                        Material mat;
                        
                    }
                }

                

            }


            if (std::strcmp(command, "circle") == 0) {
                circleResult.resize(repeats);
                for (int i = 0; i < repeats; ++i) {
                    circleResult[i] = timing(circleSampling, iter);
                }
                Stats statsCircle = statistics(circleResult);
                statsCircle.iterations = repeats;
                statsCircle.trials = iter;
                statsCircle.length = length;
                printf("Circle Function Statistics\n");
                printStats(statsCircle);
            } else if (std::strcmp(command, "buffon") == 0) {
                buffonResult.resize(repeats);
                for (int i = 0; i < repeats; ++i) {
                    buffonResult[i] = timing(buffonNeedle, iter, length);
                }
                Stats statsBuffon = statistics(buffonResult);
                statsBuffon.iterations = repeats;
                statsBuffon.trials = iter;
                statsBuffon.length = length;
                printf("Buffon Function Statistics\n");
                printStats(statsBuffon);
            } else {
                std::cout << "unknown command: " << command << '\n';
            }
        } else {
            std::cout << "Command fail: " << line << '\n';
        }
    }

}