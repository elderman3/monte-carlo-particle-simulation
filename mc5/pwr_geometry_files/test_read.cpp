
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <vector>
#include <utility>
#include <cmath>
#include "mc.h"

// Minimal implementations for missing linker symbols (not included in provided subset)
std::vector<double> logspace(double start, double end, int n){
    std::vector<double> out;
    out.reserve(n);
    double log_start = std::log10(start);
    double log_end = std::log10(end);
    double step = (log_end - log_start) / (n - 1);
    for(int i=0;i<n;i++){
        out.push_back(std::pow(10.0, log_start + step*i));
    }
    return out;
}

double valueInterp(const std::vector<std::pair<double,double>>& values, double x){
    if(values.empty()) return 0.0;
    if(x <= values.front().first) return values.front().second;
    if(x >= values.back().first) return values.back().second;
    // linear search (fine for test)
    for(size_t i=1;i<values.size();i++){
        if(x <= values[i].first){
            double x0 = values[i-1].first, y0=values[i-1].second;
            double x1 = values[i].first, y1=values[i].second;
            double t = (x - x0) / (x1 - x0);
            return y0 + t*(y1 - y0);
        }
    }
    return values.back().second;
}

int main(int argc, char** argv){
    std::string fname = "geometry/geometry.txt";
    if(argc>1) fname = argv[1];
    std::ifstream in(fname);
    if(!in){
        std::cerr << "Failed to open " << fname << "\n";
        return 1;
    }
    Universe u;
    if(!readUniverseFile(in, u)){
        std::cerr << "readUniverseFile returned false\n";
        return 2;
    }
    std::cout << "Loaded universe: " << u.name << "\n";
    std::cout << "Subuniverses: " << u.subUniverse.size() << "\n";
    std::cout << "Geometries: " << u.geometries.size() << "\n";
    std::function<void(const Universe&, int&, int&)> rec = [&](const Universe& uu, int& nU, int& nGeom){
        nU += 1;
        nGeom += (int)uu.geometries.size();
        for(const auto& su: uu.subUniverse) rec(su, nU, nGeom);
    };
    int nU=0,nGeom=0;
    rec(u,nU,nGeom);
    std::cout << "Total universes in tree: " << nU << "\n";
    std::cout << "Total geometry objects: " << nGeom << "\n";
    return 0;
}
