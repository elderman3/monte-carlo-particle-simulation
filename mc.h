#include <vector>
#include <utility>
#include <string>

typedef struct Material{
    std::string sym;
    int z;
    int a;
    double aw;
    int T;
    double rho;
    double proportion; // molar proportion of composition
    std::vector<std::pair<double, double>> neutrons;
    std::map<int, std::pair<double, std::vector<std::pair<double,double>>>> mt;
} Material;

typedef struct Result{
    double runningTime;
    float pi;
} Result;

typedef struct Stats{
    int iterations;
    int trials;
    float length; 
    float piMean;
    float piStd;
    float time;
    float timeMean;
    float timeStd;
    float fom;
    float jb;
    bool normalized;
} Stats;

double randomVal(float min = 0.f, float max = 1.f);

float circleSampling(int iter, float l = 0.f);

float buffonNeedle(int iter, float length);

Result timing(float (*f)(int, float), int iter, float length = 0.f);

Stats statistics(const std::vector<Result> inputV);

void printStats(const Stats& s);





