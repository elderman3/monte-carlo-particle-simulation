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
    std::map<int, std::vector<std::pair<double,double>>> mt;
    std::map<int, double> Qvals;
} Material;

typedef struct Collisions{
    std::vector<int> num;
    std::vector<double> sumEnergy;
} Collisions;

typedef struct SimRes{
    std::vector<std::vector<std::vector<int>>> statM;
    std::vector<Collisions> collisions;
    std::vector<int> fNeutrons;
} SimRes;

struct StatsOut {
    std::vector<std::vector<double>> mean;
    std::vector<std::vector<double>> relErr;
    std::vector<std::vector<int>> sum;
};

typedef struct Neutron{
    double energy;
    int collisions = 0;
} Neutron;

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


