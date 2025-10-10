#include <vector>
#include <utility>
#include <string>
#include <map>
#include <array>
#include <deque>

typedef struct Material{
    std::string sym{};
    int z = 0;
    int a = 0;
    double aw = 0.0;
    int T = 0;
    double rho = 0;
    double proportion = 0; // molar proportion of composition
    std::vector<std::pair<double, double>> neutrons;
    std::map<int, std::vector<std::pair<double,double>>> mt;
    std::map<int, double> Qvals;
} Material;

typedef struct Collisions{
    std::vector<int> num;
    std::vector<double> sumEnergy;
} Collisions;

typedef struct FourTally{
    int fissionBirthsTotal = 0;
    int fissionBirthsThermal = 0;
    int absThTotal = 0;
    int absThFuel = 0;
    int started = 0;
    int reachedThermal = 0;
    int resAbsBeforeThermal = 0;
} FourTally;

typedef struct FourFactors {
    double eta = 0.0;
    double eps = 0.0;
    double p = 0.0; 
    double f = 0.0;
    double keff = 0.0;
} FourFactors;

typedef struct TimeHist{
    double dt; int nbins; std::vector<long long> counts;
    TimeHist(double dt_, int nb_): dt(dt_), nbins(nb_), counts(nb_,0) {}
    inline void add(double t){ int k = int(t/dt); if (0<=k && k<nbins) ++counts[k]; }
} TimeHist;

typedef struct SimRes{
    std::vector<std::vector<std::vector<int>>> statM;
    std::vector<Collisions> collisions;
    std::vector<int> fNeutrons;
    std::vector<FourTally> fVec;
    TimeHist timeHist;
} SimRes;

struct StatsOut {
    std::vector<std::vector<double>> mean;
    std::vector<std::vector<double>> relErr;
    std::vector<std::vector<int>> sum;
};

typedef struct Neutron{
    double energy;
    int collisions = 0;
    double time = 0;
    bool reachedTh;
    bool isSource;
} Neutron;

// A collection of Universes and Geometries
typedef struct Universe{
    std::array<double, 3> pos{0.f, 0.f, 0.f};
    std::array<double, 3> rotation{0.f, 0.f, 0.f};
    std::vector<Universe> subUniverse;
    std::vector<Geometry> geometries;
} Universe;

// A single compound shape e.g., Cube, Cylinder
typedef struct Geometry{
    std::array<double, 3> pos{0.f, 0.f, 0.f};
    std::array<double, 3> rotation{0.f, 0.f, 0.f};
    std::array<double, 3> absDim{0.f, 0.f, 0.f};
    std::vector<Shape> shape;
    Node rules; // Set of operations that define the shape of Geometry
    // binary tree of intersections etc to define the shape fully
} Geometry;

// A single Plane of a Geometry (Quadratic and Torus)
typedef struct Shape{
    bool torus = 0; // For torus, A, B, C are used as a,b,R
    double A; double B; double C; double D; double E; double F; double G; double H; double I; double J;
} Shape;

typedef struct Node{
    int type; // This is 1 for not, 2 for union, 3 for intersection
    int shape; // Index number of universe/geometry that this
    Node* parent;
    std::vector<Node> children; // If empty, end of branch
} Node;

double randomVal(float min = 0.f, float max = 1.f);


