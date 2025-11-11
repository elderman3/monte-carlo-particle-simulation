#include <vector>
#include <utility>
#include <string>
#include <map>
#include <array>
#include <deque>

using std::array; using std::vector;

enum class Op : uint8_t { L=0, N=1, U=2, I=3 };
git a
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
    std::array<double, 3> pos;
    std::array<double, 3> dir;
    double vel; 
    double energy;
    int collisions = 0;
    double time = 0;
    bool reachedTh;
    bool isSource;
} Neutron;

typedef struct Node{
  Op op;
  int shape; // Index of Shape
  int left; // child index or -1
  int right; // child index or -1
} Node;

// A single Plane of a Geometry (Quadratic and Torus)
typedef struct Shape{
    bool torus = 0; // For torus, A, B, C are used as a,b,R. Additionally, DEF are center and GHI Axis
    double A,B,C,D,E,F,G,H,I,J;
} Shape;

// A single compound shape e.g., Cube, Cylinder
typedef struct Geometry{
    std::array<double, 3> pos{0.f, 0.f, 0.f};
    std::array<double, 3> rot{0.f, 0.f, 0.f}; // Always keep these 0, Universes are not rotatable. 
    std::vector<Shape> shapes;
    int shape; // Defines the possible shape, tabled elsewhere
    int nodeRoot;
    std::vector<Material> mats; // Vector for storing the material that this geometry is made of
    std::vector<Node> nodes; // Set of operations that define the shape of Geometry
    // binary tree of intersections etc to define the shape fully
} Geometry;

// A collection of Universes and Geometries
typedef struct Universe{
    std::string name;
    std::array<double, 3> pos{0.f, 0.f, 0.f};
    std::array<double, 3> rot{0.f, 0.f, 0.f}; // Always keep these 0, Universes are not rotatable. 
    std::array<double, 3> boundDim{0.f, 0.f, 0.f}; // Bounding box construction args
    Geometry boundingGeometry; // Universe bounding box
    int universeShape; // Defines how the Universe Bounding box is constructed
    std::array<int, 2> lattice{0, 0}; // Only valid for square or hex major Universes
    int latticeType = 0; // 0 For none, 1 for square, 2 for hex
    std::vector<Universe> subUniverse;
    std::vector<Geometry> geometries;
} Universe;

struct Neutron{
    std::array<double, 3> pos{0.0,0.0,0.0};
    std::array<double, 3> dir{0.0,0.0,1.0}; // unit
    double vel = 0.0; 
    double energy = 0.0; // MeV
    int collisions = 0;
    double time = 0.0;   // s
    bool reachedTh = false;
    bool isSource = true;
    double w = 1.0;      // statistical weight (for implicit crit)
} Neutron;

enum class Tracking { Surface, Delta };
enum class SourceMode { External, Criticality };

struct Mesh3D{
    std::array<double,3> pmin{0,0,0}, pmax{0,0,0};
    int nx=0, ny=0, nz=0;
    std::vector<double> cfe_density; // sum (1/v)/Sigma per (accepted) collision
    std::vector<double> tle_density; // sum track_length/v over segments (surface only)
    inline int idx(int i,int j,int k) const { return (k*ny + j)*nx + i; }
    inline bool inside(const std::array<double,3>& p) const {
        return (p[0]>=pmin[0] && p[0]<=pmax[0] &&
                p[1]>=pmin[1] && p[1]<=pmax[1] &&
                p[2]>=pmin[2] && p[2]<=pmax[2]);
    }
    inline void zero(){
        cfe_density.assign(std::max(1,nx*ny*nz), 0.0);
        tle_density.assign(std::max(1,nx*ny*nz), 0.0);
    }
} Mesh3D;

struct StatsOut{ 
    std::vector<std::vector<double>> mean; 
    relErr; 
    std::vector<std::vector<int>> sum; 
} StatsOut;

struct RxSample {
    int mt = 0;
    int matIndex = -1;
    const Material* mat = nullptr;
} RxSample;

struct FlightResult {
    bool leaked=false;
    bool collided=false;
    bool virtualCollision=false;
    array<double,3> pos{};
    const Geometry* geom=nullptr;
    double SigmaLocal=0.0;
    double traveled=0.0;
} FlightResult;

struct TallyBook {
    std::vector<std::string> matNames;
    std::map<std::string,int> mat2idx;
    std::vector<std::vector<std::vector<int>>> statM;
    Mesh3D* mesh = nullptr;
    bool useTLE = false;
    bool useCFE = true;
    bool deltaMode = false;
    double SigmaM = 0.0;

    void ensure_batch(int batch, int M, int R) {
        if ((int)statM.size()<=batch)
            statM.resize(batch+1, std::vector<std::vector<int>>(M, std::vector<int>(R,0)));
    }
    int matIndex(const Material& m) {
        auto it = mat2idx.find(m.sym);
        if (it!=mat2idx.end()) return it->second;
        int id = (int)matNames.size();
        mat2idx[m.sym]=id; matNames.push_back(m.sym);
        for (auto& B: statM) B.push_back(std::vector<int>(B.empty()?0:B[0].size(),0));
        return id;
    }
} TallyBook;

struct RunParams{
    Tracking track = Tracking::Surface;
    SourceMode src = SourceMode::External;
    int historiesPerBatch = 10000;
    int batches = 10;
    int maxSteps = 10000;
    bool inelastic = true;
    double sourceE = 1.0;
    array<double,3> sourcePos{0,0,0};
    Mesh3D* mesh = nullptr;
} RunParams;

struct RunOutputs {
    TallyBook T;
    std::vector<Collisions> collisions;
    std::vector<int> fissionChildren;
    std::vector<double> keff_history;
} RunOutputs;

double randomVal(float min = 0.f, float max = 1.f);


