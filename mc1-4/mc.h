#pragma once
#include <vector>
#include <utility>
#include <string>
#include <map>
#include <array>
#include <deque>

using std::array; using std::vector;

enum class Op : uint8_t { L=0, N=1, U=2, I=3 };

typedef struct Material{
    std::string sym{};
    int z = 0;
    int a = 0;
    double aw = 0.0;
    int T = 0;
    double rho = 0;
    double proportion = 0; // molar proportion of composition - Turns out to be useless 
    vector<std::pair<double, double>> neutrons;
    std::map<int, vector<std::pair<double,double>>> mt;
    std::map<int, double> Qvals;
} Material;

typedef struct Collisions{
    vector<int> num;
    vector<double> sumEnergy;
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
    double dt; int nbins; vector<long long> counts;
    TimeHist(double dt_, int nb_): dt(dt_), nbins(nb_), counts(nb_,0) {}
    inline void add(double t) { int k = int(t/dt); if (0<=k && k<nbins) ++counts[k]; }
} TimeHist;

typedef struct SimRes{
    vector<vector<vector<int>>> statM;
    vector<Collisions> collisions;
    vector<int> fNeutrons;
    vector<FourTally> fVec;
    TimeHist timeHist;
} SimRes;

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
    array<double, 3> pos{0.f, 0.f, 0.f};
    array<double, 3> rot{0.f, 0.f, 0.f}; // Always keep these 0, Universes are not rotatable. 
    vector<Shape> shapes;
    int shape; // Defines the possible shape, tabled elsewhere
    int nodeRoot;
    vector<Material> mats; // Vector for storing the material that this geometry is made of
    vector<Node> nodes; // Set of operations that define the shape of Geometry
    // binary tree of intersections etc to define the shape fully
} Geometry;

// A collection of Universes and Geometries
typedef struct Universe{
    std::string name;
    array<double, 3> pos{0.f, 0.f, 0.f};
    array<double, 3> rot{0.f, 0.f, 0.f}; // Always keep these 0, Universes are not rotatable. 
    array<double, 3> boundDim{0.f, 0.f, 0.f}; // Bounding box construction args
    Geometry boundingGeometry; // Universe bounding box
    int universeShape; // Defines how the Universe Bounding box is constructed
    array<int, 2> lattice{0, 0}; // Only valid for square or hex major Universes
    int latticeType = 0; // 0 For none, 1 for square, 2 for hex
    vector<Universe> subUniverse;
    vector<Geometry> geometries;
} Universe;

typedef struct Neutron{
    array<double, 3> pos{0.0,0.0,0.0};
    array<double, 3> dir{0.0,0.0,1.0}; // unit
    double vel = 0.0; 
    double energy = 0.0; // MeV
    int collisions = 0;
    double time = 0.0; // s
    bool reachedTh = false;
    bool isSource = true;
    double w = 1.0; // statistical weight (for implicit crit)
} Neutron;

enum class Tracking { Surface, Delta };
enum class SourceMode { External, Criticality };

typedef struct Mesh3D{
    array<double,3> pmin{0,0,0}, pmax{0,0,0};
    int nx=0, ny=0, nz=0;
    vector<double> cfe_density;
    vector<double> tle_density;
    std::vector<double> meshAnalogColl;
    inline int idx(int i,int j,int k) const { return (k*ny + j)*nx + i; }
    inline bool inside(const array<double,3>& p) const {
        return (p[0]>=pmin[0] && p[0]<=pmax[0] &&
                p[1]>=pmin[1] && p[1]<=pmax[1] &&
                p[2]>=pmin[2] && p[2]<=pmax[2]);
    }
    inline void zero() {
        cfe_density.assign(std::max(1,nx*ny*nz), 0.0);
        tle_density.assign(std::max(1,nx*ny*nz), 0.0);
    }
} Mesh3D;

typedef struct StatsOut{ 
    vector<vector<double>> mean; 
    vector<vector<double>> relErr; 
    vector<vector<int>> sum; 
} StatsOut;

typedef struct RxSample {
    int mt = 0;
    int matIndex = -1;
    const Material* mat = nullptr;
} RxSample;

typedef struct FlightResult {
    bool leaked=false;
    bool collided=false;
    bool virtualCollision=false;
    array<double,3> pos{};
    const Geometry* geom=nullptr;
    double SigmaLocal=0.0;
    double traveled=0.0;
} FlightResult;

typedef struct TallyBook{
    vector<std::string> matNames;
    std::map<std::string,int> mat2idx;
    vector<vector<vector<int>>> statM;
    TimeHist* timeHist = nullptr;
    Mesh3D* mesh = nullptr;
    bool useTLE = false;
    bool useCFE = true;
    bool deltaMode = false;
    double SigmaM = 0.0;
    int Rcols = 0;

    std::vector<int> leaks;
    std::vector<double> cfe_global_time;
    std::vector<std::vector<double>> cfe_Rtot;
    std::vector<std::vector<double>> cfe_Rabs;
    std::vector<std::vector<double>> tle_Rtot;
    std::vector<std::vector<double>> tle_Rabs;

    void ensureBatchAll(int batch, int M) {
        ensureBatch(batch, M, Rcols);
        if ((int)leaks.size() <= batch) leaks.resize(batch+1, 0);
        if ((int)cfe_global_time.size() <= batch) cfe_global_time.resize(batch+1, 0.0);
        auto grow2d = [&](std::vector<std::vector<double>>& A) {
            if ((int)A.size() <= batch) A.resize(batch+1);
            if ((int)A[batch].size() < M) A[batch].resize(M, 0.0);
        };
        grow2d(cfe_Rtot); grow2d(cfe_Rabs); grow2d(tle_Rtot); grow2d(tle_Rabs);
    }

    void ensureBatch(int batch, int M, int R) {
        if ((int)statM.size() <= batch) statM.resize(batch+1);
        Rcols = std::max(Rcols, R);
        auto& B = statM[batch];
        if ((int)B.size() < M) B.resize(M, std::vector<int>(Rcols, 0));
        for (auto& row : B) if ((int)row.size() < Rcols) row.resize(Rcols, 0);
    }
    int matIndex(const Material& m) {
        auto it = mat2idx.find(m.sym);
        if (it!=mat2idx.end()) return it->second;
        int id = (int)matNames.size();
        mat2idx[m.sym]=id; matNames.push_back(m.sym);
        for (auto& B: statM) {
            int cols = B.empty() ? Rcols : (int)B[0].size();
            B.push_back(std::vector<int>(cols, 0));
        }
        // extend new material column for already existing batches
        auto growCol = [&](std::vector<std::vector<double>>& A) {
            for (auto& row: A) if ((int)row.size() <= id) row.resize(id+1, 0.0);
        };
        growCol(cfe_Rtot); growCol(cfe_Rabs); growCol(tle_Rtot); growCol(tle_Rabs);
        return id;
    }
} TallyBook;

typedef struct RunParams{
    Tracking track = Tracking::Surface;
    SourceMode src = SourceMode::External;
    int historiesPerBatch = 10000;
    int batches = 10;
    int maxSteps = 100000;
    bool inelastic = true;
    double sourceE = 1.0;
    array<double,3> sourcePos{0,0,0};
    Mesh3D* mesh = nullptr;
    double sourceRate = 1.0e6;
} RunParams;

typedef struct PerfOut{
    double elapsed_s = 0.0;     // wall time for the run
    long long histories = 0;    // total source histories simulated
} PerfOut;

typedef struct RunOutputs{
    TallyBook T;
    vector<Collisions> collisions;
    vector<int> fissionChildren;
    vector<double> keff_history;
    PerfOut perf; 
} RunOutputs;


const Geometry* findGeometryAt(const Universe&, const std::array<double,3>& pWorld);
const Geometry* findGeomAtRecursive(const Universe&, const std::array<double,3>& pLocal);
bool pointInGeom(const std::array<double,3>&, const Geometry&);
std::array<double,3> boundingBox(Universe&);
std::array<double,3> squareCellCenter(const Universe&, int, int);
std::array<double,3> hexCellCenter(const Universe&, int, int);
double sigmaTot(const std::vector<Material>&, double E, const std::vector<int>& mts_total);
double macroscopicSigmaAt(const Universe&, const std::array<double,3>& pWorld, double E, const std::vector<int>& mts_total, const Geometry** gOut);
double majorSigma(const Universe&, double E, const std::vector<int>& mts_total);
double nearestCollision(const Neutron&, const Universe&);
bool readUniverseFile(std::istream&, Universe&);
bool writeVTKStructuredPoints(const Mesh3D&, const std::vector<double>&, const std::string& vtkPath, const std::string& fieldName);
bool storeDatakeff(const std::vector<double>&);
RunOutputs runExternal(const Universe&, const RunParams&);
RunOutputs runCriticality(const Universe&, const RunParams&, int cycles);
std::vector<double> linspace(double a, double b, int n);
std::vector<double> logspace(double a, double b, int n);
double valueInterp(const std::vector<std::pair<double,double>>& xy, double x);
int rxCol(const std::vector<int>& MTs, int mt);
void printStatsOut(const StatsOut&, const std::vector<std::string>& labels, const std::vector<int>& mts_total, std::ostream& os);
void volumePointMethod(Universe& u, int iter);
void volumeLineMethod(Universe& u, int iter);
void volumeLineMethodTorus(Universe& u, int iter);
