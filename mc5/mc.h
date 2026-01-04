#pragma once
#include <vector>
#include <utility>
#include <string>
#include <map>
#include <array>
#include <deque>
#include <iostream>
#include <cstdint>
#include <cmath>

using std::array; using std::vector;

#ifndef __CUDACC__
  #ifndef __host__
    #define __host__
  #endif
  #ifndef __device__
    #define __device__
  #endif
  #ifndef __global__
    #define __global__
  #endif
  #ifndef __shared__
    #define __shared__
  #endif
  #ifndef __align__
    #define __align__(n)
  #endif
  #ifndef __launch_bounds__
    #define __launch_bounds__(t, b)
  #endif
#endif


enum class Op : std::uint8_t { L=0, N=1, U=2, I=3 };

typedef struct Nuclide{
    std::string sym{};
    int z = 0;
    int a = 0;
    float aw = 0.0f;
    int T = 0;
    vector<std::pair<float, float>> neutrons;
    std::map<int, vector<std::pair<float,float>>> mt;
    std::map<int, float> Qvals;
} Nuclide;

// A mixture component: points to cached nuclide data + per-mixture properties
typedef struct Material{
    const Nuclide* nuc = nullptr; // non-owning pointer to cached nuclide data
    float rho = 0.0f;
    float proportion = 0.0f; // molar proportion of composition
} Material;

// Material caching
typedef struct MatView {
    const Material* ptr = nullptr;
    int n = 0;

    inline bool empty() const { return n <= 0 || ptr == nullptr; }
    inline int size() const { return n; }
    inline const Material& operator[](int i) const { return ptr[i]; }
    inline const Material* data() const { return ptr; }
    inline const Material* begin() const { return ptr; }
    inline const Material* end() const { return ptr + n; }
} MatView;

typedef struct Node{
  Op op;
  int shape; // Index of Shape
  int left; // child index or -1
  int right; // child index or -1
} Node;

// A single Plane of a Geometry (Quadratic and Torus)
typedef struct Shape{
    bool torus = 0; // For torus, A, B, C are used as a,b,R. Additionally, DEF are center and GHI Axis
    float A,B,C,D,E,F,G,H,I,J;
} Shape;

// A single compound shape e.g., Cube, Cylinder
typedef struct Geometry{
    array<float, 3> pos{0.f, 0.f, 0.f};
    array<float, 3> rot{0.f, 0.f, 0.f}; // Always keep these 0, Universes are not rotatable. 
    vector<Shape> shapes;
    int shape; // Defines the possible shape, tabled elsewhere
    int nodeRoot;
    MatView mats; // Vector for storing the material that this geometry is made of
    int matSetId = -1;
    vector<Node> nodes; // Set of operations that define the shape of Geometry
    // binary tree of intersections etc to define the shape fully
} Geometry;

// A collection of Universes and Geometries
typedef struct Universe{
    std::string name;
    array<float, 3> pos{0.f, 0.f, 0.f};
    array<float, 3> rot{0.f, 0.f, 0.f}; // Always keep these 0, Universes are not rotatable. 
    array<float, 3> boundDim{0.f, 0.f, 0.f}; // Bounding box construction args
    Geometry boundingGeometry; // Universe bounding box
    int universeShape; // Defines how the Universe Bounding box is constructed
    array<int, 2> lattice{0, 0}; // Only valid for square or hex major Universes
    int latticeType = 0; // 0 For none, 1 for square, 2 for hex
    vector<Universe> subUniverse;
    vector<Geometry> geometries;
} Universe;

typedef struct Neutron{
    array<float, 3> pos{0.0f,0.0f,0.0f};
    array<float, 3> dir{0.0f,0.0f,1.0}; // unit
    float vel = 0.0f; 
    float energy = 0.0f; // MeV
    int collisions = 0;
    float time = 0.0f; // s
    bool reachedTh = false;
    bool isSource = true;
    float w = 1.0; // statistical weight (for implicit crit)
} Neutron;

enum class Tracking { Surface, Delta };
enum class SourceMode { External, Criticality };

typedef struct Mesh3D{
    array<float,3> pmin{0,0,0}, pmax{0,0,0};
    int nx=0, ny=0, nz=0;
    vector<float> cfe_density;
    vector<float> tle_density;
    std::vector<float> meshAnalogColl;
    inline int idx(int i,int j,int k) const { return (k*ny + j)*nx + i; }
    inline bool inside(const array<float,3>& p) const {
        return (p[0]>=pmin[0] && p[0]<=pmax[0] &&
                p[1]>=pmin[1] && p[1]<=pmax[1] &&
                p[2]>=pmin[2] && p[2]<=pmax[2]);
    }
    inline void zero() {
        cfe_density.assign(std::max(1,nx*ny*nz), 0.0f);
        tle_density.assign(std::max(1,nx*ny*nz), 0.0f);
    }
} Mesh3D;

typedef struct StatsOut{ 
    vector<vector<float>> mean; 
    vector<vector<float>> relErr; 
    vector<vector<int>> sum; 
} StatsOut;

typedef struct TallyBook{
    vector<std::string> matNames;
    std::map<std::string,int> mat2idx;
    vector<vector<vector<int>>> statM;
    Mesh3D* mesh = nullptr;
    bool useTLE = false;
    bool useCFE = true;
    bool deltaMode = false;
    float SigmaM = 0.0f;
    int Rcols = 0;

    std::vector<int> leaks;
    std::vector<float> cfe_global_time;
    std::vector<std::vector<float>> cfe_Rtot;
    std::vector<std::vector<float>> cfe_Rabs;
    std::vector<std::vector<float>> tle_Rtot;
    std::vector<std::vector<float>> tle_Rabs;

    // Makes sure all tallies have space for batch
    void ensureBatchAll(int batch, int M) {
        ensureBatch(batch, M, Rcols);
        if ((int)leaks.size() <= batch) leaks.resize(batch+1, 0);
        if ((int)cfe_global_time.size() <= batch) cfe_global_time.resize(batch+1, 0.0f);
        auto grow2d = [&](std::vector<std::vector<float>>& A) {
            if ((int)A.size() <= batch) A.resize(batch+1);
            if ((int)A[batch].size() < M) A[batch].resize(M, 0.0f);
        };
        grow2d(cfe_Rtot); grow2d(cfe_Rabs); grow2d(tle_Rtot); grow2d(tle_Rabs);
    }

    // Makes sure all tallies have space for batch
    void ensureBatch(int batch, int M, int R) {
        if ((int)statM.size() <= batch) statM.resize(batch+1);
        Rcols = std::max(Rcols, R);
        auto& B = statM[batch];
        if ((int)B.size() < M) B.resize(M, std::vector<int>(Rcols, 0));
        for (auto& row : B) if ((int)row.size() < Rcols) row.resize(Rcols, 0);
    }

    // Get material index, adding new material if needed
    int matIndex(const Material& m) {
        const std::string key = (m.nuc ? m.nuc->sym : std::string());
        auto it = mat2idx.find(key);
        if (it!=mat2idx.end()) return it->second;
        int id = (int)matNames.size();
        mat2idx[key]=id; matNames.push_back(key);
        for (auto& B: statM) {
            int cols = B.empty() ? Rcols : (int)B[0].size();
            B.push_back(std::vector<int>(cols, 0));
        }
        // extend new material column for already existing batches
        auto growCol = [&](std::vector<std::vector<float>>& A) {
            for (auto& row: A) if ((int)row.size() <= id) row.resize(id+1, 0.0f);
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
    float sourceE = 1.0;
    array<float,3> sourcePos{0,0,0};
    Mesh3D* mesh = nullptr;
    float sourceRate = 1.0e6;
} RunParams;

typedef struct PerfOut{
    float elapsed_s = 0.0f;     // wall time for the run
    long long histories = 0;    // total source histories simulated
} PerfOut;

typedef struct RunOutputs{
    TallyBook T;
    vector<int> fissionChildren;
    vector<float> keff_history;
    PerfOut perf; 
} RunOutputs;

typedef struct ScalarStat {
    float mean=0, stddev=0, relErr=0, fom=0;
    int batches=0;
} ScalarStat;

typedef struct RegionMap { 
    const char* label; 
    std::vector<std::string> mats; 
} RegionMap;

typedef struct  BatchMetrics {
    std::vector<float> leak_nps;
    std::vector<float> abs_water_nps;
    std::vector<float> abs_steel_nps;
    std::vector<float> tle_abs_water_nps;
    std::vector<float> cfe_abs_water_nps;
    std::vector<float> cfe_tau_s;
} BatchMetrics;

// simulation.cu
array<float,3> scaleByC(const array<float,3>& a, float s);
vector<float> logspace(float ea,float eb,int n);
float neutronSpeed(float E);
float valueInterp(const vector<std::pair<float,float>>& data, float target);
RunOutputs runExternal(const Universe& U, const RunParams& P);
RunOutputs runCriticality(const Universe& U, const RunParams& P, int inactive=5);

// simulation.cu CUDA
RunOutputs runExternalCuda(const Universe& U, const RunParams& P);
RunOutputs runCriticalityCuda(const Universe& U, const RunParams& P, int inactive=5);
RunOutputs runCuda(const Universe& U, const RunParams& P, int inactive=5);
bool compiledWithCuda();

// geomops.cpp
float dot3(const array<float,3>& a, const array<float,3>& b);
array<float,3> add3(const array<float,3>& a,const array<float,3>& b);
array<float,3> sub3(const array<float,3>& a,const array<float,3>& b);
array<float,3> cross3(const array<float,3>& a, const array<float,3>& b);
array<float,3> madd(const array<float,3>& p, const array<float,3>& d, float s);
array<float,3> normed(const array<float,3>& v);
int solveQuartic(float a,float b,float c,float d,float e,float r[4]);
bool pointInGeom(const array<float, 3>& pos, const Geometry& g);
float geometryCollision(const Neutron& n, const Geometry& g);
float nearestCollision(const Neutron& n, const Universe& u);
array<float,3> cellBoxDims(const Universe& u);
array<float, 3> boundingBox(Universe& u);
array<float,3> squareCellCenter(const Universe& u, int i, int j);
array<float,3> hexCellCenter(const Universe& u, int q, int r);
const Geometry* findGeometryAt(const Universe& u, const array<float,3>& pWorld);


// Function declarations

// geominput.cpp
void createBall(float a, float b, float c, Geometry& g);
void createCube(float a, float b, float c, Geometry& g);
void createCylinder(float a, float b, float c, Geometry& g);
void createCylinderOpen(float a, float b, float c, Geometry& g);
void createPlane(float a, float b, float c, Geometry& g);
void createCuboid(float a, float b, float c, Geometry& g);
void createCuboidOpen(float a, float b, float c, Geometry& g);
void createHexPrism(float a, float b, float c, Geometry& g);
void createTorus(float a, float b, float c, Geometry& g);

// io.cpp
bool readUniverseFile(std::istream& in, Universe& u);
bool storeDatakeff(const vector<float>& data);
bool writeVTKStructuredPoints(const Mesh3D& M, const vector<float>& field, const std::string& basePath, const std::string& name);
std::vector<float> buildGeometryTagField(const Mesh3D& M, const Universe& U);

// stats.cpp
StatsOut computeStats(const vector<vector<vector<int>>>& statM);
void printStatsOut(const StatsOut& S, const vector<std::string>& matNames, const vector<int>& MTs, std::ostream& os = std::cout);
void estimateNeutronDensity1W(const RunParams& P, const RunOutputs& R, float mixtureVolume_m3, float meanNuBar = 2.43);
void printFoms(const RunParams& P, const RunOutputs& R, const std::vector<int>& MTs, const std::string& tag);
Mesh3D makeMeshFromUniverse(const Universe& U, int nx,int ny,int nz);
void printAndStoreRates(const RunParams& P, const RunOutputs& R, const StatsOut& S, const std::string& tag);
void tallyTrackToMesh(const Mesh3D& M, vector<float>& dst, const array<float,3>& p0, const array<float,3>& dHat, float segLen, float weightPerLength);
int vindex(const Mesh3D& M, float x, float y, float z);
void scoreCFE(TallyBook& T, int batch, const array<float,3>& pos, float E, float SigmaRef);
void scoreCFERxPerMaterial(TallyBook& T, int batch, int mi, const Material& m, float E, const std::vector<int>& mts_total, float SigmaRef);
void scoreTLESegmentPerGeom(TallyBook& T, int batch, const Geometry& g0, float E, float segLen, const std::vector<int>& mts_total);
void scoreCFEDensityGlobal(TallyBook& T, int batch, float E, float SigmaRef);

// volume.cpp
void volumePointMethod(Universe& u, int iter);
void volumeLineMethod(Universe& u, int iter);
void volumeLineMethodTorus(Universe& u, int iter);
