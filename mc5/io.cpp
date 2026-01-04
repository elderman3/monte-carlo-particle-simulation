#include "mc.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <cctype>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <cstdlib>


using std::array; using std::vector;
using Matrix = vector<vector<float>>;
using std::cos; using std::sin;


// --- MATERIAL READ ---

bool readNuclideBlock(std::istream& in, Nuclide& nuc) {
    if (!(in >> nuc.sym >> nuc.z >> nuc.a >> nuc.aw >> nuc.T))
        return false;

    int neu_num;
    if (!(in >> neu_num)) return false;
    nuc.neutrons.clear();
    nuc.neutrons.reserve(neu_num);
    for (int i=0; i<neu_num; ++i) {
        float E, nu;
        if (!(in >> E >> nu)) return false;
        nuc.neutrons.emplace_back(E, nu);
    }

    int mt;
    float Q;
    int nc;
    nuc.mt.clear();
    nuc.Qvals.clear();
    while (in >> mt >> Q >> nc) {
        nuc.Qvals[mt] = Q;
        auto& xs = nuc.mt[mt];
        xs.clear();
        xs.reserve(nc);
        for (int j=0; j<nc; ++j) {
            float E, sigma;
            if (!(in >> E >> sigma)) return false;
            xs.emplace_back(E, sigma);
        }
    }
    return true;
}

static void buildDerivedInelasticTotal(Nuclide& nuc, const std::vector<float>& Egrid) {
    // Build MT4 total inelastic as sum of MT51-91 channels.
    // This mirrors the previous fillData() behavior but runs once per nuclide.
    std::vector<int> inel;
    inel.reserve(64);
    for (const auto& kv : nuc.mt) {
        const int mt = kv.first;
        if (mt > 50 && mt < 92) inel.push_back(mt);
    }

    auto& out4 = nuc.mt[4];
    out4.clear();
    out4.reserve(Egrid.size());

    for (float E : Egrid) {
        float sum4 = 0.0f;
        for (int mt : inel) {
            auto it = nuc.mt.find(mt);
            if (it != nuc.mt.end()) sum4 += valueInterp(it->second, E);
        }
        out4.emplace_back(E, sum4);
    }
}


// --- MATERIAL CACHES ---

namespace {
    struct MaterialSetRec {
        int inelastic = 1;
        std::vector<Material> comps; // lightweight components (nuclide pointer + rho + proportion)
    };

    // Cache: material-set filename -> set id
    static std::vector<MaterialSetRec> g_materialSets;
    static std::unordered_map<std::string,int> g_materialSetByName;

    // Cache: nuclide token (data/<token>.dat) -> unique nuclide data
    static std::vector<std::unique_ptr<Nuclide>> g_nuclides;
    static std::unordered_map<std::string,int> g_nuclideByToken;

    static const std::vector<float>& commonEgrid() {
        static std::vector<float> grid = logspace(-11.0, std::log10(20.0f), 500);
        return grid;
    }

    static const Nuclide* loadNuclideCached(const std::string& token) {
        auto it = g_nuclideByToken.find(token);
        if (it != g_nuclideByToken.end()) return g_nuclides[it->second].get();

        const std::string path = "data/" + token + ".dat";
        std::ifstream materialdata(path);
        if (!materialdata) {
            std::cerr << "Data file opening failed: " << path << " (referenced by material file)\n";
            std::exit(EXIT_FAILURE);
        }

        auto nuc = std::make_unique<Nuclide>();
        if (!readNuclideBlock(materialdata, *nuc)) {
            std::cerr << "Block read fail in " << path << "\n";
            std::exit(EXIT_FAILURE);
        }

        // Derived MT4 (inelastic sum) on a common E-grid.
        buildDerivedInelasticTotal(*nuc, commonEgrid());

        const int id = (int)g_nuclides.size();
        g_nuclides.push_back(std::move(nuc));
        g_nuclideByToken[token] = id;
        return g_nuclides[id].get();
    }

    static inline MatView viewOfMaterialSet(int id) {
        const auto& v = g_materialSets[id].comps;
        MatView w; w.ptr = v.empty() ? nullptr : v.data(); w.n = (int)v.size();
        return w;
    }

    static int loadMaterialSetCached(const std::string& filename) {
        auto it = g_materialSetByName.find(filename);
        if (it != g_materialSetByName.end()) return it->second;

        std::ifstream file("material/" + filename + ".txt");
        if (!file) {
            std::cerr << "MatFile opening failed: material/" << filename << ".txt\n";
            std::exit(EXIT_FAILURE);
        }

        int inelastic = 1, nMaterials = 0;
        if (!(file >> inelastic >> nMaterials)) {
            std::cerr << "Failed reading header from material/" << filename << ".txt\n";
            std::exit(EXIT_FAILURE);
        }

        std::vector<Material> comps;
        comps.reserve(std::max(0, nMaterials));

        for (int i = 0; i < nMaterials; ++i) {
            std::string token; float rho = 0.0f, relativeMoles = 0.0f;
            if (!(file >> token >> rho >> relativeMoles)) {
                std::cerr << "Failed reading material entry " << i << " from material/" << filename << ".txt\n";
                std::exit(EXIT_FAILURE);
            }

            const Nuclide* nuc = loadNuclideCached(token);
            Material m;
            m.nuc = nuc;
            m.rho = rho;
            m.proportion = relativeMoles;
            comps.push_back(m);
        }

        MaterialSetRec rec;
        rec.inelastic = inelastic;
        rec.comps = std::move(comps);

        const int id = (int)g_materialSets.size();
        g_materialSets.push_back(std::move(rec));
        g_materialSetByName[filename] = id;
        return id;
    }
}

MatView readMaterial(const std::string& filename, int* outId = nullptr) {
    const int id = loadMaterialSetCached(filename);
    if (outId) *outId = id;
    return viewOfMaterialSet(id);
}

// --- GEOMETRY MATH ---

Matrix mmult(const Matrix& a, const Matrix& b) {
    if (a.empty() || b.empty()) throw std::runtime_error("empty matrix\n");
    size_t ai = a.size(); size_t aj = a[0].size(); size_t bi = b.size(); size_t bj = b[0].size();
    Matrix c(ai, vector<float>(bj, 0.f));
    if (aj != bi) throw std::runtime_error("Mismatched Matrices\n");
    for (std::size_t i = 0; i < ai; ++i) {
        for (std::size_t k = 0; k < aj; ++k) {
            const float aik = a[i][k];
            for (std::size_t j = 0; j < bj; ++j) {
                c[i][j] += aik * b[k][j];
            }
        }
    }
    return c;
}

Matrix msum(const Matrix& a, const Matrix& b) {
    if (a.empty() || b.empty()) throw std::runtime_error("empty matrix\n");
    size_t ai = a.size(); size_t aj = a[0].size(); size_t bi = b.size(); size_t bj = b[0].size();
    Matrix c(ai, vector<float>(aj, 0.f));
    if (ai != bi || aj != bj) throw std::runtime_error("Mismatched Matrices\n");
    for (size_t i = 0; i < ai; ++i) {for (size_t j = 0; j < aj; ++j) {c[i][j] = a[i][j] + b[i][j];}}
    return c;
}

Matrix T(const Matrix& a) {
    size_t ai = a.size(); size_t aj = a[0].size();
    Matrix c(aj, vector<float>(ai, 0.f));
    for (size_t i = 0; i < ai; ++i) {for (size_t j = 0; j < aj; ++j) {c[j][i] = a[i][j];}}
    return c;
}

void transformGeometry(Geometry& g, array<float, 3> pos, array<float, 3> rot) {
    float phi = rot[0]; float th = rot[1]; float psi = rot[2]; // x, y, z
    const float cPsi = cos(psi), sPsi = sin(psi);
    const float cTh = cos(th),  sTh = sin(th);
    const float cPhi = cos(phi), sPhi = sin(phi);
    Matrix R = {
            { cPsi*cTh, cPsi*sTh*sPhi - sPsi*cPhi, cPsi*sTh*cPhi + sPsi*sPhi },
            { sPsi*cTh, sPsi*sTh*sPhi + cPsi*cPhi, sPsi*sTh*cPhi - cPsi*sPhi },
            { -sTh,   cTh*sPhi,            cTh*cPhi            }
                };
    Matrix t = {
            {pos[0]},
            {pos[1]},
            {pos[2]}
                };
    
    for (Shape& s : g.shapes) {
        if (!s.torus) {
            Matrix A = {
                { s.A, s.D/2, s.F/2 },
                { s.D/2, s.B, s.E/2 },
                { s.F/2, s.E/2, s.C }
                        };
            Matrix b = {
                {s.G/2},
                {s.H/2},
                {s.I/2}
                    };
            float c = s.J;
            Matrix An = mmult(mmult(T(R), A), R);
            Matrix bn = mmult(T(R), msum(mmult(A, t), b));
            float cn = mmult(mmult(T(t), A), t)[0][0] + 2 * mmult(T(b), t)[0][0] + c;
            s.A=An[0][0]; s.B=An[1][1]; s.C=An[2][2]; s.D=2*An[0][1]; s.E=2*An[1][2]; s.F=2*An[0][2]; s.G=2*bn[0][0]; s.H=2*bn[1][0]; s.I=2*bn[2][0]; s.J=cn;
        } else {
            Matrix c = {
                {s.D},
                {s.E},
                {s.F}
                    };
            Matrix u = {
                {s.G},
                {s.H},
                {s.I}
                    };
            Matrix cn = msum(mmult(R, c), t);
            Matrix un = mmult(R, u);
            s.D = cn[0][0]; s.E = cn[1][0]; s.F = cn[2][0]; s.G = un[0][0]; s.H = un[1][0]; s.I = un[2][0];
        }
    }
}


// --- GEOMETRY READ ---

using CreateFn = void(*)(float,float,float,Geometry&);
static const vector<CreateFn> kCreateById = {
  nullptr,            // 0 = general, special handling
  &createBall,        // 1
  &createCylinder,    // 2
  &createCylinderOpen,// 3
  &createPlane,       // 4
  &createCuboid,      // 5
  &createCuboidOpen,  // 6
  &createHexPrism,    // 7
  &createTorus        // 8
};

int readNodeDef(std::string& str, vector<Node>& nodes) {
    vector<int> st;
    auto push = [&](Node n) { nodes.push_back(n); st.push_back((int)nodes.size()-1); };

    std::istringstream iss(str);
    std::string tok;
    while (iss >> tok) {
        if (tok=="n") {
            if (st.size()<1) throw std::runtime_error("underflow: n");
                int a=st.back(); st.pop_back();
                push({Op::N,-1,a,-1});
        } else if (tok=="u" || tok=="i") {
            if (st.size()<2) throw std::runtime_error("underflow: binop");
                int b=st.back(); st.pop_back();
                int a=st.back(); st.pop_back();
                push({tok=="u"?Op::U:Op::I,-1,a,b});
        } else {
            size_t pos=0;
            int id = std::stoi(tok, &pos, 10);
            if (pos!=tok.size()) throw std::runtime_error("bad token: "+tok);
                push({Op::L,id,-1,-1});
        }
    }
    if (st.size()!=1) throw std::runtime_error("invalid expression");
    return st.back();
}

bool readShape(std::istream& in, Geometry& g) {
    Shape s;
    if (!(in >> s.A >> s.B >> s.C >> s.D >> s.E >> s.F >> s.G >> s.H >> s.I >> s.J >> s.torus)) { std::cout << "Failed reading line10\n"; return false; }
    g.shapes.push_back(std::move(s));
    return true;
}

bool readGeometry(std::istream& in, Geometry& g) {
    std::string command; float a, b, c; array<float, 3> pos; array<float, 3> rot;
    if (!(in >> command >> a >> b >> c >> pos[0] >> pos[1] >> pos[2] >> rot[0] >> rot[1] >> rot[2])) { std::cout << "Failed reading line1\n"; return false; }
    if (command == "general") {
        int nShapes; 
        if (!(in >> nShapes)) { std::cout << "Failed reading line2\n"; return false; }

        std::string nodeLine;
        std::getline(in, nodeLine);
        if (!std::getline(in, nodeLine)) { std::cout << "Failed reading line3\n"; return false; }

        for (int i = 0; i < nShapes; ++i) {
            if (!readShape(in, g)) { std::cout << "Failed reading line11\n"; return false; }
        }

        g.shape = 0;
        vector<Node> nodes;
        g.nodeRoot = readNodeDef(nodeLine, nodes);
        g.nodes = std::move(nodes);
    } else if (command == "ball") {
        createBall(a, b, c, g);
    } else if (command == "cylinder") {
        createCylinder(a, b, c, g);
    } else if (command == "cylinderOpen") {
        createCylinderOpen(a, b, c, g);
    } else if (command == "plane") {
        createPlane(a, b, c, g);
    } else if (command == "cuboid") {
        createCuboid(a, b, c, g);
    } else if (command == "cuboidOpen") {
        createCuboidOpen(a, b, c, g);
    } else if (command == "hexPrism") {
        createHexPrism(a, b, c, g);
    } else if (command == "torus") {
        createTorus(a, b, c, g);
    } else {
        std::cout << "No shape found\n";
    }
    std::string fname;
    if (!(in >> fname)) { std::cout << "Failed to read material file name"; }
    g.mats = readMaterial(fname, &g.matSetId);

    transformGeometry(g, pos, rot);
    return true;
}

bool readUniverseFile(std::istream& in, Universe& u) {
    int nSubUniverse; int nGeometry;
    // Read basic universe information
    array<float, 3> uniBou;
    if (!(in >> u.name >> u.pos[0] >> u.pos[1] >> u.pos[2] >> u.rot[0] >> u.rot[1] >> u.rot[2] 
        >> u.latticeType >> nSubUniverse >> nGeometry >> uniBou[0] >> uniBou[1] >> uniBou[2] >> u.universeShape)) { std::cout << "Failed reading line6\n"; return false; }
    if (!u.universeShape) { std::cout << "Failed reading line4\n"; return false; }
    Geometry bounds;
    u.boundDim = uniBou;
    kCreateById[u.universeShape](uniBou[0], uniBou[1], uniBou[2], bounds);
    u.boundingGeometry = bounds;
    if (u.latticeType) {
        if (!(in >> u.lattice[0] >> u.lattice[1])) { std::cout << "Failed reading line5\n"; return false; }
    }
    // Read subuniverses

    array<float, 3> SUpos; array<float, 3> SUrot; std::string SUfname;
    for (int i = 0; i < nSubUniverse; i++) {
        Universe su;
        if (!(in >> SUfname >> SUpos[0] >> SUpos[1] >> SUpos[2] >> SUrot[0] >> SUrot[1] >> SUrot[2])) { std::cout << "Failed reading line7\n"; return false; }
        std::string path = "geometry/" + SUfname + ".txt";
        std::ifstream SUfile(path);
        if (!(readUniverseFile(SUfile, su))) { std::cout << "Failed reading line8\n"; return false; }
        su.pos = SUpos; su.rot = SUrot;
        u.subUniverse.push_back(std::move(su));
    }

    for (int i = 0; i < nGeometry; i++) {
        Geometry g;
        if (!(readGeometry(in, g))) { std::cout << "Failed reading line9\n"; return false; }
        u.geometries.push_back(std::move(g));
    }

    return true;    
}


// --- STORE HELPERS ---

std::string sanitizeSym(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : s) {
        if (std::isalnum(c)) out.push_back((char)std::toupper(c));
    }
    return out;
}

// Coarse geometry/material tagging for visualization
// 0 = outside/void
// 1 = water-like
// 2 = fuel-like
// 3 = zirconium-like
// 4 = steel/iron-like
// 5 = absorber/boron-like
// 6 = other
int geometryTag(const Geometry* g) {
    if (!g) return 0;
    bool hasH=false, hasO=false, hasU=false, hasPu=false, hasZr=false, hasFe=false, hasB=false, hasC=false;
    for (const auto& m : g->mats) {
        const std::string sym = sanitizeSym(m.nuc ? m.nuc->sym : std::string());
        if (sym.rfind("H",0)==0) hasH=true;
        if (sym.rfind("O",0)==0) hasO=true;
        if (sym.rfind("U",0)==0) hasU=true;
        if (sym.rfind("PU",0)==0) hasPu=true;
        if (sym.rfind("ZR",0)==0) hasZr=true;
        if (sym.rfind("FE",0)==0) hasFe=true;
        if (sym.rfind("B",0)==0) hasB=true;
        if (sym.rfind("C",0)==0) hasC=true;
    }
    if (hasU || hasPu) return 2;
    if (hasZr) return 3;
    if (hasFe) return 4;
    if (hasB) return 5;
    if (hasH || hasO) return 1;
    if (hasC) return 5;
    return 6;
}


// --- STORE DATA ---

bool storeDatakeff(const vector<float>& data) {
    using clock = std::chrono::system_clock;
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(clock::now().time_since_epoch()).count();

    std::ofstream os("output/keff_" + std::to_string(secs) + ".csv");
    if (!os) return false;
    os << std::scientific << std::setprecision(6);
    for (size_t k = 0; k < data.size(); ++k)
        os << k << "," << data[k] << "\n";
    return true;
}

bool writeVTKStructuredPoints(const Mesh3D& M, const vector<float>& field, const std::string& basePath, const std::string& name) {
    std::string path = "output/" + basePath + ".vtk";
    std::ofstream os(path);
    if (!os) return false;
    os << "# vtk DataFile Version 3.0\n" << name << "\nASCII\nDATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << M.nx << " " << M.ny << " " << M.nz << "\n";
    os << "ORIGIN " << M.pmin[0] << " " << M.pmin[1] << " " << M.pmin[2] << "\n";
    os << "SPACING "
       << (M.pmax[0]-M.pmin[0])/std::max(1,M.nx-1) << " "
       << (M.pmax[1]-M.pmin[1])/std::max(1,M.ny-1) << " "
       << (M.pmax[2]-M.pmin[2])/std::max(1,M.nz-1) << "\n";
    const int N = std::max(1, M.nx*M.ny*M.nz);
    os << "POINT_DATA " << N << "\n";
    os << "SCALARS " << name << " float 1\nLOOKUP_TABLE default\n";
    for (int k=0;k<N;++k) os << field[k] << "\n";
    return true;
}

std::vector<float> buildGeometryTagField(const Mesh3D& M, const Universe& U) {
    std::vector<float> out((size_t)M.nx * M.ny * M.nz, 0.0f);
    const float dx = (M.pmax[0] - M.pmin[0]) / (float)M.nx;
    const float dy = (M.pmax[1] - M.pmin[1]) / (float)M.ny;
    const float dz = (M.pmax[2] - M.pmin[2]) / (float)M.nz;
    for (int k=0;k<M.nz;++k) {
        const float z = M.pmin[2] + (k + 0.5) * dz;
        for (int j=0;j<M.ny;++j) {
            const float y = M.pmin[1] + (j + 0.5) * dy;
            for (int i=0;i<M.nx;++i) {
                const float x = M.pmin[0] + (i + 0.5) * dx;
                const Geometry* g = findGeometryAt(U, {x,y,z});
                out[(size_t)(k*M.ny + j)*M.nx + i] = (float)geometryTag(g);
            }
        }
    }
    return out;
}