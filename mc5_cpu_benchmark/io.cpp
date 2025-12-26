#include "mc.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <chrono>

using std::array; using std::vector;
using Matrix = vector<vector<double>>;
using std::cos; using std::sin;


// --- MATERIAL READ ---

bool readMaterialBlock(std::istream& in, Material& mat) {
    if (!(in >> mat.sym >> mat.z >> mat.a >> mat.aw >> mat.T))
        return false;

    int neu_num;
    if (!(in >> neu_num)) return false;
    mat.neutrons.reserve(neu_num);
    for (int i=0; i<neu_num; ++i) {
        double E, nu;
        if (!(in >> E >> nu)) return false;
        mat.neutrons.emplace_back(E, nu);
    }

    int mt;
    double Q;
    int nc;
    while (in >> mt >> Q >> nc) {
        mat.Qvals[mt] = Q;
        auto& xs = mat.mt[mt];
        xs.reserve(nc);
        for (int j=0; j<nc; ++j) {
            double E, sigma;
            if (!(in >> E >> sigma)) return false;
            xs.emplace_back(E, sigma);
        }
    }
    return true;
}

void fillData(vector<Material>& mats, vector<double>& x, int inelastic) {
    for (auto& m : mats) {
        auto& out1 = m.mt[1];
        auto& out4 = m.mt[4];
        out1.reserve(x.size());
        out4.reserve(x.size());
        // Get all availabe MT values
        vector<int> MTs;
        for (const auto& [k,v] : m.mt) if (k != 1 && k != 4) MTs.push_back(k);
        for (double d : x) {
            double sum1 = 0, sum4 = 0;
            for (auto mt : MTs) {
                double val = valueInterp(m.mt[mt], d);
                if (mt > 50 && mt < 92) {
                    sum4 += val;
                } else {
                    sum1 += val;
                }
            }
            if (!inelastic) {
                sum1 += sum4;
            }
        out1.emplace_back(d, sum1);
        out4.emplace_back(d, sum4);
        }
    }
}

// --- GEOMETRY MATH ---

Matrix mmult(const Matrix& a, const Matrix& b) {
    if (a.empty() || b.empty()) throw std::runtime_error("empty matrix\n");
    size_t ai = a.size(); size_t aj = a[0].size(); size_t bi = b.size(); size_t bj = b[0].size();
    Matrix c(ai, vector<double>(bj, 0.f));
    if (aj != bi) throw std::runtime_error("Mismatched Matrices\n");
    for (std::size_t i = 0; i < ai; ++i) {
        for (std::size_t k = 0; k < aj; ++k) {
            const double aik = a[i][k];
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
    Matrix c(ai, vector<double>(aj, 0.f));
    if (ai != bi || aj != bj) throw std::runtime_error("Mismatched Matrices\n");
    for (size_t i = 0; i < ai; ++i) {for (size_t j = 0; j < aj; ++j) {c[i][j] = a[i][j] + b[i][j];}}
    return c;
}

Matrix T(const Matrix& a) {
    size_t ai = a.size(); size_t aj = a[0].size();
    Matrix c(aj, vector<double>(ai, 0.f));
    for (size_t i = 0; i < ai; ++i) {for (size_t j = 0; j < aj; ++j) {c[j][i] = a[i][j];}}
    return c;
}

void transformGeometry(Geometry& g, array<double, 3> pos, array<double, 3> rot) {
    double phi = rot[0]; double th = rot[1]; double psi = rot[2]; // x, y, z
    const double cPsi = cos(psi), sPsi = sin(psi);
    const double cTh = cos(th),  sTh = sin(th);
    const double cPhi = cos(phi), sPhi = sin(phi);
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
            double c = s.J;
            Matrix An = mmult(mmult(T(R), A), R);
            Matrix bn = mmult(T(R), msum(mmult(A, t), b));
            double cn = mmult(mmult(T(t), A), t)[0][0] + 2 * mmult(T(b), t)[0][0] + c;
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

using CreateFn = void(*)(double,double,double,Geometry&);
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

vector<Material> readMaterial(std::string filename) {
    std::ifstream file("material/" + filename + ".txt");
    if (!file) { std::cout << "MatFile opening failed\n"; }
    int inelastic, nMaterials;
    if (!(file >> inelastic >> nMaterials)) { std::cout << "Failed reading line for mat properties\n"; }
    const vector<int> MTs = inelastic ? vector<int>{2,4,18,102} : vector<int>{2,18,102};
    auto x = logspace(-11.0, std::log10(20.0), 500);

    vector<Material> mats;
    mats.reserve(nMaterials);
    for (int i = 0; i < nMaterials; ++i) {
        std::string fname; double rho, relativeMoles; 
        if (!(file >> fname >> rho >> relativeMoles)) { std::cout << "Failed reading material properties\n"; }
        Material mat;
        mat.rho = rho;
        mat.proportion = relativeMoles;
        std::string path = "data/" + fname + ".dat";
        std::ifstream materialdata(path);
        if (!readMaterialBlock(materialdata, mat)) {
            std::cerr << "Block read fail " << i << "\n";
        }
        mats.push_back(std::move(mat));
    }
    fillData(mats, x, inelastic);
    return mats; 
}

bool readGeometry(std::istream& in, Geometry& g) {
    std::string command; double a, b, c; array<double, 3> pos; array<double, 3> rot;
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
    g.mats = readMaterial(fname);

    transformGeometry(g, pos, rot);
    return true;
}

bool readUniverseFile(std::istream& in, Universe& u) {
    int nSubUniverse; int nGeometry;
    // Read basic universe information
    array<double, 3> uniBou;
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

    array<double, 3> SUpos; array<double, 3> SUrot; std::string SUfname;
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


// --- STORE DATA ---

bool storeDatakeff(const vector<double>& data) {
    using clock = std::chrono::system_clock;
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(clock::now().time_since_epoch()).count();

    std::ofstream os("output/keff_" + std::to_string(secs) + ".csv");
    if (!os) return false;
    os << std::scientific << std::setprecision(6);
    for (size_t k = 0; k < data.size(); ++k)
        os << k << "," << data[k] << "\n";
    return true;
}

bool writeVTKStructuredPoints(const Mesh3D& M, const vector<double>& field, const std::string& basePath, const std::string& name) {
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
    os << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
    for (int k=0;k<N;++k) os << field[k] << "\n";
    return true;
}
