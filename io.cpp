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
// Materials

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

// Geometry

int addQuadric(Geometry& g, double A,double B,double C,double D,double E,double F, double G,double H,double I,double J) {
    Shape s{};
    s.torus = 0;
    s.A=A; s.B=B; s.C=C; s.D=D; s.E=E; s.F=F; s.G=G; s.H=H; s.I=I; s.J=J;
    g.shapes.push_back(s);
    return (int)g.shapes.size()-1;
}

int addPlane(Geometry& g, double nx, double ny, double nz, double d) {
    return addQuadric(g, 0,0,0, 0,0,0, nx,ny,nz, -d);
}

int addTorus(Geometry& g, double a, double b, double R) {
    Shape s{}; s.torus = 1; s.A=a; s.B=b; s.C=R;
    s.D=s.E=s.F=s.G=s.H=s.J=0.0;
    s.I=1; // Along z axis
    g.shapes.push_back(s);
    return (int)g.shapes.size()-1;
}

int pushLeaf(Geometry& g, int shapeIdx) {
    g.nodes.push_back(Node{Op::L, shapeIdx, -1, -1});
    return (int)g.nodes.size()-1;
}

int makeI(Geometry& g, int L, int R) {
    g.nodes.push_back(Node{Op::I, -1, L, R});
    return (int)g.nodes.size()-1;
}

int intersectAll(Geometry& g, const vector<int>& leaves) {
    if (leaves.empty()) return -1;
    int cur = leaves[0];
    for (size_t i=1;i<leaves.size();++i) cur = makeI(g, cur, leaves[i]);
    return cur;
}

// a = radius
void createBall(double a, double b, double c, Geometry& g) {
    int si = addQuadric(g, 1,1,1, 0,0,0, 0,0,0, -a*a); // x^2 + y^2 + z^2 - a^2 <= 0
    int li = pushLeaf(g, si);
    g.nodeRoot = li;
    g.shape = 1;
}
// a = radius, b = height
void createCylinder(double a, double b, double c, Geometry& g) {
    const double r = a;
    const double h2 = 0.5*b;

    int s_cyl = addQuadric(g, 1,1,0, 0,0,0, 0,0,0, -r*r); // x^2 + y^2 - r^2 <= 0
    int s_top = addPlane(g, 0,0, 1, h2);
    int s_bot = addPlane(g, 0,0,-1, h2);

    int l0 = pushLeaf(g, s_cyl);
    int l1 = pushLeaf(g, s_top);
    int l2 = pushLeaf(g, s_bot);

    g.nodeRoot = intersectAll(g, {l0,l1,l2});
    g.shape = 2;
}

// a = radius, b = height
void createCylinderOpen(double a, double b, double c, Geometry& g) {
    const double r = a;
    int si = addQuadric(g, 1,1,0, 0,0,0, 0,0,0, -r*r);
    int li = pushLeaf(g, si);
    g.nodeRoot = li;
    g.shape = 3;
}
// a < z half space
// This becomes a general plane with transformations
void createPlane(double a, double b, double c, Geometry& g) {
    int s = addPlane(g, 0,0,1, (double)a);
    int l = pushLeaf(g, s);
    g.nodeRoot = l;
    g.shape = 4;
}
// side lengths a, b, c
void createCuboid(double a, double b, double c, Geometry& g) {
    const double hx = 0.5*a, hy = 0.5*b, hz = 0.5*c;

    int sxp = addPlane(g,  1,0,0, hx);
    int sxn = addPlane(g, -1,0,0, hx);
    int syp = addPlane(g, 0, 1,0, hy);
    int syn = addPlane(g, 0,-1,0, hy);
    int szp = addPlane(g, 0,0, 1, hz);
    int szn = addPlane(g, 0,0,-1, hz);

    int l0 = pushLeaf(g, sxp);
    int l1 = pushLeaf(g, sxn);
    int l2 = pushLeaf(g, syp);
    int l3 = pushLeaf(g, syn);
    int l4 = pushLeaf(g, szp);
    int l5 = pushLeaf(g, szn);

    g.nodeRoot = intersectAll(g, {l0,l1,l2,l3,l4,l5});
    g.shape = 5;
}
// side lengths without top or bottom
void createCuboidOpen(double a, double b, double c, Geometry& g) {
    const double hx = 0.5*a, hy = 0.5*b;

    int sxp = addPlane(g,  1,0,0, hx);
    int sxn = addPlane(g, -1,0,0, hx);
    int syp = addPlane(g, 0, 1,0, hy);
    int syn = addPlane(g, 0,-1,0, hy);

    int l0 = pushLeaf(g, sxp);
    int l1 = pushLeaf(g, sxn);
    int l2 = pushLeaf(g, syp);
    int l3 = pushLeaf(g, syn);

    g.nodeRoot = intersectAll(g, {l0,l1,l2,l3});
    g.shape = 6;
}
// a = side length, b = height
void createHexPrism(double a, double b, double c, Geometry& g) {
    const double s  = a;
    const double rA = s * std::sqrt(3.0) * 0.5;
    const double h2 = 0.5*b;

    const double nx1=1.0, ny1=0.0;
    const double nx2=0.5, ny2= std::sqrt(3.0)*0.5;
    const double nx3=0.5, ny3=-std::sqrt(3.0)*0.5;

    int s1p = addPlane(g,  nx1, ny1, 0, rA);
    int s1n = addPlane(g, -nx1,-ny1, 0, rA);
    int s2p = addPlane(g,  nx2, ny2, 0, rA);
    int s2n = addPlane(g, -nx2,-ny2, 0, rA);
    int s3p = addPlane(g,  nx3, ny3, 0, rA);
    int s3n = addPlane(g, -nx3,-ny3, 0, rA);

    int szp = addPlane(g, 0,0, 1, h2);
    int szn = addPlane(g, 0,0,-1, h2);

    int l0 = pushLeaf(g, s1p);
    int l1 = pushLeaf(g, s1n);
    int l2 = pushLeaf(g, s2p);
    int l3 = pushLeaf(g, s2n);
    int l4 = pushLeaf(g, s3p);
    int l5 = pushLeaf(g, s3n);
    int l6 = pushLeaf(g, szp);
    int l7 = pushLeaf(g, szn);

    g.nodeRoot = intersectAll(g, {l0,l1,l2,l3,l4,l5,l6,l7});
    g.shape = 7;
}
// A=a, B=b, C=R
void createTorus(double a, double b, double c, Geometry& g) {
    int st = addTorus(g, a, b, c); 
    int lt = pushLeaf(g, st);
    g.nodeRoot = lt;
    g.shape = 8;
}

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

using Matrix = vector<vector<double>>;
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

using std::cos; using std::sin;
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

// Storing data

std::string timePathCol() {
    using clock = std::chrono::system_clock;
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(
                        clock::now().time_since_epoch()).count();
    return "output/col_" + std::to_string(secs) + ".csv";
}

std::string timePathTime() {
    using clock = std::chrono::system_clock;
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(
                        clock::now().time_since_epoch()).count();
    return "output/time_" + std::to_string(secs) + ".csv";
}

std::string timePathkeff() {
    using clock = std::chrono::system_clock;
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(
                        clock::now().time_since_epoch()).count();
    return "output/keff_" + std::to_string(secs) + ".csv";
}

bool storeDatakeff(const vector<double>& data) {
    std::ofstream os(timePathkeff());
    if (!os) return false;
    os << std::scientific << std::setprecision(6);
    for (size_t k = 0; k < data.size(); ++k)
        os << k << "," << data[k] << "\n";
    return true;
}

bool storeDataCol(const vector<double>& data) {
    std::ofstream os(timePathCol());
    if (!os) return false;
    os << std::scientific << std::setprecision(6);
    for (size_t k = 0; k < data.size(); ++k)
        os << k << "," << data[k] << "\n";
    return true;
}

static void storeTimeHist(const TimeHist& H) {
    std::ofstream os(timePathTime());
    os << std::scientific << std::setprecision(6);
    for (int k=0;k<H.nbins;++k) {
        const double tmid = (k + 0.5) * H.dt;
        os << tmid << "," << H.counts[k] << "\n";
    }
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
