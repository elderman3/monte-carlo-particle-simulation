#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <chrono>
#include <random>
#include <cmath>
#include <map>
#include <array>
#include <algorithm>
#include <deque>
#include <iomanip>
#include <cstring>
#include <cctype>
#include <stack>

using std::array;
using std::vector;

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
    g.nodesRoot = li;
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

    g.nodesRoot = intersectAll(g, {l0,l1,l2});
    g.shape = 2;
}

// a = radius, b = height
void createCylinderOpen(double a, double b, double c, Geometry& g) {
    int si = addQuadric(g, 1,1,0, 0,0,0, 0,0,0, -r*r);
    int li = pushLeaf(g, si);
    g.nodesRoot = li;
    g.shape = 3;
}
// a < z half space
// This becomes a general plane with transformations
void createPlane(double a, double b, double c, Geometry& g) {
    int s = addPlane(g, 0,0,1, (double)a);
    int l = pushLeaf(g, s);
    g.nodesRoot = l;
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

    g.nodesRoot = intersectAll(g, {l0,l1,l2,l3,l4,l5});
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

    g.nodesRoot = intersectAll(g, {l0,l1,l2,l3});
    g.shape = 6;
}
// a = side length, b = height
void createHexPrism(double a, double b, double c, Geometry& g) {
    const double s  = a;
    const double rA = s * std::sqrt(3.0) * 0.5; // distance to each side
    const double h2 = 0.5*b;

    // Normals for flat-top hex:
    const double nx1=1.0, ny1=0.0;
    const double nx2=0.5, ny2= std::sqrt(3.0)*0.5;
    const double nx3=0.5, ny3=-std::sqrt(3.0)*0.5;

    // Six XY half-spaces: |n·(x,y)| <= rA
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

    g.nodesRoot = intersectAll(g, {l0,l1,l2,l3,l4,l5,l6,l7});
    g.shape = 7;
}
// A=a, B=b, C=R
void createTorus(double a, double b, double c, Geometry& g) {
    int st = addTorus(g, a, b, c); 
    int lt = pushLeaf(g, st);
    g.nodesRoot = lt;
    g.shape = 8;
}

using CreateFn = void(*)(double,double,double,Geometry&);

static const vector<CreateFn> kCreateById = {
  nullptr,            // 0 = GENERAL (handled specially)
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
    auto push = [&](Node n){ nodes.push_back(n); st.push_back((int)nodes.size()-1); };

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
    if (!(in >> s.A >> s.B >> s.C >> s.D >> s.E >> s.F >> s.G >> s.H >> s.I >> s.J >> s.torus)) return false;
    g.shapes.push_back(std::move(s));
}

using Matrix = vector<vector<double>>;
Matrix mmult(Matrix& a, Matrix& b) {
    if (a.empty() || b.empty()) throw std::runtime_error("empty matrix");
    size_t ai = a.size(); size_t aj = a[0].size(); size_t bi = b.size(); size_t bj = b[0].size();
    Matrix c(ai, vector<double>(bj, 0.f));
    if (aj != bi) throw std::runtime_error("Mismatched Matrices");
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

Matrix msum(Matrix& a, Matrix& b) {
    if (a.empty() || b.empty()) throw std::runtime_error("empty matrix");
    size_t ai = a.size(); size_t aj = a[0].size(); size_t bi = b.size(); size_t bj = b[0].size();
    Matrix c(ai, vector<double>(aj, 0.f));
    if (ai != bi || aj != bj) throw std::runtime_error("Mismatched Matrices");
    for (size_t i = 0; i < ai; ++i) {for (size_t j = 0; j < aj; ++j) {c[i][j] = a[i][j] + b[i][j];}}
    return c;
}

Matrix T(Matrix& a) {
    size_t ai = a.size(); size_t aj = a[0].size();
    Matrix c(aj, vector<double>(ai, 0.f));
    for (size_t i = 0; i < ai; ++i) {for (size_t j = 0; j < aj; ++j) {c[j][i] = a[i][j];}}
    return c;
}
using std::cos; using std::sin;
void transformGeometry(Geometry& g, array<double, 3> pos, array<double, 3> rot) {
    double phi = rot[0]; double th = rot[1]; double psi = rot[2]; // x, y, z
    const double cψ = cos(psi), sψ = sin(psi);
    const double cθ = cos(th),  sθ = sin(th);
    const double cφ = cos(phi), sφ = sin(phi);

    Matrix R = {
            { cψ*cθ, cψ*sθ*sφ - sψ*cφ, cψ*sθ*cφ + sψ*sφ },
            { sψ*cθ, sψ*sθ*sφ + cψ*cφ, sψ*sθ*cφ - cψ*sφ },
            { -sθ,   cθ*sφ,            cθ*cφ            }
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

bool readGeometry(std::istream& in, Geometry& g) {
    std::string command; double a, b, c; array<double, 3> pos; array<double, 3> rot;
    if (!(in >> command >> a >> b >> c >> pos[0] >> pos[1] >> pos[2] >> rot[0] >> rot[1] >> rot[2])) return false;
    if (std::strcmp(command, "general") == 0) {
        int nShapes; std::string NodeDef;
        if (!(in >> nShapes >> NodeDef)) return false;
        for (int i = 0; i < nShapes; i++) {
            if (!(readShape(in, g))) return false; 
        }
        g.shape = 0;
        vector<Node> nodes;
        g.nodesRoot = readNodeDef(NodeDef, nodes);
        g.nodes = nodes; 
    } else if (std::strcmp(command, "ball") == 0) {
        createBall(a, b, c, g);
    } else if (std::strcmp(command, "cylinder") == 0) {
        createCylinder(a, b, c, g);
    } else if (std::strcmp(command, "cylinderOpen") == 0) {
        createCylinderOpen(a, b, c, g);
    } else if (std::strcmp(command, "plane") == 0) {
        createPlane(a, b, c, g);
    } else if (std::strcmp(command, "cuboid") == 0) {
        createCuboid(a, b, c, g);
    } else if (std::strcmp(command, "cuboidOpen") == 0) {
        createCuboidOpen(a, b, c, g);
    } else if (std::strcmp(command, "hexPrism") == 0) {
        createHexPrism(a, b, c, g);
    } else if (std::strcmp(command, "torus") == 0) {
        createTorus(a, b, c, g);
    } else {
        std:cerr << "No shape found";
    }
    transformGeometry(g, pos, rot);
}

bool readUniverseFile(std::istream& in, Universe& u) {
    int nSubUniverse; int nGeometry;
    // Read basic universe information
    array<double, 3> uniBou;
    if (!(in >> u.name >> u.pos[0] >> u.pos[1] >> u.pos[2] >> u.rot[0] >> u.rot[1] >> u.rot[2] 
        >> u.latticeType >> nSubUniverse >> nGeometry >> uniBou[0] >> uniBou[1] >> uniBou[2] >> u.universeShape)) return false;
    if (!u.universeShape) return false;
    Geometry bounds;
    u.boundDim = uniBou;
    kCreateById[u.universeShape](uniBou[0], uniBou[1], uniBou[2], bounds);
    u.boundingGeometry = bounds;
    if (u.latticeType) {
        if (!(in >> u.lattice[0] >> u.lattice[1])) return false;
    }
    // Read subuniverses
    array<double, 3> SUpos; array<double, 3> SUrot; std::string SUfname;
    for (int i = 0; i < nSubUniverse; i++) {
        Universe su;
        if (!(in >> SUfname >> SUpos[0] >> SUpos[1] >> SUpos[2] >> SUrot[0] >> SUrot[1] >> SUrot[2])) return false;
        std::string path = "data/" + SUfname + ".dat";
        std::ifstream SUfile(path);
        if (!(readUniverseFile(SUfile, su))) return false;
        su.pos = SUpos; su.rot = SUrot
        u.subUniverse.push_back(std::move(su));
    }
    for (int i = 0; i < nGeometry; i++) {
        Geometry g;
        if (!(readGeometry(in, g))) return false;
        u.geometries.push_back(std::move(g));
    }    
}

static double dot3(const array<double,3>& a, const array<double,3>& b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static inline std::array<double,3> add3(const std::array<double,3>& a,const std::array<double,3>& b){return{a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
static inline std::array<double,3> sub3(const std::array<double,3>& a,const std::array<double,3>& b){return{a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
static array<double,3> cross3(const array<double,3>& a, const array<double,3>& b){
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
static array<double,3> normed(const array<double,3>& v){
    double L = std::sqrt(dot3(v,v));
    if (L==0.0) return {0,0,0};
    return {v[0]/L, v[1]/L, v[2]/L};
}

static array<double,3> fmadd3(const array<double,3>& p, const array<double,3>& d, double s){
    return {p[0]+s*d[0], p[1]+s*d[1], p[2]+s*d[2]};
}

static int solveQuadratic(double A,double B,double C,double r[2]){
    const double EPS=1e-14;
    if (std::abs(A) < EPS){
        if (std::abs(B) < EPS) return 0;
        r[0] = -C/B; return 1;
    }
    double disc = B*B - 4*A*C;
    if (disc < 0) return 0;
    double sdisc = std::sqrt(std::max(0.0,disc));
    // Kahan: avoid catastrophic cancellation
    double q = -0.5*(B + (B>=0 ? sdisc : -sdisc));
    double t0 = q/A;
    double t1 = C/q;
    if (t0>t1) std::swap(t0,t1);
    r[0]=t0; r[1]=t1;
    return (disc==0.0)?1:2;
}

static int solveCubicMonic(double a,double b,double c,double r[3]){
    const double TWOPI=6.28318530717958647692;
    const double EPS=1e-14;
    double a3 = a/3.0;
    double p = b - a*a/3.0;
    double q = (2.0*a*a*a)/27.0 - (a*b)/3.0 + c;
    double D = 0.25*q*q + (p*p*p)/27.0;
    if (D > EPS){
        double s = std::sqrt(D);
        double u = std::cbrt(-0.5*q + s);
        double v = std::cbrt(-0.5*q - s);
        r[0] = u + v - a3;
        return 1;
    } else if (std::abs(D) <= EPS){
        double u = std::cbrt(-0.5*q);
        r[0] = 2*u - a3;
        r[1] = -u - a3;
        return 2;
    } else {
        double rho = std::sqrt(-p/3.0);
        double theta = std::acos( (-0.5*q) / (rho*rho*rho) );
        r[0] = 2*rho*std::cos(theta/3.0) - a3;
        r[1] = 2*rho*std::cos((theta+TWOPI)/3.0) - a3;
        r[2] = 2*rho*std::cos((theta+2*TWOPI)/3.0) - a3;
        // sort ascending
        if (r[0]>r[1]) std::swap(r[0],r[1]);
        if (r[1]>r[2]) std::swap(r[1],r[2]);
        if (r[0]>r[1]) std::swap(r[0],r[1]);
        return 3;
    }
}

static int solveQuartic(double a,double b,double c,double d,double e,double r[4]){
    const double EPS=1e-14;
    if (std::abs(a) < EPS){
        if (std::abs(b) < EPS){
            return solveQuadratic(c,d,e,r);
        } else {
            double roots[3];
            int cnt = solveCubicMonic(c/b, d/b, e/b, roots);
            for (int i=0;i<cnt;++i) r[i]=roots[i];
            return cnt;
        }
    }
    double A=b/a, B=c/a, C=d/a, D=e/a;

    double A2=A*A;
    double p = -3.0*A2/8.0 + B;
    double q =  A*A2/8.0 - 0.5*A*B + C;
    double r0= -3.0*A2*A2/256.0 + A2*B/16.0 - 0.25*A*C + D;

    double cr[3]; int cc = solveCubicMonic(2*p, (p*p - 4*r0), -q*q, cr);
    double z = -1e300;
    for (int i=0;i<cc;++i) if (cr[i] > z) z = cr[i];
    if (z < 0) z = 0;
    double alpha = std::sqrt(z);

    double m = p + z;
    double beta, gamma;
    if (alpha > 1e-15){
        beta  = 0.5*(m - q/alpha);
        gamma = 0.5*(m + q/alpha);
    } else {
        beta = gamma = 0.5*m;
    }

    double yroots[4]; int k=0;
    {
        double qr[2]; int qc = solveQuadratic(1.0,  alpha, beta, qr);
        for (int i=0;i<qc;++i) yroots[k++] = qr[i];
    }
    {
        double qr[2]; int qc = solveQuadratic(1.0, -alpha, gamma, qr);
        for (int i=0;i<qc;++i) yroots[k++] = qr[i];
    }

    for (int i=0;i<k;++i) r[i] = yroots[i] - A/4.0;

    std::sort(r, r+k);
    return k;
}

vector<double> surfaceDist(const Neutron& n, const Shape& s)
{
    const double EPS = 1e-10;
    const array<double,3>& P = n.pos;
    const array<double,3>& Dv = n.dir;
    double Dlen = std::sqrt(dot3(Dv,Dv));
    if (Dlen == 0.0) return {};

    vector<double> out;

    if (!s.torus){
        // General quadric: M t^2 + L t + K = 0
        const double x=P[0], y=P[1], z=P[2];
        const double u=Dv[0], v=Dv[1], w=Dv[2];

        const double K = s.A*x*x + s.B*y*y + s.C*z*z
                       + s.D*x*y + s.E*y*z + s.F*x*z
                       + s.G*x + s.H*y + s.I*z + s.J;

        const double L = 2.0*(s.A*u*x + s.B*v*y + s.C*w*z)
                       + s.D*(v*x + u*y) + s.E*(w*y + v*z) + s.F*(w*x + u*z)
                       + s.G*u + s.H*v + s.I*w;

        const double M = s.A*u*u + s.B*v*v + s.C*w*w
                       + s.D*u*v + s.E*v*w + s.F*u*w;

        double roots[2]; int rc = solveQuadratic(M,L,K,roots);
        for (int i=0;i<rc;++i){
            double t = roots[i];
            if (t > EPS){
                out.push_back(t * Dlen);
            }
        }
        std::sort(out.begin(), out.end());
        return out;
    }

    const double a = s.A;
    const double b = s.B;
    const double R = s.C;

    array<double,3> Cn = {s.D, s.E, s.F};
    array<double,3> Ax = {s.G, s.H, s.I};
    array<double,3> ez = normed(Ax);
    if (ez == array<double,3>{0,0,0}) ez = {0,0,1};

    array<double,3> tmp = (std::fabs(ez[2]) < 0.9) ? array<double,3>{0,0,1}
                                                   : array<double,3>{1,0,0};
    array<double,3> ex = normed(cross3(tmp, ez));
    array<double,3> ey = cross3(ez, ex);

    array<double,3> P0 = sub3(P, Cn);
    auto toLocal = [&](const array<double,3>& v)->array<double,3>{
        return { dot3(v,ex), dot3(v,ey), dot3(v,ez) };
    };
    array<double,3> Pl = toLocal(P0);
    array<double,3> Dl = toLocal(Dv);

    const double x0=Pl[0], y0=Pl[1], z0=Pl[2];
    const double ux=Dl[0], uy=Dl[1], uz=Dl[2];

    const double b2 = b*b, a2 = a*a, R2 = R*R;
    const double q2 = ux*ux + uy*uy;
    const double q1 = 2.0*(x0*ux + y0*uy);
    const double q0 = x0*x0 + y0*y0;

    const double s2 = b2*q2 + a2*uz*uz;
    const double s1 = b2*q1 + 2.0*a2*z0*uz;
    const double s0 = b2*q0 + a2*z0*z0 + b2*(R2 - a2);

    const double C4 = s2*s2;
    const double C3 = 2.0*s2*s1;
    const double C2 = 2.0*s2*s0 + s1*s1 - 4.0*R2*b2*b2*q2;
    const double C1 = 2.0*s1*s0 - 4.0*R2*b2*b2*q1;
    const double C0 = s0*s0 - 4.0*R2*b2*b2*q0;

    double roots[4]; int rc = solveQuartic(C4,C3,C2,C1,C0,roots);
    for (int i=0;i<rc;++i){
        double t = roots[i];
        if (t > EPS){
            out.push_back(t * Dlen);
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

static bool insideLeaf(const Shape& s, const array<double,3>& P){
    if (!s.torus){
        const double x=P[0], y=P[1], z=P[2];
        const double f = s.A*x*x + s.B*y*y + s.C*z*z
                       + s.D*x*y + s.E*y*z + s.F*x*z
                       + s.G*x + s.H*y + s.I*z + s.J;
        return f <= 0.0;
    }
    // Torus: localize by center D,E,F and axis G,H,I
    const double a=s.A, b=s.B, R=s.C;
    array<double,3> Cn{ s.D, s.E, s.F };
    array<double,3> Ax{ s.G, s.H, s.I };
    array<double,3> ez = normed(Ax);
    if (ez==array<double,3>{0,0,0}) ez = {0,0,1};
    array<double,3> tmp = (std::fabs(ez[2]) < 0.9) ? array<double,3>{0,0,1}
                                                        : array<double,3>{1,0,0};
    array<double,3> ex = normed(cross3(tmp, ez));
    array<double,3> ey = cross3(ez, ex);

    array<double,3> Pl = sub3(P, Cn);
    const double x = dot3(Pl, ex), y = dot3(Pl, ey), z = dot3(Pl, ez);

    const double b2=b*b, a2=a*a, R2=R*R;
    const double r2 = x*x + y*y;
    const double S  = b2*r2 + a2*z*z + b2*(R2 - a2);
    const double F  = S*S - 4.0*R2*b2*b2*r2;
    return F <= 0.0;
}

bool pointInGeom(const array<double, 3>& pos, const Geometry& g){
    if (g.nodes.empty()) return false;
    std::function<bool(int)> eval = [&](int idx)->bool{
        const Node& n = g.nodes[idx];
        switch (n.op){
            case Op::L: {
                const Shape& s = g.shapes[n.shape];
                return insideLeaf(s, pos);
            }
            case Op::N: return !eval(n.left);
            case Op::U: return  eval(n.left) || eval(n.right);
            case Op::I: return  eval(n.left) && eval(n.right);
            default:    return false;
        }
    };
    return eval(g.nodeRoot);
}

double geometryCollision(const Neutron& n, const Geometry& g){
    const double EPS = 1e-9;
    const double dlen = std::sqrt(dot3(n.dir,n.dir));
    if (dlen==0.0) return -1.0;
    array<double,3> dHat{ n.dir[0]/dlen, n.dir[1]/dlen, n.dir[2]/dlen };

    vector<double> cands;
    for (const Shape& s : g.shapes){
        vector<double> ds = surfaceDist(n, s);
        cands.insert(cands.end(), ds.begin(), ds.end());
    }
    if (cands.empty()) return -1.0;
    std::sort(cands.begin(), cands.end());
    vector<double> uniq;
    for (double v : cands){
        if (v <= EPS) continue;
        if (uniq.empty() || std::fabs(v-uniq.back()) > 1e-8) uniq.push_back(v);
    }

    bool inside0 = pointInGeom(n.pos, g);
    for (double sDist : uniq){
        double s0 = std::max(0.0, sDist - 1e-8);
        double s1 = sDist + 1e-8;
        bool inA = pointInGeom(madd(n.pos, dHat, s0), g);
        bool inB = pointInGeom(madd(n.pos, dHat, s1), g);
        if (inA != inB){
            if (sDist <= EPS && inB==inside0) continue;
            return sDist;
        }
    }
    return -1.0;
}

double nearestCollision(const Neutron& n, const Universe& u){
    const double INF = 1e300;
    double best = INF;

    for (const Geometry& g : u.geometries){
        double d = geometryCollision(n, g);
        if (d > 0.0 && d < best) best = d;
    }

    for (const Universe& su : u.subUniverse){
        Neutron nLocal = n;
        nLocal.pos[0] -= su.pos[0];
        nLocal.pos[1] -= su.pos[1];
        nLocal.pos[2] -= su.pos[2];
        double dBound = geometryCollision(nLocal, su.boundingGeometry);
        if (dBound > 0.0 && dBound <= best){
            double dIn = nearestCollision(nLocal, su);
            if (dIn > 0.0 && dIn < best) best = dIn;
        }
    }

    return (best<INF)? best : -1.0;
}

Result timing(float (*f)(int, float), int iter, float length) {
    auto start = std::chrono::steady_clock::now();
    float pi = f(iter, length);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    Result result{elapsed.count(), pi};
    return result;
}

array<double, 3> boundingBox(Universe& u) {
    if (u.latticeType == 1) {
        return { scale(u.subUniverse[0].boundDim[0]),
                 scale(u.subUniverse[0].boundDim[1]),
                 scale(u.subUniverse[0].boundDim[2]) };
    } else if (u.latticeType == 2) {
        double sX = 0.5 * u.subUniverse[0].boundDim[0];
        double sY = (u.subUniverse[0].boundDim[1] > 0.0) ? (u.subUniverse[0].boundDim[1] / std::sqrt(3.0)) : sX;
        const double s = std::max(1e-9, (sX > 0 && sY > 0) ? std::min(sX, sY) : std::max(sX, sY));
        const double Wc = 2.0 * s;
        const double Hc = std::sqrt(3.0) * s;
        const double Zc = u.subUniverse[0].boundDim[2];
        return { scale(Wc), scale(Hc), scale(Zc) };
    }
    // Open shapes get 1 as height parameter
    if (u.universeShape == 1) {
        const double a = u.boundDim[0];
        return { scale(2*a), scale(2*a), scale(2*a) };
    } else if (u.universeShape == 2 || u.universeShape == 3) {
        const double a = u.boundDim[0], b = u.boundDim[1] ? u.boundDim[1] : 1;
        return { scale(2*a), scale(2*a), scale(b) };
    } else if (u.universeShape == 4) {
        return { scale(u.boundDim[0]), scale(u.boundDim[1]), scale(u.boundDim[2]) };
    } else if (u.universeShape == 5 || u.universeShape == 6) {
        return { scale(u.boundDim[0]), scale(u.boundDim[1]), scale(u.boundDim[2] ? u.boundDim[2] : 1) };
    } else if (u.universeShape == 7) {
        const double a = u.boundDim[0], b = u.boundDim[1];
        return { scale(2*a), scale(std::sqrt(3.0)*a), scale(b) };
    } else if (u.universeShape == 8) {
        const double a = u.boundDim[0], b = u.boundDim[1], R = u.boundDim[2];
        return { scale(2*(R + a)), scale(2*(R + a)), scale(2*b) };
    }
    return { 0.f, 0.f, 0.f};
}

static double timeNow(){
    using clk=std::chrono::high_resolution_clock;
    static const auto t0 = clk::now();
    return std::chrono::duration<double>(clk::now()-t0).count();
}

static double relFOM(double mean, double stderr, double sec){
    if (mean<=0.0 || stderr<=0.0 || sec<=0.0) return 0.0;
    const double R = stderr/mean;
    return 1.0/(R*R*sec);
}

void printStats(std::string str, double mean, double stderr, double sec, double FOM) {
    const double rel = (mean != 0.0) ? (stderr / mean) : 0.0;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << label << "\n"
              << "  volume_mean   = " << mean   << "\n"
              << "  volume_stderr = " << stderr << "\n"
              << "  rel_error     = " << rel    << "\n"
              << "  time_sec      = " << sec    << "\n"
              << "  FOM           = " << FOM    << "\n";
}

static array<double,3> squareCellCenter(const Universe& u, int i, int j){
    const int nx=u.lattice[0], ny=u.lattice[1];
    const double pitchX=u.boundDim[0], pitchY=u.boundDim[1];
    const double cx=0.5*(nx-1)*pitchX, cy=0.5*(ny-1)*pitchY;
    return { i*pitchX - cx, j*pitchY - cy, 0.0 };
}
static array<double,3> hexCellCenter(const Universe& u, int q, int r){
    const double sX=0.5*u.boundDim[0];
    const double sY=(u.boundDim[1]>0.0)?(u.boundDim[1]/std::sqrt(3.0)):sX;
    const double s = std::max(1e-9,(sX>0&&sY>0)?std::min(sX,sY):std::max(sX,sY));
    const double pitchX=1.5*s, pitchY=std::sqrt(3.0)*s;
    const int cols=u.lattice[0], rows=u.lattice[1];
    const double cx=0.5*(cols-1)*pitchX;
    const double cy=0.5*((rows-1)*pitchY + (cols>1?0.5*pitchY:0.0));
    const double x=q*pitchX - cx;
    const double y=r*pitchY + ((q&1)?0.5*pitchY:0.0) - cy;
    return {x,y,0.0};
}

static void tallyPointRecursive(const Universe& u, const array<double,3>& pLocal, std::unordered_map<const Geometry*, long long>& hits) {
    for (const Geometry& g : u.geometries){
        if (pointInGeom(pLocal, g)) ++hits[&g];
    }
    for (const Universe& su : u.subUniverse){
        std::array<double,3> pSub = sub3(pLocal, su.pos);
        tallyPointRecursive(su, pSub, hits);
    }
}

static void collectGeometries(const Universe& u, const std::string& prefix, vector<std::pair<const Geometry*,std::string>>& out) {
    for (size_t i=0;i<u.geometries.size();++i)
        out.emplace_back(&u.geometries[i], prefix + "geom[" + std::to_string(i) + "]");
    for (size_t j=0;j<u.subUniverse.size();++j)
        
    u.subUniverse[j], prefix + "su[" + std::to_string(j) + "]/", out);
}

void volumePointMethod(Universe& u, int iter) {
    std::mt19937_64 rng(0xC0FFEE);
    const bool isSquare = (u.latticeType==1);
    const bool isHex    = (u.latticeType==2);

    vector<std::pair<const Geometry*,std::string>> roster;
    collectGeometries(u, "", roster);
    const auto boxAll = boundingBox(u);
    const double Vall = boxAll[0]*boxAll[1]*boxAll[2];
    int nCells = 1;
    if (isSquare || isHex) nCells = std::max(1, u.lattice[0]*u.lattice[1]);
    const auto cdim = cellBoxDims(u);
    const double Vcell = cdim[0]*cdim[1]*cdim[2];

    const double hx=0.5*boxAll[0], hy=0.5*boxAll[1], hz=0.5*boxAll[2];
    std::uniform_real_distribution<double> Ux(-hx,hx), Uy(-hy,hy), Uz(-hz,hz);
    const double hx_c=0.5*cdim[0], hy_c=0.5*cdim[1], hz_c=0.5*cdim[2];
    std::uniform_real_distribution<double> Ucx(-hx_c,hx_c), Ucy(-hy_c,hy_c), Ucz(-hz_c,hz_c);

    std::unordered_map<const Geometry*, long long> hits;
    hits.reserve(roster.size()*2);

    const double t0 = nowSec();

    if (!isSquare && !isHex){
        for (int k=0;k<iter;++k){
            array<double,3> p{Ux(rng),Uy(rng),Uz(rng)};
            tallyPointRecursive(u, p, hits);
        }
    } else {
        const int itPerCell = std::max(1, iter / nCells);
        if (isSquare){
            for (int j=0;j<u.lattice[1];++j){
                for (int i=0;i<u.lattice[0];++i){
                    const auto C = squareCellCenter(u,i,j);
                    for (int k=0;k<itPerCell;++k){
                        array<double,3> pLocal{Ucx(rng),Ucy(rng),Ucz(rng)};
                        array<double,3> pUniverse = add3(pLocal, C);
                        for (const Geometry& g : u.geometries){
                            if (pointInGeom(pLocal, g)) ++hits[&g];
                        }
                        for (const Universe& su : u.subUniverse){
                            array<double,3> pSub = sub3(pUniverse, su.pos);
                            tallyPointRecursive(su, pSub, hits);
                        }
                    }
                }
            }
        } else {
            for (int r=0;r<u.lattice[1];++r){
                for (int q=0;q<u.lattice[0];++q){
                    const auto C = hexCellCenter(u,q,r);
                    for (int k=0;k<itPerCell;++k){
                        array<double,3> pLocal{Ucx(rng),Ucy(rng),Ucz(rng)};
                        array<double,3> pUniverse = add3(pLocal, C);
                        for (const Geometry& g : u.geometries){
                            if (pointInGeom(pLocal, g)) ++hits[&g];
                        }
                        for (const Universe& su : u.subUniverse){
                            array<double,3> pSub = sub3(pUniverse, su.pos);
                            tallyPointRecursive(su, pSub, hits);
                        }
                    }
                }
            }
        }
    }

    const double t1 = nowSec();
    const double Veffective = (isSquare||isHex) ? ( (double)nCells * Vcell ) : Vall;
    const double Ntot = std::max(1.0, (double)iter);
    const double sec = t1 - t0;

    for (auto& pr : roster){
        const Geometry* gp = pr.first;
        const std::string& label = pr.second;
        const double cnt = (double)hits[gp];
        const double p = cnt / Ntot;
        const double mean   = Veffective * p;
        const double stderr = Veffective * std::sqrt(std::max(0.0, p*(1.0-p))) / std::sqrt(Ntot);
        const double FOM    = relFOM(mean, stderr, sec);
        printStats("PointVol " + label, mean, stderr, sec, FOM);
    }
}

static double marchLengthOneGeom(const Geometry& g, array<double,3> p0, const array<double,3>& dir, double Ltot) {
    Neutron ray; ray.pos=p0; ray.dir=dir;
    const double EPS=1e-9;
    double acc=0.0, traveled=0.0;
    bool inside = pointInGeom(ray.pos, g);
    while (traveled < Ltot - EPS){
        double d = geometryCollision(ray, g);
        if (d < EPS || traveled + d > Ltot){
            double seg = (Ltot - traveled);
            if (inside && seg>0.0) acc+=seg;
            break;
        }
        if (inside) acc += d;
        ray.pos[0]+=dir[0]*d; ray.pos[1]+=dir[1]*d; ray.pos[2]+=dir[2]*d;
        ray.pos[0]+=dir[0]*1e-8; ray.pos[1]+=dir[1]*1e-8; ray.pos[2]+=dir[2]*1e-8;
        traveled += d + 1e-8;
        inside = !inside;
    }
    return acc;
}

static void tallyLineRecursive(const Universe& u, const array<double,3>& p0_local, const array<double,3>& dir_local, double Ltot, std::unordered_map<const Geometry*, double>& lenSum) {
    for (const Geometry& g : u.geometries)
        lenSum[&g] += marchLengthOneGeom(g, p0_local, dir_local, Ltot);

    for (const Universe& su : u.subUniverse){
        std::array<double,3> p_child = sub3(p0_local, su.pos);
        tallyLineRecursive(su, p_child, dir_local, Ltot, lenSum);
    }
}

void volumeLineMethod(Universe& u, int iter){
    std::mt19937_64 rng(0xBADC0DE);
    const bool isSquare = (u.latticeType==1);
    const bool isHex = (u.latticeType==2);

    std::vector<std::pair<const Geometry*,std::string>> roster;
    collectGeometries(u, "", roster);
    const auto boxAll = boundingBox(u);
    const double Vall = boxAll[0]*boxAll[1]*boxAll[2];
    int nCells = 1;
    if (isSquare || isHex) nCells = std::max(1, u.lattice[0]*u.lattice[1]);

    const auto cdim = cellBoxDims(u);
    const double Vcell = cdim[0]*cdim[1]*cdim[2];

    const double hx=0.5*boxAll[0], hy=0.5*boxAll[1], hz=0.5*boxAll[2];
    std::uniform_real_distribution<double> Ux(-hx,hx), Uy(-hy,hy), Uz(-hz,hz);

    const double hx_c=0.5*cdim[0], hy_c=0.5*cdim[1], hz_c=0.5*cdim[2];
    std::uniform_real_distribution<double> Ucx(-hx_c,hx_c), Ucy(-hy_c,hy_c), Ucz(-hz_c,hz_c);

    const double Lx_c = cdim[0], Ly_c = cdim[1], Lz_c = cdim[2];

    std::unordered_map<const Geometry*, double> sumLen, sumLen2;
    sumLen.reserve(roster.size()*2); sumLen2.reserve(roster.size()*2);

    const double t0 = nowSec();

    if (!isSquare && !isHex){
        for (int k=0;k<iter;++k){
            int axis = k % 3;
            std::array<double,3> p0, dir; double Ltot=0.0;
            if (axis==0){ p0={-hx, Uy(rng), Uz(rng)}; dir={1,0,0}; Ltot=boxAll[0]; }
            else if (axis==1){ p0={Ux(rng), -hy, Uz(rng)}; dir={0,1,0}; Ltot=boxAll[1]; }
            else { p0={Ux(rng), Uy(rng), -hz}; dir={0,0,1}; Ltot=boxAll[2]; }

            std::unordered_map<const Geometry*, double> inc;
            tallyLineRecursive(u, p0, dir, Ltot, inc);
            for (auto& pr : inc){
                double f = (Ltot>0.0)? (pr.second / Ltot) : 0.0;
                sumLen[pr.first]  += f;
                sumLen2[pr.first] += f*f;
            }
        }
    } else {
        const int itPerCell = std::max(1, iter / nCells);

        auto do_cell = [&](const std::array<double,3>& C){
            for (int k=0;k<itPerCell;++k){
                int axis = k % 3;
                std::array<double,3> pLocal, dir; double Ltot=0.0;
                if (axis==0){ pLocal={-hx_c, Ucy(rng), Ucz(rng)}; dir={1,0,0}; Ltot=Lx_c; }
                else if (axis==1){ pLocal={Ucx(rng), -hy_c, Ucz(rng)}; dir={0,1,0}; Ltot=Ly_c; }
                else { pLocal={Ucx(rng), Ucy(rng), -hz_c}; dir={0,0,1}; Ltot=Lz_c; }

                std::unordered_map<const Geometry*, double> inc;

                for (const Geometry& g : u.geometries){
                    double len = marchLengthOneGeom(g, pLocal, dir, Ltot);
                    if (len>0) inc[&g] += len;
                }
                for (const Universe& su : u.subUniverse){
                    std::array<double,3> pParent = add3(pLocal, C);
                    std::array<double,3> pChild  = sub3(pParent, su.pos);
                    tallyLineRecursive(su, pChild, dir, Ltot, inc);
                }

                for (auto& pr : inc){
                    double f = (Ltot>0.0)? (pr.second / Ltot) : 0.0;
                    sumLen[pr.first]  += f;
                    sumLen2[pr.first] += f*f;
                }
            }
        };

        if (isSquare){
            for (int j=0;j<u.lattice[1];++j)
                for (int i=0;i<u.lattice[0];++i)
                    do_cell(squareCellCenter(u,i,j));
        } else {
            for (int r=0;r<u.lattice[1];++r)
                for (int q=0;q<u.lattice[0];++q)
                    do_cell(hexCellCenter(u,q,r));
        }
    }

    const double t1 = nowSec();
    const double N = std::max(1, iter);

    const double Veffective = (isSquare||isHex) ? ( (double)nCells * Vcell ) : Vall;

    for (auto& rp : roster){
        const Geometry* gp = rp.first;
        const std::string& label = rp.second;

        double fbar = sumLen[gp] / N;
        double m2 = sumLen2[gp] / N;
        double varf = std::max(0.0, m2 - fbar*fbar);
        double stderr = Veffective * std::sqrt(varf) / std::sqrt((double)N);
        double mean = Veffective * fbar;
        double sec = t1 - t0;
        double FOM = relFOM(mean, stderr, sec);

        printStats("LineVol " + label, mean, stderr, sec, FOM);
    }
}

static int hitGeomIndexAt(const Universe& u,
                          const array<double,3>& p_world,
                          const vector<std::pair<const Geometry*,std::string>>& roster,
                          int& cursor,
                          const array<double,3>& p_local) {
    for (size_t i=0;i<u.geometries.size();++i){
        const Geometry& g = u.geometries[i];
        int idx = cursor++;
        if (pointInGeom(p_local, g)) return idx;
    }
    for (const Universe& su : u.subUniverse){
        array<double,3> p_child = sub3(p_local, su.pos); // rot ignored by design
        int idx = hitGeomIndexAt(su, p_world, roster, cursor, p_child);
        if (idx >= 0) return idx;
    }
    return -1;
}

static int hitGeomIndex(const Universe& u,
                        const array<double,3>& p_world,
                        const vector<std::pair<const Geometry*,std::string>>& roster) {
    int cursor = 0;
    return hitGeomIndexAt(u, p_world, roster, cursor, p_world);
}

enum class SliceAxis { X, Y, Z };

static std::array<double,3> pixelToWorld(SliceAxis ax, double coord,
                                         double xmin, double xmax,
                                         double ymin, double ymax,
                                         int ix, int iy, int W, int H) {
    double u = (ix + 0.5) / W;
    double v = (iy + 0.5) / H;
    double X = xmin + u*(xmax - xmin);
    double Y = ymin + v*(ymax - ymin);
    if (ax == SliceAxis::Z) return {X, Y, coord};
    if (ax == SliceAxis::Y) return {X, coord, Y};
    return {coord, X, Y};
}

void renderSliceASCII(const Universe& u,
                      SliceAxis ax, double coord,
                      int W=120, int H=60,
                      double xspan=-1, double yspan=-1) {
    auto box = boundingBox(const_cast<Universe&>(u));
    double sx = (ax==SliceAxis::X)? box[1] : box[0];
    double sy = (ax==SliceAxis::Z)? box[1] : box[2];
    if (ax==SliceAxis::Y){ sx = box[0]; sy = box[2]; }

    if (xspan > 0) sx = xspan;
    if (yspan > 0) sy = yspan;

    double xmin = -0.5*sx, xmax = 0.5*sx;
    double ymin = -0.5*sy, ymax = 0.5*sy;

    std::vector<std::pair<const Geometry*,std::string>> roster;
    collectGeometries(u, "", roster);

    const char lut[] =
        " .:-=+*#%@ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    const int lutN = (int)sizeof(lut)-1;

    for (int iy=H-1; iy>=0; --iy){
        for (int ix=0; ix<W; ++ix){
            auto p = pixelToWorld(ax, coord, xmin,xmax,ymin,ymax, ix,iy, W,H);
            int gi = hitGeomIndex(u, p, roster);
            char ch = (gi<0) ? ' ' : lut[(gi % lutN)];
            std::putchar(ch);
        }
        std::putchar('\n');
    }
}

int main() {
    std::ifstream file("geometry.txt");
    if (!file) { std::cerr << "File opening failed\n"; return 1; }
    Universe u; // This universe is used in all further calculations, where the possible lattice will be constructed 
    Universe singleUniverse;
    std::string line;
    
    if (!readUniverseFile(&singleUniverse, file);) return 1;
    if (singleUniverse.latticeType) {
        if (singleUniverse.latticeType == 1) {
            const int cols = singleUniverse.lattice[0];
            const int rows = singleUniverse.lattice[1];
            if (cols <= 0 || rows <= 0) return 1;

            for (int i = 0; i < cols; ++i) {
                for (int j = 0; j < rows; ++j) {
                    Universe temp = singleUniverse;
                    temp.latticeType = 1;
                    temp.pos[0] = i * singeUniverse.boundDim[0];
                    temp.pos[1] = j * singeUniverse.boundDim[1];
                    u.subUniverse.push_back(std::move(temp));
                }
            }
        } else if (singleUniverse.latticeType == 2) {
            const int cols = singleUniverse.lattice[0];
            const int rows = singleUniverse.lattice[1];
            if (cols <= 0 || rows <= 0) return 1;
            
            double t = singleUniverse.boundDim[0];
            double d = t * std::sqrt(3);
            double lat = 1.5 * t; 
            for (int i = 0; i < cols; ++i) {
                for (int j = 0; j < rows; ++j) {
                    Universe temp = singleUniverse;
                    temp.latticeType = 2;
                    temp.pos[0] = lat * i;
                    temp.pos[1] =  d * j + (i & 1) * t;
                    u.subUniverse.push_back(std::move(temp));
                }
            }
        }
    } else {
        u = singleUniverse;
    }
    int iter = 100; 
    volumePointMethod(u, iter);
    volumeLineMethod(u, iter);
    renderSliceASCII(u, SliceAxis::Z, 0.0, 100, 50)
    
    // number indicates geometry file
    // volume calculation for:
    // 1     - Fuel pin slide 13
    // 2     - Hollow cylinder slide 25
    // 3     - Square lattice
    // 4     - Hex prism
    // 5     - Elliptical Torus

    // Visualization for:
    // 6     - Translations and rotations
    // 7     - Hex lattices

    
    // calculate volume 2 ways
    // surface test and distance function
    // ASCII viz



}
