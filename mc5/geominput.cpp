#include "mc.h"
#include <vector>

using std::array; using std::vector;

// --- GEOMETRY READ ---

int addQuadric(Geometry& g, float A,float B,float C,float D,float E,float F, float G,float H,float I,float J) {
    Shape s{};
    s.torus = 0;
    s.A=A; s.B=B; s.C=C; s.D=D; s.E=E; s.F=F; s.G=G; s.H=H; s.I=I; s.J=J;
    g.shapes.push_back(s);
    return (int)g.shapes.size()-1;
}

int addPlane(Geometry& g, float nx, float ny, float nz, float d) {
    return addQuadric(g, 0,0,0, 0,0,0, nx,ny,nz, -d);
}

int addTorus(Geometry& g, float a, float b, float R) {
    Shape s{}; s.torus = 1; s.A=a; s.B=b; s.C=R;
    s.D=s.E=s.F=s.G=s.H=s.J=0.0f;
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
void createBall(float a, float b, float c, Geometry& g) {
    int si = addQuadric(g, 1,1,1, 0,0,0, 0,0,0, -a*a); // x^2 + y^2 + z^2 - a^2 <= 0
    int li = pushLeaf(g, si);
    g.nodeRoot = li;
    g.shape = 1;
}
// a = radius, b = height
void createCylinder(float a, float b, float c, Geometry& g) {
    const float r = a;
    const float h2 = 0.5*b;

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
void createCylinderOpen(float a, float b, float c, Geometry& g) {
    const float r = a;
    int si = addQuadric(g, 1,1,0, 0,0,0, 0,0,0, -r*r);
    int li = pushLeaf(g, si);
    g.nodeRoot = li;
    g.shape = 3;
}
// a < z half space
// This becomes a general plane with transformations
void createPlane(float a, float b, float c, Geometry& g) {
    int s = addPlane(g, 0,0,1, (float)a);
    int l = pushLeaf(g, s);
    g.nodeRoot = l;
    g.shape = 4;
}
// side lengths a, b, c
void createCuboid(float a, float b, float c, Geometry& g) {
    const float hx = 0.5*a, hy = 0.5*b, hz = 0.5*c;

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
void createCuboidOpen(float a, float b, float c, Geometry& g) {
    const float hx = 0.5*a, hy = 0.5*b;

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
void createHexPrism(float a, float b, float c, Geometry& g) {
    const float s  = a;
    const float rA = s * std::sqrt(3.0) * 0.5;
    const float h2 = 0.5*b;

    const float nx1=1.0, ny1=0.0f;
    const float nx2=0.5, ny2= std::sqrt(3.0)*0.5;
    const float nx3=0.5, ny3=-std::sqrt(3.0)*0.5;

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
void createTorus(float a, float b, float c, Geometry& g) {
    int st = addTorus(g, a, b, c); 
    int lt = pushLeaf(g, st);
    g.nodeRoot = lt;
    g.shape = 8;
}