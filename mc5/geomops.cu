#include "mc.h"
#include <algorithm>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <functional>
#include <unordered_map>
#include <random>
#include <chrono>
#include <iostream>
#include <iomanip>

using std::array; using std::vector;


// --- MATH FUNCTIONS ---

float dot3(const array<float,3>& a, const array<float,3>& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
array<float,3> add3(const array<float,3>& a,const array<float,3>& b) {return{a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
array<float,3> sub3(const array<float,3>& a,const array<float,3>& b) {return{a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
array<float,3> cross3(const array<float,3>& a, const array<float,3>& b) {return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};}
array<float,3> madd(const array<float,3>& p, const array<float,3>& d, float s) {return {p[0]+s*d[0], p[1]+s*d[1], p[2]+s*d[2]};}
array<float,3> normed(const array<float,3>& v) {
    float L = std::sqrt(dot3(v,v));
    if (L==0.0f) return {0,0,0};
    return {v[0]/L, v[1]/L, v[2]/L};
}

int solveQuadratic(float A,float B,float C,float r[2]) {
    const float EPS=1e-14;
    if (std::abs(A) < EPS) {
        if (std::abs(B) < EPS) return 0;
        r[0] = -C/B; return 1;
    }
    float disc = B*B - 4*A*C;
    if (disc < 0) return 0;
    float sdisc = std::sqrt(std::max(0.0f,disc));
    float q = -0.5*(B + (B>=0 ? sdisc : -sdisc));
    float t0 = q/A;
    float t1 = C/q;
    if (t0>t1) std::swap(t0,t1);
    r[0]=t0; r[1]=t1;
    return (disc==0.0f)?1:2;
}

int solveCubicMonic(float a,float b,float c,float r[3]) {
    const float TWOPI=6.28318530717958647692;
    const float EPS=1e-14;
    float a3 = a/3.0;
    float p = b - a*a/3.0;
    float q = (2.0*a*a*a)/27.0 - (a*b)/3.0 + c;
    float D = 0.25*q*q + (p*p*p)/27.0;
    if (D > EPS) {
        float s = std::sqrt(D);
        float u = std::cbrt(-0.5*q + s);
        float v = std::cbrt(-0.5*q - s);
        r[0] = u + v - a3;
        return 1;
    } else if (std::abs(D) <= EPS) {
        float u = std::cbrt(-0.5*q);
        r[0] = 2*u - a3;
        r[1] = -u - a3;
        return 2;
    } else {
        float rho = std::sqrt(-p/3.0);
        float theta = std::acos( (-0.5*q) / (rho*rho*rho) );
        r[0] = 2*rho*std::cos(theta/3.0) - a3;
        r[1] = 2*rho*std::cos((theta+TWOPI)/3.0) - a3;
        r[2] = 2*rho*std::cos((theta+2*TWOPI)/3.0) - a3;
        if (r[0]>r[1]) std::swap(r[0],r[1]);
        if (r[1]>r[2]) std::swap(r[1],r[2]);
        if (r[0]>r[1]) std::swap(r[0],r[1]);
        return 3;
    }
}

int solveQuartic(float a,float b,float c,float d,float e,float r[4]) {
    const float EPS=1e-14;
    if (std::abs(a) < EPS) {
        if (std::abs(b) < EPS) {
            return solveQuadratic(c,d,e,r);
        } else {
            float roots[3];
            int cnt = solveCubicMonic(c/b, d/b, e/b, roots);
            for (int i=0;i<cnt;++i) r[i]=roots[i];
            return cnt;
        }
    }
    float A=b/a, B=c/a, C=d/a, D=e/a;

    float A2=A*A;
    float p = -3.0*A2/8.0 + B;
    float q =  A*A2/8.0 - 0.5*A*B + C;
    float r0= -3.0*A2*A2/256.0 + A2*B/16.0 - 0.25*A*C + D;

    float cr[3]; int cc = solveCubicMonic(2*p, (p*p - 4*r0), -q*q, cr);
    float z = -1e300;
    for (int i=0;i<cc;++i) if (cr[i] > z) z = cr[i];
    if (z < 0) z = 0;
    float alpha = std::sqrt(z);

    float m = p + z;
    float beta, gamma;
    if (alpha > 1e-15) {
        beta  = 0.5*(m - q/alpha);
        gamma = 0.5*(m + q/alpha);
    } else {
        beta = gamma = 0.5*m;
    }

    float yroots[4]; int k=0;
    {
        float qr[2]; int qc = solveQuadratic(1.0,  alpha, beta, qr);
        for (int i=0;i<qc;++i) yroots[k++] = qr[i];
    }
    {
        float qr[2]; int qc = solveQuadratic(1.0, -alpha, gamma, qr);
        for (int i=0;i<qc;++i) yroots[k++] = qr[i];
    }

    for (int i=0;i<k;++i) r[i] = yroots[i] - A/4.0;

    std::sort(r, r+k);
    return k;
}


// --- GEOMETRY INFERENCE ---

vector<float> surfaceDist(const Neutron& n, const Shape& s) {
    const float EPS = 1e-10;
    const array<float,3>& P = n.pos;
    const array<float,3>& Dv = n.dir;
    float Dlen = std::sqrt(dot3(Dv,Dv));
    if (Dlen == 0.0f) return {};

    vector<float> out;

    if (!s.torus) {
        const float x=P[0], y=P[1], z=P[2];
        const float u=Dv[0], v=Dv[1], w=Dv[2];

        const float K = s.A*x*x + s.B*y*y + s.C*z*z
                       + s.D*x*y + s.E*y*z + s.F*x*z
                       + s.G*x + s.H*y + s.I*z + s.J;

        const float L = 2.0*(s.A*u*x + s.B*v*y + s.C*w*z)
                       + s.D*(v*x + u*y) + s.E*(w*y + v*z) + s.F*(w*x + u*z)
                       + s.G*u + s.H*v + s.I*w;

        const float M = s.A*u*u + s.B*v*v + s.C*w*w
                       + s.D*u*v + s.E*v*w + s.F*u*w;

        float roots[2]; int rc = solveQuadratic(M,L,K,roots);
        for (int i=0;i<rc;++i) {
            float t = roots[i];
            if (t > EPS) {
                out.push_back(t * Dlen);
            }
        }
        std::sort(out.begin(), out.end());
        return out;
    }

    const float a = s.A;
    const float b = s.B;
    const float R = s.C;

    array<float,3> Cn = {s.D, s.E, s.F};
    array<float,3> Ax = {s.G, s.H, s.I};
    array<float,3> ez = normed(Ax);
    if (ez == array<float,3>{0,0,0}) ez = {0,0,1};

    array<float,3> tmp = (std::fabs(ez[2]) < 0.9) ? array<float,3>{0,0,1}
                                                   : array<float,3>{1,0,0};
    array<float,3> ex = normed(cross3(tmp, ez));
    array<float,3> ey = cross3(ez, ex);

    array<float,3> P0 = sub3(P, Cn);
    auto toLocal = [&](const array<float,3>& v)->array<float,3> {
        return { dot3(v,ex), dot3(v,ey), dot3(v,ez) };
    };
    array<float,3> Pl = toLocal(P0);
    array<float,3> Dl = toLocal(Dv);

    const float x0=Pl[0], y0=Pl[1], z0=Pl[2];
    const float ux=Dl[0], uy=Dl[1], uz=Dl[2];

    const float b2 = b*b, a2 = a*a, R2 = R*R;
    const float q2 = ux*ux + uy*uy;
    const float q1 = 2.0*(x0*ux + y0*uy);
    const float q0 = x0*x0 + y0*y0;

    const float s2 = b2*q2 + a2*uz*uz;
    const float s1 = b2*q1 + 2.0*a2*z0*uz;
    const float s0 = b2*q0 + a2*z0*z0 + b2*(R2 - a2);

    const float C4 = s2*s2;
    const float C3 = 2.0*s2*s1;
    const float C2 = 2.0*s2*s0 + s1*s1 - 4.0*R2*b2*b2*q2;
    const float C1 = 2.0*s1*s0 - 4.0*R2*b2*b2*q1;
    const float C0 = s0*s0 - 4.0*R2*b2*b2*q0;

    float roots[4]; int rc = solveQuartic(C4,C3,C2,C1,C0,roots);
    for (int i=0;i<rc;++i) {
        float t = roots[i];
        if (t > EPS) {
            out.push_back(t * Dlen);
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

bool insideLeaf(const Shape& s, const array<float,3>& P) {
    if (!s.torus) {
        const float x=P[0], y=P[1], z=P[2];
        const float f = s.A*x*x + s.B*y*y + s.C*z*z
                       + s.D*x*y + s.E*y*z + s.F*x*z
                       + s.G*x + s.H*y + s.I*z + s.J;
        return f <= 0.0f;
    }
    const float a=s.A, b=s.B, R=s.C;
    array<float,3> Cn{ s.D, s.E, s.F };
    array<float,3> Ax{ s.G, s.H, s.I };
    array<float,3> ez = normed(Ax);
    if (ez==array<float,3>{0,0,0}) ez = {0,0,1};
    array<float,3> tmp = (std::fabs(ez[2]) < 0.9) ? array<float,3>{0,0,1}
                                                        : array<float,3>{1,0,0};
    array<float,3> ex = normed(cross3(tmp, ez));
    array<float,3> ey = cross3(ez, ex);

    array<float,3> Pl = sub3(P, Cn);
    const float x = dot3(Pl, ex), y = dot3(Pl, ey), z = dot3(Pl, ez);

    const float b2=b*b, a2=a*a, R2=R*R;
    const float r2 = x*x + y*y;
    const float S  = b2*r2 + a2*z*z + b2*(R2 - a2);
    const float F  = S*S - 4.0*R2*b2*b2*r2;
    return F <= 0.0f;
}

bool pointInGeom(const array<float, 3>& pos, const Geometry& g) {
    if (g.nodes.empty()) return false;
    std::function<bool(int)> eval = [&](int idx)->bool{
        const Node& n = g.nodes[idx];
        switch (n.op) {
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

float geometryCollision(const Neutron& n, const Geometry& g) {
    const float EPS = 1e-9;
    const float dlen = std::sqrt(dot3(n.dir,n.dir));
    if (dlen==0.0f) return -1.0;
    array<float,3> dHat{ n.dir[0]/dlen, n.dir[1]/dlen, n.dir[2]/dlen };

    vector<float> cands;
    for (const Shape& s : g.shapes) {
        vector<float> ds = surfaceDist(n, s);
        cands.insert(cands.end(), ds.begin(), ds.end());
    }
    if (cands.empty()) return -1.0;
    std::sort(cands.begin(), cands.end());
    vector<float> uniq;
    for (float v : cands) {
        if (v <= EPS) continue;
        if (uniq.empty() || std::fabs(v-uniq.back()) > 1e-8) uniq.push_back(v);
    }

    bool inside0 = pointInGeom(n.pos, g);
    for (float sDist : uniq) {
        float s0 = std::max(0.0, sDist - 1e-8);
        float s1 = sDist + 1e-8;
        bool inA = pointInGeom(madd(n.pos, dHat, s0), g);
        bool inB = pointInGeom(madd(n.pos, dHat, s1), g);
        if (inA != inB) {
            if (sDist <= EPS && inB==inside0) continue;
            return sDist;
        }
    }
    return -1.0;
}

float nearestCollision(const Neutron& n, const Universe& u) {
    const float INF = 1e300;
    float best = INF;

    for (const Geometry& g : u.geometries) {
        float d = geometryCollision(n, g);
        if (d > 0.0f && d < best) best = d;
    }

    for (const Universe& su : u.subUniverse) {
        Neutron nLocal = n;
        nLocal.pos[0] -= su.pos[0];
        nLocal.pos[1] -= su.pos[1];
        nLocal.pos[2] -= su.pos[2];
        float dBound = geometryCollision(nLocal, su.boundingGeometry);
        if (dBound > 0.0f && dBound <= best) {
            float dIn = nearestCollision(nLocal, su);
            if (dIn > 0.0f && dIn < best) best = dIn;
        }
    }

    return (best<INF)? best : -1.0;
}


// --- GEOMETRY HELPERS ---

array<float,3> cellBoxDims(const Universe& u) {
    if (u.latticeType==1) return {u.boundDim[0], u.boundDim[1], u.boundDim[2]};
    if (u.latticeType==2) return {u.boundDim[0], u.boundDim[1], u.boundDim[2]};
    return {u.boundDim[0], u.boundDim[1], u.boundDim[2]};
}

array<float, 3> boundingBox(Universe& u) {
    auto scale = [](float v) -> float { return 1.1 * std::max(0.0f, v); };
    if (!u.subUniverse.empty() && (u.latticeType==1 || u.latticeType==2)) {
        const int cols = std::max(1, u.lattice[0]);
        const int rows = std::max(1, u.lattice[1]);

        if (u.latticeType==1) {
            const float W = cols * u.boundDim[0];
            const float H = rows * u.boundDim[1];
            const float Z = u.boundDim[2];
            return { scale(W), scale(H), scale(Z) };
        } else {
            const float t = u.boundDim[0];
            const float s  = t;
            const float pitchX = 1.5*s;
            const float pitchY = std::sqrt(3.0)*s;
            const float W = (cols>1) ? ((cols-1)*pitchX + 2*s) : (2*s);
            const float H = rows*pitchY + ((cols>1)? 0.5*pitchY : 0.0f);
            const float Z = u.boundDim[2];
            return { scale(W), scale(H), scale(Z) };
        }
    }
    if (u.universeShape == 1) {
        const float a = u.boundDim[0];
        return { scale(2*a), scale(2*a), scale(2*a) };
    } else if (u.universeShape == 2 || u.universeShape == 3) {
        const float a = u.boundDim[0], b = u.boundDim[1] ? u.boundDim[1] : 1;
        return { scale(2*a), scale(2*a), scale(b) };
    } else if (u.universeShape == 4) {
        return { scale(u.boundDim[0]), scale(u.boundDim[1]), scale(u.boundDim[2]) };
    } else if (u.universeShape == 5 || u.universeShape == 6) {
        return { scale(u.boundDim[0]), scale(u.boundDim[1]), scale(u.boundDim[2] ? u.boundDim[2] : 1) };
    } else if (u.universeShape == 7) {
        const float a = u.boundDim[0], b = u.boundDim[1];
        return { scale(2*a), scale(std::sqrt(3.0)*a), scale(b) };
    } else if (u.universeShape == 8) {
        const float a = u.boundDim[0], b = u.boundDim[1], R = u.boundDim[2];
        return { scale(2*(R + a)), scale(2*(R + a)), scale(2*b) };
    }
    return { 0.f, 0.f, 0.f};
}

array<float,3> squareCellCenter(const Universe& u, int i, int j) {
    const int nx=u.lattice[0], ny=u.lattice[1];
    const float pitchX=u.boundDim[0], pitchY=u.boundDim[1];
    const float cx=0.5*(nx-1)*pitchX, cy=0.5*(ny-1)*pitchY;
    return { i*pitchX - cx, j*pitchY - cy, 0.0f };
}

array<float,3> hexCellCenter(const Universe& u, int q, int r) {
    const double sX=0.5*u.boundDim[0];
    const double sY=(u.boundDim[1]>0.0f)?(u.boundDim[1]/std::sqrt(3.0)):sX;
    const double s = std::max(1e-9,(sX>0&&sY>0)?std::min(sX,sY):std::max(sX,sY));
    const float pitchX=1.5*s, pitchY=std::sqrt(3.0)*s;
    const int cols=u.lattice[0], rows=u.lattice[1];
    const float cx=0.5*(cols-1)*pitchX;
    const float cy=0.5*((rows-1)*pitchY + (cols>1?0.5*pitchY:0.0f));
    const float x=q*pitchX - cx;
    const float y=r*pitchY + ((q&1)?0.5*pitchY:0.0f) - cy;
    return {x,y,0.0f};
}


// --- SIMULATION GEOMETRY ---

static bool pointInUniverseLocal(const Universe& u, const array<float,3>& pLocal) {
    return pointInGeom(pLocal, u.boundingGeometry);
}

const Geometry* findGeomAtRecursive(const Universe& u, const array<float,3>& pLocal) {
    // Preference for child Universes -> Faster
    for (const Universe& su : u.subUniverse) {
        array<float,3> ps = sub3(pLocal, su.pos);
        if (!pointInUniverseLocal(su, ps)) continue;
        if (auto g = findGeomAtRecursive(su, ps)) return g;
    }
    for (auto it = u.geometries.rbegin(); it != u.geometries.rend(); ++it) {
        if (pointInGeom(pLocal, *it)) return &(*it);
    }
    return nullptr;
}

const Geometry* findGeometryAt(const Universe& u, const array<float,3>& pWorld) {
    array<float,3> pLocal = sub3(pWorld, u.pos);
    if (!pointInUniverseLocal(u, pLocal)) return nullptr;
    if (auto g = findGeomAtRecursive(u, pLocal)) return g;
    return nullptr;
}

