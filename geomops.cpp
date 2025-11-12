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

double dot3(const array<double,3>& a, const array<double,3>& b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
inline array<double,3> add3(const array<double,3>& a,const array<double,3>& b) {return{a[0]+b[0],a[1]+b[1],a[2]+b[2]};}
inline array<double,3> sub3(const array<double,3>& a,const array<double,3>& b) {return{a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
array<double,3> cross3(const array<double,3>& a, const array<double,3>& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
array<double,3> normed(const array<double,3>& v) {
    double L = std::sqrt(dot3(v,v));
    if (L==0.0) return {0,0,0};
    return {v[0]/L, v[1]/L, v[2]/L};
}

array<double,3> madd(const array<double,3>& p, const array<double,3>& d, double s) {
    return {p[0]+s*d[0], p[1]+s*d[1], p[2]+s*d[2]};
}

int solveQuadratic(double A,double B,double C,double r[2]) {
    const double EPS=1e-14;
    if (std::abs(A) < EPS) {
        if (std::abs(B) < EPS) return 0;
        r[0] = -C/B; return 1;
    }
    double disc = B*B - 4*A*C;
    if (disc < 0) return 0;
    double sdisc = std::sqrt(std::max(0.0,disc));
    double q = -0.5*(B + (B>=0 ? sdisc : -sdisc));
    double t0 = q/A;
    double t1 = C/q;
    if (t0>t1) std::swap(t0,t1);
    r[0]=t0; r[1]=t1;
    return (disc==0.0)?1:2;
}

int solveCubicMonic(double a,double b,double c,double r[3]) {
    const double TWOPI=6.28318530717958647692;
    const double EPS=1e-14;
    double a3 = a/3.0;
    double p = b - a*a/3.0;
    double q = (2.0*a*a*a)/27.0 - (a*b)/3.0 + c;
    double D = 0.25*q*q + (p*p*p)/27.0;
    if (D > EPS) {
        double s = std::sqrt(D);
        double u = std::cbrt(-0.5*q + s);
        double v = std::cbrt(-0.5*q - s);
        r[0] = u + v - a3;
        return 1;
    } else if (std::abs(D) <= EPS) {
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
        if (r[0]>r[1]) std::swap(r[0],r[1]);
        if (r[1]>r[2]) std::swap(r[1],r[2]);
        if (r[0]>r[1]) std::swap(r[0],r[1]);
        return 3;
    }
}

int solveQuartic(double a,double b,double c,double d,double e,double r[4]) {
    const double EPS=1e-14;
    if (std::abs(a) < EPS) {
        if (std::abs(b) < EPS) {
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
    if (alpha > 1e-15) {
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

vector<double> surfaceDist(const Neutron& n, const Shape& s) {
    const double EPS = 1e-10;
    const array<double,3>& P = n.pos;
    const array<double,3>& Dv = n.dir;
    double Dlen = std::sqrt(dot3(Dv,Dv));
    if (Dlen == 0.0) return {};

    vector<double> out;

    if (!s.torus) {
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
        for (int i=0;i<rc;++i) {
            double t = roots[i];
            if (t > EPS) {
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
    auto toLocal = [&](const array<double,3>& v)->array<double,3> {
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
    for (int i=0;i<rc;++i) {
        double t = roots[i];
        if (t > EPS) {
            out.push_back(t * Dlen);
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

bool insideLeaf(const Shape& s, const array<double,3>& P) {
    if (!s.torus) {
        const double x=P[0], y=P[1], z=P[2];
        const double f = s.A*x*x + s.B*y*y + s.C*z*z
                       + s.D*x*y + s.E*y*z + s.F*x*z
                       + s.G*x + s.H*y + s.I*z + s.J;
        return f <= 0.0;
    }
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

bool pointInGeom(const array<double, 3>& pos, const Geometry& g) {
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

double geometryCollision(const Neutron& n, const Geometry& g) {
    const double EPS = 1e-9;
    const double dlen = std::sqrt(dot3(n.dir,n.dir));
    if (dlen==0.0) return -1.0;
    array<double,3> dHat{ n.dir[0]/dlen, n.dir[1]/dlen, n.dir[2]/dlen };

    vector<double> cands;
    for (const Shape& s : g.shapes) {
        vector<double> ds = surfaceDist(n, s);
        cands.insert(cands.end(), ds.begin(), ds.end());
    }
    if (cands.empty()) return -1.0;
    std::sort(cands.begin(), cands.end());
    vector<double> uniq;
    for (double v : cands) {
        if (v <= EPS) continue;
        if (uniq.empty() || std::fabs(v-uniq.back()) > 1e-8) uniq.push_back(v);
    }

    bool inside0 = pointInGeom(n.pos, g);
    for (double sDist : uniq) {
        double s0 = std::max(0.0, sDist - 1e-8);
        double s1 = sDist + 1e-8;
        bool inA = pointInGeom(madd(n.pos, dHat, s0), g);
        bool inB = pointInGeom(madd(n.pos, dHat, s1), g);
        if (inA != inB) {
            if (sDist <= EPS && inB==inside0) continue;
            return sDist;
        }
    }
    return -1.0;
}

double nearestCollision(const Neutron& n, const Universe& u) {
    const double INF = 1e300;
    double best = INF;

    for (const Geometry& g : u.geometries) {
        double d = geometryCollision(n, g);
        if (d > 0.0 && d < best) best = d;
    }

    for (const Universe& su : u.subUniverse) {
        Neutron nLocal = n;
        nLocal.pos[0] -= su.pos[0];
        nLocal.pos[1] -= su.pos[1];
        nLocal.pos[2] -= su.pos[2];
        double dBound = geometryCollision(nLocal, su.boundingGeometry);
        if (dBound > 0.0 && dBound <= best) {
            double dIn = nearestCollision(nLocal, su);
            if (dIn > 0.0 && dIn < best) best = dIn;
        }
    }

    return (best<INF)? best : -1.0;
}

array<double,3> cellBoxDims(const Universe& u) {
    if (u.latticeType==1) return {u.boundDim[0], u.boundDim[1], u.boundDim[2]};
    if (u.latticeType==2) return {u.boundDim[0], u.boundDim[1], u.boundDim[2]};
    return {u.boundDim[0], u.boundDim[1], u.boundDim[2]};
}

array<double, 3> boundingBox(Universe& u) {
    auto scale = [](double v) -> double { return 1.1 * std::max(0.0, v); };
    if (!u.subUniverse.empty() && (u.latticeType==1 || u.latticeType==2)) {
        const int cols = std::max(1, u.lattice[0]);
        const int rows = std::max(1, u.lattice[1]);

        if (u.latticeType==1) {
            const double W = cols * u.boundDim[0];
            const double H = rows * u.boundDim[1];
            const double Z = u.boundDim[2];
            return { scale(W), scale(H), scale(Z) };
        } else {
            const double t = u.boundDim[0];
            const double s  = t;
            const double pitchX = 1.5*s;
            const double pitchY = std::sqrt(3.0)*s;
            const double W = (cols>1) ? ((cols-1)*pitchX + 2*s) : (2*s);
            const double H = rows*pitchY + ((cols>1)? 0.5*pitchY : 0.0);
            const double Z = u.boundDim[2];
            return { scale(W), scale(H), scale(Z) };
        }
    }
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

array<double,3> squareCellCenter(const Universe& u, int i, int j) {
    const int nx=u.lattice[0], ny=u.lattice[1];
    const double pitchX=u.boundDim[0], pitchY=u.boundDim[1];
    const double cx=0.5*(nx-1)*pitchX, cy=0.5*(ny-1)*pitchY;
    return { i*pitchX - cx, j*pitchY - cy, 0.0 };
}
array<double,3> hexCellCenter(const Universe& u, int q, int r) {
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

void tallyPointRecursive(const Universe& u, const array<double,3>& pLocal, std::unordered_map<const Geometry*, long long>& hits) {
    for (const Geometry& g : u.geometries) {
        if (pointInGeom(pLocal, g)) ++hits[&g];
    }
    for (const Universe& su : u.subUniverse) {
        array<double,3> pSub = sub3(pLocal, su.pos);
        tallyPointRecursive(su, pSub, hits);
    }
}

void collectGeometries(const Universe& u, const std::string& prefix, vector<std::pair<const Geometry*, std::string>>& out) {
    for (size_t i = 0; i < u.geometries.size(); ++i)
        out.emplace_back(&u.geometries[i], prefix + "geom[" + std::to_string(i) + "]");

    for (size_t j = 0; j < u.subUniverse.size(); ++j)
        collectGeometries(u.subUniverse[j], prefix + "su[" + std::to_string(j) + "]/", out);
}

inline double nowSec() {
    using clock = std::chrono::steady_clock;
    static const auto t0 = clock::now();
    return std::chrono::duration<double>(clock::now() - t0).count();
}

double relFOM(double mean, double stdErr, double seconds) {
    if (mean<=0.0 || stdErr<=0.0 || seconds<=0.0) return 0.0;
    const double R = stdErr/mean;
    return 1.0/(R*R*seconds);
}

void printStats(std::string str, double mean, double stdErr, double seconds, double FOM) {
    const double rel = (mean != 0.0) ? (stdErr / mean) : 0.0;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << str << "\n"

              << "  volume_mean   = " << mean   << "\n"
              << "  volume_stdErr = " << stdErr << "\n"
              << "  rel_error     = " << rel    << "\n"
              << "  time_sec      = " << seconds    << "\n"
              << "  FOM           = " << FOM    << "\n";
}

void volumePointMethod(Universe& u, int iter) {
    std::mt19937_64 rng(0x294823178913);
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

    if (!isSquare && !isHex) {
        for (int k=0;k<iter;++k) {
            array<double,3> p{Ux(rng),Uy(rng),Uz(rng)};
            tallyPointRecursive(u, p, hits);
        }
    } else {
        const int itPerCell = std::max(1, iter / nCells);
        if (isSquare) {
            for (int j=0;j<u.lattice[1];++j) {
                for (int i=0;i<u.lattice[0];++i) {
                    const auto C = squareCellCenter(u,i,j);
                    for (int k=0;k<itPerCell;++k) {
                        array<double,3> pLocal{Ucx(rng),Ucy(rng),Ucz(rng)};
                        array<double,3> pUniverse = add3(pLocal, C);
                        for (const Geometry& g : u.geometries) {
                            if (pointInGeom(pUniverse, g)) ++hits[&g];
                        }
                        for (const Universe& su : u.subUniverse) {
                            array<double,3> pSub = sub3(pUniverse, su.pos);
                            tallyPointRecursive(su, pSub, hits);
                        }
                    }
                }
            }
        } else {
            for (int r=0;r<u.lattice[1];++r) {
                for (int q=0;q<u.lattice[0];++q) {
                    const auto C = hexCellCenter(u,q,r);
                    for (int k=0;k<itPerCell;++k) {
                        array<double,3> pLocal{Ucx(rng),Ucy(rng),Ucz(rng)};
                        array<double,3> pUniverse = add3(pLocal, C);
                        for (const Geometry& g : u.geometries) {
                            if (pointInGeom(pUniverse, g)) ++hits[&g];
                        }
                        for (const Universe& su : u.subUniverse) {
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

    for (auto& pr : roster) {
        const Geometry* gp = pr.first;
        const std::string& label = pr.second;
        const double cnt = (double)hits[gp];
        const double p = cnt / Ntot;
        const double mean   = Veffective * p;
        const double stdErr = std::sqrt(std::max(0.0, p*(1.0-p))) * Veffective / std::sqrt(Ntot);
        const double FOM    = relFOM(mean, stdErr, sec);
        printStats("PointVol " + label, mean, stdErr, sec, FOM);
    }
}

double marchLengthOneGeom(const Geometry& g, array<double,3> p0, const array<double,3>& dir, double Ltot) {
    Neutron ray; ray.pos=p0; ray.dir=dir;
    const double EPS=1e-9;
    double acc=0.0, traveled=0.0;
    bool inside = pointInGeom(ray.pos, g);
    while (traveled < Ltot - EPS) {
        double d = geometryCollision(ray, g);
        if (d < EPS || traveled + d > Ltot) {
            double seg = (Ltot - traveled);
            if (inside && seg>0.0) acc+=seg;
            break;
        }
        if (inside) acc += d;
        ray.pos[0]+=dir[0]*d; ray.pos[1]+=dir[1]*d; ray.pos[2]+=dir[2]*d;
        ray.pos[0]+=dir[0]*1e-8; ray.pos[1]+=dir[1]*1e-8; ray.pos[2]+=dir[2]*1e-8;
        traveled += d + 1e-8;
        inside = pointInGeom(ray.pos, g);
    }
    return acc;
}

void tallyLineRecursive(const Universe& u, const array<double,3>& p0_local, const array<double,3>& dir_local, double Ltot, std::unordered_map<const Geometry*, double>& lenSum) {
    for (const Geometry& g : u.geometries)
        lenSum[&g] += marchLengthOneGeom(g, p0_local, dir_local, Ltot);

    for (const Universe& su : u.subUniverse) {
        array<double,3> p_child = sub3(p0_local, su.pos);
        tallyLineRecursive(su, p_child, dir_local, Ltot, lenSum);
    }
}

void volumeLineMethod(Universe& u, int iter) {
    std::mt19937_64 rng(0x81723561289736);
    const bool isSquare = (u.latticeType==1);
    const bool isHex = (u.latticeType==2);

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

    const double Lx_c = cdim[0], Ly_c = cdim[1], Lz_c = cdim[2];

    std::unordered_map<const Geometry*, double> sumLen, sumLen2;
    sumLen.reserve(roster.size()*2); sumLen2.reserve(roster.size()*2);

    const double t0 = nowSec();

    if (!isSquare && !isHex) {
        for (int k=0;k<iter;++k) {
            int axis = k % 3;
            array<double,3> p0, dir; double Ltot=0.0;
            if (axis==0) { p0={-hx, Uy(rng), Uz(rng)}; dir={1,0,0}; Ltot=boxAll[0]; }
            else if (axis==1) { p0={Ux(rng), -hy, Uz(rng)}; dir={0,1,0}; Ltot=boxAll[1]; }
            else { p0={Ux(rng), Uy(rng), -hz}; dir={0,0,1}; Ltot=boxAll[2]; }

            std::unordered_map<const Geometry*, double> inc;
            tallyLineRecursive(u, p0, dir, Ltot, inc);
            for (auto& pr : inc) {
                double f = (Ltot>0.0)? (pr.second / Ltot) : 0.0;
                sumLen[pr.first]  += f;
                sumLen2[pr.first] += f*f;
            }
        }
    } else {
        const int itPerCell = std::max(1, iter / nCells);

        auto do_cell = [&](const array<double,3>& C) {
            for (int k=0;k<itPerCell;++k) {
                int axis = k % 3;
                array<double,3> pLocal, dir; double Ltot=0.0;
                if (axis==0) { pLocal={-hx_c, Ucy(rng), Ucz(rng)}; dir={1,0,0}; Ltot=Lx_c; }
                else if (axis==1) { pLocal={Ucx(rng), -hy_c, Ucz(rng)}; dir={0,1,0}; Ltot=Ly_c; }
                else { pLocal={Ucx(rng), Ucy(rng), -hz_c}; dir={0,0,1}; Ltot=Lz_c; }

                std::unordered_map<const Geometry*, double> inc;

                for (const Geometry& g : u.geometries) {
                    double len = marchLengthOneGeom(g, pLocal, dir, Ltot);
                    if (len>0) inc[&g] += len;
                }
                for (const Universe& su : u.subUniverse) {
                    array<double,3> pParent = add3(pLocal, C);
                    array<double,3> pChild  = sub3(pParent, su.pos);
                    tallyLineRecursive(su, pChild, dir, Ltot, inc);
                }

                for (auto& pr : inc) {
                    double f = (Ltot>0.0)? (pr.second / Ltot) : 0.0;
                    sumLen[pr.first]  += f;
                    sumLen2[pr.first] += f*f;
                }
            }
        };

        if (isSquare) {
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

    for (auto& rp : roster) {
        const Geometry* gp = rp.first;
        const std::string& label = rp.second;

        double fbar = sumLen[gp] / N;
        double m2 = sumLen2[gp] / N;
        double varf = std::max(0.0, m2 - fbar*fbar);
        double stdErr = std::sqrt(varf) * Veffective / std::sqrt((double)N);
        double mean = Veffective * fbar;
        double sec = t1 - t0;
        double FOM = relFOM(mean, stdErr, sec);

        printStats("LineVol " + label, mean, stdErr, sec, FOM);
    }
}

vector<double> intersectTorusDistances(const Shape& s, const array<double,3>& P, const array<double,3>& D, double tmin = 0.0, double tmax = 1e300) {
    const double EPS_DIR = 1e-16;
    const double EPS_KEEP = 1e-10;
    const double EPS_MERGE = 1e-7;
    vector<double> out;

    double d2 = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
    if (d2 < EPS_DIR) return out;
    const double invL = 1.0 / std::sqrt(d2);
    array<double,3> Du{ D[0]*invL, D[1]*invL, D[2]*invL };

    const double a = s.A;
    const double b = s.B;
    const double R = s.C;
    array<double,3> Cn{ s.D, s.E, s.F };
    array<double,3> Ax{ s.G, s.H, s.I };

    array<double,3> ez = normed(Ax);
    if (ez[0]==0.0 && ez[1]==0.0 && ez[2]==0.0) ez = {0.0,0.0,1.0};
    const array<double,3> tmp = (std::fabs(ez[2]) < 0.9) ? array<double,3>{0,0,1}
                                                              : array<double,3>{1,0,0};
    array<double,3> ex = normed(cross3(tmp, ez));
    array<double,3> ey = cross3(ez, ex);

    auto toLocal = [&](const array<double,3>& v)->array<double,3>{
        return { dot3(v,ex), dot3(v,ey), dot3(v,ez) };
    };

    array<double,3> P0{ P[0]-Cn[0], P[1]-Cn[1], P[2]-Cn[2] };
    array<double,3> Pl = toLocal(P0);
    array<double,3> Dl = toLocal(Du);

    const double x0 = Pl[0], y0 = Pl[1], z0 = Pl[2];
    const double ux = Dl[0], uy = Dl[1], uz = Dl[2];

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

    if (!(std::isfinite(C4) && std::isfinite(C3) && std::isfinite(C2) &&
          std::isfinite(C1) && std::isfinite(C0))) {
        return out;
    }

    double roots[4];
    const int rc = solveQuartic(C4, C3, C2, C1, C0, roots);

    vector<double> tmpRoots; tmpRoots.reserve(rc);
    for (int i=0;i<rc;++i) {
        const double t = roots[i];
        if (!std::isfinite(t)) continue;
        if (t <= EPS_KEEP) continue;
        if (t < tmin - EPS_MERGE) continue;
        if (t > tmax + EPS_MERGE) continue;
        tmpRoots.push_back(t);
    }
    if (tmpRoots.empty()) return out;

    auto torusF = [&](double t)->double {
        const double x = x0 + ux*t;
        const double y = y0 + uy*t;
        const double z = z0 + uz*t;
        const double r2 = x*x + y*y;
        const double S  = b2*r2 + a2*z*z + b2*(R2 - a2);
        return S*S - 4.0*R2*b2*b2*r2;
    };

    vector<double> filtered; filtered.reserve(tmpRoots.size());
    for (double t : tmpRoots) {
        const double x = x0 + ux*t, y = y0 + uy*t, z = z0 + uz*t;
        const double r2 = x*x + y*y;
        const double S  = b2*r2 + a2*z*z + b2*(R2 - a2);
        const double F  = (S*S - 4.0*R2*b2*b2*r2);
        if (std::isfinite(F) && std::abs(F) <= 1e-8 * (1.0 + S*S)) filtered.push_back(t);
    }
    tmpRoots.swap(filtered);
    if (tmpRoots.empty()) return out;

    std::sort(tmpRoots.begin(), tmpRoots.end());
    out.reserve(tmpRoots.size());
    for (double t : tmpRoots) {
        if (out.empty() || std::fabs(t - out.back()) > EPS_MERGE) out.push_back(t);
    }
    if (out.size() & 1U) {
        vector<double> fixed; fixed.reserve(out.size());
        for (size_t i=0;i<out.size();) {
            if (i+1<out.size() && std::fabs(out[i+1]-out[i]) <= 10*EPS_MERGE) {
                fixed.push_back(0.5*(out[i]+out[i+1]));
                i += 2;
            } else {
                fixed.push_back(out[i]);
                i += 1;
            }
        }
        out.swap(fixed);
        if (out.size() & 1U) {
            out.erase(out.begin());
        }
    }
    if (out.size() > 4) out.resize(4);
    return out;
}

void volumeLineMethodTorus(Universe& u, int iter) {
    std::mt19937_64 rng(0x919823789);

    vector<std::pair<const Geometry*,std::string>> roster;
    collectGeometries(u, "", roster);
    if (roster.empty()) { printStats("LineVol (none)",0,0,0,0); return; }
    const Geometry* gp = roster[0].first;
    const std::string label = roster[0].second;
    const Shape* T = nullptr;
    for (const Shape& s : gp->shapes) if (s.torus) { T = &s; break; }
    if (!T) { printStats("LineVol " + label, 0,0,0,0); return; }
    const double a = T->A, b = T->B, R = T->C;
    const array<double,3> Cn{T->D, T->E, T->F};
    const array<double,3> Ax{T->G, T->H, T->I};
    array<double,3> ez = normed(Ax);
    if (ez == array<double,3>{0,0,0}) ez = {0,0,1};
    array<double,3> tmp = (std::fabs(ez[2]) < 0.9) ? array<double,3>{0,0,1} : array<double,3>{1,0,0};
    array<double,3> ex = normed(cross3(tmp, ez));
    array<double,3> ey = cross3(ez, ex);

    auto toWorld = [&](const array<double,3>& v)->array<double,3>{
        return add3(Cn, { v[0]*ex[0] + v[1]*ey[0] + v[2]*ez[0],
                          v[0]*ex[1] + v[1]*ey[1] + v[2]*ez[1],
                          v[0]*ex[2] + v[1]*ey[2] + v[2]*ez[2] });
    };
    const double Lx = 2.0*(R + b);
    const double Ly = 2.0*(R + b);
    const double Lz = 2.0*a;
    const array<double,3> dirW = { ey[0], ey[1], ey[2] };
    const double tmin = 0.0;
    const double tmax = Ly;
    const double EPS0 = 1e-9;
    std::uniform_real_distribution<double> Ux(-0.5*Lx, 0.5*Lx);
    std::uniform_real_distribution<double> Uz(-0.5*Lz, 0.5*Lz);
    double sumChord=0.0, sumChord2=0.0;
    const double t0 = nowSec();
    for (int n=0; n<std::max(1,iter); ++n) {
        const double x = Ux(rng);
        const double z = Uz(rng);
        const array<double,3> p0W = toWorld({x, -0.5*Ly - EPS0, z});
        auto ts = intersectTorusDistances(*T, p0W, dirW, tmin, tmax + 2*EPS0);
        double chord = 0.0;
        if (ts.size() >= 2) chord += std::max(0.0, std::min(ts[1], tmax) - std::max(ts[0], tmin));
        if (ts.size() >= 4) chord += std::max(0.0, std::min(ts[3], tmax) - std::max(ts[2], tmin));
        sumChord  += chord;
        sumChord2 += chord*chord;
    }

    const double t1 = nowSec();
    const int N = std::max(1, iter);
    const double W = Lx * Lz;
    const double meanChord = sumChord / double(N);
    const double m2 = sumChord2 / double(N);
    const double var = std::max(0.0, m2 - meanChord*meanChord);
    const double mean   = W * meanChord;
    const double stdErr = W * std::sqrt(var) / std::sqrt(double(N));
    const double sec    = t1 - t0;
    const double FOM    = relFOM(mean, stdErr, sec);

    printStats("LineVol " + label, mean, stdErr, sec, FOM);
}

static bool pointInUniverseLocal(const Universe& u, const array<double,3>& pLocal) {
    return pointInGeom(pLocal, u.boundingGeometry);
}

const Geometry* findGeomAtRecursive(const Universe& u, const array<double,3>& pLocal) {
    // Preference for child Universes -> Faster
    for (const Universe& su : u.subUniverse) {
        array<double,3> ps = sub3(pLocal, su.pos);
        if (!pointInUniverseLocal(su, ps)) continue;
        if (auto g = findGeomAtRecursive(su, ps)) return g;
    }
    for (auto it = u.geometries.rbegin(); it != u.geometries.rend(); ++it) {
        if (pointInGeom(pLocal, *it)) return &(*it);
    }
    return nullptr;
}

const Geometry* findGeometryAt(const Universe& u, const array<double,3>& pWorld) {
    array<double,3> pLocal = sub3(pWorld, u.pos);
    if (!pointInUniverseLocal(u, pLocal)) return nullptr;
    if (auto g = findGeomAtRecursive(u, pLocal)) return g;
    return nullptr;
}

double macroscopicSigmaAt(const Universe& u, const array<double,3>& pWorld, double E, const vector<int>& mts_total, const Geometry** gOut) {
    if (gOut) *gOut = nullptr;
    const Geometry* g = findGeometryAt(u, pWorld);
    if (!g) return 0.0;
    if (gOut) *gOut = g;
    return sigmaTot(g->mats, E, mts_total);
}

double majorSigma(const Universe& u, double E, const vector<int>& mts_total) {
    double m = 0.0;
    for (const Geometry& g : u.geometries)
        m = std::max(m, sigmaTot(g.mats, E, mts_total));
    for (const Universe& su : u.subUniverse)
        m = std::max(m, majorSigma(su, E, mts_total));
    return m;
}
