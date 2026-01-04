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


// --- TIMING FUNCTIONS ---

float nowSec() {
    using clock = std::chrono::steady_clock;
    static const auto t0 = clock::now();
    return std::chrono::duration<float>(clock::now() - t0).count();
}

float relFOM(float mean, float stdErr, float seconds) {
    if (mean<=0.0f || stdErr<=0.0f || seconds<=0.0f) return 0.0f;
    const float R = stdErr/mean;
    return 1.0f/(R*R*seconds);
}


// --- PRINTING FUNCTIONS ---

void printVolumeStats(std::string str, float mean, float stdErr, float seconds, float FOM) {
    const float rel = (mean != 0.0f) ? (stdErr / mean) : 0.0f;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << str << "\n"

              << "  volume_mean   = " << mean   << "\n"
              << "  volume_stdErr = " << stdErr << "\n"
              << "  rel_error     = " << rel    << "\n"
              << "  time_sec      = " << seconds    << "\n"
              << "  FOM           = " << FOM    << "\n";
}


// --- GEOMETRY HELPERS ---

void tallyPointRecursive(const Universe& u, const array<float,3>& pLocal, std::unordered_map<const Geometry*, long long>& hits) {
    for (const Geometry& g : u.geometries) {
        if (pointInGeom(pLocal, g)) ++hits[&g];
    }
    for (const Universe& su : u.subUniverse) {
        array<float,3> pSub = sub3(pLocal, su.pos);
        tallyPointRecursive(su, pSub, hits);
    }
}

void collectGeometries(const Universe& u, const std::string& prefix, vector<std::pair<const Geometry*, std::string>>& out) {
    for (size_t i = 0; i < u.geometries.size(); ++i)
        out.emplace_back(&u.geometries[i], prefix + "geom[" + std::to_string(i) + "]");

    for (size_t j = 0; j < u.subUniverse.size(); ++j)
        collectGeometries(u.subUniverse[j], prefix + "su[" + std::to_string(j) + "]/", out);
}

// --- VOLUME POINT ---

void volumePointMethod(Universe& u, int iter) {
    std::mt19937_64 rng(0x294823178913);
    const bool isSquare = (u.latticeType==1);
    const bool isHex    = (u.latticeType==2);

    vector<std::pair<const Geometry*,std::string>> roster;
    collectGeometries(u, "", roster);
    const auto boxAll = boundingBox(u);
    const float Vall = boxAll[0]*boxAll[1]*boxAll[2];
    int nCells = 1;
    if (isSquare || isHex) nCells = std::max(1, u.lattice[0]*u.lattice[1]);
    const auto cdim = cellBoxDims(u);
    const float Vcell = cdim[0]*cdim[1]*cdim[2];

    const float hx=0.5*boxAll[0], hy=0.5*boxAll[1], hz=0.5*boxAll[2];
    std::uniform_real_distribution<float> Ux(-hx,hx), Uy(-hy,hy), Uz(-hz,hz);
    const float hx_c=0.5*cdim[0], hy_c=0.5*cdim[1], hz_c=0.5*cdim[2];
    std::uniform_real_distribution<float> Ucx(-hx_c,hx_c), Ucy(-hy_c,hy_c), Ucz(-hz_c,hz_c);

    std::unordered_map<const Geometry*, long long> hits;
    hits.reserve(roster.size()*2);

    const float t0 = nowSec();

    if (!isSquare && !isHex) {
        for (int k=0;k<iter;++k) {
            array<float,3> p{Ux(rng),Uy(rng),Uz(rng)};
            tallyPointRecursive(u, p, hits);
        }
    } else {
        const int itPerCell = std::max(1, iter / nCells);
        if (isSquare) {
            for (int j=0;j<u.lattice[1];++j) {
                for (int i=0;i<u.lattice[0];++i) {
                    const auto C = squareCellCenter(u,i,j);
                    for (int k=0;k<itPerCell;++k) {
                        array<float,3> pLocal{Ucx(rng),Ucy(rng),Ucz(rng)};
                        array<float,3> pUniverse = add3(pLocal, C);
                        for (const Geometry& g : u.geometries) {
                            if (pointInGeom(pUniverse, g)) ++hits[&g];
                        }
                        for (const Universe& su : u.subUniverse) {
                            array<float,3> pSub = sub3(pUniverse, su.pos);
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
                        array<float,3> pLocal{Ucx(rng),Ucy(rng),Ucz(rng)};
                        array<float,3> pUniverse = add3(pLocal, C);
                        for (const Geometry& g : u.geometries) {
                            if (pointInGeom(pUniverse, g)) ++hits[&g];
                        }
                        for (const Universe& su : u.subUniverse) {
                            array<float,3> pSub = sub3(pUniverse, su.pos);
                            tallyPointRecursive(su, pSub, hits);
                        }
                    }
                }
            }
        }
    }

    const float t1 = nowSec();
    const float Veffective = (isSquare||isHex) ? ( (float)nCells * Vcell ) : Vall;
    const float Ntot = std::max(1.0f, (float)iter);
    const float sec = t1 - t0;

    for (auto& pr : roster) {
        const Geometry* gp = pr.first;
        const std::string& label = pr.second;
        const float cnt = (float)hits[gp];
        const float p = cnt / Ntot;
        const float mean   = Veffective * p;
        const float stdErr = std::sqrt(std::max(0.0f, p*(1.0f-p))) * Veffective / std::sqrt(Ntot);
        const float FOM    = relFOM(mean, stdErr, sec);
        printVolumeStats("PointVol " + label, mean, stdErr, sec, FOM);
    }
}


// --- VOLUME LINE ---

float marchLengthOneGeom(const Geometry& g, array<float,3> p0, const array<float,3>& dir, float Ltot) {
    Neutron ray; ray.pos=p0; ray.dir=dir;
    const float EPS=1e-9;
    float acc=0.0f, traveled=0.0f;
    bool inside = pointInGeom(ray.pos, g);
    while (traveled < Ltot - EPS) {
        float d = geometryCollision(ray, g);
        if (d < EPS || traveled + d > Ltot) {
            float seg = (Ltot - traveled);
            if (inside && seg>0.0f) acc+=seg;
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

void tallyLineRecursive(const Universe& u, const array<float,3>& p0_local, const array<float,3>& dir_local, float Ltot, std::unordered_map<const Geometry*, float>& lenSum) {
    for (const Geometry& g : u.geometries)
        lenSum[&g] += marchLengthOneGeom(g, p0_local, dir_local, Ltot);

    for (const Universe& su : u.subUniverse) {
        array<float,3> p_child = sub3(p0_local, su.pos);
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
    const float Vall = boxAll[0]*boxAll[1]*boxAll[2];
    int nCells = 1;
    if (isSquare || isHex) nCells = std::max(1, u.lattice[0]*u.lattice[1]);

    const auto cdim = cellBoxDims(u);
    const float Vcell = cdim[0]*cdim[1]*cdim[2];

    const float hx=0.5*boxAll[0], hy=0.5*boxAll[1], hz=0.5*boxAll[2];
    std::uniform_real_distribution<float> Ux(-hx,hx), Uy(-hy,hy), Uz(-hz,hz);

    const float hx_c=0.5*cdim[0], hy_c=0.5*cdim[1], hz_c=0.5*cdim[2];
    std::uniform_real_distribution<float> Ucx(-hx_c,hx_c), Ucy(-hy_c,hy_c), Ucz(-hz_c,hz_c);

    const float Lx_c = cdim[0], Ly_c = cdim[1], Lz_c = cdim[2];

    std::unordered_map<const Geometry*, float> sumLen, sumLen2;
    sumLen.reserve(roster.size()*2); sumLen2.reserve(roster.size()*2);

    const float t0 = nowSec();

    if (!isSquare && !isHex) {
        for (int k=0;k<iter;++k) {
            int axis = k % 3;
            array<float,3> p0, dir; float Ltot=0.0f;
            if (axis==0) { p0={-hx, Uy(rng), Uz(rng)}; dir={1,0,0}; Ltot=boxAll[0]; }
            else if (axis==1) { p0={Ux(rng), -hy, Uz(rng)}; dir={0,1,0}; Ltot=boxAll[1]; }
            else { p0={Ux(rng), Uy(rng), -hz}; dir={0,0,1}; Ltot=boxAll[2]; }

            std::unordered_map<const Geometry*, float> inc;
            tallyLineRecursive(u, p0, dir, Ltot, inc);
            for (auto& pr : inc) {
                float f = (Ltot>0.0f)? (pr.second / Ltot) : 0.0f;
                sumLen[pr.first]  += f;
                sumLen2[pr.first] += f*f;
            }
        }
    } else {
        const int itPerCell = std::max(1, iter / nCells);

        auto do_cell = [&](const array<float,3>& C) {
            for (int k=0;k<itPerCell;++k) {
                int axis = k % 3;
                array<float,3> pLocal, dir; float Ltot=0.0f;
                if (axis==0) { pLocal={-hx_c, Ucy(rng), Ucz(rng)}; dir={1,0,0}; Ltot=Lx_c; }
                else if (axis==1) { pLocal={Ucx(rng), -hy_c, Ucz(rng)}; dir={0,1,0}; Ltot=Ly_c; }
                else { pLocal={Ucx(rng), Ucy(rng), -hz_c}; dir={0,0,1}; Ltot=Lz_c; }

                std::unordered_map<const Geometry*, float> inc;

                for (const Geometry& g : u.geometries) {
                    float len = marchLengthOneGeom(g, pLocal, dir, Ltot);
                    if (len>0) inc[&g] += len;
                }
                for (const Universe& su : u.subUniverse) {
                    array<float,3> pParent = add3(pLocal, C);
                    array<float,3> pChild  = sub3(pParent, su.pos);
                    tallyLineRecursive(su, pChild, dir, Ltot, inc);
                }

                for (auto& pr : inc) {
                    float f = (Ltot>0.0f)? (pr.second / Ltot) : 0.0f;
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

    const float t1 = nowSec();
    const float N = std::max(1, iter);

    const float Veffective = (isSquare||isHex) ? ( (float)nCells * Vcell ) : Vall;

    for (auto& rp : roster) {
        const Geometry* gp = rp.first;
        const std::string& label = rp.second;

        float fbar = sumLen[gp] / N;
        float m2 = sumLen2[gp] / N;
        float varf = std::max(0.0f, m2 - fbar*fbar);
        float stdErr = std::sqrt(varf) * Veffective / std::sqrt((float)N);
        float mean = Veffective * fbar;
        float sec = t1 - t0;
        float FOM = relFOM(mean, stdErr, sec);

        printVolumeStats("LineVol " + label, mean, stdErr, sec, FOM);
    }
}


// --- VOLUME TORUS

vector<float> intersectTorusDistances(const Shape& s, const array<float,3>& P, const array<float,3>& D, float tmin = 0.0f, float tmax = 1e300) {
    const float EPS_DIR = 1e-16;
    const float EPS_KEEP = 1e-10;
    const float EPS_MERGE = 1e-7;
    vector<float> out;

    float d2 = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
    if (d2 < EPS_DIR) return out;
    const float invL = 1.0f / std::sqrt(d2);
    array<float,3> Du{ D[0]*invL, D[1]*invL, D[2]*invL };

    const float a = s.A;
    const float b = s.B;
    const float R = s.C;
    array<float,3> Cn{ s.D, s.E, s.F };
    array<float,3> Ax{ s.G, s.H, s.I };

    array<float,3> ez = normed(Ax);
    if (ez[0]==0.0f && ez[1]==0.0f && ez[2]==0.0f) ez = {0.0f,0.0f,1.0f};
    const array<float,3> tmp = (std::fabs(ez[2]) < 0.9) ? array<float,3>{0,0,1}
                                                              : array<float,3>{1,0,0};
    array<float,3> ex = normed(cross3(tmp, ez));
    array<float,3> ey = cross3(ez, ex);

    auto toLocal = [&](const array<float,3>& v)->array<float,3>{
        return { dot3(v,ex), dot3(v,ey), dot3(v,ez) };
    };

    array<float,3> P0{ P[0]-Cn[0], P[1]-Cn[1], P[2]-Cn[2] };
    array<float,3> Pl = toLocal(P0);
    array<float,3> Dl = toLocal(Du);

    const float x0 = Pl[0], y0 = Pl[1], z0 = Pl[2];
    const float ux = Dl[0], uy = Dl[1], uz = Dl[2];

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

    if (!(std::isfinite(C4) && std::isfinite(C3) && std::isfinite(C2) &&
          std::isfinite(C1) && std::isfinite(C0))) {
        return out;
    }

    float roots[4];
    const int rc = solveQuartic(C4, C3, C2, C1, C0, roots);

    vector<float> tmpRoots; tmpRoots.reserve(rc);
    for (int i=0;i<rc;++i) {
        const float t = roots[i];
        if (!std::isfinite(t)) continue;
        if (t <= EPS_KEEP) continue;
        if (t < tmin - EPS_MERGE) continue;
        if (t > tmax + EPS_MERGE) continue;
        tmpRoots.push_back(t);
    }
    if (tmpRoots.empty()) return out;

    auto torusF = [&](float t)->float {
        const float x = x0 + ux*t;
        const float y = y0 + uy*t;
        const float z = z0 + uz*t;
        const float r2 = x*x + y*y;
        const float S  = b2*r2 + a2*z*z + b2*(R2 - a2);
        return S*S - 4.0*R2*b2*b2*r2;
    };

    vector<float> filtered; filtered.reserve(tmpRoots.size());
    for (float t : tmpRoots) {
        const float x = x0 + ux*t, y = y0 + uy*t, z = z0 + uz*t;
        const float r2 = x*x + y*y;
        const float S  = b2*r2 + a2*z*z + b2*(R2 - a2);
        const float F  = (S*S - 4.0*R2*b2*b2*r2);
        if (std::isfinite(F) && std::abs(F) <= 1e-8 * (1.0f + S*S)) filtered.push_back(t);
    }
    tmpRoots.swap(filtered);
    if (tmpRoots.empty()) return out;

    std::sort(tmpRoots.begin(), tmpRoots.end());
    out.reserve(tmpRoots.size());
    for (float t : tmpRoots) {
        if (out.empty() || std::fabs(t - out.back()) > EPS_MERGE) out.push_back(t);
    }
    if (out.size() & 1U) {
        vector<float> fixed; fixed.reserve(out.size());
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
    if (roster.empty()) { printVolumeStats("LineVol (none)",0,0,0,0); return; }
    const Geometry* gp = roster[0].first;
    const std::string label = roster[0].second;
    const Shape* T = nullptr;
    for (const Shape& s : gp->shapes) if (s.torus) { T = &s; break; }
    if (!T) { printVolumeStats("LineVol " + label, 0,0,0,0); return; }
    const float a = T->A, b = T->B, R = T->C;
    const array<float,3> Cn{T->D, T->E, T->F};
    const array<float,3> Ax{T->G, T->H, T->I};
    array<float,3> ez = normed(Ax);
    if (ez == array<float,3>{0,0,0}) ez = {0,0,1};
    array<float,3> tmp = (std::fabs(ez[2]) < 0.9) ? array<float,3>{0,0,1} : array<float,3>{1,0,0};
    array<float,3> ex = normed(cross3(tmp, ez));
    array<float,3> ey = cross3(ez, ex);

    auto toWorld = [&](const array<float,3>& v)->array<float,3>{
        return add3(Cn, { v[0]*ex[0] + v[1]*ey[0] + v[2]*ez[0],
                          v[0]*ex[1] + v[1]*ey[1] + v[2]*ez[1],
                          v[0]*ex[2] + v[1]*ey[2] + v[2]*ez[2] });
    };
    const float Lx = 2.0*(R + b);
    const float Ly = 2.0*(R + b);
    const float Lz = 2.0*a;
    const array<float,3> dirW = { ey[0], ey[1], ey[2] };
    const float tmin = 0.0f;
    const float tmax = Ly;
    const float EPS0 = 1e-9;
    std::uniform_real_distribution<float> Ux(-0.5*Lx, 0.5*Lx);
    std::uniform_real_distribution<float> Uz(-0.5*Lz, 0.5*Lz);
    float sumChord=0.0f, sumChord2=0.0f;
    const float t0 = nowSec();
    for (int n=0; n<std::max(1,iter); ++n) {
        const float x = Ux(rng);
        const float z = Uz(rng);
        const array<float,3> p0W = toWorld({x, -0.5*Ly - EPS0, z});
        auto ts = intersectTorusDistances(*T, p0W, dirW, tmin, tmax + 2*EPS0);
        float chord = 0.0f;
        if (ts.size() >= 2) chord += std::max(0.0f, std::min(ts[1], tmax) - std::max(ts[0], tmin));
        if (ts.size() >= 4) chord += std::max(0.0f, std::min(ts[3], tmax) - std::max(ts[2], tmin));
        sumChord  += chord;
        sumChord2 += chord*chord;
    }

    const float t1 = nowSec();
    const int N = std::max(1, iter);
    const float W = Lx * Lz;
    const float meanChord = sumChord / float(N);
    const float m2 = sumChord2 / float(N);
    const float var = std::max(0.0f, m2 - meanChord*meanChord);
    const float mean   = W * meanChord;
    const float stdErr = W * std::sqrt(var) / std::sqrt(float(N));
    const float sec    = t1 - t0;
    const float FOM    = relFOM(mean, stdErr, sec);

    printVolumeStats("LineVol " + label, mean, stdErr, sec, FOM);
}

