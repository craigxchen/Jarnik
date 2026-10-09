"""Sanity check of Prop C.2 (lattice-line incidence bound) on actual clusters.

For a primitive v in Z[i] and the minimal arc Gamma containing a cluster (angular width D,
length L = R D), with eta = angular distance from Gamma to the line R v:
   M <= |v| L^2 / R + 2                 if eta = 0,
   M <= |v| L (eta + L/(2R)) + 1        if eta > 0.
Floating point is used for angles (double precision is ample for N <= 1e10); the bound is
asserted with a 1e-9 relative slack.  Also records, for each cluster, the smallest bound over
|v| <= 40 (how well the incidence method does on actual data).
"""
import random
from math import atan2, pi, sqrt, gcd, cos
from gx import reps, clusters

random.seed(2)


def angdist(a, b):
    d = (a - b) % (2 * pi)
    return min(d, 2 * pi - d)


def bound_for(cl, N, v):
    R = sqrt(N)
    angs = [atan2(z[1], z[0]) for z in cl]
    a0 = angs[0]
    rel = [((a - a0 + pi) % (2 * pi)) - pi for a in angs]
    lo, hi = min(rel), max(rel)
    D = hi - lo
    L = R * D
    av = atan2(v[1], v[0])
    # eta: distance from arc [a0+lo, a0+hi] to av + pi Z
    best = None
    for t in (av, av + pi):
        c = ((t - a0 + pi) % (2 * pi)) - pi
        if lo <= c <= hi:
            e = 0.0
        else:
            e = min(abs(c - lo), abs(c - hi))
        best = e if best is None else min(best, e)
    eta = best
    nv = sqrt(v[0] ** 2 + v[1] ** 2)
    if eta == 0.0:
        return nv * L * L / R + 2, eta, L
    return nv * L * (eta + L / (2 * R)) + 1, eta, L


def main():
    split = [5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113]
    Ns = [1176852625]
    for _ in range(600):
        N = 1
        for p in random.sample(split, random.randint(2, 5)):
            N *= p ** random.randint(1, 3)
        Ns.append(N)
    vs = [(a, b) for a in range(-40, 41) for b in range(0, 41)
          if (a, b) != (0, 0) and gcd(abs(a), b) == 1 and a * a + b * b <= 1600]
    nchk = 0
    ncl = 0
    rows = []
    for N in Ns:
        for cl in clusters(N, 8.3, 3):
            ncl += 1
            M = len(cl)
            bmin = None
            for v in vs:
                b, eta, L = bound_for(cl, N, v)
                assert M <= b * (1 + 1e-9), (N, cl, v, b)
                nchk += 1
                bmin = b if bmin is None else min(bmin, b)
            rows.append((M, bmin, N))
    print("Prop C.2: %d (cluster, v) checks on %d clusters: OK" % (nchk, ncl))
    rows.sort(key=lambda t: (-t[0], t[1]))
    print("largest clusters: (M, best incidence bound over |v|<=40, N)")
    for r in rows[:8]:
        print("  M=%d  bound=%.2f  N=%d" % r)
    ratio = sum(r[1] / r[0] for r in rows) / len(rows)
    print("mean (best bound / M) over clusters: %.2f" % ratio)


if __name__ == "__main__":
    main()
