#!/usr/bin/env python3
"""R6. Consistency check of the referee's derivation (report §2.4) of a gcd lemma from the
Cilleruelo-Cordoba theorem (CC: an arc of length sqrt2 * rho^(1/2 - 1/(4 floor(m/2) + 2)) of a circle of radius
rho, rho^2 in N, contains at most m lattice points; Proc. AMS 115 (1992) 899-905, quoted, not re-proved here).

Derived statement (dividing s points on an arc of length L = C_arc sqrt R by their gcd g, n = Norm g):
   s even:  log n < W/s + (4(s-1)/s) log C_arc - (2(s-1)/s) log 2,
   s odd:   log n < (W + 4 s log C_arc - 2 s log 2)/(s+1).
Compared with Prop 2.2: log n <= W/s + (4(s-1)/s) log C_B - 8F*(s)/s^2 (C_B chord-based <= C_arc).
Tested on all subsets of angular windows (arc = the shorter arc spanned by the subset), mpmath 50 digits.
Prints which bound is smaller, per s. This is a sanity check of the derivation, not a proof of CC.
"""
import sys
from itertools import combinations
from math import atan2, pi
import mpmath
from rlib import factorize, lattice_points, gnorm, gsub, ggcd_list, F8vec

mpmath.mp.dps = 50


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60000
    CMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 6.0
    fails = 0
    tests = 0
    better = {}
    for N in range(5, NMAX + 1):
        fac = factorize(N)
        if any(q % 4 == 3 and e % 2 for q, e in fac.items()):
            continue
        pts = sorted(lattice_points(fac), key=lambda z: atan2(z[1], z[0]) % (2 * pi))
        n = len(pts)
        if n < 6:
            continue
        R = mpmath.sqrt(N)
        W = mpmath.log(N)
        ang = [mpmath.atan2(z[1], z[0]) for z in pts]
        for i in range(0, n, 2):
            win = [(i + t) % n for t in range(min(7, n // 2))]
            for s in range(2, len(win) + 1):
                for B in combinations(win, s):
                    a0 = ang[B[0]]
                    rel = [((ang[k] - a0 + mpmath.pi) % (2 * mpmath.pi)) - mpmath.pi for k in B]
                    width = max(rel) - min(rel)
                    if width >= mpmath.pi / 2:
                        continue
                    Carc = width * R / mpmath.sqrt(R)
                    if Carc > CMAX:
                        continue
                    Z = [pts[k] for k in B]
                    lnn = mpmath.log(gnorm(ggcd_list(Z)))
                    lc = mpmath.log(Carc)
                    if s % 2 == 0:
                        cc = W / s + 4 * mpmath.mpf(s - 1) / s * lc - 2 * mpmath.mpf(s - 1) / s * mpmath.log(2)
                    else:
                        cc = (W + 4 * s * lc - 2 * s * mpmath.log(2)) / (s + 1)
                    X = max(gnorm(gsub(Z[a], Z[b])) for a in range(s) for b in range(a))
                    CB = mpmath.sqrt(X) / mpmath.mpf(N) ** mpmath.mpf('0.25')
                    f8 = sum(c * mpmath.log(q) for q, c in F8vec(s).items())
                    p22 = W / s + 4 * mpmath.mpf(s - 1) / s * mpmath.log(CB) - f8 / s ** 2
                    tests += 1
                    if not lnn < cc:
                        fails += 1
                        print("CC-derived bound fails:", N, s, float(lnn), float(cc))
                    key = s
                    d = better.setdefault(key, [0, 0])
                    d[0 if cc < p22 else 1] += 1
    print(f"N <= {NMAX}, arc constant <= {CMAX}: {tests} subsets; CC-derived bound violated: {fails}")
    for s in sorted(better):
        print(f"   s = {s}: CC-derived bound smaller in {better[s][0]} cases, Prop 2.2 smaller in {better[s][1]}")
    print("ALL CHECKS PASSED:", fails == 0)


if __name__ == "__main__":
    main()
