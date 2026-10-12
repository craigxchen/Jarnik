#!/usr/bin/env python3
"""R5. Exhaustive EXACT check of Theorem 3.1 (smoothness) on every circle N <= NMAX.

For every N <= NMAX and every window of M >= 4 angularly consecutive lattice points with chord constant
C = maxchord / N^(1/4) <= CMAX, and every prime q | N (content included), test in integers
     exists s in [ceil(M/2), M]:  q^(s^2) * exp(8F*(s)) * N^(s(s-2)) <= X^(2 s (s-1)),
X = max squared chord (this is log q <= max_s b_C(s), the refined form (3.1)).
Also a crude form, exactly: the test  q^(2M) * N^M < N^4 * max(N, X^2)^M  is equivalent to
log q < 2W/M + 2 log^+ C (C^4 = X^2/N), which is STRONGER than the stated crude form
log q < 2W/M + 4 log^+ C, so passing it implies the stated form.
Reports the tightest ratio log q / max_s b_C(s) per M.
Usage: python3 r5_smooth_exhaustive.py NMAX CMAX
"""
import sys
from math import atan2, pi, log
from rlib import lattice_points, gsub, gnorm, F8int, F8float


def spf_table(n):
    spf = list(range(n + 1))
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2 * 10 ** 6
    CMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
    spf = spf_table(NMAX)
    F8 = {s: F8int(s) for s in range(2, 40)}
    F8f = {s: F8float(s) for s in range(2, 40)}
    nwin = 0
    fails = 0
    tight = {}
    for N in range(5, NMAX + 1):
        m, fac = N, {}
        while m > 1:
            p = spf[m]
            fac[p] = fac.get(p, 0) + 1
            m //= p
        if any(q % 4 == 3 and e % 2 for q, e in fac.items()):
            continue
        r = 1
        for q, e in fac.items():
            if q % 4 == 1:
                r *= e + 1
        if 4 * r < 4:
            continue
        q4 = N ** 0.25
        # quick reject: need >= 4 points within chord CMAX*q4; the 4r points are spread on the circle
        pts = sorted(lattice_points(fac), key=lambda z: atan2(z[1], z[0]) % (2 * pi))
        n = len(pts)
        if n < 5:
            continue
        lim2 = (CMAX * q4) ** 2 * (1 + 1e-12)
        W = log(N)
        for i in range(n):
            for M in range(4, n):
                Z = [pts[(i + t) % n] for t in range(M)]
                X = max(gnorm(gsub(Z[a], Z[b])) for a in range(M) for b in range(a))
                if X > lim2:
                    break
                C = (X ** 0.5) / q4
                if C > CMAX:
                    break
                nwin += 1
                ss = range((M + 1) // 2, M + 1)
                for q in fac:
                    if not any(q ** (s * s) * F8[s] * N ** (s * (s - 2)) <= X ** (2 * s * (s - 1)) for s in ss):
                        fails += 1
                        print("FAIL refined", N, M, q, C)
                    if not q ** (2 * M) * N ** M < N ** 4 * max(N, X * X) ** M:
                        fails += 1
                        print("FAIL crude", N, M, q, C)
                    bnd = max(W / s + 4 * (s - 1) / s * log(C) - F8f[s] / s ** 2 for s in ss)
                    rat = log(q) / bnd if bnd > 0 else float('inf')
                    if M not in tight or rat > tight[M][0]:
                        tight[M] = (round(rat, 4), N, round(C, 4), q)
    print(f"N <= {NMAX}, C <= {CMAX}: {nwin} windows with M >= 4 tested exactly; failures {fails}")
    for M in sorted(tight):
        print(f"   M = {M}: tightest log q / max_s b_C(s) = {tight[M]}")
    print("ALL CHECKS PASSED:", fails == 0)


if __name__ == "__main__":
    main()
