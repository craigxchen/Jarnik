#!/usr/bin/env python3
"""Larger scan of the smoothness theorem (Thm 3.1) and the omega bound (Thm 3.3).

For every N <= NMAX that is a sum of two squares, every window of M >= 4 consecutive lattice
points with C = maxchord/N^(1/4) <= CMAX is tested (double precision; any case within 1e-9
of the bound would be re-tested exactly, none occurs):
  log P_+(N) <= max_{ceil(M/2)<=s<=M} b(s),      b(s) = W/s + (4(s-1)/s) log C - 8F*(s)/s^2,
  omega(N') > W'/(4W'/M + 8 log^+ C)            (N' = normalised norm).
Prints the tightest instances per M.  Usage: python3 check_scan_fast.py NMAX CMAX
"""
import sys
from math import log, sqrt, atan2, pi
from stab_lib import points_of_norm, ggcd, norm, Fstar


def spf_sieve(n):
    spf = list(range(n + 1))
    i = 2
    while i * i <= n:
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    return spf


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    CMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
    spf = spf_sieve(NMAX)
    F = {s: float(Fstar(s)) for s in range(2, 40)}
    tight = {}
    tightO = {}
    nwin = 0
    fails = 0
    for N in range(5, NMAX + 1):
        m, fac = N, {}
        while m > 1:
            p = spf[m]
            fac[p] = fac.get(p, 0) + 1
            m //= p
        if any(q % 4 == 3 and e % 2 for q, e in fac.items()):
            continue
        tau = 1
        for q, e in fac.items():
            if q % 4 == 1:
                tau *= e + 1
        if tau < 2:
            continue
        pts = points_of_norm(N, fac)
        n = len(pts)
        W = log(N)
        q4 = N ** 0.25
        Pplus = max(fac)
        lP = log(Pplus)
        for i in range(n):
            for M in range(4, min(30, n - 1) + 1):
                Z = [pts[(i + t) % n] for t in range(M)]
                ch2 = max(norm((Z[a][0] - Z[b][0], Z[a][1] - Z[b][1]))
                          for a in range(M) for b in range(a))
                C = sqrt(ch2) / q4
                if C > CMAX:
                    break
                nwin += 1
                lC = log(C)
                bound = max(W / s + 4 * (s - 1) / s * lC - 8 * F[s] / (s * s)
                            for s in range((M + 1) // 2, M + 1))
                r = lP / bound
                if r > 1 - 1e-9:
                    fails += 1
                if M not in tight or r > tight[M][0]:
                    tight[M] = (round(r, 4), N, round(C, 4))
                g = Z[0]
                for z in Z[1:]:
                    g = ggcd(g, z)
                Np = N // norm(g)
                if Np > 1:
                    om = 0
                    m = Np
                    while m > 1:
                        p = spf[m]
                        om += 1
                        while m % p == 0:
                            m //= p
                    Wp = log(Np)
                    lb = Wp / (4 * Wp / M + 8 * max(lC, 0))
                    if om <= lb + 1e-9:
                        fails += 1
                    ro = om / lb
                    if M not in tightO or ro < tightO[M][0]:
                        tightO[M] = (round(ro, 4), N, round(C, 4), om)
    print("scan N <=", NMAX, " CMAX =", CMAX, " windows with M >= 4:", nwin)
    print("tightest smoothness ratio log P_+(N)/bound per M (must be < 1):")
    for M in sorted(tight):
        print("   M =", M, tight[M])
    print("tightest omega(N')/lower bound per M (must be > 1):")
    for M in sorted(tightO):
        print("   M =", M, tightO[M])
    print("ALL CHECKS PASSED:", fails == 0)


if __name__ == "__main__":
    main()
