#!/usr/bin/env python3
"""Exact checks of the sub-cluster identity (Lemma 2.1), the gcd lemma (Prop. 2.2), the
smoothness theorem (Thm 3.1) and the omega bound (Thm 3.3) of stability.md on actual clusters.

For every N in the scan range that is a sum of two squares with enough points, every window of
M >= 2 consecutive lattice points (angular order) whose chord constant C = maxchord/N^(1/4) is
<= CMAX is tested:

  (I)  exact integer identity, for every subset B of the window with |B| = s >= 2:
         prod_{B pairs} |z_i - z_j|^8 = N^(s(s-2)) * Norm(G_B)^(s^2)
                                       * prod_{non-constant split cuts on B} p^((2 S_tau - s)^2)
                                       * prod_{B pairs} |c_ij|^8 ,
       with G_B = gcd_B z and c_ij = (z_i - z_j)/gcd(z_i, z_j)   [only for windows with M <= MSUB];
  (G)  gcd lemma: log Norm(G_B) <= b(s) := W/s + (4(s-1)/s) log C_B - 8 F*(s)/s^2,
       C_B = max_B |z_i - z_j| / N^(1/4)                       [only for windows with M <= MSUB];
  (A)  smoothness: for M >= 4, every prime q | N has  log q <= max_{ceil(M/2)<=s<=M} b(s)
       (with the window's C), and the cruder log q < (2/M) log N + 4 log^+ C;
  (O)  for M >= 4: omega(N') > W'/(4W'/M + 8 log^+ C) where N' is the normalised norm,
       e'_p log p < 4 W'/M + 8 log^+ C for every p^e' || N', and the number of distinct cut
       partitions of the normalised profile exceeds W'/(2W'/M + 4 log^+ C).
Usage: python3 check_subcluster.py NMAX CMAX MSUB MMAX
"""
import sys
from itertools import combinations
import mpmath
TOL = mpmath.mpf(10) ** -40
from stab_lib import (factor, points_of_norm, ggcd, norm, sub, exact_div, gauss_prime,
                      conj, bfun, Fstar)


def levels(z, fac):
    """dict p -> v_pi(z) for split p | N."""
    out = {}
    for p, e in fac.items():
        if p % 4 != 1:
            continue
        pq = gauss_prime(p)
        a = 0
        w = z
        while True:
            t = (w[0] * pq[0] + w[1] * pq[1], w[1] * pq[0] - w[0] * pq[1])  # w * conj(pq)
            if t[0] % p == 0 and t[1] % p == 0:
                w = (t[0] // p, t[1] // p)
                a += 1
            else:
                break
        out[p] = a
    return out


def gcd_list(zs):
    g = zs[0]
    for z in zs[1:]:
        g = ggcd(g, z)
    return g


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    CMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
    MSUB = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    MMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    stats = dict(windows=0, subsets=0, ident_fail=0, gcd_fail=0, A_fail=0, A2_fail=0, O_fail=0,
                 pp_fail=0, part_fail=0, Atests=0)
    worstG = None   # largest gamma_B - b(s)
    worstA = None   # largest log q / bound
    worstO = None
    best_per_M = {}
    for N in range(2, NMAX + 1):
        fac = factor(N)
        if any(q % 4 == 3 and e % 2 for q, e in fac.items()):
            continue
        prod = 1
        for q, e in fac.items():
            if q % 4 == 1:
                prod *= e + 1
        if 4 * prod < 4:
            continue
        pts = points_of_norm(N, fac)
        n = len(pts)
        if n < 4:
            continue
        W = mpmath.log(N)
        q4 = mpmath.mpf(N) ** mpmath.mpf('0.25')
        lev = [levels(z, fac) for z in pts]
        for i in range(n):
            for M in range(2, min(MMAX, n - 1) + 1):
                idx = [(i + t) % n for t in range(M)]
                Z = [pts[k] for k in idx]
                ch2 = max(norm(sub(Z[a], Z[b])) for a in range(M) for b in range(a))
                Cw = mpmath.sqrt(ch2) / q4
                if Cw > CMAX:
                    break
                stats['windows'] += 1
                logC = mpmath.log(Cw)
                logCp = max(logC, 0)
                # ---- (A) smoothness, (O) omega
                if M >= 4:
                    stats['Atests'] += 1
                    bound = max(bfun(s, W, logC) for s in range((M + 1) // 2, M + 1))
                    crude = 2 * W / M + 4 * logCp
                    for q in fac:
                        lq = mpmath.log(q)
                        if lq > bound + TOL:
                            stats['A_fail'] += 1
                        if lq >= crude:
                            stats['A2_fail'] += 1
                        r = lq / bound
                        if worstA is None or r > worstA[0]:
                            worstA = (float(r), N, M, float(Cw), q, float(bound))
                    G = gcd_list(Z)
                    Np = N // norm(G)
                    facp = factor(Np)
                    Wp = mpmath.log(Np)
                    om = len(facp)
                    if om <= Wp / (4 * Wp / M + 8 * logCp):
                        stats['O_fail'] += 1
                    ratio = om / (M / 4.0)
                    key = M
                    if Cw <= 1:
                        if key not in best_per_M or ratio < best_per_M[key][0]:
                            best_per_M[key] = (ratio, N, float(Cw), om, Np)
                    for p, e in facp.items():
                        if e * mpmath.log(p) >= 4 * Wp / M + 8 * logCp:
                            stats['pp_fail'] += 1
                    # distinct cut partitions of the normalised profile
                    U = [exact_div(z, G) for z in Z]
                    parts = set()
                    for p, e in facp.items():
                        pq = gauss_prime(p)
                        la = [levels(u, {p: e})[p] for u in U]
                        for tau in range(1, e + 1):
                            parts.add(frozenset(k for k in range(M) if la[k] < tau))
                    if len(parts) <= Wp / (2 * Wp / M + 4 * logCp):
                        stats['part_fail'] += 1
                # ---- (I), (G) on subsets
                if M <= MSUB and M == min(MSUB, n - 1) or (M <= MSUB and i % 3 == 0):
                    for s in range(2, M + 1):
                        for B in combinations(range(M), s):
                            zs = [Z[k] for k in B]
                            GB = gcd_list(zs)
                            nG = norm(GB)
                            maxch = max(norm(sub(zs[a], zs[b])) for a in range(s) for b in range(a))
                            CB = mpmath.sqrt(maxch) / q4
                            gam = mpmath.log(nG)
                            slack = gam - bfun(s, W, mpmath.log(CB))
                            stats['subsets'] += 1
                            if slack > TOL:
                                stats['gcd_fail'] += 1
                            if worstG is None or slack > worstG[0]:
                                worstG = (float(slack), N, s, float(CB))
                            # exact identity
                            lhs = 1
                            rhs = N ** (s * (s - 2)) * nG ** (s * s)
                            for a in range(s):
                                for b in range(a):
                                    d = sub(zs[a], zs[b])
                                    lhs *= norm(d) ** 4
                                    g = ggcd(zs[a], zs[b])
                                    c = exact_div(d, g)
                                    rhs *= norm(c) ** 4
                            for p, e in fac.items():
                                if p % 4 != 1:
                                    continue
                                la = [lev[idx[k]][p] for k in B]
                                for tau in range(1, e + 1):
                                    S = sum(1 for x in la if x < tau)
                                    if 0 < S < s:
                                        rhs *= p ** ((2 * S - s) ** 2)
                            if lhs != rhs:
                                stats['ident_fail'] += 1
    print("scan N <=", NMAX, " CMAX =", CMAX, " MSUB =", MSUB, " MMAX =", MMAX)
    for k, v in stats.items():
        print(f"  {k:12s} {v}")
    print("  worst gcd-lemma slack gamma_B - b(s) (must be <= 0):", worstG)
    print("  worst smoothness ratio log q / bound (must be <= 1):", worstA)
    print("  tightest omega/(M/4) with C <= 1, per M:")
    for M in sorted(best_per_M):
        print("    M =", M, best_per_M[M])
    ok = all(stats[k] == 0 for k in ('ident_fail', 'gcd_fail', 'A_fail', 'A2_fail', 'O_fail',
                                      'pp_fail', 'part_fail'))
    print("ALL CHECKS PASSED:", ok)


if __name__ == "__main__":
    main()
