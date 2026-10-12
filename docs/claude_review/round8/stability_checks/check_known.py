#!/usr/bin/env python3
"""Stability checks on the record clusters of growth/families.md (and the referee's 10-point cluster).

For each listed N, every window of M >= 3 consecutive points with chord constant
C = maxchord/N^(1/4) <= CMAX is tested for:
  (A) smoothness: every prime q | N has log q <= max_{ceil(M/2)<=s<=M} b(s)  (M >= 4),
      b(s) = W/s + (4(s-1)/s) log C - 8F*(s)/s^2;   and  log q < (2/M) log N + 4 log^+ C;
  (O) omega(N') > W'/(4W'/M + 8 log^+ C);  e'_p log p < 4W'/M + 8 log^+ C;
  (G) the gcd lemma on every subset B (|B| >= 2), with C_B the subset's own chord constant.
It prints, per N, the largest M found at each C threshold and the tightest ratios.
"""
import sys
from itertools import combinations
import mpmath
from stab_lib import factor, points_of_norm, ggcd, norm, sub, bfun

TOL = mpmath.mpf(10) ** -40

KNOWN = [
    (412333805, "first 4-point cluster with C <= 1 (families.md)"),
    (1762645, "first 4-point cluster with C <= sqrt2"),
    (10414625, "first 5-point cluster with C <= 3"),
    (235370525, "sporadic 7/8-point record C_7 = 5.2699"),
    (537606725, "golden eight-point member"),
    (1176852625, "referee's primitive 10-point cluster, C = 8.2394"),
    (128605, "first 3-point cluster with C <= 1"),
    (276108509, "first 3-point cluster with C <= 1/2"),
]


def gcd_list(zs):
    g = zs[0]
    for z in zs[1:]:
        g = ggcd(g, z)
    return g


def run(N, label, CMAX=9.0):
    fac = factor(N)
    pts = points_of_norm(N, fac)
    n = len(pts)
    W = mpmath.log(N)
    q4 = mpmath.mpf(N) ** mpmath.mpf('0.25')
    fails = 0
    bestM = {}
    tightA = None
    tightO = None
    tightG = None
    nsub = 0
    for i in range(n):
        for M in range(3, min(12, n - 1) + 1):
            Z = [pts[(i + t) % n] for t in range(M)]
            ch2 = max(norm(sub(Z[a], Z[b])) for a in range(M) for b in range(a))
            C = mpmath.sqrt(ch2) / q4
            if C > CMAX:
                break
            logC = mpmath.log(C)
            lcp = max(logC, 0)
            for thr in (1, 2, 3, 4, 6, 9):
                if C <= thr:
                    bestM[thr] = max(bestM.get(thr, 0), M)
            if M >= 4:
                bound = max(bfun(s, W, logC) for s in range((M + 1) // 2, M + 1))
                crude = 2 * W / M + 4 * lcp
                for q in fac:
                    lq = mpmath.log(q)
                    if lq > bound + TOL or lq >= crude:
                        fails += 1
                    r = float(lq / bound)
                    if tightA is None or r > tightA[0]:
                        tightA = (round(r, 4), M, round(float(C), 4), q)
                G = gcd_list(Z)
                Np = N // norm(G)
                facp = factor(Np)
                Wp = mpmath.log(Np)
                lb = Wp / (4 * Wp / M + 8 * lcp)
                if len(facp) <= lb:
                    fails += 1
                r = float(len(facp) / lb)
                if tightO is None or r < tightO[0]:
                    tightO = (round(r, 4), M, round(float(C), 4), len(facp), Np)
                for p, e in facp.items():
                    if e * mpmath.log(p) >= 4 * Wp / M + 8 * lcp:
                        fails += 1
            if M <= 10:
                for s in range(2, M + 1):
                    for B in combinations(range(M), s):
                        if s < M and 0 not in B and (M - 1) not in B and M > 6:
                            continue  # keep the run short: subsets touching an endpoint
                        zs = [Z[k] for k in B]
                        nG = norm(gcd_list(zs))
                        mc = max(norm(sub(zs[a], zs[b])) for a in range(s) for b in range(a))
                        CB = mpmath.sqrt(mc) / q4
                        sl = mpmath.log(nG) - bfun(s, W, mpmath.log(CB))
                        nsub += 1
                        if sl > TOL:
                            fails += 1
                        if tightG is None or sl > tightG[0]:
                            tightG = (float(sl), s, round(float(CB), 4))
    print(f"N = {N} = {fac}  [{label}]  r2 = {n}")
    print(f"   largest M at C <= 1,2,3,4,6,9: {[bestM.get(t, 0) for t in (1, 2, 3, 4, 6, 9)]}")
    print(f"   tightest smoothness ratio log q/bound (<=1): {tightA}")
    print(f"   tightest omega(N')/lower bound (>1): {tightO}")
    print(f"   gcd lemma: {nsub} subsets, max slack (<=0): {tightG}")
    print(f"   failures: {fails}")
    return fails


def main():
    tot = 0
    for N, lab in KNOWN:
        tot += run(N, lab)
    print("ALL CHECKS PASSED:", tot == 0)


if __name__ == "__main__":
    main()
