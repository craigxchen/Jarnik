"""Referee statistics on the re-scanned clusters (own code).  Exact checks: Thm 1.1 (beta = unit*gamma_s*B,
B = D_s(A_v) != 0), gcd formula Norm g = prod p^(e-|a_x-a_y|), exact integer form of (1.1), bisector
parallelism (u_x + u_y parallel to g or g(1+-i)), Lemma 1.3 (float, C irrational).  Uniformity statistics with
two definitions of "K without complementary pairs".
Usage: python3 stats_ref.py ref_1e9.pkl
"""
import sys, pickle, itertools
from math import log, sqrt
import numpy as np
from refgauss import *
from records import profile, Dsel, GAMMA

data = pickle.load(open(sys.argv[1], 'rb'))
rows_out = []
nfail = 0; nbis = 0; nl13 = 0; npairs = 0
l13_slack_min = 1e9
for (N, fac, cl, C) in data:
    M = len(cl)
    u, G = normalise(cl)
    assert gnorm(G) == 1
    Np, fc, pis, rows, ks = profile(u)
    W = log(Np)
    mus = []; comp = []
    for x, y in itertools.combinations(range(M), 2):
        npairs += 1
        g = ggcd(u[x], u[y])
        ng = 1
        for p, e in fc.items():
            ng *= p ** (e - abs(rows[x][p] - rows[y][p]))
        assert gnorm(g) == ng
        beta = gexactdiv(gsub(u[x], u[y]), g)
        A = (1, 0)
        for p in fc:
            v = rows[x][p] - rows[y][p]
            if v > 0:
                A = gmul(A, gpow(pis[p], v))
            elif v < 0:
                A = gmul(A, gpow(gconj(pis[p]), -v))
        s = (ks[y] - ks[x]) % 4
        B = Dsel(A, s)
        cand = gmul(GAMMA[s], (B, 0))
        if B == 0 or not any(gmul(e_, cand) == beta for e_ in UNITS):
            nfail += 1
        assert gnorm(gsub(u[x], u[y])) == ng * gnorm(GAMMA[s]) * B * B
        S = gadd(u[x], u[y])
        dirs = [g, gmul(g, (0, 1))] if s % 2 == 0 else [gmul(g, (1, 1)), gmul(g, (1, -1))]
        if not any(S[0] * dd[1] - S[1] * dd[0] == 0 for dd in dirs):
            nbis += 1
        slack = sqrt(2) * C ** 2 * sqrt(ng) + 2 - M
        l13_slack_min = min(l13_slack_min, slack)
        if slack < -1e-9:
            nl13 += 1
        mus.append(log(Np) - log(ng) - W / 2); comp.append(ng == 1)
    mus = np.array(mus); comp = np.array(comp)
    mbar = mus.mean()
    K = mus.max() / mbar if mbar > 0 else np.inf
    if (~comp).any():
        nc = mus[~comp]
        K2auth = nc.max() / mbar if mbar > 0 else np.inf
        K2nat = nc.max() / nc.mean() if nc.mean() > 0 else np.inf
    else:
        K2auth = K2nat = np.nan
    rows_out.append((M, C, W / (M * log(M)), mbar, mus.max(), K, comp.any(), K2auth, K2nat))

print("clusters", len(data), "pairs", npairs, "Thm1.1 failures", nfail, "bisector failures", nbis,
      "Lemma1.3 violations", nl13, "min slack of Lemma 1.3 (sqrt2 C^2|g|+2-M)", round(l13_slack_min, 4))
for M in (4, 5):
    R = [r for r in rows_out if r[0] == M]
    K = np.array([r[5] for r in R]); fin = np.isfinite(K)
    print(f"M={M}: n={len(R)} minC={min(r[1] for r in R):.3f} minA={min(r[2] for r in R):.2f} "
          f"#mubar<=0: {int((~fin).sum())}  median K={np.median(K):.2f} pct25/75/95={np.round(np.percentile(K[fin],[25,75,95]),2)} "
          f"median mubar={np.median([r[3] for r in R]):.2f} median max mu={np.median([r[4] for r in R]):.2f} "
          f"frac complementary={np.mean([r[6] for r in R]):.3f}")
    K2a = np.array([r[7] for r in R if r[3] > 0], dtype=float)
    K2n = np.array([r[8] for r in R if r[3] > 0], dtype=float)
    print(f"      author-def K2 (max over non-comp / mean over ALL pairs), mubar>0 only: n={len(K2a)} median={np.nanmedian(K2a):.2f}")
    print(f"      natural K2 (max and mean both over non-comp pairs), mubar>0: median={np.nanmedian(K2n[np.isfinite(K2n)]):.2f}"
          f"  #inf={int(np.isinf(K2n).sum())}")
    # restrict to clusters with a complementary pair
    Kc = np.array([r[5] for r in R if r[6] and r[3] > 0]); Kn = np.array([r[5] for r in R if (not r[6]) and r[3] > 0])
    print(f"      K median among clusters WITH comp pair: {np.median(Kc):.2f}; WITHOUT comp pair: {np.median(Kn):.2f} (n={len(Kn)})")
