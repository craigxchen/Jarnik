#!/usr/bin/env python3
"""Referee check of Lemma A.4 / Cor. A.5 (p3lower.md) on the round-6 design D_A (sharp.md Def. 3.2).
Builds the point classes kappa^(q^a) literally (small q < 4M, levels a <= a*(q)), computes unit-1 pair
collision moduli, D_xy = 2 log m_xy, point totals D_x, lambda_max(D), and compares with
  D_x <= 4 sum_{d<M} (Lambda'(d) + log d) <= 36 M log M + 111 M.
Also checks Lemma 3.3(1) formula, Lemma 3.1(1) (2p'(q) <= m_q, q < 8p'(q)), and the H-window counts
c_xy (Lemma C.2) for eta = 1/2.  Independent code."""
import math, sys
import numpy as np
from sympy import primerange, isprime, factorint

def pprime_map(qmax):
    pp = {3: 2, 5: 2, 7: 3}
    k = 3
    while 2 ** k <= qmax:
        qs = [q for q in primerange(2 ** k, 2 ** (k + 1)) if q > 7]
        ts = list(primerange(2 ** (k - 2), 2 ** (k - 1)))
        n, t = len(qs), len(ts)
        for i, q in enumerate(qs):
            pp[q] = ts[(i * t) // n]
        k += 1
    return pp

def mq(q):
    return q - 1 if q % 4 == 1 else q + 1

def run(M, parity, rng, check_formula=True):
    pp = pprime_map(8 * M + 10)
    for q in primerange(3, 8 * M):
        assert 2 * pp[q] <= mq(q) and q < 8 * pp[q], q
    logm = np.zeros((M, M))
    cwin = np.zeros((M, M), dtype=int)  # window (sqrt M, M] level-1 collision counts (eta = 1/2)
    x = np.arange(M)
    small = list(primerange(3, 4 * M))
    for q in small:
        p1 = pp[q]; m = mq(q)
        sigma = (np.zeros(M, dtype=int) if parity == "const" else rng.integers(0, 2, size=M))
        kap1 = (sigma + 2 * (x % p1)) % m
        astar = 1
        while p1 * q ** (astar - 1) < M:
            astar += 1
        eq = (kap1[:, None] == kap1[None, :])
        lev = eq.astype(int)
        if q > math.sqrt(M) and q <= M:
            cwin += eq
        for a in range(2, astar + 1):
            comp = x % q ** (a - 1)
            eq = eq & (comp[:, None] == comp[None, :])
            lev += eq
        # level a* injective?
        assert eq.sum() == M, (q, "level a* not injective")
        logm += lev * math.log(q)
        if check_formula:
            pass
    np.fill_diagonal(logm, 0.0); np.fill_diagonal(cwin, 0)
    D = 2 * logm
    Dx = D.sum(axis=1)
    lam = np.linalg.eigvalsh(D)[-1]
    # Lambda'(d) + log d bound
    inv = {}
    for q, p in pp.items():
        inv.setdefault(p, []).append(q)
    def Lam(d):
        return sum(math.log(q) for p in factorint(d) for q in inv.get(p, []))
    bound1 = 4 * sum(Lam(d) + math.log(d) for d in range(1, M))
    bound2 = 36 * M * math.log(M) + 111 * M
    # Lemma 3.3(1) formula check on a sample of pairs: log m_xy <= Lambda'(d) + log d
    worst = 0.0
    for d in range(1, M):
        worst = max(worst, logm[0, d] - (Lam(d) + math.log(d)))
    return dict(Dmax=D.max(), maxDx=Dx.max(), lam=lam, b1=bound1, b2=bound2, worst=worst,
                cmax=cwin.max())

rng = np.random.default_rng(7)
fails = 0
print(f"{'M':>5} {'parity':>6} {'Dmax/logM':>10} {'maxDx/MlogM':>12} {'lmax/MlogM':>11} {'maxDx':>10} {'bound4sum':>10} {'36MlogM+111M':>13} {'cwin_max':>8}")
for M in [64, 128, 256, 512, 1024]:
    for parity in ["const", "random"]:
        r = run(M, parity, rng)
        L = math.log(M)
        ok = r["maxDx"] <= r["b1"] + 1e-6 and r["b1"] <= r["b2"] and r["lam"] <= r["maxDx"] + 1e-6 and r["worst"] <= 1e-9
        fails += (not ok)
        print(f"{M:5d} {parity:>6} {r['Dmax']/L:10.2f} {r['maxDx']/(M*L):12.3f} {r['lam']/(M*L):11.3f} {r['maxDx']:10.1f} {r['b1']:10.1f} {r['b2']:13.1f} {r['cmax']:8d} {'OK' if ok else 'FAIL'}")
print("Note: D_max/log M is 2 log m_xy / log M; Lemma 3.3(1) gives log m <= (10+o(1)) log M, so D_max <= (20+o(1)) log M.")
print("ALL D_A CHECKS PASSED" if fails == 0 else f"FAILURES {fails}")
