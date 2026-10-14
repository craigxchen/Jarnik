"""Prop R / R' certificates (round 14).  Usage: python3 anticode_cert.py K "p1,..." 
For every q0 in F_K^ex: for anticode shapes D of l1-diameter <= K-1 (balls S_h, double balls S_h u (e_j+S_h)):
  cover  : genuine c = L(q0) mod n maps D onto Z/n            (=> Prop R: no F2 survivor; every non-genuine label
                                                                  has an extra collision in (D-D) u (D+S_1-D))
  cert   : additionally every v in ((D-D) u (D+S_1-D)) \\ Lambda_gen(q0), v != 0, is unaffordable at q0 for every unit
           (=> Prop R': q0 is one-modulus rigid; checked against the same allowed-table as deform.py)
Prints, per q0, the best shape found."""
import sys, math, itertools
import numpy as np
from rlib14 import *
from forced import forced_set
def ball(r, h):
    return [v for v in VK(r, h, lo=0)]
def shapes(r, d):
    h = d // 2
    S = ball(r, h)
    if d % 2 == 0: return [("S%d" % h, S)]
    out = []
    for j in range(r):
        e = [0] * r; e[j] = 1
        D = sorted(set(S) | {tuple(a + b for a, b in zip(v, e)) for v in S})
        out.append(("S%d+e%d" % (h, j + 1), D))
    return out
def main(K, ps, Fmax=60000):
    r = len(ps); pis = gaussian_primes(ps)
    F = forced_set(ps, K, Fmax)
    shp = []
    for d in range(1, K):  # all diameters up to K-1, largest first
        shp = shapes(r, d) + shp
    E1 = [tuple(x) for x in VK(r, 1, lo=0)]
    res = []
    for q0 in F:
        m, L = genuine_labels(ps, q0); n = m // 4
        others = [q for q in F if q != q0]
        memo = {}
        def unaff(v):
            if v in memo: return memo[v]
            z = Zv(pis, v); nz = z[0] ** 2 + z[1] ** 2; ok = True
            for s in range(4):
                D = Ds(z, s); pr = q0
                for q in others: pr *= q ** vq(D, q)
                if pr * pr <= nz * (2 if s % 2 else 1): ok = False; break
            memo[v] = ok; return ok
        best = None
        for name, D in shp:
            if len(D) < n: continue
            img = {sum(a * b for a, b in zip(v, L)) % n for v in D}
            if len(img) < n: continue
            diffs = set()
            for a in D:
                for b in D:
                    diffs.add(tuple(x - y for x, y in zip(a, b)))
                    for e in E1:
                        diffs.add(tuple(x + z - y for x, z, y in zip(a, e, b)))
            diffs.discard(tuple([0] * r))
            extra = [v for v in diffs if sum(x * y for x, y in zip(v, L)) % n != 0]
            bad = [v for v in extra if not unaff(v)]
            cand = (name, len(D), True, len(bad) == 0, len(bad))
            if best is None or (cand[3] and not best[3]): best = cand
            if cand[3]: break
        res.append((q0, n, best))
        if best is None: print(f"K={K} q0={q0} n={n}: no covering anticode", flush=True)
        else: print(f"K={K} q0={q0} n={n}: cover by {best[0]} (|D|={best[1]}); Prop R' certificate={'YES' if best[3] else 'no'} (affordable extras={best[4]})", flush=True)
    return res
if __name__ == "__main__":
    main(int(sys.argv[1]), [int(x) for x in sys.argv[2].split(",")])
