"""Additive (Walsh-type) structure of actual clusters versus random subsets of the same sign cube.
For a cluster X of sign vectors on k primes (weights log p), over the VARYING primes:
  r(Q) = (W - T_Q)/W = 2 * weight{j : s1j s2j s3j s4j = -1} / W   (relative quartic defect),
  exact parallelogram  <=> r(Q) = 0 (Hadamard-4 core / affine 2-cube),
  flip-like quadruple  <=> odd support has <= 4 primes (Walsh 2-frame with flips).
Also the pair slack d_xy - logN/2 + 2 log C (>= 0 is the necessary pair condition)."""
import json, math, random, sys, os, itertools, glob
import numpy as np

def stats(S, w):
    M, k = S.shape
    var = [j for j in range(k) if len(set(S[:, j])) > 1]
    if not var: return None
    Sv = S[:, var]; wv = w[var]; W = wv.sum()
    rmin = 9; rsum = 0; npar = 0; nflip = 0; nq = 0
    for Qd in itertools.combinations(range(M), 4):
        prod = Sv[Qd[0]] * Sv[Qd[1]] * Sv[Qd[2]] * Sv[Qd[3]]
        odd = prod < 0
        r = 2 * wv[odd].sum() / W
        rmin = min(rmin, r); rsum += r; nq += 1
        if odd.sum() == 0: npar += 1
        if odd.sum() <= 4: nflip += 1
    return dict(rmin=rmin, rmean=rsum / nq, npar=npar, nflip=nflip, nq=nq, kvar=len(var), W=W)

def pair_slack(S, w, logN, C):
    M = S.shape[0]
    out = []
    for x in range(M):
        for y in range(x + 1, M):
            d = w[S[x] != S[y]].sum()
            out.append(d - logN / 2 + 2 * math.log(C))
    return min(out), float(np.mean(out))

rows = []
files = sorted(glob.glob("large_clusters_*.json")) + sorted(glob.glob("best_clusters_*.json"))
random.seed(5)
for fn in files:
    data = json.load(open(fn))
    for Mk, rec in data.items():
        M = int(Mk)
        if M < 4 or M > 16: continue
        P = rec["P"]; k = len(P)
        S = np.array([[(-1) ** b for b in s] for s in rec["signs"]])
        w = np.array([math.log(p) for p in P])
        logN = w.sum()
        st = stats(S, w)
        ps = pair_slack(S, w, logN, rec["C"])
        # random baseline: random M-subsets of {+-1}^k
        base = []
        for t in range(60):
            Sr = np.array([[random.choice((1, -1)) for _ in range(k)] for _ in range(M)])
            sr = stats(Sr, w)
            if sr: base.append(sr)
        bmin = np.mean([b["rmin"] for b in base]); bmean = np.mean([b["rmean"] for b in base])
        bpar = np.mean([b["npar"] for b in base]); bflip = np.mean([b["nflip"] for b in base])
        rows.append((fn, M, rec["C"], rec.get("lambda", float("nan")), k, st, ps, bmin, bmean, bpar, bflip))
        print("%-28s M=%2d C=%8.2f lam=%5.2f k=%2d | rmin=%.3f rmean=%.3f par=%d flip<=4=%d (of %d) | base rmin=%.3f rmean=%.3f par=%.2f flip=%.2f | pairslack min=%.2f mean=%.2f"
              % (fn[:28], M, rec["C"], rec.get("lambda", float("nan")), k, st["rmin"], st["rmean"], st["npar"], st["nflip"], st["nq"], bmin, bmean, bpar, bflip, ps[0], ps[1]))
