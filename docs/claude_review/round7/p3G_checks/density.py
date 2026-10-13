"""Evidence: tight density of actual Gaussian-prime relations vs the 'random odd integer' model.

Model prediction for a relation with value D: dens(|D|) = #{d | L_X : d in [|D|/2, 2|D|]} / #{odd n in [|D|/2, 2|D|]}.
Divisors of L_X are enumerated exactly (DFS bounded by Ymax).
"""
import random, bisect, sys, math
import numpy as np
from gtight import primes_upto, two_squares, Ds, smooth_mask
from trend import rand_primes

def divisors_upto(X, Ymax):
    pp = []
    for q in primes_upto(X):
        if q == 2: continue
        Q = q; pw = []
        while Q <= X:
            pw.append(Q); Q *= q
        pp.append(pw)
    out = []
    def dfs(k, d):
        if k == len(pp):
            out.append(d); return
        dfs(k + 1, d)
        for Q in pp[k]:
            if d * Q > Ymax: break
            dfs(k + 1, d * Q)
    dfs(0, 1)
    out.sort()
    return out

def run(M, K, nsample=200000, seed=0):
    rng = random.Random(seed)
    X = 2 * M; P = int(M ** K)
    r = 400
    ps = rand_primes(r, P, 2 * P, rng)
    gp = [two_squares(p) for p in ps]
    gp = [(a, b if rng.random() < .5 else -b) for a, b in gp]
    a = np.array([g[0] for g in gp], dtype=np.int64); b = np.array([g[1] for g in gp], dtype=np.int64)
    divs = divisors_upto(X, 8 * P)
    def dens(Y):
        lo, hi = Y / 2, 2 * Y
        c = bisect.bisect_right(divs, hi) - bisect.bisect_left(divs, lo)
        return c / ((hi - lo) / 2)
    res = {}
    # singles
    for s in range(4):
        D = Ds(a, b, s)
        obs = smooth_mask(D, X).sum()
        pred = sum(dens(abs(int(x))) for x in D)
        res[('S', s)] = (len(D), int(obs), pred)
    # pairs (random sample)
    I = np.array([rng.randrange(r) for _ in range(nsample)]); J = np.array([rng.randrange(r) for _ in range(nsample)])
    m = I != J; I, J = I[m], J[m]
    for tau in (1, -1):
        Xp = a[I] * a[J] - tau * b[I] * b[J]
        Yp = a[I] * b[J] * tau + b[I] * a[J]
        for s in range(4):
            D = Ds(Xp, Yp, s)
            obs = smooth_mask(D, X).sum()
            Dl = np.abs(D).astype(np.float64)
            # vectorised density via bisect on divs
            dv = np.array(divs, dtype=np.float64)
            c = np.searchsorted(dv, 2 * Dl, side='right') - np.searchsorted(dv, Dl / 2, side='left')
            pred = float((c / (0.75 * Dl)).sum())
            res[('P', tau, s)] = (len(D), int(obs), pred)
    return res

if __name__ == '__main__':
    for M, K in [(44, 3), (44, 4), (84, 3), (84, 4), (164, 3), (164, 4)]:
        res = run(M, K)
        tS = [sum(v[i] for k, v in res.items() if k[0] == 'S') for i in range(3)]
        tP = [sum(v[i] for k, v in res.items() if k[0] == 'P') for i in range(3)]
        print(f'M={M} K={K} X={2*M}: singles obs/pred = {tS[1]}/{tS[2]:.1f} of {tS[0]} (rate {tS[1]/tS[0]:.4f}); '
              f'pairs obs/pred = {tP[1]}/{tP[2]:.1f} of {tP[0]} (rate {tP[1]/tP[0]:.5f})', flush=True)
        for k in sorted(res):
            n, o, p = res[k]
            print('   ', k, n, o, f'{p:.1f}')
