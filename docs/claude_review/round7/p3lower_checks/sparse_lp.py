"""Sparse cut LP f_D(r) (Section D): least total weight of a cut system with at most r cuts
whose margins dominate D (real weights).

Why it matters: columns are distinct primes p_j = 1 mod 4 (> X in the round-5/6 classes), so r
columns cost W >= sum_j log p_j >= r log X, and also W >= (1+o(1)) r log r.  The real-weight
optimum of Theorem A typically has support binom(M,2) (an LP vertex) with weights ~ W/binom(M,2),
far below log 5.  So the pair channel with prime weights is governed by
    min_r max( f_D(r), r log X ).
Part 1 (M = 8, exact MILP): f_D(r) for r = 7..14, with the valid bound lam_s <= U = 2(M-1)D_max
        (a Sylvester Hadamard matrix of order 8 gives a feasible 7-cut solution of that cost).
Part 2 (M = 10, 12, 16; heuristic UPPER bounds): prune the LP support one cut at a time
        (remove the smallest weight, re-solve the LP on the remaining support).
Demand model: random independent residues at odd prime powers <= X = 3M.
Evidence only (floating point LP/MILP).
"""
import math, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds, linprog

rng = np.random.default_rng(2718)

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def cuts(M):
    bits = np.arange(1, 2 ** (M - 1))
    S = np.ones((len(bits), M), dtype=np.int8)
    for i in range(M - 1):
        S[:, i + 1] = np.where((bits >> i) & 1, -1, 1)
    return S

def demand_random(M, X):
    D = np.zeros((M, M))
    for q in primes_upto(X):
        if q == 2:
            continue
        m = q - 1 if q % 4 == 1 else q + 1
        cls = rng.integers(0, m, size=M)
        Q = q
        while Q <= X:
            D += 2 * math.log(q) * (cls[:, None] == cls[None, :])
            cls = cls * q + rng.integers(0, q, size=M)
            Q *= q
    np.fill_diagonal(D, 0)
    return D

def setup(D, S):
    M = D.shape[0]
    iu = np.triu_indices(M, 1)
    A = 0.5 * (S[:, iu[0]] != S[:, iu[1]]).T.astype(float) - 0.25 * 2 * (S[:, iu[0]] == S[:, iu[1]]).T.astype(float)
    # A[(xy), s] = +1/2 if separated, -1/2 otherwise
    b = D[iu]
    return A, b

def lp(A, b, cols):
    res = linprog(np.ones(len(cols)), A_ub=-A[:, cols], b_ub=-b, bounds=(0, None), method='highs')
    if res.status != 0:
        return None, None
    return res.fun, res.x

def f_milp(A, b, r, U, tlim):
    n = A.shape[1]
    c = np.concatenate([np.ones(n), np.zeros(n)])
    cons = [LinearConstraint(np.hstack([A, np.zeros_like(A)]), lb=b, ub=np.inf),
            LinearConstraint(np.hstack([np.eye(n), -U * np.eye(n)]), lb=-np.inf, ub=0),
            LinearConstraint(np.concatenate([np.zeros(n), np.ones(n)])[None, :], lb=0, ub=r)]
    integrality = np.concatenate([np.zeros(n), np.ones(n)])
    bounds = Bounds(np.zeros(2 * n), np.concatenate([np.full(n, U), np.ones(n)]))
    res = milp(c, constraints=cons, integrality=integrality, bounds=bounds,
               options=dict(time_limit=tlim, mip_rel_gap=1e-6, disp=False))
    return (None if res.x is None else res.fun), res.status, getattr(res, 'mip_dual_bound', None)

def prune(A, b):
    cols = np.arange(A.shape[1])
    W, x = lp(A, b, cols)
    keep = cols[x > 1e-9]
    W, x = lp(A, b, keep)
    best = {len(keep): W}
    while True:
        order = np.argsort(x)
        progressed = False
        for idx in order[:12]:
            trial = np.delete(keep, idx)
            W2, x2 = lp(A, b, trial)
            if W2 is not None:
                k2 = trial[x2 > 1e-9]
                W3, x3 = lp(A, b, k2)
                keep, x, W = k2, x3, W3
                best[len(keep)] = min(best.get(len(keep), math.inf), W)
                progressed = True
                break
        if not progressed:
            break
    # monotone envelope: f(r) <= min_{r' <= r} best[r']
    return best

def main():
    import sys as _s
    part = _s.argv[1] if len(_s.argv) > 1 else 'both'
    if part in ('both', '1'):
        part1()
    if part in ('both', '2'):
        part2()
    print("done (evidence only)")

def part1():
    print("Part 1: exact MILP, M = 8")
    M = 8
    S = cuts(M)
    for t in range(2):
        D = demand_random(M, 3 * M)
        A, b = setup(D, S)
        Wc, x = lp(A, b, np.arange(len(S)))
        lam = max(np.linalg.eigvalsh(D)); unif = 2 * (M - 1) * D.max()
        print(f" trial {t}: 2lam={2*lam:.2f} W_cut={Wc:.2f} (support {(x>1e-9).sum()}) uniform={unif:.2f}")
        for r in (7, 8, 10, 14):
            t0 = time.time()
            f, st, lb = f_milp(A, b, r, unif, 240)
            tag = "optimal" if st == 0 else f"time limit; dual bound {lb}"
            print(f"   r={r:2d}: f={f:8.2f}  f/W_cut={f/Wc:6.3f}  f/uniform={f/unif:6.3f}  [{tag}, {time.time()-t0:.0f}s]")
            sys.stdout.flush()

def part2():
    print("Part 2: support pruning (upper bounds on f(r))")
    for M, ntr in ((10, 3), (12, 3), (16, 2)):
        S = cuts(M)
        for t in range(ntr):
            D = demand_random(M, 3 * M)
            A, b = setup(D, S)
            lam = max(np.linalg.eigvalsh(D)); unif = 2 * (M - 1) * D.max()
            t0 = time.time()
            best = prune(A, b)
            Wc = best[max(best)]
            rs = sorted(best)
            env, cur = {}, math.inf
            for r in rs:
                cur = min(cur, best[r]); env[r] = cur
            def at(r):
                ok = [rr for rr in rs if rr <= r]
                return env[max(ok)] if ok else None
            pts = []
            for r in (M - 1, M, M + 2, 2 * M, 3 * M, 4 * M):
                v = at(r)
                pts.append(f"r={r}:{'-' if v is None else f'{v/Wc:.2f}'}")
            print(f" M={M} trial {t}: 2lam={2*lam:.1f} W_cut={Wc:.1f} (support {max(best)}) uniform={unif:.1f} "
                  f"[unif/W_cut={unif/Wc:.2f}] min support reached {min(best)} (f/W_cut={best[min(best)]/Wc:.2f}); "
                  f"f(r)/W_cut: {' '.join(pts)} [{time.time()-t0:.0f}s]")
            sys.stdout.flush()

if __name__ == '__main__':
    main()
