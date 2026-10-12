"""Equal-weight discrete pair problem (Question Q, Section 5), solved exactly by integer programming.

Columns are distinct primes in [P, P(1 + 1/r)], so every column has weight tau (+ O(1/r)), and a
column is a cut s (cuts may repeat: distinct primes may share a sign pattern).  With n_s columns on
cut s, the pair margin is tau * sum_s n_s ([s_x != s_y] - 1/2).  The least number of columns is
    r*(D, tau) = min { sum n_s : n in Z_{>=0}^{cuts},  sum_s n_s ([s_x != s_y] - 1/2) >= D_xy / tau }.
We compare r* with:
    Perron lower bound      2 lam_max(D)/tau            (Theorem A(2), any weights),
    real LP                  W_cut(D)/tau,
    uniform margins          r_unif = min{ sum n_s : margins all >= D_max/tau }  (also an IP),
for demand matrices from the independent random residue model at X = 3M, and tau = kappa log M.
Exact MILP (HiGHS, integer variables, gap 0) for M = 8, 10; M = 12 with a time limit (reported).
Evidence only.
"""
import math, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds, linprog

rng = np.random.default_rng(99)

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

def Amat(S, M):
    iu = np.triu_indices(M, 1)
    return np.where(S[:, iu[0]] != S[:, iu[1]], 0.5, -0.5).T, iu

def ip(A, b, tlim):
    n = A.shape[1]
    res = milp(np.ones(n), constraints=[LinearConstraint(A, lb=b, ub=np.inf)],
               integrality=np.ones(n), bounds=Bounds(np.zeros(n), np.full(n, np.inf)),
               options=dict(time_limit=tlim, mip_rel_gap=0.0, disp=False))
    lb = getattr(res, 'mip_dual_bound', None)
    return (None if res.x is None else round(res.fun)), res.status, lb, res.x

def main():
    import sys as _s
    plan = ((8, (1.0, 2.0), 3, 120), (10, (1.0, 2.0), 2, 600)) if len(_s.argv) < 2 else ((12, (2.0,), 1, int(_s.argv[1])),)
    for M, kappas, ntr, tlim in plan:
        S = cuts(M)
        A, iu = Amat(S, M)
        for kappa in kappas:
            for t in range(ntr):
                D = demand_random(M, 3 * M)
                tau = kappa * math.log(M)
                b = D[iu] / tau
                lam = max(np.linalg.eigvalsh(D))
                lpres = linprog(np.ones(A.shape[1]), A_ub=-A, b_ub=-b, bounds=(0, None), method='highs')
                t0 = time.time()
                r_star, st, lb, x = ip(A, b, tlim)
                r_unif, st2, lb2, _ = ip(A, np.full_like(b, b.max()), tlim)
                kk = np.maximum(1, np.ceil(2 * D / tau - 1e-12)); np.fill_diagonal(kk, 0)
                tag = "opt" if st == 0 else f"time limit, lower bound {lb:.1f}"
                tag2 = "opt" if st2 == 0 else f"time limit, lower bound {lb2:.1f}"
                print(f"M={M:2d} tau={kappa}logM #{t}: D_max/tau={D.max()/tau:5.2f} Perron 2lam/tau={2*lam/tau:6.2f} "
                      f"LP={lpres.fun:6.2f} r*={r_star} [{tag}] r_unif={r_unif} [{tag2}] "
                      f"r*/LP={r_star/lpres.fun:5.2f} r_unif/r*={r_unif/r_star:5.2f} distinct cuts used={(x>0.5).sum()} "
                      f"[{time.time()-t0:.0f}s]")
                sys.stdout.flush()
    print("done (evidence only)")

if __name__ == '__main__':
    main()
