"""Theorem A checks (evidence + exact LP values).

For a nonnegative symmetric demand matrix D (zero diagonal) on M points, the cut LP is
    W_cut(D) = min sum_s lam_s  s.t.  sum_s lam_s ([s_x != s_y] - 1/2) >= D_xy  for all x<y,
over nonconstant cuts s (mod +-).  Theorem A: 2 lam_max(D) <= W_cut(D) <= pi lam_max(D),
and W_cut(D) = max_nu <nu,D>/surplus(nu).

This script solves the LP exactly (HiGHS, double precision) for M <= 11 and several demand
models, checks the sandwich, checks LP duality numerically, and compares with the cost of
uniform margins, 2 lam_max(D_max (J-I)) = 2 (M-1) D_max.
Also verifies, on the LP optimum, the cut identity sum_{x<y} mu_xy = (M W - sum_s lam_s (sum s)^2)/4.
Evidence only (floating point); no proof uses it.
"""
import itertools, math, sys
import numpy as np
from scipy.optimize import linprog

rng = np.random.default_rng(20261009)

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def cuts(M):
    # nonconstant cuts mod +-: fix s_0 = +1, exclude all-ones
    out = []
    for bits in range(1, 2 ** (M - 1)):
        s = np.ones(M, dtype=int)
        for i in range(M - 1):
            if (bits >> i) & 1:
                s[i + 1] = -1
        out.append(s)
    return np.array(out)

def solve_lp(D, S):
    M = D.shape[0]
    pairs = [(x, y) for x in range(M) for y in range(x + 1, M)]
    # A[(xy), s] = [s_x != s_y] - 1/2
    A = np.array([[(0.5 if S[k, x] != S[k, y] else -0.5) for k in range(len(S))] for (x, y) in pairs])
    b = np.array([D[x, y] for (x, y) in pairs])
    res = linprog(c=np.ones(len(S)), A_ub=-A, b_ub=-b, bounds=(0, None), method='highs')
    assert res.status == 0, res.message
    dual = -res.ineqlin.marginals  # nu >= 0
    return res.fun, res.x, dual, A, b, pairs

def demand_random_collision(M, X, rng):
    """Random independent residue model: each point gets a uniform class mod m_Q at each odd
    prime power Q = q^a <= X (lifts uniform), unit-1 collision demand 2 sum log q."""
    D = np.zeros((M, M))
    for q in primes_upto(X):
        if q == 2:
            continue
        m = q - 1 if q % 4 == 1 else q + 1
        cls = rng.integers(0, m, size=M)
        a, Q = 1, q
        while Q <= X:
            same = (cls[:, None] == cls[None, :])
            D += 2 * math.log(q) * same
            # lift
            cls = cls * q + rng.integers(0, q, size=M)
            a += 1; Q *= q
    np.fill_diagonal(D, 0)
    return D

def report(name, D, S):
    M = D.shape[0]
    lam = max(np.linalg.eigvalsh(D))
    W, x, nu, A, b, pairs = solve_lp(D, S)
    # duality check: <nu, D> / max_s (cut_nu(s) - |nu|/2) == W
    surplus = max(A.T @ nu)  # per s: sum nu_xy([.]-1/2)
    dual_val = (nu @ b) / surplus if surplus > 0 else float('nan')
    Dmax = D.max()
    unif = 2 * (M - 1) * Dmax
    # cut identity on optimum
    mu = A @ x
    lhs = mu.sum()
    rhs = (M * x.sum() - sum(x[k] * (S[k].sum()) ** 2 for k in range(len(S)))) / 4
    ok = (2 * lam - 1e-7 <= W <= math.pi * lam + 1e-7) and abs(dual_val - W) < 1e-6 * max(1, W) and abs(lhs - rhs) < 1e-6 * max(1, W)
    print(f"{name:28s} M={M:2d} 2lam={2*lam:9.3f} W_cut={W:9.3f} pi*lam={math.pi*lam:9.3f} "
          f"W/2lam={W/(2*lam):6.4f} dual={dual_val:9.3f} unif={unif:9.3f} unif/W={unif/W:6.3f} "
          f"#cuts>0={int((x>1e-9).sum()):3d} {'OK' if ok else 'FAIL'}")
    return ok

def main():
    allok = True
    for M in (6, 8, 10, 11):
        S = cuts(M)
        # uniform
        D = np.ones((M, M)) - np.eye(M)
        allok &= report("uniform J-I", D, S)
        # single heavy pair
        D = np.ones((M, M)) - np.eye(M); D[0, 1] = D[1, 0] = 10.0
        allok &= report("J-I + heavy pair (10)", D, S)
        # matching pattern, L = M
        D = np.ones((M, M)) - np.eye(M)
        for i in range(0, M - 1, 2):
            D[i, i + 1] = D[i + 1, i] = float(M)
        allok &= report("J-I + matching (L=M)", D, S)
        # random collision model, X = 3M and X = M^2
        for X in (3 * M, M * M):
            for t in range(3):
                D = demand_random_collision(M, X, rng)
                if D.max() == 0:
                    continue
                allok &= report(f"random coll. X={X} #{t}", D, S)
        # iid heavy-tailed (Poisson-type collision counts at large primes)
        for t in range(3):
            cnt = rng.poisson(1.0, size=(M, M)); cnt = np.triu(cnt, 1); cnt = cnt + cnt.T
            D = 2.0 * cnt * math.log(M) + 1.0 * (1 - np.eye(M))
            allok &= report(f"Poisson-tail iid #{t}", D, S)
        sys.stdout.flush()
    print("ALL LP/GW SANDWICH AND DUALITY CHECKS PASSED" if allok else "SOME CHECK FAILED")

if __name__ == '__main__':
    main()
