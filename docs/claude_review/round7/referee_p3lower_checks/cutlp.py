#!/usr/bin/env python3
"""Referee check of Theorem A (p3lower.md): exact cut LP W_cut(B), sandwich
2 lambda_max(B) <= W_cut(B) <= pi lambda_max(B) for B >= 0, LP duality,
Perron/hypermetric dual certificates; adversarial search for the ratio.
Independent code (does not import the author's scripts)."""
import itertools, math, sys
import numpy as np
from scipy.optimize import linprog

rng = np.random.default_rng(12345)

def cuts(M):
    # nonconstant s in {+-1}^M modulo s ~ -s : fix s_0 = +1, exclude all-ones
    out = []
    for bits in itertools.product([1, -1], repeat=M - 1):
        s = (1,) + bits
        if all(b == 1 for b in bits):
            continue
        out.append(np.array(s))
    return out

_cache = {}
def cutmat(M):
    if M in _cache: return _cache[M]
    S = cuts(M)
    pairs = [(x, y) for x in range(M) for y in range(x + 1, M)]
    A = np.array([[0.5 if s[x] != s[y] else -0.5 for s in S] for (x, y) in pairs])  # rows pairs, cols cuts
    _cache[M] = (S, pairs, A)
    return _cache[M]

def wcut(B):
    M = B.shape[0]
    S, pairs, A = cutmat(M)
    b = np.array([B[x, y] for (x, y) in pairs])
    # primal: min 1.lam  s.t. A lam >= b, lam >= 0
    res = linprog(np.ones(len(S)), A_ub=-A, b_ub=-b, bounds=(0, None), method="highs")
    assert res.status == 0, res.message
    # dual: max b.nu s.t. A^T nu <= 1, nu >= 0
    resd = linprog(-b, A_ub=A.T, b_ub=np.ones(len(S)), bounds=(0, None), method="highs")
    assert resd.status == 0, resd.message
    return res.fun, -resd.fun, res.x, resd.x

def lmax(B):
    return float(np.linalg.eigvalsh(B)[-1])

def rand_sym(M, gen):
    B = np.zeros((M, M))
    for x in range(M):
        for y in range(x + 1, M):
            B[x, y] = B[y, x] = gen()
    return B

def residue_demand(M, X, rng):
    # independent uniform classes at odd prime powers Q <= X; D_xy = 2 sum log q [collision mod Q]
    D = np.zeros((M, M))
    for q in range(3, X + 1, 2):
        if any(q % p == 0 for p in range(3, int(q ** 0.5) + 1, 2)):
            continue
        chi = 1 if q % 4 == 1 else -1
        Q = q; kap = None
        while Q <= X:
            m = (q - chi) * (Q // q)
            if kap is None:
                kap = rng.integers(0, m, size=M)
            else:
                kap = kap + (m // q) * rng.integers(0, q, size=M)  # random lift
            for x in range(M):
                for y in range(x + 1, M):
                    if kap[x] == kap[y]:
                        D[x, y] += 2 * math.log(q); D[y, x] = D[x, y]
            Q *= q
    return D

fails = 0
ratios = []
print("== sandwich on test families ==")
tests = []
for M in [4, 5, 6, 7, 8, 9]:
    tests.append((f"uniform M={M}", np.ones((M, M)) - np.eye(M)))
    H = np.zeros((M, M)); H[0, 1] = H[1, 0] = 10.0
    tests.append((f"heavy pair M={M}", H + 0.5 * (np.ones((M, M)) - np.eye(M))))
    tests.append((f"rand exp M={M}", rand_sym(M, lambda: rng.exponential())))
    tests.append((f"rand sparse M={M}", rand_sym(M, lambda: rng.exponential() * (rng.random() < 0.3))))
    tests.append((f"residues X=3M M={M}", residue_demand(M, 3 * M, rng)))
    tests.append((f"residues X=M^2 M={M}", residue_demand(M, M * M, rng)))
for name, B in tests:
    if not B.any():
        continue
    w, wd, lam, nu = wcut(B)
    L = lmax(B)
    r = w / (2 * L)
    ratios.append(r)
    ok = (2 * L <= w + 1e-7) and (w <= math.pi * L + 1e-7) and abs(w - wd) < 1e-6 * max(1, w)
    # cut identity at optimum: sum_{x<y} mu = M W/4 - (1/4) sum lam sigma^2
    M = B.shape[0]; S, pairs, A = cutmat(M)
    mu = A @ lam
    lhs = mu.sum(); rhs = M * lam.sum() / 4 - 0.25 * sum(l * (s.sum()) ** 2 for l, s in zip(lam, S))
    ok = ok and abs(lhs - rhs) < 1e-6 * max(1, abs(lhs))
    supp = int((lam > 1e-9).sum())
    if not ok:
        fails += 1
    print(f"{name:24s} Wcut={w:9.4f} dual={wd:9.4f} 2lmax={2*L:9.4f} pi*lmax={math.pi*L:9.4f} ratio={r:.4f} supp={supp} (binom={M*(M-1)//2}) {'OK' if ok else 'FAIL'}")

print("ratio W_cut/(2 lmax) range on test families: [%.4f, %.4f]" % (min(ratios), max(ratios)))

# Perron and hypermetric dual certificates are feasible: check |u|^2-(u.s)^2 <= 1 and (b.s)^2>=1
print("== dual-certificate feasibility (exact over all cuts) ==")
for M in [5, 6, 7]:
    S, pairs, A = cutmat(M)
    for trial in range(200):
        u = rng.random(M); u /= np.linalg.norm(u)
        nu = np.array([4 * u[x] * u[y] for (x, y) in pairs])
        if (A.T @ nu).max() > 1 + 1e-12: fails += 1; print("Perron cert infeasible", M)
        b = rng.integers(0, 4, size=M)
        if b.sum() % 2 == 1 and (b * b).sum() >= 2:
            nub = np.array([4 * b[x] * b[y] / ((b * b).sum() - 1) for (x, y) in pairs])
            if (A.T @ nub).max() > 1 + 1e-12: fails += 1; print("hypermetric cert infeasible", M, b)
print("dual certificates feasible" if fails == 0 else f"FAILS={fails}")

# adversarial search: maximise W_cut/(2 lmax) over nonnegative B, M=5..7 (hill climbing)
print("== adversarial search for max W_cut/(2 lmax) (evidence) ==")
for M in [5, 6, 7]:
    best = 0; bestB = None
    for start in range(6):
        B = rand_sym(M, lambda: rng.exponential() * (rng.random() < 0.5))
        if not B.any(): continue
        cur = wcut(B)[0] / (2 * lmax(B))
        for it in range(150):
            B2 = B.copy(); x, y = rng.choice(M, 2, replace=False)
            B2[x, y] = B2[y, x] = max(0.0, B2[x, y] + rng.normal() * 0.5)
            if not B2.any(): continue
            r2 = wcut(B2)[0] / (2 * lmax(B2))
            if r2 > cur: B, cur = B2, r2
        if cur > best: best, bestB = cur, B
    print(f"M={M}: best ratio found {best:.4f}  (pi/2 = {math.pi/2:.4f})")
    if best > math.pi / 2 + 1e-7: fails += 1
print("ALL THEOREM-A CHECKS PASSED" if fails == 0 else f"FAILURES: {fails}")
