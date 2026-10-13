#!/usr/bin/env python3
"""Checks for L5.md Theorem 4.1 (volume form of the certification tension).

(1) Builds the sharp.md section 4 Paley twin profile (normalised Paley-I Hadamard [1|H'], b columns per label,
    canonical flips/twins), verifies the twin relations and the Hadamard property, and computes EXACTLY
    log det Gram(a_x - a_x0) (Bareiss), compared with ((M-1)/2) log(bM/4) + (1/2) log M.
(2) Tests the singular-value inequality behind Theorem 4.1 on random integer matrices G with full column
    rank:  prod_x ||G b_x|| >= sqrt(det Gram_b) * binom(r, d)^(-1/2) * B^(-(r-d))  with d = rank of the b_x,
    B = max column norm.  (Floating point; the inequality has large slack, reported.)
Run: python3 check_volume.py  (output check_volume_out.txt)
"""
import math, random
from fractions import Fraction
random.seed(7)

def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def paley_profile(q0, b):
    M = q0 + 1
    # normalised Paley I: rows inf, 0..q0-1 ; columns: [1 | H'], H' indexed by a in F_q0
    rows = ['inf'] + list(range(q0))
    Hp = {}
    for x in rows:
        for a in range(q0):
            if x == 'inf':
                Hp[(x, a)] = 1
            else:
                Hp[(x, a)] = -(legendre(a - x, q0) + (1 if a == x else 0))
    # Hadamard check of [1|H']
    full = [[1] + [Hp[(x, a)] for a in range(q0)] for x in rows]
    for i in range(M):
        for k in range(M):
            ip = sum(full[i][c] * full[k][c] for c in range(M))
            assert ip == (M if i == k else 0), 'not Hadamard'
    astar = 0
    alpha = {x: (x if x != 'inf' else astar) for x in rows}
    cols = [(a, i) for a in range(q0) for i in range(b)]
    fcol = {x: ((x, 0) if x != 'inf' else (astar, 2)) for x in rows}
    scol = {x: ((x, 1) if x != 'inf' else (astar, 3)) for x in rows}
    S = {}
    for x in rows:
        for (a, i) in cols:
            S[(x, (a, i))] = Hp[(x, a)]
    for x in rows:
        S[(x, fcol[x])] = -S[(x, fcol[x])]
    A = [[(1 + S[(x, c)]) // 2 for c in cols] for x in rows]
    # twin check
    ci = {c: k for k, c in enumerate(cols)}
    for xi, x in enumerate(rows):
        fj, sj = ci[fcol[x]], ci[scol[x]]
        diff = [A[y][fj] - A[y][sj] for y in range(M)]
        assert all(diff[y] == 0 for y in range(M) if y != xi) and abs(diff[xi]) == 1
    return A

def bareiss_det(Mx):
    n = len(Mx)
    a = [row[:] for row in Mx]
    sign = 1; prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if a[i][k] != 0), None)
            if sw is None:
                return 0
            a[k], a[sw] = a[sw], a[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) // prev
        prev = a[k][k]
    return sign * a[n - 1][n - 1]

def part1():
    out = []
    for q0, b in [(43, 33), (47, 33), (59, 33), (43, 4), (47, 8)]:
        A = paley_profile(q0, b)
        M = len(A); r = len(A[0])
        x0 = 0
        B = [[A[x][j] - A[x0][j] for j in range(r)] for x in range(1, M)]
        G = [[sum(u * v for u, v in zip(B[i], B[k])) for k in range(M - 1)] for i in range(M - 1)]
        d = bareiss_det(G)
        ld = math.log(d)
        pred = (M - 1) / 2 * math.log(b * M / 4) + 0.5 * math.log(M)
        lb = 0.5 * math.lgamma(r + 1) - 0.5 * math.lgamma(M) - 0.5 * math.lgamma(r - M + 2)  # (1/2) log binom(r, M-1)
        out.append('q0=%d b=%d M=%d r=%d: (1/2)log det Gram = %.2f ; ((M-1)/2)log(bM/4)+(1/2)log M = %.2f ; '
                   '(1/2)log binom(r,M-1) = %.2f ; margin = %.2f' % (q0, b, M, r, ld / 2, pred / 2 * 2 / 2 if False else pred / 1, lb, ld / 2 - lb))
    return out

def part2():
    import numpy as np
    out = []; worst = 1e9
    for trial in range(200):
        r = random.randint(4, 9); d = random.randint(1, r - 1); P = r + random.randint(0, 3)
        while True:
            G = np.array([[random.randint(-2, 2) for _ in range(r)] for _ in range(P)], dtype=float)
            if np.linalg.matrix_rank(G) == r:
                break
        Bv = np.array([[random.randint(-1, 1) for _ in range(r)] for _ in range(d)], dtype=float)
        if np.linalg.matrix_rank(Bv) < d:
            continue
        gram = Bv @ Bv.T
        lhs = sum(math.log(np.linalg.norm(G @ Bv[i])) for i in range(d))
        Bn = max(np.linalg.norm(G[:, j]) for j in range(r))
        rhs = 0.5 * math.log(np.linalg.det(gram)) - 0.5 * math.log(math.comb(r, d)) - (r - d) * math.log(Bn)
        worst = min(worst, lhs - rhs)
        if lhs < rhs - 1e-9:
            out.append('VIOLATION trial %d' % trial)
    out.append('part2: 200 random instances, min slack (log) = %.3f (must be >= 0)' % worst)
    return out

if __name__ == '__main__':
    res = part1() + part2()
    txt = '\n'.join(res)
    print(txt)
    open(__file__.replace('.py', '_out.txt'), 'w').write(txt + '\n')
