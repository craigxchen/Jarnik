#!/usr/bin/env python3
"""Referee check of L5.md Theorem 4.1 and Corollary 4.2(1) with independent code.

(1) Paley-I twin profile of sharp.md section 4 built from scratch (normalised Paley-I Hadamard matrix,
    b copies of each nonconstant column, one flipped entry per row in a distinct column), exact
    det Gram_A (Bareiss, integers), lambda_min(Gram_A) vs (b-24)M/4, and max |E_xy| vs 6.
(2) Inequality (4.1) on random and on adversarially perturbed integer matrices G (injective), with
    b_x from random 0/1 profiles; reports the minimum slack found.
"""
import random, math
import numpy as np
random.seed(3)

def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))

def paley1(q):
    chi = [0] * q
    for x in range(1, q): chi[(x * x) % q] = 1
    chi = [0 if a == 0 else (1 if chi[a] else -1) for a in range(q)]
    M = q + 1
    S = [[0] * M for _ in range(M)]
    for j in range(1, M): S[0][j] = 1; S[j][0] = -1
    for i in range(q):
        for j in range(q): S[1 + i][1 + j] = chi[(j - i) % q]
    H = [[S[i][j] + (1 if i == j else 0) for j in range(M)] for i in range(M)]
    # normalise: first row and first column all +1
    for j in range(M):
        if H[0][j] == -1:
            for i in range(M): H[i][j] = -H[i][j]
    for i in range(M):
        if H[i][0] == -1:
            H[i] = [-v for v in H[i]]
    for i in range(M):
        for k in range(M):
            assert sum(H[i][c] * H[k][c] for c in range(M)) == (M if i == k else 0)
    return H

def bareiss_det(Mx):
    a = [row[:] for row in Mx]; n = len(a); sign = 1; prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if a[i][k] != 0), None)
            if sw is None: return 0
            a[k], a[sw] = a[sw], a[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) // prev
        prev = a[k][k]
    return sign * a[n - 1][n - 1]

def twin_profile(q0, b):
    H = paley1(q0); M = q0 + 1
    cols = [(a, c) for a in range(q0) for c in range(b)]     # label a = column 1+a of H
    ci = {cc: k for k, cc in enumerate(cols)}
    S = [[H[x][1 + a] for (a, c) in cols] for x in range(M)]
    flips = []
    for x in range(M):           # row 0 = infinity, row 1+u = u
        fc = (0, 2) if x == 0 else (x - 1, 0)
        S[x][ci[fc]] = -S[x][ci[fc]]
        flips.append(ci[fc])
    assert len(set(flips)) == M
    A = [[(1 + v) // 2 for v in row] for row in S]
    return A, S, H, M, len(cols)

print('(1) Paley twin profile: exact det Gram_A, lambda_min, |E|')
for (q0, b) in [(43, 33), (47, 33), (59, 33), (43, 25), (47, 26), (11, 30)]:
    if not (is_prime(q0) and q0 % 4 == 3): continue
    A, S, H, M, r = twin_profile(q0, b)
    x0 = 0
    B = [[A[x][j] - A[x0][j] for j in range(r)] for x in range(M) if x != x0]
    G = [[sum(u * v for u, v in zip(bx, by)) for by in B] for bx in B]
    d = bareiss_det(G)
    lam = np.linalg.eigvalsh(np.array(G, dtype=float)).min()
    E = [[4 * G[i][j] - b * M * ((1 if i == j else 0) + 1) for j in range(M - 1)] for i in range(M - 1)]
    maxE = max(abs(v) for row in E for v in row) / 4
    lhs = 0.5 * math.log(d)
    print('  q0=%2d b=%2d M=%2d r=%4d: (1/2)log det=%.4f, ((M-1)/2)log(bM/4)+(1/2)log M=%.4f, '
          'bound ((M-1)/2)log((b-24)M/4)=%.4f; lambda_min=%.2f vs (b-24)M/4=%.2f; max|E|=%.2f (<=6)'
          % (q0, b, M, r, lhs, (M - 1) / 2 * math.log(b * M / 4) + 0.5 * math.log(M),
             (M - 1) / 2 * math.log(max(b - 24, 1e-9) * M / 4), lam, (b - 24) * M / 4, maxE))

print('\n(2) inequality (4.1): sum_x log||G b_x|| >= (1/2)log det Gram - (1/2)log C(r,M-1) - (r-M+1)log B')
def slack(G, Bv):
    G = np.array(G, dtype=float); Bv = np.array(Bv, dtype=float)
    r = G.shape[1]; Mm1 = Bv.shape[0]
    if np.linalg.matrix_rank(G) < r: return 1e9, False, True
    lhs = sum(math.log(np.linalg.norm(G @ bx)) for bx in Bv)
    Gram = Bv @ Bv.T
    sgn, ld = np.linalg.slogdet(Gram)
    Bmax = max(np.linalg.norm(G[:, j]) for j in range(r))
    rhs = 0.5 * ld - 0.5 * (math.lgamma(r + 1) - math.lgamma(Mm1 + 1) - math.lgamma(r - Mm1 + 1)) \
          - (r - Mm1) * math.log(Bmax)
    return lhs - rhs, np.linalg.matrix_rank(G) == r, sgn > 0

minsl = 1e9; tested = 0
for trial in range(600):
    M = random.randint(3, 7); r = random.randint(M, 12); P = random.randint(r, r + 4)
    A = [[random.randint(0, 1) for _ in range(r)] for _ in range(M)]
    Bv = [[A[x][j] - A[0][j] for j in range(r)] for x in range(1, M)]
    if np.linalg.matrix_rank(np.array(Bv, dtype=float)) < M - 1: continue
    mode = trial % 3
    if mode == 0:
        G = [[random.randint(-2, 2) for _ in range(r)] for _ in range(P)]
    elif mode == 1:   # nearly dependent columns: G = I + rank-one integer perturbation
        u = [random.randint(-3, 3) for _ in range(P)]; v = [random.randint(-3, 3) for _ in range(r)]
        G = [[(1 if i == j else 0) + u[i] * v[j] for j in range(r)] for i in range(P)]
    else:             # adversarial: G maps the b_x to short vectors (columns combine to cancel)
        G = [[0] * r for _ in range(P)]
        for i in range(min(P, r)): G[i][i] = 1
        for i in range(P):
            for j in range(r): G[i][j] += random.choice([0, 0, 0, 1, -1])
    sl, inj, pos = slack(G, Bv)
    if not inj or not pos: continue
    tested += 1
    # local search to minimise slack
    for it in range(60):
        i = random.randrange(P); j = random.randrange(r)
        G2 = [row[:] for row in G]; G2[i][j] += random.choice([1, -1])
        s2, inj2, _ = slack(G2, Bv)
        if inj2 and s2 < sl: G, sl = G2, s2
    minsl = min(minsl, sl)
print('  tested %d injective instances (with 60-step slack-minimising search each); minimum slack %.4f (must be >= 0)'
      % (tested, minsl))
