#!/usr/bin/env python3
"""Referee checks of L5.md section 5 (independent code).

(1) Prop 5.3(2) proof: 'kappa' = d_M u with u injective mod n/gcd(d_M,n) has no pair collision'.
    The pair data are kappa = A nut + 2 kappa' (mod m = 2n), not kappa'.  Search for (A, n, nut, u)
    where kappa' is injective but kappa has a unit-1 pair collision; then verify by brute force that a
    collision-free kappa in A nut + 2 H does exist when n/gcd(d_M, n) >= M (the corrected argument).
(2) Level: Prop 5.3 uses n = n_Q (top level).  'No collision at q' must hold at level 1 (mod m_q).
    Pigeonhole example: q = 3, Q = 81, n_Q = 54, M = 5, d_M = 2: n_Q/gcd = 27 >= 5, yet every kappa
    collides mod m_3 = 4.
(3) Hadamard cores (J+H)/2: own Smith normal form; first column constant (=> not a valid profile);
    the valid core profile (nonconstant columns, any number of copies) has rank M-1, so Delta_M = 0;
    the identity B^(-1) = 2 H^T/M - e_1 e_1^T (=> d_M = M/2 for EVERY normalised Hadamard H, M >= 4).
(4) Random 0/1 profiles: Delta_M = 1 frequency.
"""
import random, itertools, math
from fractions import Fraction
random.seed(5)

def snf_diag(Ain):
    """invariant factors of an integer matrix (list of nonzero d_i), plain Python ints."""
    A = [row[:] for row in Ain]
    m, n = len(A), len(A[0]); diag = []
    t = 0
    while t < min(m, n):
        # find nonzero pivot with minimal abs value in submatrix
        piv = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] != 0 and (piv is None or abs(A[i][j]) < abs(A[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None: break
        i, j = piv
        A[t], A[i] = A[i], A[t]
        for row in A: row[t], row[j] = row[j], row[t]
        while True:
            done = True
            for i in range(t + 1, m):
                if A[i][t]:
                    qq = A[i][t] // A[t][t]
                    A[i] = [a - qq * b for a, b in zip(A[i], A[t])]
                    if A[i][t]:
                        done = False
                        if abs(A[i][t]) < abs(A[t][t]): A[t], A[i] = A[i], A[t]
            for j in range(t + 1, n):
                if A[t][j]:
                    qq = A[t][j] // A[t][t]
                    for row in A: row[j] -= qq * row[t]
                    if A[t][j]:
                        done = False
                        if abs(A[t][j]) < abs(A[t][t]):
                            for row in A: row[t], row[j] = row[j], row[t]
            if done:
                # divisibility condition
                bad = None
                for i in range(t + 1, m):
                    for j in range(t + 1, n):
                        if A[i][j] % A[t][t]: bad = (i, j); break
                    if bad: break
                if bad:
                    A[t] = [a + b for a, b in zip(A[t], A[bad[0]])]
                    continue
                break
        diag.append(abs(A[t][t])); t += 1
    return diag

def Delta(A):
    d = snf_diag(A)
    if len(d) < len(A): return 0
    return math.prod(d)

# ---------------- (1) kappa versus kappa' ---------------------------------------------------------
print('(1) Prop 5.3(2): injective kappart vs collision-free kappa = A nut + 2 kappart (mod 2n)')
found_gap = None; checked = 0; corrected_ok = True
for trial in range(4000):
    M = random.randint(2, 3); r = random.randint(M, 4)
    A = [[random.randint(0, 2) for _ in range(r)] for _ in range(M)]
    d = snf_diag(A)
    if len(d) < M: continue
    dM = d[-1]
    if dM == 1: continue
    for n in range(2, 13):
        g = math.gcd(dM, n)
        if g == 1 or n // g < M: continue
        checked += 1
        H = set()
        for lp in itertools.product(range(n), repeat=r):
            H.add(tuple(sum(A[x][j] * lp[j] for j in range(r)) % n for x in range(M)))
        for nut in itertools.product((0, 1), repeat=r):
            base = [sum(A[x][j] * nut[j] for j in range(r)) for x in range(M)]
            m = 2 * n
            # literal recipe: kappart = dM*u, u injective mod n/g
            for u in itertools.permutations(range(n // g), M):
                kp = [dM * ux % n for ux in u]
                kap = [(base[x] + 2 * kp[x]) % m for x in range(M)]
                if len(set(kap)) < M and found_gap is None:
                    found_gap = (A, n, nut, kp, kap)
            # corrected claim: some kappart in H gives collision-free kappa
            if not any(len(set((base[x] + 2 * kp[x]) % m for x in range(M))) == M for kp in H):
                corrected_ok = False
print('  instances (A, n) checked: %d' % checked)
if found_gap:
    A, n, nut, kp, kap = found_gap
    print('  GAP: A=%s, n=%d, nut=%s: kappart=%s injective mod n, but kappa=%s (mod %d) has a collision'
          % (A, n, nut, kp, kap, 2 * n))
print('  corrected statement (collision-free kappa in A nut + 2H exists when n/gcd(d_M,n) >= M): %s'
      % corrected_ok)

# ---------------- (2) level 1 versus top level ----------------------------------------------------
print('\n(2) level: q = 3, X = 81 (A_q = 4), m_Q = 108, n_Q = 54, M = 5')
while True:
    A = [[random.randint(0, 1) for _ in range(7)] for _ in range(5)]
    d = snf_diag(A)
    if len(d) == 5 and d[-1] == 2 and math.prod(d) == 2: break
print('  profile with invariant factors %s: n_Q/gcd(d_M,n_Q) = %d >= M = 5' % (d, 54 // math.gcd(2, 54)))
# any kappa: mod m_3 = 4 the 5 values lie in 2 parity classes of 2 residues each -> collision
cnt = 0
for trial in range(2000):
    nut = [random.randint(0, 1) for _ in range(7)]
    lp = [random.randrange(54) for _ in range(7)]
    kap = [sum(A[x][j] * (nut[j] + 2 * lp[j]) for j in range(7)) % 108 for x in range(5)]
    if len(set(k % 4 for k in kap)) < 5: cnt += 1
print('  random parity-respecting labels with a level-1 (mod 4) pair collision: %d / 2000' % cnt)

# ---------------- (3) Hadamard cores --------------------------------------------------------------
def sylvester(k):
    H = [[1]]
    for _ in range(k):
        H = [row + row for row in H] + [row + [-v for v in row] for row in H]
    return H

def paley1(q):
    sq = set((x * x) % q for x in range(1, q))
    chi = lambda a: 0 if a % q == 0 else (1 if a % q in sq else -1)
    M = q + 1
    S = [[0] * M for _ in range(M)]
    for j in range(1, M): S[0][j] = 1; S[j][0] = -1
    for i in range(q):
        for j in range(q): S[1 + i][1 + j] = chi(j - i)
    H = [[S[i][j] + (1 if i == j else 0) for j in range(M)] for i in range(M)]
    for j in range(M):
        if H[0][j] == -1:
            for i in range(M): H[i][j] = -H[i][j]
    for i in range(M):
        if H[i][0] == -1: H[i] = [-v for v in H[i]]
    return H

print('\n(3) Hadamard cores')
for name, H in [('Sylvester 8', sylvester(3)), ('Sylvester 16', sylvester(4)), ('Paley 12', paley1(11)),
                ('Paley 20', paley1(19)), ('Paley 24', paley1(23))]:
    M = len(H)
    B = [[(1 + v) // 2 for v in row] for row in H]
    d = snf_diag(B)
    const_cols = [j for j in range(M) if len(set(B[i][j] for i in range(M))) == 1]
    core = [[B[i][j] for j in range(M) if j not in const_cols] for i in range(M)]
    core3 = [row * 3 for row in core]          # three copies of every nonconstant column
    dc = snf_diag(core3)
    # identity B^{-1} = 2H^T/M - e1 e1^T
    Binv = [[Fraction(2 * H[j][i], M) - (1 if i == j == 0 else 0) for j in range(M)] for i in range(M)]
    ident = all(sum(B[i][k] * Binv[k][j] for k in range(M)) == (1 if i == j else 0)
                for i in range(M) for j in range(M))
    print('  %-12s M=%2d: invariant factors of (J+H)/2 = %s; constant columns %s; '
          'valid core profile (nonconstant cols x3 copies): rank %d (Delta_M = %d); B^-1 identity: %s'
          % (name, M, d, const_cols, len(dc), 0 if len(dc) < M else math.prod(dc), ident))

# ---------------- (4) random 0/1 profiles ---------------------------------------------------------
print('\n(4) random 0/1 profiles with distinct rows and no constant column: Delta_M = 1 frequency')
for (M, r) in [(8, 16), (12, 24), (16, 32), (16, 64)]:
    ones = tot = 0
    while tot < 20:
        A = [[random.randint(0, 1) for _ in range(r)] for _ in range(M)]
        if len(set(map(tuple, A))) < M: continue
        if any(len(set(A[i][j] for i in range(M))) == 1 for j in range(r)): continue
        tot += 1
        if Delta(A) == 1: ones += 1
    print('  M=%2d r=%2d: Delta_M = 1 in %d/%d' % (M, r, ones, tot))
