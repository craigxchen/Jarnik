"""Part 2: exact rank / dependency structure of the full quadruple (Ptolemy = Pluecker) system.

(R1) At random rational rank-two configurations, the Jacobian of all C(k,4) three-term relations
     F_abcd = p_ac p_bd - p_ab p_cd - p_ad p_bc  (in the C(k,2) bracket coordinates) has rank
     (k-2)(k-3)/2; the (k-2)(k-3)/2 relations through the fixed pair {0,1} already have full rank.
(R2) Linear form: for fixed nonzero residual data (Y_i, t_ij) the anchored relations
     Y_i t_jk G_jk - Y_j t_ik G_ik + Y_k t_ij G_ij = 0 are linear in G=(G_ij)_{1<=i<j<=m};
     the solution space is exactly {((X_i Y_j - X_j Y_i)/t_ij)} of dimension m-1.
(R3) Generation: every 4-term relation among nonanchor rows lies in the span of the
     anchored linear relations once the G's are substituted (checked by rank).
(R4) Print the counts of unknowns/equations and the heuristic exponents (see ptolemy.md, Sec. 4).
"""
from fractions import Fraction as F
import random, itertools

random.seed(7)

def rank(rows):
    M = [list(r) for r in rows]
    if not M:
        return 0
    ncol = len(M[0]); rk = 0
    for c in range(ncol):
        piv = None
        for r in range(rk, len(M)):
            if M[r][c] != 0:
                piv = r; break
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        pv = M[rk][c]
        for r in range(len(M)):
            if r != rk and M[r][c] != 0:
                f = M[r][c] / pv
                M[r] = [x - f * y for x, y in zip(M[r], M[rk])]
        rk += 1
    return rk

def plucker_jacobian(k):
    vecs = [(F(random.randint(-50, 50)), F(random.randint(-50, 50))) for _ in range(k)]
    pairs = list(itertools.combinations(range(k), 2))
    idx = {pq: n for n, pq in enumerate(pairs)}
    p = {(a, b): vecs[a][0]*vecs[b][1] - vecs[b][0]*vecs[a][1] for a, b in pairs}
    rows = []; rows01 = []
    for a, b, c, d in itertools.combinations(range(k), 4):
        assert p[a, c]*p[b, d] - p[a, b]*p[c, d] - p[a, d]*p[b, c] == 0
        row = [F(0)] * len(pairs)
        row[idx[a, c]] += p[b, d]; row[idx[b, d]] += p[a, c]
        row[idx[a, b]] -= p[c, d]; row[idx[c, d]] -= p[a, b]
        row[idx[a, d]] -= p[b, c]; row[idx[b, c]] -= p[a, d]
        rows.append(row)
        if a == 0 and b == 1:
            rows01.append(row)
    return rank(rows), rank(rows01), len(rows), len(rows01)

print("(R1) Jacobian ranks of the full Ptolemy/Pluecker system")
for k in range(4, 10):
    r, r01, n, n01 = plucker_jacobian(k)
    expect = (k - 2) * (k - 3) // 2
    assert r == expect and r01 == expect and n01 == expect
    print(f"  k={k}: C(k,4)={n} relations, rank {r} = (k-2)(k-3)/2; relations through {{0,1}}: {n01}, rank {r01}")

print("(R2,R3) anchored linear system in G for fixed small residues")
for m in range(3, 8):
    for trial in range(3):
        Y = [random.choice([-3, -2, -1, 1, 2, 3]) for _ in range(m + 1)]
        t = {}
        for i, j in itertools.combinations(range(1, m + 1), 2):
            t[i, j] = random.choice([-4, -3, -2, -1, 1, 2, 3, 4])
        pairs = list(itertools.combinations(range(1, m + 1), 2))
        idx = {pq: n for n, pq in enumerate(pairs)}
        rows = []
        for i, j, k in itertools.combinations(range(1, m + 1), 3):
            row = [F(0)] * len(pairs)
            row[idx[j, k]] += Y[i] * t[j, k]
            row[idx[i, k]] -= Y[j] * t[i, k]
            row[idx[i, j]] += Y[k] * t[i, j]
            rows.append(row)
        rk = rank(rows)
        assert rk == len(pairs) - (m - 1), (m, rk)
        # an explicit solution from rows (X_i, Y_i): G_ij = (X_i Y_j - X_j Y_i)/t_ij
        X = [F(random.randint(-10**6, 10**6)) for _ in range(m + 1)]
        G = [(X[i]*Y[j] - X[j]*Y[i]) / t[i, j] for i, j in pairs]
        for row in rows:
            assert sum(a * b for a, b in zip(row, G)) == 0
        # (R3) every nonanchor 4-term relation (quadratic in G) holds on the solution space
        for a, b, c, d in itertools.combinations(range(1, m + 1), 4):
            pr = lambda u, v: t[u, v] * G[idx[u, v]]
            assert pr(a, c)*pr(b, d) == pr(a, b)*pr(c, d) + pr(a, d)*pr(b, c)
    print(f"  m={m} (k={m+1}): {len(pairs)} unknown G_ij, {len(rows)} anchored relations, rank {rk} = C(m,2)-(m-1)")

print("(R4) counts: blocks relevant to matchings, independent relations, heuristic exponents (units of w)")
print("  k   m  blocks(|T|>=2)  indep.rel  angular-only  with-radial  lead-naive")
for k in range(5, 11):
    m = k - 1
    blocks = 2**m - 1 - m
    indep = (k - 2) * (k - 3) // 2
    angular = 2**(m - 2) * (7 - m) - m - 2
    radial = 2**(m - 2) * (4 - m) - 1
    naive = blocks - (k*(k-1)*(k-2)*(k-3)//24) * 2**(k - 4)
    print(f"  {k:2d} {m:2d}  {blocks:6d}        {indep:4d}       {angular:6d}       {radial:6d}      {naive:6d}")
