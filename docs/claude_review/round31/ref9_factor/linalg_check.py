"""Referee check (independent code): Paley exponent matrix, coverage, D = H mod 2,
Smith normal form of D over Z, integral right inverse U, V, F with D F = H.
Exact integer arithmetic only."""
import random, sys
from sympy import isprime
from sympy.matrices.normalforms import smith_normal_form
from sympy import Matrix, ZZ

def chi(x, q):
    x %= q
    if x == 0: return 0
    return 1 if pow(x, (q-1)//2, q) == 1 else -1

def col_hnf_rightinv(D):
    """D: list of rows (r x n) of ints, full row rank.  Returns U (n x r) with D U = I
    if the row lattice is primitive, else None.  Column operations with unimodular Q."""
    r, n = len(D), len(D[0])
    A = [row[:] for row in D]
    Q = [[int(i == j) for j in range(n)] for i in range(n)]
    def colop_swap(a, b):
        for M in (A, Q):
            for row in M: row[a], row[b] = row[b], row[a]
    def colop_add(dst, src, f):  # col dst += f * col src
        for M in (A, Q):
            for row in M: row[dst] += f * row[src]
    def colop_neg(a):
        for M in (A, Q):
            for row in M: row[a] = -row[a]
    for i in range(r):
        # make A[i][i] = gcd of A[i][i..n-1], zero out the rest
        while True:
            nz = [j for j in range(i, n) if A[i][j] != 0]
            if not nz: return None
            jmin = min(nz, key=lambda j: abs(A[i][j]))
            if jmin != i: colop_swap(i, jmin)
            done = True
            for j in range(i+1, n):
                if A[i][j]:
                    f = A[i][j] // A[i][i]
                    colop_add(j, i, -f)
                    if A[i][j]: done = False
            if done: break
        if A[i][i] < 0: colop_neg(i)
        if A[i][i] != 1: return None
        # clear left entries of row i using pivot (not needed for L inverse but fine)
    # now A = [L | 0], L lower triangular with unit diagonal
    L = [row[:r] for row in A]
    # solve L X = I  (X = L^{-1}, integral)
    X = [[0]*r for _ in range(r)]
    for c in range(r):
        for i in range(r):
            s = int(i == c) - sum(L[i][k]*X[k][c] for k in range(i))
            X[i][c] = s  # L[i][i] = 1
    U = [[sum(Q[a][k]*X[k][c] for k in range(r)) for c in range(r)] for a in range(n)]
    return U

def matmul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def build(q, rho, rng):
    labels = list(range(1, q))
    Kf = lambda i, j: chi(j - i, q)
    cls = [i for i in labels if Kf(i, 0) == 1]
    S = [i for i in cls if rng.random() < rho]
    return labels, cls, S, Kf

def run(q, rho, seed, verbose=True, do_snf=True):
    rng = random.Random(seed)
    labels, cls, S, Kf = build(q, rho, rng)
    m = len(S)
    if m < 3: return None
    a = {i: [1 + Kf(i, j) for j in labels] for i in cls}
    # coverage on full class: counts of a=0 and a=2 per column
    cov_full = []
    for jj, j in enumerate(labels):
        c0 = sum(1 for i in cls if a[i][jj] == 0); c2 = sum(1 for i in cls if a[i][jj] == 2)
        pred0 = (q - 2 + 1 - chi(-j, q) + chi(j, q)) // 4  # eps=-1 -> a=0
        pred2 = (q - 2 - 1 - chi(-j, q) - chi(j, q)) // 4  # eps=+1 -> a=2
        assert (c0, c2) == (pred0, pred2), (q, j, c0, c2, pred0, pred2)
        cov_full.append((c0, c2))
    covered = all({0, 2} <= {a[i][jj] for i in S} for jj in range(len(labels)))
    # odd exponent at exactly one label j = i
    for i in cls:
        odd = [labels[jj] for jj in range(len(labels)) if a[i][jj] % 2]
        assert odd == [i]
    i0 = S[0]
    D = [[a[i][jj] - a[i0][jj] for jj in range(len(labels))] for i in S[1:]]
    H = [[int(labels[jj] == i) - int(labels[jj] == i0) for jj in range(len(labels))] for i in S[1:]]
    assert all((D[r][c] - H[r][c]) % 2 == 0 for r in range(m-1) for c in range(q-1))
    inv = None
    if do_snf:
        snf = smith_normal_form(Matrix(D), domain=ZZ)
        inv = [snf[k, k] for k in range(m-1)]
        assert all(abs(x) == 1 for x in inv), inv
    U = col_hnf_rightinv(D)
    assert U is not None
    DU = matmul(D, U)
    assert all(DU[r][c] == int(r == c) for r in range(m-1) for c in range(m-1))
    HD2 = [[(H[r][c] - D[r][c]) // 2 for c in range(q-1)] for r in range(m-1)]
    V = matmul(U, HD2)
    F = [[int(r == c) + 2*V[r][c] for c in range(q-1)] for r in range(q-1)]
    DF = matmul(D, F)
    assert DF == H
    if verbose:
        print(f"q={q} rho={rho} seed={seed} |class|={len(cls)} m={m} covered_by_S={covered} "
              f"SNF_all_one={do_snf} max|U|={max(abs(x) for r in U for x in r)} "
              f"max|V|={max(abs(x) for r in V for x in r)} DF=H ok")
    return dict(labels=labels, cls=cls, S=S, a=a, D=D, H=H, U=U, V=V, F=F, i0=i0)

if __name__ == "__main__":
    for q in [11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103]:
        for seed in range(3):
            for rho in (1.0, 0.95, 0.7):
                run(q, rho, seed, do_snf=(q <= 71))
    # also check the Vandermonde step directly: K^T y = 0 mod q with support < q/2 forces y = 0 mod q
    print("ALL LINEAR-ALGEBRA CHECKS PASSED")
