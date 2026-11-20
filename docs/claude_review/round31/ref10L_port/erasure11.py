"""Paley order 11 with up to two missing rows: kernel of T = [ones; missing rows of S], images a = S x."""
from fractions import Fraction
from itertools import combinations
m = 11
QR = {x*x % m for x in range(1, m)}
S = [[(1 if (v-i) % m in QR else -1) if v != i else -1 for v in range(m)] for i in range(m)]

def kernel(T):
    A = [[Fraction(x) for x in row] for row in T]
    piv, r = [], 0
    for c in range(m):
        k = next((i for i in range(r, len(A)) if A[i][c]), None)
        if k is None: continue
        A[r], A[k] = A[k], A[r]
        A[r] = [x/A[r][c] for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]; A[i] = [x - f*y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    basis = []
    for f in range(m):
        if f in piv: continue
        x = [Fraction(0)]*m; x[f] = Fraction(1)
        for i, c in enumerate(piv): x[c] = -A[i][f]
        basis.append(x)
    return r, basis

sets = vecs = 0
for k in range(3):
    for om in combinations(range(m), k):
        T = [[1]*m] + [S[i] for i in om]
        r, B = kernel(T)
        assert r == k+1
        sets += 1
        for x in B:
            a = [sum(S[i][v]*x[v] for v in range(m)) for i in range(m)]
            assert any(a) and sum(a) == 0 and all(a[i] == 0 for i in om)
            for v in range(m):
                if x[v] == 0:
                    assert sum(a[i]*S[i][v] for i in range(m)) == 0
            vecs += 1
print(sets, "row sets,", vecs, "kernel basis vectors (all checks pass)")
