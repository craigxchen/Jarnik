"""Referee checks for Part X.3 (Paley cuts). Exact integer arithmetic, own code."""
from itertools import combinations
from fractions import Fraction
import sys, math

def legendre(a, m):
    a %= m
    if a == 0:
        return 0
    return 1 if pow(a, (m - 1) // 2, m) == 1 else -1

def paley_S(m):
    return [[(legendre(v - i, m) if v != i else -1) for v in range(m)] for i in range(m)]

def rank_Q(M):
    M = [[Fraction(x) for x in row] for row in M]
    r = 0
    rows, cols = len(M), (len(M[0]) if M else 0)
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r

def det_int(M):
    # Bareiss
    n = len(M)
    if n == 0:
        return 1
    A = [row[:] for row in M]
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]

def check_structure(m):
    S = paley_S(m)
    n = m
    # Gram identities
    for i in range(n):
        for j in range(n):
            rr = sum(S[i][v] * S[j][v] for v in range(n))
            cc = sum(S[v][i] * S[v][j] for v in range(n))
            want = (m + 1 - 1) if i == j else -1
            assert rr == want and cc == want, (m, i, j)
    assert all(sum(row) == -1 for row in S)
    assert all(sum(S[i][v] for i in range(n)) == -1 for v in range(n))
    # pair separations
    seps = {sum(1 for v in range(n) if S[i][v] != S[j][v]) for i, j in combinations(range(n), 2)}
    # triple constant columns
    consts = set()
    ntrip = 0
    for i, j, k in combinations(range(n), 3):
        c = sum(1 for v in range(n) if S[i][v] == S[j][v] == S[k][v])
        consts.add(c)
        ntrip += 1
    # column-pair difference vectors a = S_v - S_w: support and sums
    supp = set()
    for v, w in combinations(range(n), 2):
        a = [S[i][v] - S[i][w] for i in range(n)]
        supp.add(sum(1 for x in a if x))
        assert sum(a) == 0
        for h in range(n):
            if h not in (v, w):
                plus = sum(a[i] for i in range(n) if S[i][h] == 1)
                minus = sum(a[i] for i in range(n) if S[i][h] == -1)
                assert plus == 0 and minus == 0
    d = det_int(S)
    return seps, consts, ntrip, supp, d

def deletion_check(m, kmax):
    S = paley_S(m)
    sets = ann = 0
    worst_sup = {}
    for k in range(kmax + 1):
        for M in combinations(range(m), k):
            T = [[1] * m] + [S[i][:] for i in M]
            r = rank_Q(T)
            # choose r independent rows of T (all rows here) and pivot columns
            rows = T  # rank checked equals len(T) below
            assert r == len(T), (m, M, r)
            pivots = next(c for c in combinations(range(m), r)
                          if det_int([[row[j] for j in c] for row in rows]) != 0)
            D = det_int([[row[j] for j in pivots] for row in rows])
            images = []
            bound = r ** (r / 2)
            for f in range(m):
                if f in pivots:
                    continue
                x = [0] * m
                x[f] = D
                # Cramer: sum_j T[:,pivot_j] x_pj = -D*T[:,f]
                for idx, pj in enumerate(pivots):
                    cols = list(pivots)
                    cols[idx] = f
                    x[pj] = -det_int([[row[j] for j in cols] for row in rows])
                assert all(sum(row[j] * x[j] for j in range(m)) == 0 for row in T)
                sup = max(abs(t) for t in x)
                assert sup <= bound + 1e-9
                assert sum(1 for t in x if t) <= r + 1
                worst_sup[r] = max(worst_sup.get(r, 0), sup)
                a = [sum(S[i][v] * x[v] for v in range(m)) for i in range(m)]
                assert sum(a) == 0 and all(a[i] == 0 for i in M)
                for v in range(m):
                    if x[v] == 0:
                        assert sum(a[i] * S[i][v] for i in range(m)) == 0
                    else:
                        assert sum(a[i] * S[i][v] for i in range(m)) == (m + 1) * x[v]
                images.append(a)
                ann += 1
            # span check: images span {a : a_M = 0, sum a = 0}, dim m-1-k
            assert rank_Q(images) == m - 1 - k
            sets += 1
    return sets, ann, worst_sup

if __name__ == "__main__":
    for m in [7, 11, 19, 23, 43]:
        seps, consts, nt, supp, d = check_structure(m)
        print(f"m={m}: pair separations {seps} (want {(m+1)//2}); triples {nt}, "
              f"constant-column counts {consts} (want {(m-3)//4}); "
              f"supp(S_v-S_w) {supp}; det S = {d} vs (m+1)^((m-1)/2) = {(m+1)**((m-1)//2)}")
        sys.stdout.flush()
    sets, ann, ws = deletion_check(11, 2)
    print("m=11 deletion, k<=2:", sets, "row sets,", ann, "annihilators, max sup by rank", ws)
    sets, ann, ws = deletion_check(19, 2)
    print("m=19 deletion, k<=2:", sets, "row sets,", ann, "annihilators, max sup by rank", ws)
    sets, ann, ws = deletion_check(11, 4)
    print("m=11 deletion, k<=4:", sets, "row sets,", ann, "annihilators, max sup by rank", ws)
