# Exact check: in the Paley skew double core of order M = q+1, every pair of rows is
# separated by exactly M of the 2q cut types (so all pair margins are equal when the
# 2q primes have equal weight). Used for the "exactly uniform margins" remark in README §3.
def chi(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1

def core(q):
    idx = [None] + list(range(q))          # None = border row/column 0
    M = q + 1
    def S(x, y):
        if x is None and y is None: return 0
        if x is None: return -1
        if y is None: return 1
        return chi(y - x, q)
    K = [[(1 if i == j else 0) + S(idx[i], idx[j]) for j in range(M)] for i in range(M)]
    # checks: K + K^T = 2I, K K^T = M I
    for i in range(M):
        for j in range(M):
            assert K[i][j] + K[j][i] == (2 if i == j else 0)
            assert sum(K[i][t] * K[j][t] for t in range(M)) == (M if i == j else 0)
    types = []
    for a in range(1, M):                  # column a (finite labels)
        s1 = [K[x][a] for x in range(M)]
        s2 = [K[a][x] * K[0][x] for x in range(M)]   # K^T(x,a) K^T(x,0)
        types += [s1, s2]
    return K, types

for q in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
    K, types = core(q)
    M = q + 1
    seps = set()
    for x in range(M):
        for y in range(x + 1, M):
            seps.add(sum(1 for s in types if s[x] != s[y]))
    print(q, M, len(types), sorted(seps))
    assert seps == {M}
print("ALL CHECKS PASSED: True")
