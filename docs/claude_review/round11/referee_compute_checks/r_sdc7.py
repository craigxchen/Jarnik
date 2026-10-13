# Skew double core at q=7 (round5/exact.md Thm 4.1), built from scratch; slab-box T=1, L=2; all base points.
import sys, itertools, json
import flint
from rmodels import slab_box, moment_matrix, rank_mod
q = 7
def chi(a):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1
n = q + 1
# index 0 = infinity, 1..7 = F_7 elements 0..6
S = [[0] * n for _ in range(n)]
for y in range(q):
    S[0][1 + y] = 1; S[1 + y][0] = -1
for x in range(q):
    for y in range(q):
        S[1 + x][1 + y] = chi(y - x)
K = [[(1 if i == j else 0) + S[i][j] for j in range(n)] for i in range(n)]
# checks: skew-Hadamard
for i in range(n):
    for j in range(n):
        assert K[i][j] + K[j][i] == (2 if i == j else 0)
        assert sum(K[i][t] * K[j][t] for t in range(n)) == (n if i == j else 0)
D = [K[x][0] for x in range(n)]
K = [[D[i] * K[i][j] * D[j] for j in range(n)] for i in range(n)]
assert all(K[x][0] == 1 for x in range(n))
cuts = []
for a in range(1, n):
    Sa = frozenset(x for x in range(n) if K[x][a] == -1)
    Spa = Sa ^ frozenset([0, a])
    cuts += [Sa, Spa]
assert len(set(cuts)) == 14 and all(len(c) == 4 for c in cuts)
# pair separation: every pair separated by exactly n cuts
for x, y in itertools.combinations(range(n), 2):
    assert sum((x in c) != (y in c) for c in cuts) == n
print("skew double core q=7: 14 distinct balanced cuts, every pair separated by exactly 8 (Plotkin-tight)")
res = {}
for base in range(n):
    A, others = slab_box(n, cuts, 1, 2, base=base)
    ranks = {}
    for m in range(6, 11):
        rows, _ = moment_matrix(A, others, m, 1)
        ranks[m] = rank_mod(rows)
    print(f"base={base}: |A|={len(A)}  dimV_p(A,m) for m=6..10: {[len(A)-ranks[m] for m in range(6,11)]}", flush=True)
json.dump([sorted(c) for c in cuts], open('r_sdc7_cuts.json', 'w'))
