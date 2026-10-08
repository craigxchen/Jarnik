# Rows b_j = b0 + M*Q_j, M = 1+3X+X^2 (min poly of q=-phi^-2 up to reversal: X^2+3X+1 is palindromic),
# Q_j(1) even, so every pair difference vanishes at q and has even coefficient sum: contact order
# nu>=3 for the Fibonacci blocks (eta purely imaginary).  Bounded C needs total width E<=12.
import itertools
from golden import template_best
M = [1, 3, 1]
def mul(Q):
    out = [0]*(len(Q)+2)
    for i, a in enumerate(Q):
        for j, b in enumerate(M):
            out[i+j] += a*b
    return out
L = 4
Qs = [Q for Q in itertools.product(range(-1, 2), repeat=L) if sum(Q) % 2 == 0]
vec = {Q: mul(Q) for Q in Qs}
zero = tuple([0]*L)
others = [Q for Q in Qs if Q != zero]
best = []
def width(rowset):
    cols = list(zip(*[vec[Q] for Q in rowset]))
    return sum(max(c) - min(c) for c in cols)
found = {}
def dfs(cur, start, E):
    k = len(cur)
    if k >= 5:
        found.setdefault(k, []).append((E, tuple(cur)))
    if k == 6: return
    for idx in range(start, len(others)):
        nxt = cur + [others[idx]]
        e = width(nxt)
        if e <= 12:
            dfs(nxt, idx+1, e)
dfs([zero], 0, 0)
for k in sorted(found):
    print("k=%d sets=%d minE=%d" % (k, len(found[k]), min(e for e, _ in found[k])))
import pickle
pickle.dump(found, open('nu3_sets.pkl', 'wb'))
