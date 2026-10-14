"""Referee: enumerate all multisets of 7 heavy forms (e_a +- e_b, 42 forms) in Y_3(2) coordinates with
axis counts <= 4 and line counts <= 5; classify by width-count feasibility in the 43 cheap directions."""
import itertools, numpy as np, pickle, sys
from collections import Counter
Y, D, T = pickle.load(open(sys.argv[1], 'rb'))
cheapD = D[T <= 6]; cheapT = T[T <= 6]
axes = list(range(1, 8))
heavy = []
for i, j in itertools.combinations(range(7), 2):
    for s in (1, -1):
        g = np.zeros(7, dtype=np.int64); g[i] = 1; g[j] = s; heavy.append(g)
heavy = np.array(heavy)
lines = [(a, c, a ^ c) for a in range(1, 8) for c in range(a + 1, 8) if (a ^ c) > c]
Lmask = np.array([[1 if (k + 1) in L else 0 for k in range(7)] for L in lines])
touchL = (np.abs(heavy) @ Lmask.T > 0).astype(int)      # 42 x 7
touchA = (heavy != 0).astype(int)                        # 42 x 7
touchD = ((heavy @ cheapD.T) != 0).astype(int)           # 42 x 43
out = []
def rec(start, ch, ca, cl):
    if len(ch) == 7:
        out.append(tuple(ch)); return
    for i in range(start, 42):
        na = ca + touchA[i]; nl = cl + touchL[i]
        if na.max() <= 4 and nl.max() <= 5:
            ch.append(i); rec(i, ch, na, nl); ch.pop()
rec(0, [], np.zeros(7, int), np.zeros(7, int))
print('multisets (axis<=4, line<=5):', len(out))
stats = Counter()
feas = []
for ms in out:
    idx = list(ms)
    ax = touchA[idx].sum(0).max()
    cd = touchD[idx].sum(0)
    viol = (cd > cheapT).any()
    # which kind of direction is violated
    kinds = set()
    for k in np.nonzero(cd > cheapT)[0]:
        sup = (cheapD[k] != 0).sum(); kinds.add({1: 'axis', 3: 'line', 7: 'ord4'}[sup])
    stats[(int(ax), bool(viol), tuple(sorted(kinds)))] += 1
    if not viol: feas.append(ms)
for k, v in sorted(stats.items()): print('  maxaxis=%d width-violated=%s kinds=%s : %d' % (k[0], k[1], k[2], v))
print('width-count feasible multisets:', len(feas), Counter(int(touchA[list(m)].sum(0).max()) for m in feas))
pickle.dump((heavy, out, feas), open(sys.argv[2], 'wb'))
