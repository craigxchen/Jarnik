"""Search for product (zonotope) certificates of order m at k* in the widened slab box (UM_p at k*):
multisets {g_1..g_m}, l(g_i) >= 1, with every subset sum k* + sum_{i in I} g_i in A_p."""
import sys, itertools
import numpy as np
from had_sl import build

def search(A, lev, kstar, m):
    Aset = set(map(tuple, A.tolist()))
    ks = np.array(kstar)
    lstar = lev[[tuple(a) for a in A.tolist()].index(tuple(kstar))]
    G = [tuple(int(t) for t in (a - ks)) for a, l in zip(A, lev) if l > lstar]
    G = sorted(set(G))
    def ok(sums, g):
        new = [tuple(s[j] + g[j] for j in range(len(g))) for s in sums]
        return all(x in Aset for x in new), new
    found = []
    def dfs(start, sums, chosen):
        if len(chosen) == m:
            found.append(list(chosen)); return True
        for i in range(start, len(G)):
            good, new = ok(sums, G[i])
            if good:
                if dfs(i, sums + new, chosen + [G[i]]):
                    return True
        return False
    dfs(0, [tuple(kstar)], [])
    return found

if __name__ == '__main__':
    kdim, T, w, m = (int(a) for a in sys.argv[1:5])
    A, lev, d = build(kdim, T, w, m, 0)
    for kstar in [(0, 0, 0, 1, -1, 1, -1), (-1, 0, 0, 0, 0, 1, -1), (-1, 0, 0, 0, 0, 1, 0), (0, 1, 0, 0, 0, 0, 0)]:
        f = search(A, lev, kstar, m)
        print(kstar, "product certificate:", f[0] if f else None)
