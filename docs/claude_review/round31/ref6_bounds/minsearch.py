# Pruned search for the least reduced degree D' of a c=0 cluster formed by ALL s=Q/2 marks of a
# hyperplane design with translation B (Y = T + theta B), V = F_2[t]_{<d}, over all side choices.
# Pruning is exact: adding marks can only lower deg G and raise max diff, so the criticality
# condition 2*maxdiff - degF <= deg G, once violated, stays violated; D' only grows.
import sys, itertools
from f4lib import *
from clusters import span
from test_vi3 import coords, parity

def search(d, B, T, best_cap=10**9):
    basis = [1 << i for i in range(d)]
    V = span(basis); co = coords(basis); Q = 1 << d
    forms = [l for l in range(1, Q) if parity(l & co[B])]
    s = len(forms)
    Y = (T, B)
    lin = {v: add4(Y, (v, 0)) for v in V}
    Fp = {}
    for l in forms:
        for c in (0, 1):
            Fp[(l, c)] = prod4([lin[v] for v in V if parity(l & co[v]) == c])
    degF = d4(next(iter(Fp.values())))
    mon = {k: monic4(x) for k, x in Fp.items()}
    keys = list(Fp)
    dd = {}
    for a in keys:
        for b in keys:
            if a != b:
                dd[(a, b)] = d4(add4(mon[a], mon[b]))
    best = [best_cap, None]
    # symmetry: fix the side of the first form? Not a symmetry in general (conjugation maps the
    # family to conjugates), so enumerate all.
    def dfs(i, G, chosen, maxd):
        g = d4(G) if G is not None else degF
        if 2 * maxd - degF > g:
            return
        if degF - g >= best[0]:
            return
        if i == s:
            best[0] = degF - g; best[1] = tuple(c for _, c in chosen)
            return
        l = forms[i]
        for c in (0, 1):
            k = (l, c)
            md = maxd; ok = True
            for k2 in chosen:
                x = dd[(k, k2)]
                if x < 0:
                    ok = False; break
                md = max(md, x)
            if not ok:
                continue
            G2 = Fp[k] if G is None else gcd4(G, Fp[k])
            dfs(i + 1, G2, chosen + [k], md)
    dfs(0, None, [], -1)
    return best[0] if best[1] is not None else None, best[1], degF

if __name__ == '__main__':
    d = int(sys.argv[1]); Bs = [int(x) for x in sys.argv[2].split(',')]
    Lmin, Lmax, low = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    for B in Bs:
        bestB = None
        for L in range(Lmin, Lmax + 1):
            for x in range(1 << low):
                T = (1 << L) ^ x
                c = T
                for b in [1 << i for i in range(d)]:
                    c = gcd2(c, b)
                if c != 1:
                    continue
                r = search(d, B, T, bestB[0] + 1 if bestB else 10**9)
                if r[0] is not None and (bestB is None or r[0] < bestB[0]):
                    bestB = (r[0], T, r[1], r[2])
        print(f"d={d} B={bin(B)}: least D' found = {bestB[0] if bestB else None}  T={bin(bestB[1]) if bestB else ''} sides={bestB[2] if bestB else ''} raw deg={bestB[3] if bestB else ''}", flush=True)
