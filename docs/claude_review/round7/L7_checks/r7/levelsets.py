"""NUMERICAL contact analysis for a curve c(t) = (c_1..c_5) (any degree) in the Kapranov psi_7 model:
configuration (c_1,...,c_5, 0, inf).  Returns Counter {frozenset T (subset of 1..6): number of contacts},
the list of points with >= 2 nontrivial level sets (codim-2 strata), and min distance between contacts."""
import numpy as np, itertools
from collections import Counter
P = np.polynomial.polynomial
def remob(cs, M=(1.3, -0.7, 0.45, 1.1)):
    a, b, g, d = M
    deg = max(len(c) for c in cs) - 1
    out = []
    for c in cs:
        c = list(c) + [0] * (deg + 1 - len(c))
        y = np.zeros(deg + 1, complex)
        for k in range(deg + 1):
            term = P.polymul(P.polypow([b, a], k), P.polypow([d, g], deg - k))
            y[:len(term)] += c[k] * term
        out.append(y)
    return out
def analyse(cs, tol=1e-6):
    cs = remob([np.array(c, complex) for c in cs])
    cs6 = cs + [np.zeros_like(cs[0])]
    scale = max(np.max(np.abs(c)) for c in cs)
    roots = []
    for i, j in itertools.combinations(range(6), 2):
        dpol = np.trim_zeros(cs6[i] - cs6[j], 'b')
        if len(dpol) > 1 and np.max(np.abs(dpol)) > 1e-12 * scale:
            roots += list(P.polyroots(dpol))
    pts = []
    for r in roots:
        if all(abs(r - q) > 1e-5 * (1 + abs(r)) for q in pts): pts.append(r)
    cnt = Counter(); bad = []
    for r in pts:
        vals = [P.polyval(r, c) for c in cs6]
        nrm = max(abs(v) for v in vals)
        groups = []
        for i in range(6):
            for g in groups:
                if abs(vals[i] - vals[g[0]]) < tol * nrm: g.append(i); break
            else: groups.append([i])
        big = [frozenset(x + 1 for x in g) for g in groups if len(g) >= 2]
        if len(big) >= 2: bad.append((r, big))
        for T in big: cnt[T] += 1
    md = min((abs(p - q) for p, q in itertools.combinations(pts, 2)), default=np.inf)
    lead = [c[-1] for c in cs6]
    inf_bad = any(abs(lead[i] - lead[j]) < 1e-9 * scale for i, j in itertools.combinations(range(6), 2))
    return cnt, bad, md, inf_bad
