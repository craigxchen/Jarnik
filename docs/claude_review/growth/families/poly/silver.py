# Pell-number blocks H_r = P_{r+1} + i P_r (norm P_{2r+1}), unit 1+sqrt2, limiting phase pi/8.
# Rows of one parity class need no prefactor (offsets are multiples of pi/4... see note): we use
# prefactor 1 and let unit rotations align (offset (2|a|-E)*pi/8 must be = 0 mod pi/2 -> |a| = E/2 mod 2).
import itertools, math, sys
from gpoly import gmul
from golden import subspans
_p = [0, 1]
def pell(r):
    while len(_p) <= r: _p.append(2*_p[-1] + _p[-2])
    return _p[r]
def H(r): return (pell(r+1), pell(r))
def conj(z): return (z[0], -z[1])
def gpow(z, e):
    p = (1, 0)
    for _ in range(e): p = gmul(p, z)
    return p
def evaluate(n, shifts, widths, rows):
    pts = []
    for a in rows:
        z = (1, 0)
        for s, w, x in zip(shifts, widths, a):
            b = H(n+s)
            z = gmul(z, gmul(gpow(b, x), gpow(conj(b), w-x)))
        pts.append(z)
    return pts
best = {}
S = 8
pats = [([0,s2,s3,s4],[1,1,1,1]) for s2,s3,s4 in itertools.combinations(range(1,S+1),3)]
pats += [([0,s2,s3],w) for s2,s3 in itertools.combinations(range(1,S+1),2) for w in ([2,1,1],[1,2,1],[1,1,2])]
for shifts, widths in pats:
    E = sum(widths)
    box = list(itertools.product(*[range(w+1) for w in widths]))
    for cls in (0, 1):
        rows = [a for a in box if (sum(a) - E//2) % 2 == cls]
        res = {}
        for n in range(40, 64):
            b, N, g = subspans(evaluate(n, shifts, widths, rows))
            if b is None: continue
            for k, v in b.items():
                if k not in best or v < best[k][0]:
                    best[k] = (v, shifts, widths, cls, n)
for k in sorted(best): print(k, best[k])
