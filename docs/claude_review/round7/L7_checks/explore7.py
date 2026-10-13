"""NUMERICAL: genuine-solution counts per variable for parametrisation param2 (b, tau, mu, sp)."""
import numpy as np, itertools
from fractions import Fraction as Fr
from collections import Counter
import param2 as P2
import cub13f as F
from modp import M, ratdeg
Pn = np.polynomial.polynomial

def curve_float(d):
    a, b, c = d['a'], d['b'], d['c']
    l2, l3, l5 = d['l2'], d['l3'], d['l5']
    L = lambda r: np.array([-r, 1], complex)
    y5 = l5 * Pn.polymul(L(0), L(a)); y4 = Pn.polymul(Pn.polymul(L(0), L(1)), L(b))
    y3 = l3 * Pn.polymul(L(1), L(c)); y2 = l2 * Pn.polymul(Pn.polymul(L(a), L(b)), L(c))
    base = Pn.polymul(Pn.polymul(L(d['p25']), L(d['p24'])), L(d['sigma']))
    kap = (Pn.polyval(d['p35'], y3) - Pn.polyval(d['p35'], y2)) / Pn.polyval(d['p35'], base)
    y1 = Pn.polyadd(y2, kap * base)
    return [y1, y2, y3, y4, y5]

def genuine(d):
    c = curve_float(d)
    cnt = Counter()
    for p, big in F.events(c):
        for g in big: cnt[g] += 1
    trip = {g for g in cnt if len(g) == 3}
    pairs = {g: cnt[g] for g in cnt if len(g) == 2}
    return trip == F.TRIPLES13 and all(cnt[g] == 1 for g in trip) and pairs == {(1,4):1,(1,5):1,(2,3):1,(1,6):3} and not any(len(g) >= 4 for g in cnt)

base = dict(b=Fr(3), tau=Fr(-2), mu=Fr(5, 2), sp=Fr(7, 3))
def fexact(var, x):
    a = dict(base); a[var] = x
    d = P2.stage1(a['b'], a['tau'], a['mu'], a['sp'])
    return P2.stage2(d)

from ratrec import nullspace
import qpoly as Q
def exact_num(var, dn, dd):
    pts = []; k = 0
    while len(pts) < dn + dd + 6:
        x = Fr(3 * k + 1, 11) + 2; k += 1
        try: pts.append((x, fexact(var, x)))
        except ZeroDivisionError: pass
    rows = [[x ** i for i in range(dn + 1)] + [-y * x ** j for j in range(dd + 1)] for x, y in pts]
    ns = nullspace(rows); assert len(ns) == 1
    v = ns[0]; return Q.norm(v[:dn + 1]), Q.norm(v[dn + 1:])

for var in ['sp', 'b', 'tau', 'mu']:
    def fm(x):
        a = {k: M(v) for k, v in base.items()}; a[var] = x
        return P2.stage2(P2.stage1(a['b'], a['tau'], a['mu'], a['sp']))
    dn, dd, _, _ = ratdeg(fm, 40)
    num, den = exact_num(var, dn, dd)
    rr = sorted(Q.rational_roots(num)) if max(abs(c.numerator) for c in num) < 10**30 else 'big'
    coeffs = np.array([float(c) for c in num])
    roots = Pn.polyroots(coeffs)
    gen = []
    for r in roots:
        a = {k: complex(v) for k, v in base.items()}; a[var] = complex(r)
        try:
            d = P2.stage1(a['b'], a['tau'], a['mu'], a['sp']); P2.stage2(d)
            if genuine(d): gen.append(r)
        except Exception: pass
    print(var, 'num deg', dn, 'den deg', dd, '| genuine roots:', len(gen), np.round(gen, 4), flush=True)
