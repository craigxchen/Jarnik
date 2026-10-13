"""curve from param3 data (works with complex floats via numpy, or exactly via qpoly)."""
import numpy as np
from collections import Counter
import cub13f as F
Pn = np.polynomial.polynomial
def L(r): return np.array([-r, 1], complex)
def curve_float(d):
    a, b, c, l2, l3, l5 = [complex(d[k]) for k in ('a', 'b', 'c', 'l2', 'l3', 'l5')]
    y5 = l5 * Pn.polymul(L(0), L(a)); y4 = Pn.polymul(Pn.polymul(L(0), L(1)), L(b))
    y3 = l3 * Pn.polymul(L(1), L(c)); y2 = l2 * Pn.polymul(Pn.polymul(L(a), L(b)), L(c))
    base = Pn.polymul(Pn.polymul(L(complex(d['p25'])), L(complex(d['p24']))), L(complex(d['sigma'])))
    p35 = complex(d['p35'])
    kap = (Pn.polyval(p35, y3) - Pn.polyval(p35, y2)) / Pn.polyval(p35, base)
    y1 = Pn.polyadd(y2, kap * base)
    return [y1, y2, y3, y4, y5]
def classify(c):
    cnt = Counter()
    for p, big in F.events(c):
        for g in big: cnt[g] += 1
    return cnt
def genuine(d):
    try:
        cnt = classify(curve_float(d))
    except Exception:
        return False
    trip = {g for g in cnt if len(g) == 3}
    pairs = {g: cnt[g] for g in cnt if len(g) == 2}
    return trip == F.TRIPLES13 and all(cnt[g] == 1 for g in trip) and pairs == {(1,4):1,(1,5):1,(2,3):1,(1,6):3} and not any(len(g) >= 4 for g in cnt)
