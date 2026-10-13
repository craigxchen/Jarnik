import numpy as np
from collections import Counter
import cub13f as F
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

