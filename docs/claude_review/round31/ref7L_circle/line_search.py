"""Line analogue: all roots of N on the imaginary axis (a vertical line, i.e. a 'circle with real centre at infinity').
Mark = multiset of roots i*s, s in {+beta_k,-beta_k}. Contact h <=> sum s^m equal for odd m<h (even m automatic).
Search integer beta, find triples (and 5-tuples) with 2h > D after primitive reduction."""
import itertools, sys
from collections import Counter, defaultdict
from fractions import Fraction
def contact(marks):
    m = 1
    while True:
        ps = {sum(Fraction(x)**m for x in mk) for mk in marks}
        if len(ps) > 1: return m
        m += 1
        if m > 60: return 99
def reduce_common(marks):
    C = Counter(marks[0])
    for mk in marks[1:]: C &= Counter(mk)
    out = []
    for mk in marks:
        c = Counter(mk); c.subtract(C); out.append(tuple(sorted(c.elements())))
    return out
found = defaultdict(list)
B = int(sys.argv[1]); Dmax = int(sys.argv[2]); k = int(sys.argv[3])
for D in range(3, Dmax+1):
    best = None
    for beta in itertools.combinations_with_replacement(range(1, B+1), D):
        polys = set()
        for sg in itertools.product([1, -1], repeat=D):
            polys.add(tuple(sorted(s*b for s, b in zip(sg, beta))))
        # group by p1 (necessary for h>=2)
        groups = defaultdict(list)
        for p in polys: groups[sum(p)].append(p)
        for g in groups.values():
            if len(g) < k: continue
            for T in itertools.combinations(g, k):
                R = reduce_common(list(T))
                Dr = len(R[0])
                if Dr < D: continue      # require primitive at this D
                hh = contact(R)
                if (k == 3 and 2*hh > Dr) or (k == 5 and 20*hh > 9*Dr):
                    if best is None or hh/Dr > best[0]: best = (hh/Dr, Dr, hh, T)
    print('D', D, 'best', best, flush=True)
