"""Exhaustive search for small-height Q-points of W: residual(r25,r34,r35,w) == 0 exactly (Fractions),
then exact genuineness check (class == cub13) over Q."""
import sys, itertools, time
from fractions import Fraction as Fr
import cub13exact as C
H = int(sys.argv[1]); part = int(sys.argv[2]); nparts = int(sys.argv[3])
rats = sorted({Fr(a, b) for b in range(1, H + 1) for a in range(-H, H + 1)} - {Fr(0), Fr(1)})
one = Fr(1)
hits = []
t0 = time.time(); n = 0
for i, r25 in enumerate(rats):
    if i % nparts != part: continue
    for r34, r35, w in itertools.product(rats, repeat=3):
        try:
            d = C.build(r25, r34, r35, w, one)
            res = d['residual']
        except ZeroDivisionError:
            continue
        n += 1
        if res == 0:
            hits.append((r25, r34, r35, w))
            print('ZERO', r25, r34, r35, w, flush=True)
print('done', n, 'evaluations', time.time() - t0, 's; zeros', len(hits), flush=True)
