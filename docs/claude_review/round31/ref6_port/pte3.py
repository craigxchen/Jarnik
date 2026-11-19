"""Vertical-line templates with ternary-style exponents (characteristic 0 analogue of VI.4 at the
polynomial level): z_s = prod_{x in S} (X - i x)^{1+s_x} (X + i x)^{1-s_x}, s in {-1,0,1}^d, common norm
prod (X^2+x^2)^2, divided by the common gcd. Per coordinate type: A (s in {-1,0,1}), B (s in {0,1}),
C (s in {-1,0}); g = #B + #C (forced coordinates can be deleted). A family inside the box is a c=0 cluster
iff all its odd moments sum s_x x^r agree for odd r < d - g/2 (char 0, Newton). Max class size, exhaustive."""
import itertools
import numpy as np
from collections import Counter

def run(d, Xb):
    best = (0, None)
    for S in itertools.combinations(range(1, Xb + 1), d):
        for tau in itertools.product('ABC', repeat=d):
            g = sum(1 for c in tau if c != 'A')
            thr = d - g / 2.0           # need agreement for odd r < thr
            odd = [r for r in range(1, d + 1) if r % 2 and r < thr]
            ranges = [(-1, 0, 1) if c == 'A' else ((0, 1) if c == 'B' else (-1, 0)) for c in tau]
            box = np.array(list(itertools.product(*ranges)), dtype=np.int64)
            if not odd:
                c = len(box)
            else:
                Mx = np.array([[x**r for x in S] for r in odd], dtype=np.int64).T
                _, cnt = np.unique(box @ Mx, axis=0, return_counts=True)
                c = int(cnt.max())
            if c > best[0]:
                best = (c, (S, ''.join(tau), g))
    return best

for d, Xb in ((2, 8), (3, 9), (4, 9), (5, 9)):
    b = run(d, Xb)
    print("d=%d (X-degree 2d=%d), S in [1,%d]: max cluster %d at %s" % (d, 2 * d, Xb, b[0], b[1]))
