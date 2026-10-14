# exact Hilbert functions of small Hadamard T=1 boxes at L=1 (Paley 12, Sylvester 16) and twin check for Z6
import itertools, numpy as np
from tclib import box_points
from certify import hilbert_exact
def cyc_pts(q, D, L):
    cols = [sorted((i + x) % q for x in D) for i in range(q)]
    return box_points(q, [(c, L) for c in cols]), cols
# Paley 12
pts, _ = cyc_pts(11, [0,1,3,4,5,9], 1)
H = hilbert_exact(pts); print('Paley12 L=1 |A|=%d H=%s rho=%d' % (len(pts), [H[i] for i in sorted(H)], max(H)))
# Sylvester 16 in Singer coordinates, D = trace-one exponents
D = [3, 6, 7, 9, 11, 12, 13, 14]
cols = [sorted((i + t) % 15 for i in D) for t in range(15)]
# enumerate points with |support| small (ball effect): ||k||^2 <= 3 incl. base
pts = []
for s in range(0, 4):
    for X in itertools.combinations(range(15), s):
        for signs in itertools.product((-1, 1), repeat=s):
            k = [0]*15
            for x, e in zip(X, signs): k[x] = e
            if all(abs(sum(k[x] for x in c)) <= 1 for c in cols): pts.append(tuple(k))
H = hilbert_exact(pts); print('Syl16 L=1 |A|=%d H=%s rho=%d' % (len(pts), [H[i] for i in sorted(H)], max(H)))
# Z6 twin check
cols6 = [[1,2,3,4,5],[0,1,2,3,4,5],[0,3,4],[1,4,5],[1,2,3,4],[0,1,2,3,4,5],[0,2,4,5]]
pats = [tuple(1 if x in c else 0 for c in cols6) for x in range(6)] + [tuple(0 for c in cols6)]
print('Z6 patterns distinct:', len(set(pats)) == 7, pats)
