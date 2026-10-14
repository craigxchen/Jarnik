# sections of the Paley-12 L=2 box on point subsets Y (base in Y): rho_p and verified zeta lower bound
import itertools, collections
from tclib import box_points, hilbert, cube_monomials
from zono import zeta_lower
q = 11; D = [0,1,3,4,5,9]; L = 2
cols = [sorted((i + x) % q for x in D) for i in range(q)]
res = collections.Counter()
for size in (7,):   # |Y \ base| = 7  (8 points incl. base)
    mons = cube_monomials(size)
    seen = set()
    for Y in itertools.combinations(range(q), size):
        # canonical up to translation (cyclic symmetry fixes base)
        can = min(tuple(sorted((y + t) % q for y in Y)) for t in range(q))
        if can in seen: continue
        seen.add(can)
        idx = {y: i for i, y in enumerate(Y)}
        rc = [[idx[x] for x in c if x in idx] for c in cols]
        rc = [c for c in rc if c]
        pts = box_points(size, [(c, L) for c in rc])
        H = hilbert(pts, mons=mons); rho = max(H)
        z, _ = zeta_lower(pts, rc, [L]*len(rc), size)
        print('Y=%s |A_Y|=%d rho_p=%d zeta>=%d H=%s' % (Y, len(pts), rho, z, [H.get(i,0) for i in range(rho+1)]), flush=True)
