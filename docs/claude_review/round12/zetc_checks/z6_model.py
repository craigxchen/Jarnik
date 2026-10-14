# the d=6 zonotope-free candidate from the random search: exact certification
import itertools
from tclib import box_points
from certify import hilbert_exact, zeta_certificate
from zono import zeta_lower
cols = [[1,2,3,4,5],[0,1,2,3,4,5],[0,3,4],[1,4,5],[1,2,3,4],[0,1,2,3,4,5],[0,2,4,5]]
Ls = [2,2,1,2,2,1,1]
d = 6
pts = box_points(d, list(zip(cols, Ls)))
H = hilbert_exact(pts)
print('|A|=%d exact H=%s rho=%d' % (len(pts), [H[i] for i in sorted(H)], max(H)))
funcs = [tuple(1 if x == y else 0 for x in range(d)) for y in range(d)] + [tuple(1 if x in S else 0 for x in range(d)) for S in cols]
for x, y in itertools.combinations(range(d), 2):
    f = [0]*d; f[x] = 1; f[y] = 1; funcs.append(tuple(f))
    f = [0]*d; f[x] = 1; f[y] = -1; funcs.append(tuple(f))
lpval, cert, nd = zeta_certificate(pts, funcs)
print('LP value %.6f, exact certificate value %s (= %.6f), directions %d' % (lpval, cert[0], float(cert[0]), nd))
z, edges = zeta_lower(pts, cols, Ls, d)
print('verified zonotope of length', z, edges)
