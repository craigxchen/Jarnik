"""For fixed m, vary L >= d: check all norms equal, and maxdiff - D/2 = Smax - L/2."""
import sys, itertools
sys.path.insert(0, __import__('os').path.dirname(__file__))
from f4arith import *
def par(x): return bin(x).count('1') & 1
for r in (2, 3, 4):
    m = 2**r; d = r+1; Q = 2**d
    forms = [w for w in range(1, Q) if w & 1]
    H = {w: [v for v in range(Q) if par(v & w) == 0 and v] for w in forms}
    Smax = (d-2)*2**(d-2)
    for L in range(d, 2*Smax + 12):
        T = 1 << L
        pts = []
        for w in forms:
            z = (1, 0)
            for v in H[w]: z = mul4(z, (T ^ v, 1))
            pts.append(z)
        D = (m-1)*L
        N = norm4(pts[0]); assert all(norm4(z) == N for z in pts) and len(set(pts)) == m
        md = max(deg4(add4(a, b)) for a, b in itertools.combinations(pts, 2))
        assert 2*md - D == 2*Smax - L, (m, L)
    print(f"m={m}: for d<=L<{2*Smax+12}: 2*maxdiff - D = 2*Smax - L  (Smax={Smax})")
