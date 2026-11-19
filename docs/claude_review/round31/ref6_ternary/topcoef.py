"""delta_ij = deg_t (c_i^((q-1)/2) - c_j^((q-1)/2)), the t-degree of the top X-coefficient of F_i - F_j.
For L large, deg_t (F_i - F_j)(t,t^L) = ((n+h)/2) L + delta_ij."""
import sys, itertools, json
from marks import raw_mark
from f9lib import F9P, pw
def deltas(forms, signs, q):
    cs = []
    for l, e in zip(forms, signs):
        _, A, c = raw_mark(l, e, q, F9P([0, 0, 0, 0, 0, 1]))   # c does not depend on T
        cs.append(pw(c, (q-1)//2))
    return {(a, b): (cs[a] - cs[b]).deg() for a, b in itertools.combinations(range(len(forms)), 2)}, [c.deg() for c in cs]
if __name__ == '__main__':
    for forms, signs, g in [([(1,0),(0,1),(1,1),(1,2)], [1]*4, 92),
                            ([(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1)], [1]*5, 1301),
                            ([(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,1,2)], [1,1,-1,1,1], 1379),
                            ([(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,1,2),(1,2,1)], [1,1,-1,1,1,1], 833)]:
        k = len(forms[0]); q = 27; n = q**k; h = n//q
        d, cdeg = deltas(forms, signs, q)
        dm = max(d.values())
        Lmin = -(-2*dm // (g-h))
        print('k=%d m=%d g=%d h=%d deg c^13=%s max delta=%d -> L >= %d, D = %d' % (k, len(forms), g, h, cdeg, dm, Lmin, (n-g)*Lmin))
        print('   deltas', sorted(d.values()))
