"""Transfer check (section 8): maximal order at e of sections of a T-invariant divisor sum_A d_A E_A on the
Losev-Manin 4-fold (rays 1_A, A nonempty proper subset of [5]); sections = Laurent polynomials with exponents m in
P = {m in Z^5, sum m = 0 : sum_{i in A} m_i >= -d_A}. A section vanishes to order >= r at e=(1,...,1) iff
sum_m c_m Q(m) = 0 for all polynomials Q of degree < r. Calibration: d = 1 (=-K), k=3 must give order 10.
Test: L' = 2(-K) - B_{1,4}  (d_A = 2 for |A| in {2,3}, 1 for |A| in {1,4}) = K + sum_{a<b} Z_ab.
Log general type of (Bl_e X, sum Z_ab~ + E) at the first blow-up needs tau_e(L') > 6 (= 10 - 4)."""
import itertools, sys
import numpy as np
from fractions import Fraction
import sympy as sp

subsets = [A for r in range(1, 5) for A in itertools.combinations(range(5), r)]

def lattice_points(dfun, k):
    pts = []
    R = 3*k + 2
    for m4 in itertools.product(range(-R, R + 1), repeat=4):
        m = list(m4) + [-sum(m4)]
        if all(sum(m[i] for i in A) >= -k*dfun(A) for A in subsets):
            pts.append(m[:4])
    return pts

def max_order(pts, rmax=40):
    pts = np.array(pts, dtype=object)
    n = len(pts)
    best = 0
    for r in range(1, rmax):
        monos = [al for d in range(r) for al in itertools.product(range(d + 1), repeat=4) if sum(al) == d]
        if len(monos) > 9000:
            print("cap reached at r=%d" % r)
            break
        import flint
        M = flint.fmpz_mat([[int(np.prod([int(p[i])**al[i] for i in range(4)])) for p in pts] for al in monos])
        rk = M.rank()
        if rk < n:
            best = r          # nonzero c annihilating all degree < r
        else:
            break
    return best

if __name__ == '__main__':
    which = sys.argv[1]; ks = [int(a) for a in sys.argv[2:]]
    if which == 'K':
        dfun = lambda A: 1
    else:
        dfun = lambda A: 2 if len(A) in (2, 3) else 1
    for k in ks:
        pts = lattice_points(dfun, k)
        o = max_order(pts)
        print(which, 'k=%d: %d lattice points, max order at e = %d, ratio %.4f' % (k, len(pts), o, o/k), flush=True)
