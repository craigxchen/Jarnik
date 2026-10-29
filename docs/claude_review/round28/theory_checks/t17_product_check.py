"""T17: sanity check of Theorem 2.1 (product additivity of reach and H) on small boxes, mod p and exactly (fmpz ranks).
hexagon = {k in {-1,0,1}^2 : |k_0 + k_1| <= 1}; products hex x hex, hex x [-1,1], and a 1-cut d=3 box x hexagon."""
import sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
def reach_H(A, d):
    Lp = Lpair_half(d); H = []; reach = None
    for e in range(0, 2 * d + 1):
        hs = H_modp(A, e, P1, d); H.append(hs)
        if hs == 0: break
        _, ds = shift_dims(A, Lp, e, P1, d)
        if reach is None and ds < hs: reach = e
    return reach, H
boxes = {
 'hex': (3, [((0, 1), 1)]),
 'seg': (2, []),
 'cut3': (4, [((0, 1, 2), 1)]),
 'hex x hex': (5, [((0, 1), 1), ((2, 3), 1)]),
 'hex x seg': (4, [((0, 1), 1)]),
 'cut3 x hex': (6, [((0, 1, 2), 1), ((3, 4), 1)]),
}
for name, (n, slabs) in boxes.items():
    A = build_box(n, slabs)
    r, H = reach_H(A, n - 1)
    print("%-10s d=%d |A|=%d reach(L_pair)=%s H=%s" % (name, n - 1, len(A), r, H), flush=True)
