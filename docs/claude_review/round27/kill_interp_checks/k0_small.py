"""K0: cube d=2,3,4 (Prop 3.2), hexagon; exact over Q."""
import itertools, sys
sys.path.insert(0, '.')
from klib import *
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
for d in (2, 3, 4):
    G = list(itertools.product((-1, 0, 1), repeat=d))
    H = hilbert(G, 2 * d, exact=True)
    pr("cube d=%d |G|=%d H=%s" % (d, len(G), H))
    for nm, L in (("coord", Lcoord(d)), ("pair", Lpair(d))):
        res = []
        for e in range(2 * d + 1):
            hs, sm = shift_dims_Q(G, L, e)
            res.append((e, hs, sm, hs == sm))
        firstfail = [r[0] for r in res if not r[3]]
        pr("  L=%s: (e,H,dimsum,ok) %s ; first failure e=%s" % (nm, res, firstfail[:1]))
# hexagon: n=3, slab S={0,1}, L=1  -> |k1+k2|<=1
A = build_box(3, [((0, 1), 1)])
pr("hexagon |A|=%d pts=%s H=%s" % (len(A), A, hilbert(A, 4, exact=True)))
for nm, L in (("coord", Lcoord(2)), ("pair", Lpair(2))):
    pr("  L=%s:" % nm, [(e,) + shift_dims_Q(A, L, e) for e in range(5)])
# V(A,3) for hexagon: exact kernel
lt = monos_lt(2, 3)
import flint
M = flint.fmpz_mat([[mval(k, a) for k in A] for a in lt])
X, nul = M.nullspace()
pr("  dim V(hex,3) =", nul, "vector:", [(A[i], int(X[i, 0])) for i in range(len(A)) if X[i, 0] != 0])
open("k0_small_out.txt", "w").write("\n".join(out) + "\n")
