"""t1: zeta (classes 'all' and 'pair') for the round-27 models: LP bound (exact rational) + MILP witness (exact check)."""
import sys, time
sys.path.insert(0, '.')
import slib, zono, models
out = []
def pr(s):
    print(s, flush=True); out.append(s)
M8, syl8 = models.syl(3)
M7, sdc7 = models.sdc(7)
for name, n, pats, L in [('sdc7', M7, sdc7, 2), ('syl8', M8, syl8, 2), ('syl8', M8, syl8, 3), ('sdc7', M7, sdc7, 3)]:
    slabs = [(S, L) for S in pats]
    A = slib.build_box(n, slabs)
    for cls in ('all', 'pair'):
        t0 = time.time()
        V, lams, bs = zono.candidates(n, slabs, cls)
        b, f, y = zono.lp_bound(V, lams, bs)
        w, res = zono.milp_zeta(V, lams, bs)
        chk = zono.verify_exact(A, w[0], w[1], n - 1)
        pr("%s L=%d |A|=%d class=%s #cand=%d LPbound=%s (%.4f) MILP opt=%s status=%s witness k0=%s gens=%s check=%s (%.0fs)" % (
            name, L, len(A), cls, len(V), b, float(b), -res.fun, res.status, w[0], w[1], chk, time.time() - t0))
open('t1_zeta_small_out.txt', 'w').write("\n".join(out) + "\n")
