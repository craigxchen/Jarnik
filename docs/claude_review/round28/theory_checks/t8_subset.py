"""T8: subset monotonicity.  reach(A', L_pair) for sub-boxes A' of a model box (any A' subset A certifies
Bs V(A,m) subset Bs V(A',m)).  A' = A cap {|phi| <= h_phi - t} for a linear form phi.  Usage: t8_subset.py model L phi t"""
import sys, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
model, L, phis, t = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
p = P1
n, pats = syl(3) if model == 'syl8' else sdc(7)
d = n - 1
A = build_box(n, [(S, L) for S in pats])
phi = tuple(int(x) for x in phis.split(','))
h = max(sum(a * b for a, b in zip(phi, k)) for k in A)
Ap = [k for k in A if abs(sum(a * b for a, b in zip(phi, k))) <= h - t]
print(model, "L=%d |A|=%d phi=%s h=%d  A' = |phi|<=%d : |A'|=%d" % (L, len(A), phi, h, h - t, len(Ap)), flush=True)
Lp = Lpair_half(d)
t0 = time.time()
H = []; delta = []
for e in range(0, 2 * d + 1):
    hs = H_modp(Ap, e, p, d); H.append(hs)
    if hs == 0: break
    acc = shift_sum_modp(Ap, Lp, e, p, d, stop_at=hs)
    s = 0 if acc is None else acc.ncols()
    delta.append(hs - s)
    print("  e=%d H_A'=%d dim sum=%d delta=%d (%.0fs)" % (e, hs, s, hs - s, time.time() - t0), flush=True)
rho = len([x for x in H if x]) - 1
reach = next((e for e in range(len(delta)) if delta[e] > 0), rho + 1)
print("H_A' =", H, " rho(A') =", rho, " reach(A', L_pair) =", reach)
