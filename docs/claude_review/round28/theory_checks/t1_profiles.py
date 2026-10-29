"""T1: full defect profiles delta(e) = H_A(e) - dim sum_{l in L_pair} M_e(A_l) for syl8 / sdc7, L = 1,2,3 (mod p);
also the linkage identity H_A(e) = T_d(e) - H_B(2d-e), B = cube minus A (mod p).  Usage: t1_profiles.py model L p"""
import sys, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
model, L, p = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
n, pats = syl(3) if model == 'syl8' else sdc(7)
d = n - 1
A = build_box(n, [(S, L) for S in pats])
Aset = set(A)
B = [k for k in itertools.product((-1, 0, 1), repeat=d) if k not in Aset]
t0 = time.time()
HA = []; HB = []; delta = []
for e in range(0, 2 * d + 1):
    HA.append(H_modp(A, e, p, d))
for e in range(0, 2 * d + 1):
    HB.append(H_modp(B, e, p, d) if B else 0)
rho = max(e for e in range(2 * d + 1) if HA[e] > 0)
print(model, "L=%d |A|=%d |B|=%d rho=%d" % (L, len(A), len(B), rho), flush=True)
print("H_A =", HA, flush=True)
print("H_B =", HB, flush=True)
print("linkage H_A(e) == T(e) - H_B(2d-e) for all e:", all(HA[e] == T(d, e) - HB[2 * d - e] for e in range(2 * d + 1)), flush=True)
Lp = Lpair_half(d)
for e in range(0, rho + 1):
    acc = shift_sum_modp(A, Lp, e, p, d, stop_at=HA[e])
    s = 0 if acc is None else acc.ncols()
    delta.append(HA[e] - s)
    print("  e=%d H_A=%d dim sum=%d delta=%d (%.0fs)" % (e, HA[e], s, HA[e] - s, time.time() - t0), flush=True)
print("delta profile:", delta)
reach = next((e for e in range(rho + 1) if delta[e] > 0), rho + 1)
print("reach(L_pair) =", reach, " rho =", rho)
