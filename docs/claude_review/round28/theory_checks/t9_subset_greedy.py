"""T9: greedy subset search on syl8 L=2: starting from A, repeatedly (i) find the first SHIFT failure e of the current
set S w.r.t. L_pair, (ii) find power obstructions phi^e there, (iii) remove the extreme level sets of one obstruction
form phi (|phi| = max) and recompute.  Records reach(S) along the way.  Evidence (mod p)."""
import sys, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
model, L = sys.argv[1], int(sys.argv[2]); steps = int(sys.argv[3]) if len(sys.argv) > 3 else 4
p = P1
n, pats = syl(3) if model == 'syl8' else sdc(7)
d = n - 1
S = build_box(n, [(P_, L) for P_ in pats])
Lp = Lpair_half(d)
t0 = time.time()
def reach_of(S):
    H = []
    for e in range(0, 2 * d + 1):
        hs = H_modp(S, e, p, d); H.append(hs)
        if hs == 0: return e, H, e   # reached rho+1
        hs2, ds = shift_dims(S, Lp, e, p, d)
        if ds < hs: return e, H, None
    return 2 * d + 1, H, None
for step in range(steps):
    r, H, top = reach_of(S)
    rho = max(i for i, x in enumerate(H) if x) if any(H) else -1
    print("step %d: |S|=%d reach=%d (H so far %s) (%.0fs)" % (step, len(S), r, H, time.time() - t0), flush=True)
    if top is not None: print("  reached rho+1"); break
    HA, dS, obs, dim = power_obstructions(S, r, d, Lp, p)
    print("  e=%d: H=%d dim Sigma=%d, %d power obstructions spanning %d" % (r, HA, dS, len(obs), dim), flush=True)
    if not obs: print("  no power obstruction; stop"); break
    # choose the obstruction phi whose extreme levels are smallest
    best = None
    for ph in obs:
        vals = [sum(a * b for a, b in zip(ph, k)) for k in S]
        hmax = max(abs(v) for v in vals); cnt = sum(1 for v in vals if abs(v) == hmax)
        if best is None or cnt < best[1]: best = (ph, cnt, hmax)
    ph, cnt, hmax = best
    S = [k for k in S if abs(sum(a * b for a, b in zip(ph, k))) < hmax]
    print("  removing |phi|=%d for phi=%s (%d points)" % (hmax, ph, cnt), flush=True)
