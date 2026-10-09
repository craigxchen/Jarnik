# E1: recompute shell values reg(S_{s,t}) and ball values reg(A_{s,D}) mod two primes.
# reg_p >= reg_Q, so reg_p/t is a rigorous UPPER bound on the rational value; agreement of
# two primes is reported. Lower bounds are certified separately (e2).
import sys, time
from taulib import *
M = int(sys.argv[1]); Dmax = int(sys.argv[2]); mode = sys.argv[3]  # 'shell' or 'ball'
best = (0,)
for D in range(1, Dmax + 1):
    for s in range(0, D + 1):
        if (D - s) % 2: continue
        pts = shell(M, s, D) if mode == 'shell' else ball(M, s, D)
        if not pts: continue
        t0 = time.time()
        r1, h1 = reg_mod(pts, P1)
        r2, h2 = reg_mod(pts, P2)
        flag = '' if (r1 == r2 and h1 == h2) else ' PRIME-MISMATCH'
        print(f"M={M} {mode} D={D} s={s} |A|={len(pts)} reg={r1} ratio={r1/D:.4f} h={h1}{flag} ({time.time()-t0:.1f}s)", flush=True)
        if r1 / D > best[0]: best = (r1 / D, D, s, r1)
print("BEST", M, mode, best, flush=True)
