# ball A_{s,D} for M=5, D=8 with a fixed starting degree window (one matrix), two primes.
import sys, time
from taulib import *
M, D = int(sys.argv[1]), int(sys.argv[2]); dm = int(sys.argv[3])
for s in [int(x) for x in sys.argv[4].split(',')]:
    pts = ball(M, s, D)
    for p in (P1, P2):
        t0 = time.time()
        T, degs = eval_matrix_mod(pts, dm, p)
        piv = pivot_columns_mod(T, p)
        h = [sum(1 for c in piv if degs[c] <= d) for d in range(dm + 1)]
        reg = next((d for d in range(dm + 1) if h[d] == len(pts)), None)
        print(f"M={M} ball D={D} s={s} |A|={len(pts)} p={p} reg={reg} ratio={(reg/D) if reg is not None else None} h={h} ({time.time()-t0:.0f}s)", flush=True)
