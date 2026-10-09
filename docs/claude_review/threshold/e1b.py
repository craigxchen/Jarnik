import sys, time
from taulib import *
M = int(sys.argv[1]); D = int(sys.argv[2]); mode = sys.argv[3]
for s in range(D % 2, D + 1, 2):
    pts = shell(M, s, D) if mode == 'shell' else ball(M, s, D)
    t0 = time.time()
    r1, h1 = reg_mod(pts, P1)
    r2, h2 = reg_mod(pts, P2)
    flag = '' if (r1 == r2 and h1 == h2) else ' PRIME-MISMATCH'
    print(f"M={M} {mode} D={D} s={s} |A|={len(pts)} reg={r1} ratio={r1/D:.4f}{flag} ({time.time()-t0:.1f}s)", flush=True)
