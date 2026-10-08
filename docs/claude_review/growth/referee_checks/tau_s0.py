"""dstar(A_s) = max vanishing order at w=1 of a Laurent polynomial supported on
A_s = {e' in Z^n : |s - sum e'| + |e'|_1 <= N}  (n = M-1), i.e. the max order along Y
of a weight-s section of O(N) on X_M.  tau_N(M) = max_s dstar(A_s)/N <= tau(H,E)."""
import sys, time
from beta_vanishing import points, hilbert_T
for arg in sys.argv[1:]:
    n, N, s = map(int, arg.split(','))
    t = time.time()
    pts = points(n, N, s)
    _, hist = hilbert_T(pts, 1000003)
    d1 = len(hist) - 1
    print(f"M={n+1} N={N} s={s}: |A_s|={len(pts)} dstar={d1} ratio={d1/N:.4f}  hist={hist}  ({time.time()-t:.1f}s)", flush=True)
