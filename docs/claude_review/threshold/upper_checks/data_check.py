"""Compare exact-modular reg(Q^M(p,n)) and reg(S^M(p,n)) with the proven upper bound floor(phi_M)
and the proven (constructive) lower bound L_M.  reg mod q is a rigorous UPPER bound for reg over Q;
L_M is a rigorous LOWER bound (products of antisymmetrised orbits and linear factors)."""
import sys, time, math
sys.setrecursionlimit(100000)
from hilb import hilbert, ball, shell, Q1, Q2
from peel import B, L
from fractions import Fraction as F
def phi(M, p, n):
    return min(F(M - 1) * (F(p, M - j) + F(n, j)) for j in range(1, M))
M = int(sys.argv[1]); tmax = int(sys.argv[2]); kinds = sys.argv[3].split(',')
for t in range(1, tmax + 1):
    for n in range(0, t // 2 + 1):
        p = t - n
        for kind in kinds:
            pts = (ball if kind == 'ball' else shell)(M, p, n)
            t0 = time.time()
            h1 = hilbert(pts, Q1); h2 = hilbert(pts, Q2)
            r1, r2 = len(h1) - 1, len(h2) - 1
            ph = phi(M, p, n); fl = math.floor(ph)
            lo = L(M, p, n) if kind == 'ball' else None
            ok = (r1 <= fl and r2 <= fl and (lo is None or lo <= min(r1, r2)))
            print(f"M={M} {kind:5s} (p,n)=({p},{n}) |A|={len(pts):6d} reg_q1={r1} reg_q2={r2} "
                  f"L={lo} B={B(M,p,n)} floor(phi)={fl} phi={float(ph):.3f} {'OK' if ok else 'FAIL'} "
                  f"[{time.time()-t0:.1f}s]", flush=True)
