"""Referee recount of Prop D.1 data (outside.md Section 5), independent code.

For every N <= X, sort the lattice points of x^2+y^2=N by angle; N is counted in E_k(C) if some k
points, consecutive in angle, have all pairwise chords <= C N^(1/4) -- equivalently the extreme pair
(chord of the spanned arc) has squared chord <= C^2 N^(1/2).  Exact integer test: chord^2 <= C^2 sqrt N
<=> chord^4 <= C^4 N (C rational: C = 1, 2).  Also evaluates the bound of Prop D.1.
"""
import numpy as np
from math import log, isqrt


def run(X, Cs):
    R = isqrt(X)
    xs = np.arange(-R, R + 1, dtype=np.int64)
    gx, gy = np.meshgrid(xs, xs, indexing="ij")
    gx = gx.ravel()
    gy = gy.ravel()
    N = gx * gx + gy * gy
    keep = (N <= X) & (N > 0)
    gx, gy, N = gx[keep], gy[keep], N[keep]
    ang = np.arctan2(gy.astype(float), gx.astype(float))
    order = np.lexsort((ang, N))
    gx, gy, N = gx[order], gy[order], N[order]
    starts = np.flatnonzero(np.r_[True, N[1:] != N[:-1]])
    ends = np.r_[starts[1:], len(N)]
    for C in Cs:
        C4 = C ** 4
        counts = {2: 0, 3: 0, 4: 0, 5: 0}
        for s, e in zip(starts, ends):
            m = e - s
            if m < 2:
                continue
            n = int(N[s])
            px = gx[s:e]
            py = gy[s:e]
            best = 1
            # for each k, sliding window of k consecutive points (cyclic)
            for k in (2, 3, 4, 5):
                if k > m:
                    break
                found = False
                for i in range(m):
                    j = (i + k - 1) % m
                    dx = int(px[i] - px[j])
                    dy = int(py[i] - py[j])
                    d2 = dx * dx + dy * dy
                    if d2 * d2 <= C4 * n:
                        # check the arc i..j is the short way (k-1 steps): for k < m always the span
                        found = True
                        break
                if found:
                    counts[k] += 1
                else:
                    break
        bound = 16 * C * X ** 0.75 * (1 + log(C * X ** 0.25)) + 14 * C * C * X ** 0.5
        print("X=%.0e C=%g: E2=%d E3=%d E4=%d E5=%d   (Prop D.1 bound %.3e; E2/X^(3/4)=%.3f)"
              % (X, C, counts[2], counts[3], counts[4], counts[5], bound, counts[2] / X ** 0.75))


run(10 ** 6, (1, 2))
