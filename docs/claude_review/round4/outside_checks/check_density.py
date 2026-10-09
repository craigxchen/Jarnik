"""Exact count for Prop D.1 of outside.md.

E_k(X, C) = #{ N <= X : some arc of length C N^(1/4) of x^2+y^2=N contains >= k lattice points }.
Prop D.1 (proved in the text): E_2(X, C) <= 16 C X^(3/4) (1 + log(C X^(1/4))) + 20 C^2 X^(1/2).
We compute exactly (integer arithmetic) the CHORD version E'_k: some k consecutive points whose
extreme chord satisfies |z-w|^4 <= C^4 N.  Since chord <= arc, E_k <= E'_k, and the proof of
Prop D.1 bounds E'_2 as well, so the assertion E'_2 <= bound is a valid test of the proof.
Points are generated in the first quadrant (x > 0, y >= 0); rotation by i closes the circle.
"""
import numpy as np
from math import log, pi
import sys


def run(X, Cs=(1.0, 2.0)):
    xs = []
    ys = []
    lim = int(X ** 0.5)
    for x in range(1, lim + 1):
        ymax = int((X - x * x) ** 0.5)
        while ymax * ymax > X - x * x:
            ymax -= 1
        while (ymax + 1) ** 2 <= X - x * x:
            ymax += 1
        y = np.arange(0, ymax + 1, dtype=np.int64)
        xs.append(np.full(len(y), x, dtype=np.int64))
        ys.append(y)
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    n = x * x + y * y
    ang = np.arctan2(y.astype(float), x.astype(float))
    order = np.lexsort((ang, n))
    x, y, n, ang = x[order], y[order], n[order], ang[order]
    # group boundaries
    starts = np.flatnonzero(np.r_[True, n[1:] != n[:-1]])
    ends = np.r_[starts[1:], len(n)]
    res = {}
    for C in Cs:
        E = {2: 0, 3: 0, 4: 0, 5: 0}
        E2arc = 0
        for s, e in zip(starts, ends):
            m = e - s
            if m < 1:
                continue
            N = int(n[s])
            # full circle: quadrant points then their rotations by i (angles + pi/2)
            px = list(x[s:e]) + [-int(t) for t in y[s:e]]
            py = list(y[s:e]) + [int(t) for t in x[s:e]]
            # angles in [0, pi): consecutive window, then wrap with rotation by i^2 (pi)
            pts = list(zip(px, py))
            pts = pts + [(-a, -b) for a, b in pts]
            L = len(pts)
            if L < 2:
                continue
            best = 1
            # sliding window over the cyclic sequence (angles increasing)
            j = 0
            C4N = (C ** 4) * N
            for i in range(L):
                if j < i:
                    j = i
                while j + 1 < i + L:
                    a = pts[i]
                    b = pts[(j + 1) % L]
                    d2 = (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                    if d2 * d2 <= C4N:
                        j += 1
                    else:
                        break
                best = max(best, j - i + 1)
            for k in E:
                if best >= k:
                    E[k] += 1
        res[C] = E
    return res


def main():
    X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 6
    res = run(X)
    for C, E in res.items():
        bound = 16 * C * X ** 0.75 * (1 + log(C * X ** 0.25)) + 20 * C * C * X ** 0.5
        print("X=%.0e C=%.1f: E2=%d (Prop D.1 bound %.3e; E2/X^(3/4)=%.3f)  E3=%d (E3/X^(1/2)=%.3f)  E4=%d  E5=%d"
              % (X, C, E[2], bound, E[2] / X ** 0.75, E[3], E[3] / X ** 0.5, E[4], E[5]))
        assert E[2] <= bound


if __name__ == "__main__":
    main()
