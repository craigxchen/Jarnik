"""Sanity check of the height normalisation in hypothesis (H_kappa) (Section 4 of flipped.md).

For conj-primitive Gaussian integers G, P (odd norm, nonunit) with Norm(G) <= X1, Norm(P) <= X2 and an
even M, let d(G,P) = dist(M arg G - 2 arg P, (pi/2)Z).  The counting heuristic predicts
   min_{G,P} d(G,P) * N(G) * N(P)  ~  1/(log X1 log X2)   (exponent 1 on naive heights, up to logs),
whereas an exponent-1/2 normalisation (absolute heights) would require min d*sqrt(N(G)N(P)) >~ 1,
which the data refute.  Floating point (double) angles; d is far above the rounding error.
This is evidence about the normalisation only, not about the truth of (H_kappa).
"""
import math
import numpy as np


def conj_primitive_reps(X):
    out = []
    for a in range(1, math.isqrt(X) + 1):
        for b in range(0, a + 1):
            n = a * a + b * b
            if n > X or n == 1 or n % 2 == 0:
                continue
            if math.gcd(a, b) != 1:
                continue
            # conj-primitive & odd norm: gcd(a,b)=1 and n odd suffices (no 1+i, no rational prime factor)
            out.append((a, b, n))
            if b != 0 and b != a:
                out.append((b, a, n))      # the other orientation (conjugate up to unit)
    return out


def run(M, X1, X2):
    G = conj_primitive_reps(X1)
    P = conj_primitive_reps(X2)
    argG = np.array([math.atan2(b, a) for a, b, n in G])
    nG = np.array([n for a, b, n in G], dtype=float)
    argP = np.array([math.atan2(b, a) for a, b, n in P])
    nP = np.array([n for a, b, n in P], dtype=float)
    q = math.pi / 2
    best1 = best2 = float("inf")
    excluded = [0]
    for i in range(0, len(P), 256):
        t = M * argG[None, :] - 2 * argP[i:i + 256, None]
        d = np.abs(t - q * np.round(t / q))
        # Lambda = 0 cases (G^M conj(P)^2 exactly on an axis, e.g. P = unit*G^(M/2)) are excluded by
        # hypothesis; detect them exactly and remove them
        for (r, c) in zip(*np.nonzero(d < 1e-9)):
            a1, b1, _ = G[c]
            a2, b2, _ = P[i + r]
            z = (1, 0)
            for _ in range(M):
                z = (z[0] * a1 - z[1] * b1, z[0] * b1 + z[1] * a1)
            for _ in range(2):
                z = (z[0] * a2 + z[1] * b2, -z[0] * b2 + z[1] * a2)
            assert z[0] == 0 or z[1] == 0, "tiny but nonzero form: increase precision"
            d[r, c] = np.inf
            excluded[0] += 1
        w = nG[None, :] * nP[i:i + 256, None]
        best1 = min(best1, float((d * w).min()))
        best2 = min(best2, float((d * np.sqrt(w)).min()))
    return len(G), len(P), best1, best2, excluded[0]


def main():
    for M in [12, 20]:
        for X in [500, 2000, 8000, 20000]:
            nG, nP, b1, b2, ex = run(M, X, X)
            print(f"M={M:3d} X1=X2={X:6d}  #G={nG:6d} #P={nP:6d}  min d*N(G)N(P)={b1:9.4f}"
                  f"  [x log(X)^2 = {b1*math.log(X)**2:7.3f}]   min d*sqrt(N(G)N(P))={b2:.2e}  (exact Lambda=0 pairs removed: {ex})")
    print("DONE")


if __name__ == "__main__":
    main()
