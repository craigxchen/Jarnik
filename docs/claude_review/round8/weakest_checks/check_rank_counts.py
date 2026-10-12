"""Evidence: how often do clusters of rank k (k+1 points or more) with normalised constant C <= CMAX occur,
per decade of N?  Same enumeration and normalisation as check_slw_data.py (all N <= X with prime factors
= 1 mod 4 and >= MINPTS lattice points; windows of consecutive points; cluster normalised by its own gcd).
For every circle we record, for each rank k, whether some window of rank >= k has C <= CMAX.
Usage: python3 check_rank_counts.py X
"""
import sys
import math
import numpy as np
from check_slw_data import enumerate_N, points_of

X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 7
MINPTS = 24
CMAX = 6.0


def main():
    counts = {}
    nN = {}
    for N, fac in enumerate_N(X):
        npts = 4
        for e in fac.values():
            npts *= (e + 1)
        if npts < MINPTS:
            continue
        dec = int(math.log10(N))
        nN[dec] = nN.get(dec, 0) + 1
        primes, pts, logH = points_of(fac)
        ang = sorted(((math.atan2(z[1], z[0]) % (2 * math.pi), z, vec) for z, vec in pts))
        L = len(ang)
        best_rank = 0
        for s in range(L):
            for t in range(2, min(L, 40)):
                win = [ang[(s + u) % L] for u in range(t)]
                zs = [w[1] for w in win]
                vecs = [w[2] for w in win]
                Rn2 = N
                for j, p in enumerate(primes):
                    col = [v[j] for v in vecs]
                    e = fac[p]
                    Rn2 //= p ** (min(col) + min(e - c for c in col))
                chord = max(math.hypot(z[0] - w[0], z[1] - w[1]) for z in zs for w in zs)
                Cn = (chord / math.sqrt(N / Rn2)) / math.sqrt(math.sqrt(Rn2))
                if Cn > CMAX:
                    break
                if t - 1 > best_rank:
                    D = np.array([[v[j] - vecs[0][j] for j in range(len(primes))] for v in vecs[1:]], dtype=float)
                    best_rank = max(best_rank, int(np.linalg.matrix_rank(D)))
        for k in range(1, best_rank + 1):
            counts[(dec, k)] = counts.get((dec, k), 0) + 1
    print(f"X = {X:.3g}, circles with >= {MINPTS} points, clusters with normalised C <= {CMAX}")
    print("decade  circles   rank>=3   rank>=4   rank>=5   rank>=6")
    for dec in sorted(nN):
        row = [counts.get((dec, k), 0) for k in (3, 4, 5, 6)]
        print(f"10^{dec:<3} {nN[dec]:>8}  " + "  ".join(f"{c:>8}" for c in row))


if __name__ == '__main__':
    main()
