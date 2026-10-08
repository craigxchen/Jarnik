"""Census of tight 4-point sets (angular span <= Delta = C N^{-1/4}) on circles x^2+y^2=N, N a squarefree
product of k split primes, as k grows.  For each tight 4-set we record whether it is
  * an integer rectangle  x1+x2 = x3+x4 in Z^k (some pairing)  <=>  z1 z2 = unit * z3 z4 exactly
    (pure Hadamard-4 core with one matching absent),
  * an F_2 parallelogram (all quartic signs +1)  <=> pure Hadamard-4 core (affine 2-cube),
  * a Walsh 2-frame with flips (odd quartic support <= 4 primes, cf. item 333's model at t=2).
Exact verification of the multiplicative relation is done with Gaussian integers for every rectangle found."""
import numpy as np, math, random, sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gauss_common import *
from cluster_search import angles_for, bits_of, point, SP

Q = math.pi / 2

def census(P, C):
    k = len(P)
    pis, ang = angles_for(P)
    a = np.mod(ang, Q); order = np.argsort(a); srt = a[order]; n = len(srt)
    logN = sum(math.log(p) for p in P)
    Delta = C * math.exp(-logN / 4)
    ext = np.concatenate([srt, srt + Q]); oext = np.concatenate([order, order])
    T4 = R4 = P4 = F4 = 0
    # for each left endpoint i, points j>i within Delta; count 4-sets whose minimal element is i
    right = np.searchsorted(ext, srt + Delta, side="right")
    for i in range(n):
        cnt = right[i] - i - 1
        if cnt < 3: continue
        members = [int(oext[t]) for t in range(i, right[i])]
        first = members[0]
        for trip in itertools.combinations(members[1:], 3):
            Qd = (first,) + trip
            T4 += 1
            X = [np.array(bits_of(ix, k)) for ix in Qd]
            odd = (X[0] ^ X[1] ^ X[2] ^ X[3])
            nodd = int(odd.sum())
            if nodd <= 4: F4 += 1
            if nodd == 0:
                P4 += 1
                rect = any(np.array_equal(X[p] + X[q], X[r] + X[s]) for (p, q, r, s) in [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)])
                if rect:
                    R4 += 1
                    # exact check z_p z_q = unit z_r z_s
                    zs = [point(pis, bits_of(ix, k), 0) for ix in Qd]
                    ok = False
                    for (p, q, r, s) in [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]:
                        lhs = gmul(zs[p], zs[q]); rhs = gmul(zs[r], zs[s])
                        if any(gmul(u, rhs) == lhs for u in UNITS): ok = True
                    assert ok
    return T4, R4, P4, F4, n, logN

if __name__ == "__main__":
    random.seed(int(sys.argv[1]))
    C = float(sys.argv[2]); ks = list(range(int(sys.argv[3]), int(sys.argv[4]) + 1)); reps = int(sys.argv[5])
    print("C=%g" % C)
    for k in ks:
        tot = [0, 0, 0, 0]; logs = []
        for r in range(reps):
            P = sorted(random.sample(SP[:k + int(sys.argv[6])], k))
            T4, R4, P4, F4, n, logN = census(P, C)
            tot[0] += T4; tot[1] += R4; tot[2] += P4; tot[3] += F4; logs.append(logN / math.log(10))
        T4, R4, P4, F4 = tot
        print("k=%2d  log10N~%5.1f  tight4=%7d  rectangles=%5d (%.4f)  parallelograms=%5d (%.4f)  2-frames(odd<=4)=%6d (%.4f)"
              % (k, np.mean(logs), T4, R4, R4 / max(T4, 1), P4, P4 / max(T4, 1), F4, F4 / max(T4, 1)))
