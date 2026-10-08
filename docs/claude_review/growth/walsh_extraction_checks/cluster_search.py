"""Search actual lattice-point clusters on circles x^2+y^2=N (N squarefree product of split primes)
and record, for each cluster size M, the smallest normalised arc constant C = (arc length)/sqrt(R).
Clusters are verified exactly (integer chords) and saved with their sign vectors for profile analysis.
"""
import numpy as np, math, random, json, sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gauss_common import *

Q = math.pi / 2
SP = split_primes(60)

def angles_for(P):
    pis = [gauss_prime(p) for p in P]
    phis = np.array([math.atan2(b, a) for (a, b) in pis])
    ang = np.zeros(1)
    for ph in phis:              # bit j of index = 1 means conjugate factor at prime j
        ang = np.concatenate([ang + ph, ang - ph])
    # NB: concatenation puts the new prime as the most significant bit: index = sum s_j 2^j with j reversed below
    return pis, ang

def decode(idx, k):
    # concatenation order: after processing primes 0..k-1, index bit (k-1-j)... we rebuild explicitly
    return [(idx >> (k - 1 - j)) & 1 for j in range(k)][::-1]

def point(pis, s, unitk):
    z = UNITS[unitk % 4]
    for pi, sj in zip(pis, s):
        z = gmul(z, gconj(pi) if sj else pi)
    return z

def best_windows(P, Ms):
    k = len(P)
    pis, ang = angles_for(P)
    # verify decoding convention on a few indices
    a = np.mod(ang, Q)
    order = np.argsort(a)
    srt = a[order]
    n = len(srt)
    ext = np.concatenate([srt, srt + Q])
    logN = sum(math.log(p) for p in P)
    N14 = math.exp(logN / 4)
    res = {}
    for M in Ms:
        if M > n: continue
        widths = ext[M - 1:M - 1 + n] - ext[:n]
        i = int(np.argmin(widths))
        res[M] = (float(widths[i]) * N14, [int(order[(i + t) % n]) for t in range(M)], float(widths[i]))
    return pis, ang, res

def bits_of(index, k):
    # index built by repeated concatenation: first prime processed is the least significant "block"
    # ang[idx] = sum_j (-1)^{b_j} phi_j where b_j = bit j of idx when primes are processed 0..k-1 and
    # concatenation [ang+ph, ang-ph] makes the NEW prime the most significant bit.
    return [(index >> j) & 1 for j in range(k)]

def cluster_points(P, pis, ang, idxs):
    k = len(P)
    pts = []
    th0 = None
    for idx in idxs:
        s = bits_of(idx, k)
        th = sum((-1) ** s[j] * math.atan2(pis[j][1], pis[j][0]) for j in range(k))
        assert abs(th - ang[idx]) < 1e-9
        pts.append((s, th))
    # choose units so all angles lie in one short arc: rotate each by multiple of pi/2 toward the first
    ref = pts[0][1]
    out = []
    for s, th in pts:
        kk = round((ref - th) / Q)
        out.append((s, kk % 4, th + kk * Q))
    return out

def verify(P, pis, cl, C):
    zs = [point(pis, s, u) for (s, u, th) in cl]
    N = 1
    for p in P: N *= p
    assert all(gnorm(z) == N for z in zs)
    assert len(set(zs)) == len(zs)
    # chord^2 <= (arc)^2 <= C^2 sqrt(N); check max chord^2 <= C^2 * sqrt(N) * (1+1e-9)
    mx = max(gnorm((a[0] - b[0], a[1] - b[1])) for a in zs for b in zs)
    return mx <= C * C * math.sqrt(N) * (1 + 1e-9) + 1, zs

if __name__ == "__main__":
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    Ms = [4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 32]
    best = {M: None for M in Ms}
    sets = []
    for k in range(6, 21):
        sets.append(SP[:k])
    for trial in range(400):
        k = random.randint(8, 19)
        pool = SP[:random.randint(k, min(60, k + 25))]
        sets.append(sorted(random.sample(pool, k)))
    for P in sets:
        pis, ang, res = best_windows(P, Ms)
        for M, (C, idxs, width) in res.items():
            if best[M] is None or C < best[M][0]:
                best[M] = (C, P, idxs)
    out = {}
    for M in Ms:
        C, P, idxs = best[M]
        pis, ang, _ = best_windows(P, [M])
        cl = cluster_points(P, pis, ang, idxs)
        ok, zs = verify(P, pis, cl, C)
        out[M] = {"C": C, "P": P, "signs": [s for (s, u, th) in cl], "units": [u for (s, u, th) in cl],
                  "points": [list(z) for z in zs], "verified": ok}
        print("M=%3d  best C=%10.4f  k=%2d  N~10^%.1f  verified=%s  P=%s" % (M, C, len(P), sum(math.log10(p) for p in P), ok, P))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "best_clusters_seed%s.json" % (sys.argv[1] if len(sys.argv) > 1 else "1")), "w"))
