"""Referee's independent re-scan of the author's Section-4 data set (own code, own method).

Definition reproduced from L1L2.md Section 4: every N <= Nmax all of whose prime factors are 1 mod 4 and
tau(N) >= T; every run of >= Mmin angularly consecutive lattice points of x^2+y^2=N with
(theta_last - theta_first) N^(1/4) <= Cmax; keep runs whose points have Gaussian gcd 1; deduplicate under
the 8 symmetries of Z[i].

Method differs from the author's: points generated with numpy int64 from Gaussian primes; runs are
detected vectorially; deduplication keeps only runs whose (circular) mean angle lies in [0, pi/4] and then
removes exact duplicates by a canonical key (only needed for runs symmetric about an axis/diagonal).

Usage: python3 scan_ref.py Nmax T Cmax Mmin nproc out.pkl
"""
import sys, pickle, math
import numpy as np
from multiprocessing import Pool
from refgauss import gauss_prime, ggcd, gnorm

Nmax = int(float(sys.argv[1])); T = int(sys.argv[2]); Cmax = float(sys.argv[3]); Mmin = int(sys.argv[4])
NP = int(sys.argv[5]); OUT = sys.argv[6]


def primes1mod4(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    ps = np.nonzero(s)[0]
    return [int(p) for p in ps if p % 4 == 1]


# minimal N with tau >= 8 built from primes 1 mod 4 is 5*13*17 = 1105 (or 5^3*13=1625); so the largest
# prime (exponent 1) of an admissible N is <= Nmax/1105.
PR = primes1mod4(max(Nmax // 1105 + 1, 100))


def gen_N():
    """All (N, [(p, e), ...]) with N <= Nmax, primes 1 mod 4, tau >= T. Iterative DFS."""
    out = []
    stack = [(0, 1, (), 1)]
    while stack:
        idx, cur, fac, tau = stack.pop()
        if tau >= T:
            out.append((cur, fac))
        for k in range(idx, len(PR)):
            p = PR[k]
            if cur * p > Nmax:
                break
            # upper bound on the further tau multiplier: primes >= p, product <= Nmax/cur
            rem = Nmax // cur
            L = int(math.log(rem) / math.log(p) + 1e-12)
            if tau * (2 ** L) < T:
                break
            q = cur; e = 0
            while q * p <= Nmax:
                q *= p; e += 1
                stack.append((k + 1, q, fac + ((p, e),), tau * (e + 1)))
    return out


GP = {}


def gp(p):
    if p not in GP:
        GP[p] = gauss_prime(p)
    return GP[p]


def points(fac):
    re = np.array([1], dtype=np.int64); im = np.array([0], dtype=np.int64)
    for p, e in fac:
        a, b = gp(p)
        # powers pi^k, conj(pi)^k
        pw = [(1, 0)]
        for _ in range(e):
            x, y = pw[-1]; pw.append((x * a - y * b, x * b + y * a))
        cands = []
        for k in range(e + 1):
            x1, y1 = pw[k]; x2, y2 = pw[e - k]; y2 = -y2  # conj(pi)^(e-k) = conj(pi^(e-k))
            cands.append((x1 * x2 - y1 * y2, x1 * y2 + y1 * x2))
        nre = []; nim = []
        for (x, y) in cands:
            nre.append(re * x - im * y); nim.append(re * y + im * x)
        re = np.concatenate(nre); im = np.concatenate(nim)
    # four units
    re4 = np.concatenate([re, -im, -re, im]); im4 = np.concatenate([im, re, -im, -re])
    return re4, im4


def canon(cl):
    best = None
    for c in (1, -1):
        for r in range(4):
            Q = []
            for (x, y) in cl:
                y = c * y
                for _ in range(r):
                    x, y = -y, x
                Q.append((x, y))
            k = tuple(sorted(Q))
            if best is None or k < best:
                best = k
    return best


def work(chunk):
    res = []
    for N, fac in chunk:
        re, im = points(fac)
        ang = np.mod(np.arctan2(im, re), 2 * np.pi)
        o = np.argsort(ang, kind='stable')
        ang = ang[o]; re = re[o]; im = im[o]
        n = len(ang)
        s = N ** 0.25
        # circular differences ang[(i+Mmin-1)%n] - ang[i]
        ext = np.concatenate([ang, ang[:Mmin] + 2 * np.pi])
        span = (ext[Mmin - 1:Mmin - 1 + n] - ext[:n]) * s
        starts = np.nonzero(span <= Cmax)[0]
        if len(starts) == 0:
            continue
        seen = set()
        for i in starts:
            i = int(i)
            # maximal run from i
            m = Mmin
            while m < n and ((ang[(i + m) % n] - ang[i]) % (2 * np.pi)) * s <= Cmax:
                m += 1
            for mm in range(Mmin, m + 1):
                idx = [(i + t) % n for t in range(mm)]
                cl = [(int(re[j]), int(im[j])) for j in idx]
                # circular mean angle of the run
                a0 = ang[i]
                rel = [((ang[j] - a0) % (2 * np.pi)) for j in idx]
                mid = (a0 + (rel[0] + rel[-1]) / 2) % (2 * np.pi)
                if not (-1e-12 <= mid <= np.pi / 4 + 1e-12):
                    continue
                g = cl[0]
                for z in cl[1:]:
                    g = ggcd(g, z)
                if gnorm(g) != 1:
                    continue
                key = canon(cl)
                if key in seen:
                    continue
                seen.add(key)
                C = (rel[-1]) * s
                res.append((N, fac, cl, float(C)))
    return res


if __name__ == "__main__":
    Ns = gen_N()
    print("admissible N:", len(Ns), flush=True)
    chunks = [Ns[i::NP * 50] for i in range(NP * 50)]
    allres = []
    with Pool(NP) as pool:
        for r in pool.imap_unordered(work, chunks):
            allres.extend(r)
    # global dedup (keys are per N anyway)
    d = {}
    for t in allres:
        d[canon(t[2])] = t
    allres = list(d.values())
    from collections import Counter
    print("clusters:", len(allres), "by M:", sorted(Counter(len(t[2]) for t in allres).items()), flush=True)
    pickle.dump(allres, open(OUT, "wb"))
