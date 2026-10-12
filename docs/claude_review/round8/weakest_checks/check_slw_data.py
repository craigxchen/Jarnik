"""Evidence (not proof): SLW exponents gamma_k of actual lattice-point configurations.

For every N <= X whose prime factors are all = 1 mod 4 and which has at least MINPTS lattice
points on x^2 + y^2 = N, sort the lattice points by angle and, for every window of k+1
consecutive points (k = 3, 4, 5), compute
   rank  = rank of the exponent differences (prime coordinates),
   lam   = min over base points of max_i |angle(z_i) - angle(z_base)|,
   gamma = log(1/lam) / sum_{p in U} log H(pi_p/conj pi_p),  U = primes where exponents vary.
SLW(k, g) predicts gamma <= g + o(1) for windows of rank k.  The Dirichlet heuristic predicts
typical gamma ~ 1/k at best.  We report the largest gamma per (k, rank = k), and, for every
N, the size and rank of the largest cluster on an arc of normalised constant <= CMAX, where the
cluster is normalised by its own Gaussian gcd (C = chord-diameter / sqrt(R'), R' = |z/g|).

Usage: python3 check_slw_data.py X
"""
import sys
import math
import numpy as np
import sympy

X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 8
MINPTS = 48
CMAX = 6.0

split = [p for p in sympy.primerange(5, X // 5 + 1) if p % 4 == 1] if X <= 10 ** 7 else None


def two_sq(p):
    for c in range(2, p):
        t = pow(c, (p - 1) // 4, p)
        if (t * t) % p == p - 1:
            break
    a, b = p, t
    lim = math.isqrt(p)
    while b > lim:
        a, b = b, a % b
    x = b
    y = math.isqrt(p - x * x)
    assert x * x + y * y == p
    return (x, y)


def gmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def points_of(fac):
    """All lattice points with exponent vectors, for N = prod p^e (all p = 1 mod 4)."""
    primes = sorted(fac)
    pis = [two_sq(p) for p in primes]
    pts = [((1, 0), ())]
    for p, pi in zip(primes, pis):
        e = fac[p]
        conj = (pi[0], -pi[1])
        powers = []
        for a in range(e + 1):
            z = (1, 0)
            for _ in range(a):
                z = gmul(z, pi)
            for _ in range(e - a):
                z = gmul(z, conj)
            powers.append((a, z))
        pts = [(gmul(z, w), vec + (a,)) for (z, vec) in pts for (a, w) in powers]
    out = []
    for z, vec in pts:
        u = z
        for _ in range(4):
            out.append((u, vec))
            u = (-u[1], u[0])
    logH = []
    for p, pi in zip(primes, pis):
        re2 = pi[0] * pi[0] - pi[1] * pi[1]
        logH.append(math.log(max(p, 2 * abs(re2))))
    return primes, out, logH


def enumerate_N(X):
    ps = [p for p in sympy.primerange(5, X // 5 + 1) if p % 4 == 1]
    res = []

    def rec(start, n, fac):
        res.append((n, dict(fac)))
        for i in range(start, len(ps)):
            p = ps[i]
            if n * p > X:
                break
            m, e = n * p, 1
            while m <= X:
                fac[p] = e
                rec(i + 1, m, fac)
                del fac[p]
                m *= p
                e += 1
    rec(0, 1, {})
    return res


def main():
    allN = enumerate_N(X)
    best = {k: (0.0, None) for k in (3, 4, 5)}
    hist = {k: [] for k in (3, 4, 5)}
    rank_of_big = {}
    nN = 0
    for N, fac in allN:
        npts = 4
        for e in fac.values():
            npts *= (e + 1)
        if npts < MINPTS:
            continue
        nN += 1
        primes, pts, logH = points_of(fac)
        ang = sorted(((math.atan2(z[1], z[0]) % (2 * math.pi), z, vec) for z, vec in pts))
        L = len(ang)
        R = math.sqrt(math.sqrt(N))  # = sqrt(radius)
        for k in (3, 4, 5):
            for s in range(L):
                win = [ang[(s + t) % L] for t in range(k + 1)]
                th = [w[0] for w in win]
                th = [th[0]] + [t if t >= th[0] else t + 2 * math.pi for t in th[1:]]
                vecs = [w[2] for w in win]
                U = [j for j in range(len(primes)) if len({v[j] for v in vecs}) > 1]
                if len(U) < k:
                    continue
                D = np.array([[v[j] - vecs[0][j] for j in U] for v in vecs[1:]], dtype=float)
                if np.linalg.matrix_rank(D) < k:
                    continue
                lam = min(max(abs(t - tb) for t in th) for tb in th)
                if lam <= 0:
                    continue
                g = math.log(1 / lam) / sum(logH[j] for j in U)
                hist[k].append(g)
                if g > best[k][0]:
                    best[k] = (g, (N, k, [w[1] for w in win]))
        # largest cluster with normalised constant <= CMAX (window search over consecutive points)
        bestM = 0
        for s in range(L):
            for t in range(2, min(L, 40)):
                win = [ang[(s + u) % L] for u in range(t)]
                zs = [w[1] for w in win]
                # normalise by common factor: drop primes where all exponents agree
                vecs = [w[2] for w in win]
                Rn2 = N
                for j, p in enumerate(primes):
                    col = [v[j] for v in vecs]
                    e = fac[p]
                    common = min(col) + min(e - c for c in col)
                    Rn2 //= p ** common
                chord = max(math.hypot(z[0] - w[0], z[1] - w[1]) for z in zs for w in zs)
                g = math.sqrt(N / Rn2)
                Cn = (chord / g) / math.sqrt(math.sqrt(Rn2))
                if Cn > CMAX:
                    break
                if t > bestM:
                    bestM = t
                    D = np.array([[v[j] - vecs[0][j] for j in range(len(primes))] for v in vecs[1:]], dtype=float)
                    rk = int(np.linalg.matrix_rank(D)) if t > 1 else 0
                    rank_of_big[N] = (t, rk, Cn, math.log(Rn2))
    print(f"X = {X:.3g}: circles with >= {MINPTS} points: {nN}")
    for k in (3, 4, 5):
        h = sorted(hist[k])
        if not h:
            print(f"  k={k}: no rank-{k} windows")
            continue
        q = lambda f: h[min(len(h) - 1, int(f * len(h)))]
        g, info = best[k]
        print(f"  k={k}: rank-{k} windows={len(h):>8}  median gamma={q(0.5):.4f}  99.9%={q(0.999):.4f}  "
              f"max gamma={g:.4f} at N={info[0]}  (Dirichlet 1/k={1/k:.4f}; uniform needs SLW with g<1/4)")
    by = {}
    hi = []
    for N, (t, rk, Cn, W) in rank_of_big.items():
        by.setdefault(t, []).append(rk)
        if rk >= 5:
            hi.append((N, t, rk, Cn, W))
    for t in sorted(by):
        rks = by[t]
        print(f"  largest cluster with C<={CMAX}: size {t:>2}: circles={len(rks):>6}, max rank={max(rks)}, "
              f"rank distribution={dict(sorted({r: rks.count(r) for r in set(rks)}.items()))}")
    print(f"  clusters (C<={CMAX}) of rank >= 5: {len(hi)}; largest normalised W' = log N' among them: "
          f"{max((h[4] for h in hi), default=0):.2f}")
    for h in sorted(hi, key=lambda h: -h[4])[:10]:
        print(f"    N={h[0]} size={h[1]} rank={h[2]} C={h[3]:.3f} W'={h[4]:.2f}")


if __name__ == '__main__':
    main()
