"""Referee check of Theorem 1.3 (rank reduction) of weakest.md on actual lattice points.

For every N <= X with all prime factors = 1 mod 4 and r_2(N) >= 16, and every window of m = 2..8
angularly consecutive lattice points, normalise (divide by the Gaussian gcd), compute exponent vectors at
fixed Gaussian primes, W = log N', the angular width Delta and the exact arc constant C := Delta e^(W/4).
Check exactly (integers) / at 50 digits (angles):
  (1) the exponent vectors are distinct (whenever Delta < pi/2, as step 1 of the proof requires);
  (2) separation: for all u != u' in the difference set D, w(u - u') >= rho := W/2 - 2 log C - log 2;
  (3) if rho > 0: M <= ((1 + 2W/rho)^k + 1)/2, k = rank of the exponent lattice.
"""
import sys
from fractions import Fraction
import mpmath
import sympy

mpmath.mp.dps = 50


def gmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def gdivexact(u, v):
    n = v[0] * v[0] + v[1] * v[1]
    w = gmul(u, (v[0], -v[1]))
    if w[0] % n or w[1] % n:
        return None
    return (w[0] // n, w[1] // n)


def ggcd(u, v):
    while v != (0, 0):
        n = v[0] * v[0] + v[1] * v[1]
        w = gmul(u, (v[0], -v[1]))
        q = ((2 * w[0] + n) // (2 * n), (2 * w[1] + n) // (2 * n))
        t = gmul(q, v)
        u, v = v, (u[0] - t[0], u[1] - t[1])
    return u


def gval(z, pi):
    v = 0
    while True:
        w = gdivexact(z, pi)
        if w is None:
            return v
        z, v = w, v + 1


def two_squares(p):
    for c in range(2, p):
        t = pow(c, (p - 1) // 4, p)
        if (t * t) % p == p - 1:
            break
    a, b = p, t
    lim = int(p ** 0.5)
    while b > lim:
        a, b = b, a % b
    x = b
    y = sympy.integer_nthroot(p - x * x, 2)[0]
    return (x, y)


def rank(rows):
    if not rows:
        return 0
    m = [[Fraction(x) for x in r] for r in rows]
    rk, ncol = 0, len(m[0])
    for c in range(ncol):
        piv = next((i for i in range(rk, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[rk], m[piv] = m[piv], m[rk]
        for i in range(len(m)):
            if i != rk and m[i][c] != 0:
                f = m[i][c] / m[rk][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[rk])]
        rk += 1
    return rk


def points(N):
    pts = []
    x = 0
    while x * x <= N:
        y2 = N - x * x
        y = sympy.integer_nthroot(y2, 2)[0]
        if y * y == y2:
            for sx in (1, -1):
                for sy in (1, -1):
                    pts.append((sx * x, sy * y))
        x += 1
    pts = sorted(set(pts), key=lambda z: float(mpmath.atan2(z[1], z[0])))
    return pts


def main(X):
    n_win = n_pairs = n_pack = 0
    worst_slack = None
    max_rank_seen = {}
    for N in range(5, X + 1, 4):
        fac = sympy.factorint(N)
        if any(p % 4 != 1 for p in fac):
            continue
        r2 = 4 * sympy.prod([e + 1 for e in fac.values()])
        if r2 < 16:
            continue
        pts = points(N)
        L = len(pts)
        for m in range(2, 9):
            for s in range(L):
                win = [pts[(s + j) % L] for j in range(m)]
                g = win[0]
                for z in win[1:]:
                    g = ggcd(g, z)
                zs = [gdivexact(z, g) for z in win]
                Np = zs[0][0] ** 2 + zs[0][1] ** 2
                if Np == 1:
                    continue
                f2 = sympy.factorint(Np)
                primes = sorted(f2)
                pis = {p: two_squares(p) for p in primes}
                vecs = [tuple(gval(z, pis[p]) for p in primes) for z in zs]
                for z, v in zip(zs, vecs):
                    for p, a in zip(primes, v):
                        assert a + gval(z, (pis[p][0], -pis[p][1])) == f2[p]
                W = mpmath.log(Np)
                th = [mpmath.atan2(z[1], z[0]) for z in win]
                t0 = th[0]
                rel = [mpmath.atan2(mpmath.sin(t - t0), mpmath.cos(t - t0)) for t in th]
                Delta = max(rel) - min(rel)
                if Delta >= mpmath.pi / 2 - mpmath.mpf(10) ** -40:
                    continue
                assert len(set(vecs)) == m, ("(1) fails", N, win)
                C = Delta * mpmath.e ** (W / 4)
                rho = W / 2 - 2 * mpmath.log(C) - mpmath.log(2)
                logp = [mpmath.log(p) for p in primes]
                Dset = set()
                for a in vecs:
                    for b in vecs:
                        Dset.add(tuple(x - y for x, y in zip(a, b)))
                Dl = list(Dset)
                for i in range(len(Dl)):
                    for j in range(i + 1, len(Dl)):
                        w = sum(abs(x - y) * lp for x, y, lp in zip(Dl[i], Dl[j], logp))
                        sl = w - rho
                        n_pairs += 1
                        assert sl > -mpmath.mpf(10) ** -30, ("(2) fails", N, win)
                        if worst_slack is None or sl < worst_slack[0]:
                            worst_slack = (sl, N, m)
                k = rank([tuple(x - y for x, y in zip(v, vecs[0])) for v in vecs[1:]])
                max_rank_seen[m] = max(max_rank_seen.get(m, 0), k)
                if rho > 0:
                    bound = ((1 + 2 * W / rho) ** k + 1) / 2
                    assert m <= bound, ("(3) fails", N, win)
                    n_pack += 1
                n_win += 1
    print(f"X={X}: windows checked={n_win}, separation pairs checked={n_pairs}, packing bound applicable={n_pack}")
    print(f"  smallest separation slack w(u-u') - rho = {mpmath.nstr(worst_slack[0], 8)} (N={worst_slack[1]}, m={worst_slack[2]})")
    print(f"  max rank per window size: {max_rank_seen}")
    print("  (1) distinct vectors, (2) separation, (3) packing: ALL PASS")


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 30000)
