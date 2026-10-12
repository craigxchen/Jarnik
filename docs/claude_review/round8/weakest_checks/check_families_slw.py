"""Rank and SLW exponents of the known infinite families (exact arithmetic).

(G) golden eight-point family of growth/families.md Prop. G:
    Z_w(n) = F^[p=1] conj(F)^[p=3] prod_{s=0}^3 J_{n+s}^{w_s} conj(J_{n+s})^{1-w_s},
    w in {0,1}^4 of odd weight p, J_r = f_{r+1} + i f_r, F = 1+2i, n = 1 mod 3.
(CG) Cilleruelo-Granville six points prod_{j=1}^4 (a+j+i sigma_j), sigma in {+-1}^4, sum sigma = 0.
(CGP) the CG family at a with a+1+i = (2+i)(x+iy)^2 (Pell X^2-5Y^2=1, x = X-2Y, y = Y),
      so that (a+1)^2+1 = 5 (x^2+y^2)^2 has a large square factor.

For every member: normalised N, rank of exponent differences (prime coordinates),
and gamma_k = log(1/max|Lambda_i|) / sum_{j in U} log H(alpha_j) over all (k+1)-subsets of
rank k, where Lambda_i = Log(z_i/z_0).  SLW(k, g) says gamma_k <= g + o(1).
Also reported: log N / log rad(N) and the normalised arc constant C = spread*sqrt(R)
with spread measured as the chord length max|z_i - z_j|.
"""
import itertools
import mpmath
import sympy
from g8common import (gmul, gconj, gnorm, profile, int_rank, slw_exponent)


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def golden_points(n):
    F = (1, 2)
    pts = []
    for w in itertools.product((0, 1), repeat=4):
        p = sum(w)
        if p % 2 == 0:
            continue
        z = F if p == 1 else gconj(F)
        for s in range(4):
            J = (fib(n + s + 1), fib(n + s))
            z = gmul(z, J if w[s] else gconj(J))
        pts.append(z)
    return pts


def cg_points(a):
    pts = []
    for s in itertools.product((1, -1), repeat=4):
        if sum(s) != 0:
            continue
        z = (1, 0)
        for j, sj in enumerate(s, start=1):
            z = gmul(z, (a + j, sj))
        pts.append(z)
    return pts


def pell_a_values(count):
    X, Y = -9, 4
    out = []
    for _ in range(count):
        x, y = X - 2 * Y, Y
        assert x * x - y * y + 4 * x * y == 1
        a1 = 2 * (x * x - y * y) - 2 * x * y
        assert a1 * a1 + 1 == 5 * (x * x + y * y) ** 2
        out.append(a1 - 1)
        X, Y = 9 * X - 20 * Y, -4 * X + 9 * Y
    return out


def analyse(name, pts, k):
    N, fac, primes, pis, vecs, npts = profile(pts)
    M = len(npts)
    diffs = [[a - b for a, b in zip(v, vecs[0])] for v in vecs[1:]]
    rk = int_rank(diffs)
    lograd = sum(mpmath.log(p) for p in primes)
    logN = mpmath.log(N)
    R = mpmath.sqrt(mpmath.sqrt(N))  # R = N^(1/2); here sqrt(R) = N^(1/4)
    chord = max(mpmath.sqrt(gnorm((z[0] - w[0], z[1] - w[1]))) for z in npts for w in npts)
    C = chord / R
    best = None
    nsub = 0
    for sub in itertools.combinations(range(M), k + 1):
        sv = [vecs[i] for i in sub]
        d = [[a - b for a, b in zip(v, sv[0])] for v in sv[1:]]
        if int_rank(d) < k:
            continue
        nsub += 1
        # order: try every base point (the max over forms depends on the base)
        for b in range(k + 1):
            order = [sub[b]] + [sub[i] for i in range(k + 1) if i != b]
            g, lam, U = slw_exponent([npts[i] for i in order], [vecs[i] for i in order], primes, pis)
            if best is None or g > best[0]:
                best = (g, lam, len(U))
    print(f"{name}: M={M} logN={float(logN):.2f} omega={len(primes)} logN/lograd={float(logN/lograd):.4f} "
          f"rank={rk} C={float(C):.4f} rank-{k} subsets={nsub} max gamma_{k}={best[0] if best else None:.5f}")
    return rk, best


if __name__ == '__main__':
    print("(G) golden family, k = 4 (expect rank 4, gamma_4 -> 1/4 from above as logN/lograd -> 1)")
    for n in (4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34):
        analyse(f"  n={n:>2}", golden_points(n), 4)
    print("(CG) generic Cilleruelo-Granville members, k = 3 (expect rank 3, gamma_3 ~ 1/4)")
    for a in (10**3, 10**4 + 7, 10**5 + 3, 10**6 + 1):
        analyse(f"  a={a}", cg_points(a), 3)
    print("(CGP) Pell members a+1+i = (2+i)(x+iy)^2, k = 3 (expect gamma_3 -> 2/7 = 0.2857)")
    for a in pell_a_values(4):
        analyse(f"  a={a}", cg_points(a), 3)
