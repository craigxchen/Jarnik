"""Referee check R1: vertical polynomial templates refute SLW(k, g, kappa) of weakest.md Sec. 2.1.

Construction (not a circle configuration).  Fix distinct positive integers b_1..b_D and integer vectors
sigma_1..sigma_k in Z^D, linearly independent, with sum_l sigma_il * b_l^(2j-1) = 0 for j = 1..s.
For t divisible by Q (so that b_l | t and (t/b_l)^2 + 1 is coprime to every prime <= 2 max b), put
u_l = t/b_l + i and beta_i = prod_l u_l^sigma_il (conj for negative exponents).  Then
    arg beta_i = sum_l sigma_il * arctan(b_l / t) = O(t^-(2s+1)),
the Gaussian primes of different u_l are disjoint (and no u_l meets conj u_m), so the coefficient
rows of the forms Lambda_i = Log(beta_i / conj beta_i) = 2 i arg beta_i in the alpha_p = pi_p/conj pi_p
are linearly independent, and the union height is prod_p H(alpha_p) <= 2^n prod_l N(u_l) ~ t^(2D).
So gamma := log(1/max|Lambda_i|) / sum_U log H(alpha_p) -> (2s+1)/(2D).

Everything below is exact (integers, sympy factorisation) except the angles (mpmath, 300 digits).
We also report the weighted exponent gamma_w (heights weighted by B_p, Mahler measure p),
i.e. log(1/max|Lambda|) / sum_p B_p log p, relevant for the weighted repair discussed in the report.
"""
import sys
from fractions import Fraction
from math import gcd
import mpmath
import sympy

mpmath.mp.dps = 300


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
    assert x * x + y * y == p
    return (x, y)


def gmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def gdivexact(u, v):
    n = v[0] * v[0] + v[1] * v[1]
    w = gmul(u, (v[0], -v[1]))
    if w[0] % n or w[1] % n:
        return None
    return (w[0] // n, w[1] // n)


def gval(z, pi):
    v = 0
    while True:
        w = gdivexact(z, pi)
        if w is None:
            return v
        z = w
        v += 1


def int_rank(rows):
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


def kernel_basis(b, s):
    """Integer basis of {sigma : sum sigma_l b_l^(2j-1) = 0, j=1..s}, LLL-reduced for small entries."""
    A = sympy.Matrix([[bl ** (2 * j - 1) for bl in b] for j in range(1, s + 1)])
    ns = A.nullspace()
    vecs = []
    for v in ns:
        den = sympy.ilcm(*[sympy.fraction(x)[1] for x in v])
        w = [int(x * den) for x in v]
        g = 0
        for x in w:
            g = gcd(g, abs(x))
        vecs.append([x // g for x in w])
    try:
        from sympy.matrices.normalforms import hermite_normal_form  # noqa
        L = sympy.Matrix(vecs).lll()
        vecs = [list(map(int, L.row(i))) for i in range(L.rows)]
    except Exception:
        pass
    return vecs


def modulus_Q(b):
    """Q with: b_l | t and q | t/b_l for all primes q <= 2 max b, whenever Q | t."""
    B = 2 * max(b)
    Q = 1
    for l in b:
        Q = sympy.ilcm(Q, l)
    for q in sympy.primerange(2, B + 1):
        Q *= q
    return int(Q)


def analyse(b, sigmas, t, verbose=False):
    D = len(b)
    k = len(sigmas)
    assert all(t % bl == 0 for bl in b)
    us = [(t // bl, 1) for bl in b]
    prime_owner = {}
    pis = {}
    expo = {}  # p -> (l, signed exponent of pi_p in u_l)
    for l, u in enumerate(us):
        n = u[0] * u[0] + 1
        fac = sympy.factorint(n)
        for p, e in fac.items():
            assert p > 2 * max(b), ("small prime", p)
            assert p % 4 == 1
            assert p not in prime_owner, ("shared prime", p)
            prime_owner[p] = l
            pi = two_squares(p)
            pis[p] = pi
            a = gval(u, pi)
            c = gval(u, (pi[0], -pi[1]))
            assert a * c == 0 and a + c == e
            expo[p] = a - c
    primes = sorted(prime_owner)
    # coefficient rows
    rows = []
    for sg in sigmas:
        rows.append([sg[prime_owner[p]] * expo[p] for p in primes])
    rk = int_rank(rows)
    # Lambda_i = 2 i arg beta_i, arg beta_i = sum sigma_il atan(b_l/t)
    lams = []
    for sg in sigmas:
        a = mpmath.mpf(0)
        for l in range(D):
            a += sg[l] * mpmath.atan(mpmath.mpf(b[l]) / t)
        lams.append(abs(2 * a))
    lam = max(lams)
    # heights: H(alpha_p) = max(p, 2|Re pi^2|), Mahler measure = p
    U = [p for p in primes if any(r[j] != 0 for r in rows for j in [primes.index(p)])]
    logH = sum(mpmath.log(max(p, 2 * abs(pis[p][0] ** 2 - pis[p][1] ** 2))) for p in U)
    Bp = {p: max(abs(r[primes.index(p)]) for r in rows) for p in primes}
    logW = sum(Bp[p] * mpmath.log(p) for p in U)
    gamma = mpmath.log(1 / lam) / logH
    gamma_w = mpmath.log(1 / lam) / logW
    # SLW right side with c = 1/2, kappa = 1 (coefficient part incl. B_0 <= 1 + 2 sum|b|)
    B0 = 1 + 2 * max(sum(abs(x) for x in r) for r in rows)
    coef = sum(mpmath.log(1 + Bp[p]) for p in U) + mpmath.log(1 + B0)
    n = len(U)
    return dict(n=n, rank=rk, k=k, lam=lam, gamma=gamma, gamma_w=gamma_w, logH=logH, logW=logW,
                coef=coef, logt=mpmath.log(t))


def run(label, b, s, sigmas, us, target):
    Q = modulus_Q(b)
    print(f"== {label}: b={b}, s={s}, k={len(sigmas)} forms, D={len(b)}, Q={Q}; "
          f"predicted gamma -> (2s+1)/(2D) = {Fraction(2*s+1, 2*len(b))} = {float(Fraction(2*s+1, 2*len(b))):.4f}")
    for sg in sigmas:
        for j in range(1, s + 1):
            assert sum(x * bl ** (2 * j - 1) for x, bl in zip(sg, b)) == 0
    print("   sigmas:", sigmas)
    assert int_rank(sigmas) == len(sigmas)
    for u in us:
        t = Q * u
        r = analyse(b, sigmas, t)
        assert r['rank'] == r['k'], r
        # effective exponent with c=1/2, kappa=1: (log 1/lam - (n+1) log 2 - coef)/logH
        geff = (mpmath.log(1 / r['lam']) - (r['n'] + 1) * mpmath.log(2) - r['coef']) / r['logH']
        print(f"   t=Q*{u:<14d} log t={float(r['logt']):7.2f} n={r['n']:3d} rank={r['rank']} "
              f"max|Lambda|={mpmath.nstr(r['lam'], 5):>11} gamma={float(r['gamma']):.4f} "
              f"gamma_w={float(r['gamma_w']):.4f} gamma_eff(c=1/2,kappa=1)={float(geff):.4f}")
    print(f"   -> SLW({len(sigmas)}, g, kappa) violated for every g < {target} as t -> infinity (Prop. R1)")


if __name__ == '__main__':
    us = [10 ** 3 + 7, 10 ** 6 + 3, 10 ** 9 + 7]
    # k = 5, D = 6, s = 1, coefficients in {0, +-1}: gamma -> 1/4
    run("k=5, {0,+-1}", [1, 2, 3, 4, 5, 6], 1,
        [[1, 1, -1, 0, 0, 0], [1, 0, 1, -1, 0, 0], [1, 0, 0, 1, -1, 0], [1, 0, 0, 0, 1, -1], [0, 1, 1, 0, -1, 0]],
        us, "1/4")
    # k = 4, D = 5, s = 1, {0,+-1}: gamma -> 3/10 (empties the window [1/4, 4/15) of Prop. 2.7)
    run("k=4, {0,+-1}", [1, 2, 3, 4, 5], 1,
        [[1, 1, -1, 0, 0], [1, 0, 1, -1, 0], [1, 0, 0, 1, -1], [0, 1, 1, 0, -1]], us, "3/10")
    # k = 5, D = 7, s = 2, integer coefficients: gamma -> 5/14 > 1/4
    b7 = [1, 2, 3, 4, 5, 6, 7]
    run("k=5, integer", b7, 2, kernel_basis(b7, 2), us, "5/14")
    # k = 5, D = 9, s = 4, integer coefficients: gamma -> 1/2
    b9 = list(range(1, 10))
    run("k=5, integer", b9, 4, kernel_basis(b9, 4), us[:3], "1/2")
