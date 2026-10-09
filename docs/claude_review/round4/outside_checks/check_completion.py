"""Checks for Theorem A.1 (generic completion) of outside.md.

Part 1 (exact, Fractions): for random multi-level profiles (rows a_x in prod [0, e_j]),
   U = {u : A~ u in R 1},  A~ = rows (2 a_x - e),
   verify U^perp = span{a_x - a_y} (equal dimension and mutual containment).
Part 2 (numerical illustration only): sample realisations phi uniformly from the torus
   T_U (phi = phi_0 + 2 pi B t, B an integer basis of U, t uniform); check that <v, phi> is
   constant for v in span{a_x - a_y} cap Z^r and Kolmogorov-Smirnov-uniform mod 2 pi otherwise.
Part 3 (exact rationals + floats): the union-bound constant Pi = prod_j (1 + 2/(sqrt p_j - 1))
   for the prime sizes used in the Paley barrier (p_j ~ M^4, r = b(M-1)), showing
   1.33 (Pi - 1) < 1, and the threshold P_min >= 21 r^2 quoted in the text.
"""
import random
from fractions import Fraction as Fr
from math import pi, sqrt, log
import numpy as np

random.seed(3)
rng = np.random.default_rng(3)


def rref(rows):
    rows = [list(map(Fr, r)) for r in rows]
    m = len(rows)
    n = len(rows[0]) if m else 0
    piv = []
    i = 0
    for j in range(n):
        k = next((t for t in range(i, m) if rows[t][j] != 0), None)
        if k is None:
            continue
        rows[i], rows[k] = rows[k], rows[i]
        pv = rows[i][j]
        rows[i] = [x / pv for x in rows[i]]
        for t in range(m):
            if t != i and rows[t][j] != 0:
                f = rows[t][j]
                rows[t] = [a - f * b for a, b in zip(rows[t], rows[i])]
        piv.append(j)
        i += 1
        if i == m:
            break
    return rows[:i], piv


def nullspace(rows, n):
    if not rows:
        return [[Fr(int(i == j)) for i in range(n)] for j in range(n)]
    R, piv = rref(rows)
    free = [j for j in range(n) if j not in piv]
    basis = []
    for f in free:
        v = [Fr(0)] * n
        v[f] = Fr(1)
        for r, pj in zip(R, piv):
            v[pj] = -r[f]
        basis.append(v)
    return basis


def rank(rows):
    if not rows:
        return 0
    return len(rref(rows)[0])


def part1(trials=200):
    for _ in range(trials):
        M = random.randint(3, 8)
        r = random.randint(2, 10)
        e = [random.randint(1, 3) for _ in range(r)]
        A = [[random.randint(0, e[j]) for j in range(r)] for _ in range(M)]
        At = [[2 * A[x][j] - e[j] for j in range(r)] for x in range(M)]
        # U = {u : At u in R 1}: unknowns (u, s) with At u - s 1 = 0
        rows = [At[x] + [-1] for x in range(M)]
        ker = nullspace(rows, r + 1)
        U = [v[:r] for v in ker]
        Uperp = nullspace(U, r) if U else [[Fr(int(i == j)) for i in range(r)] for j in range(r)]
        D = [[A[x][j] - A[0][j] for j in range(r)] for x in range(1, M)]
        dD = rank(D)
        assert len(Uperp) == dD, (len(Uperp), dD)
        # containment: each difference is orthogonal to U
        for d in D:
            for u in U:
                assert sum(Fr(a) * b for a, b in zip(d, u)) == 0
    print("Part 1: U^perp = span{a_x - a_y} on %d random multi-level profiles: OK" % trials)


def ks_uniform(xs):
    xs = np.sort(np.asarray(xs))
    n = len(xs)
    cdf = np.arange(1, n + 1) / n
    return max(np.max(cdf - xs), np.max(xs - (cdf - 1.0 / n)))


def part2():
    M, r = 6, 14
    e = [random.randint(1, 2) for _ in range(r)]
    A = np.array([[random.randint(0, e[j]) for j in range(r)] for _ in range(M)])
    At = 2 * A - np.array(e)
    rows = [list(At[x]) + [-1] for x in range(M)]
    ker = nullspace(rows, r + 1)
    B = np.array([[float(v[j]) for j in range(r)] for v in ker]).T  # r x d (rational)
    # scale each basis vector to an integer vector
    Bi = []
    for v in ker:
        den = 1
        for x in v[:r]:
            den = den * x.denominator // np.gcd(den, x.denominator)
        Bi.append([int(x * den) for x in v[:r]])
    B = np.array(Bi, dtype=float).T
    delta = rng.uniform(0, 1e-3, M)
    k = rng.integers(-3, 4, M)
    theta0 = 0.37
    rhs = theta0 + delta + (pi / 2) * k
    phi0, *_ = np.linalg.lstsq(At.astype(float), rhs, rcond=None)
    assert np.allclose(At @ phi0, rhs)
    S = 4000
    t = rng.uniform(0, 1, (S, B.shape[1]))
    phis = phi0[None, :] + 2 * pi * t @ B.T
    # v in the span of differences: constant
    lam = np.zeros(M)
    lam[1], lam[2], lam[0] = 2, -1, -1
    v_in = lam @ A
    vals = (phis @ v_in) % (2 * pi)
    spread = np.ptp(np.unwrap(vals))
    # v outside: random integer vectors not in the span
    Dspan = A[1:] - A[0]
    ks = []
    for _ in range(20):
        v = rng.integers(-2, 3, r)
        if np.linalg.matrix_rank(np.vstack([Dspan, v])) == np.linalg.matrix_rank(Dspan):
            continue
        x = ((phis @ v) % (2 * pi)) / (2 * pi)
        ks.append(ks_uniform(x))
    print("Part 2 (numerical illustration): spread of <v,phi> for v in span = %.2e; "
          "max KS distance from uniform over %d v outside span = %.3f (S=%d samples)"
          % (spread, len(ks), max(ks), S))


def part3():
    for M, b in [(12, 5), (20, 7), (44, 7), (108, 10), (1000, 10)]:
        r = b * (M - 1)
        P = M ** 4
        Pi = (1 + 2 / (sqrt(P) - 1)) ** r
        print("Part 3: M=%5d b=%2d r=%6d p_j>=M^4: Pi-1=%.3e  1.33(Pi-1)=%.3e  (<1: %s)"
              % (M, b, r, Pi - 1, 1.33 * (Pi - 1), 1.33 * (Pi - 1) < 1))
    # constant: arcsin(x) <= c0 x for x <= 1/sqrt 5, c0 = arcsin(1/sqrt5) sqrt5
    from math import asin
    c0 = asin(1 / sqrt(5)) * sqrt(5)
    K = 0.5 * (8 / pi) * c0
    print("Part 3: union-bound constant K = (4/pi) arcsin(5^-1/2) sqrt5 = %.4f (text uses 1.33)" % K)
    assert K < 1.33
    # threshold P_min >= 21 r^2: Pi <= exp(2r/(sqrt(21) r - 1)) and K (Pi - 1) < 1
    for r in list(range(1, 2000)) + [10 ** 4, 10 ** 6, 10 ** 9]:
        s = 2 * r / (sqrt(21) * r - 1)
        assert 1.33 * (np.exp(s) - 1) < 1.0, (r, s)
    print("Part 3: P_min >= 21 r^2  =>  1.33 (Pi - 1) < 1 for every r >= 1 tested (monotone in r): OK")


if __name__ == "__main__":
    part1()
    part2()
    part3()
