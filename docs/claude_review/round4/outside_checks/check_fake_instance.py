"""Numerical ILLUSTRATION (not a proof) of Theorem A.3: an explicit fake circle for M = 8.

Profile: Sylvester H of order 8, b = 6M copies per label, two flips per label (as in
check_fake_barrier.py), r = 6M(M-1) = 336 distinct primes p = 1 mod 4 taken consecutively above
P = 10^9 (so all weights log p agree to ~1e-5; this replaces the pigeonhole balancing of the
text, which is only needed to avoid short-interval prime theorems in the proof).
Steps:
 1. W = sum log p_j, Delta = C exp(-W/4) with C = 1/2; delta = Delta * u, u uniform in [0,1]^M.
 2. Every character c (zero-sum, ||c||_1 <= 4, exhaustive) and 3000 random larger ones: check the
    residue-strengthened gap  |c.delta|/2 >= arcsin(Q e^{-V_c/2}),  Q = lcm(1..2M), in log form.
 3. Pick a realisation phi (theta_0 random; delta is below double precision, so the realisation
    equations are solved with delta = 0, which changes <v,phi> by < 1e-700 for every v) and test the
    plain monomial gaps at 200000 random v with entries in {-2..2} on random supports (all v
    outside L, the character lattice, with probability 1).
"""
import random
from math import log, asin, sqrt, pi, exp, lcm
from itertools import combinations, product
import numpy as np
from check_fake_barrier import sylvester, normalize, build

random.seed(11)
rng = np.random.default_rng(11)


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def main():
    M, C = 8, 0.5
    H = normalize(sylvester(M))
    b = 6 * M
    S, flip_row, designated, idx = build(H, b)
    r = S.shape[1]
    primes = []
    n = 10 ** 9 + 1
    while len(primes) < r:
        if n % 4 == 1 and is_prime(n):
            primes.append(n)
        n += 4
    w = np.array([log(p) for p in primes])
    W = w.sum()
    Q = 1
    for k in range(1, 2 * M + 1):
        Q = lcm(Q, k)
    logQ = log(Q)
    Pi = np.prod(1 + 2 / (np.sqrt(np.array(primes, dtype=float)) - 1))
    print("M=%d b=%d r=%d  primes in [%d, %d]  W=%.2f  log(R)=%.2f  Pi-1=%.4f (1.33(Pi-1)=%.4f < 1)"
          % (M, b, r, primes[0], primes[-1], W, W / 2, Pi - 1, 1.33 * (Pi - 1)))
    u = rng.uniform(0, 1, M)
    logDelta = log(C) - W / 4
    # characters
    cs = []
    for k in (2, 3, 4):
        for supp in combinations(range(M), k):
            for vals in product(range(-4, 5), repeat=k):
                if 0 in vals or sum(vals) != 0 or sum(abs(v) for v in vals) > 4:
                    continue
                c = np.zeros(M, dtype=np.int64)
                c[list(supp)] = vals
                cs.append(c)
    for _ in range(3000):
        c = np.zeros(M, dtype=np.int64)
        for _ in range(random.randint(1, 3 * M)):
            i, j = random.sample(range(M), 2)
            c[i] += 1
            c[j] -= 1
        if c.any():
            cs.append(c)
    worst = None
    for c in cs:
        V = 0.5 * float(np.abs(c @ S) @ w)
        lhs = log(abs(float(c @ u))) + logDelta - log(2)  # log(|c.delta|/2)
        # need |c.delta|/2 >= arcsin(Q e^{-V/2}); arcsin(y) <= 1.04 y for y <= 5^-1/2
        rhs = log(1.04) + logQ - V / 2
        assert logQ - V / 2 < log(1 / sqrt(5))
        margin = lhs - rhs
        assert margin > 0, (c, margin)
        worst = margin if worst is None else min(worst, margin)
    print("step 2: %d characters satisfy the residue-strengthened gaps (Q = lcm(1..%d)); "
          "min log-margin %.1f" % (len(cs), 2 * M, worst))
    # realisation with delta ~ 0
    theta0 = rng.uniform(0, 2 * pi)
    rhs = np.full(M, theta0)
    phi0, *_ = np.linalg.lstsq(S.astype(float), rhs, rcond=None)
    # add a random element of ker S (generic point of the realisation torus)
    _, _, Vt = np.linalg.svd(S.astype(float))
    K = Vt[M:].T
    phi = phi0 + K @ rng.normal(0, 50, K.shape[1])
    assert np.allclose(S @ phi, theta0)
    bad = 0
    mind = None
    for _ in range(200000):
        k = random.randint(1, 12)
        v = np.zeros(r, dtype=np.int64)
        supp = rng.choice(r, k, replace=False)
        v[supp] = rng.choice([-2, -1, 1, 2], k)
        x = float(v @ phi)
        d = abs(((x + pi / 8) % (pi / 4)) - pi / 8)
        wv = float(np.abs(v) @ w)
        need = asin(exp(-wv / 2))
        if d < need:
            bad += 1
        ratio = d / need
        mind = ratio if mind is None else min(mind, ratio)
    print("step 3: 200000 random monomials v: %d gap violations; min dist/required = %.3g" % (bad, mind))


if __name__ == "__main__":
    main()
