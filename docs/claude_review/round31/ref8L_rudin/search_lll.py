"""Lattice search for large-diameter two-class sets.

For B points (y_0 = 0), pick bipartitions pi of {0..B-1} and pairwise coprime
moduli N_pi (random primes), form the lattice
  Lambda = { y in Z^{B-1} : y_i = y_j mod N_pi when i,j are on the same side of pi },
LLL-reduce it (own exact LLL), enumerate small combinations of the reduced basis,
keep vectors with pairwise distinct coordinates and compute a lower bound for
Qmax(X) (exact e_p at every prime dividing some N_pi and at every prime < 200).
Exponent E = log Qmax / log diam.
Usage: search_lll.py B trials seed
"""
import sys
import random
from fractions import Fraction
from math import log, prod, gcd
from itertools import product
from sympy import randprime, primerange
from sympy.matrices.normalforms import hermite_normal_form
from sympy import Matrix

B = int(sys.argv[1])
TRIALS = int(sys.argv[2])
rng = random.Random(int(sys.argv[3]))
SMALL = list(primerange(2, 200))


def lll(basis, delta=Fraction(3, 4)):
    b = [list(v) for v in basis]
    n = len(b)

    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))

    def gs():
        bs, mu = [], [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in b[i]]
            for j in range(i):
                mu[i][j] = Fraction(dot(b[i], bs[j])) / dot(bs[j], bs[j])
                v = [x - mu[i][j] * y for x, y in zip(v, bs[j])]
            bs.append(v)
        return bs, mu

    bs, mu = gs()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                b[k] = [x - q * y for x, y in zip(b[k], b[j])]
                bs, mu = gs()
        if dot(bs[k], bs[k]) >= (delta - mu[k][k - 1] ** 2) * dot(bs[k - 1], bs[k - 1]):
            k += 1
        else:
            b[k], b[k - 1] = b[k - 1], b[k]
            bs, mu = gs()
            k = max(k - 1, 1)
    return b


def emax(X, p):
    e, pe = 0, p
    while len({x % pe for x in X}) <= 2:
        e += 1
        pe *= p
        if e > 200:
            break
    return e


def qmax_lb(X, primes):
    Q = 1
    for p in primes:
        e = emax(X, p)
        Q *= p ** e
    return Q


def trial():
    n = B - 1
    allparts = []
    for mask in range(1, 2 ** n):  # side not containing 0, nonempty
        allparts.append([i for i in range(n) if mask >> i & 1])
    m = rng.randint(1, min(len(allparts), 3 * B))
    parts = rng.sample(allparts, m)
    bits = rng.randint(8, 22)
    Ns = set()
    while len(Ns) < m:
        Ns.add(randprime(2 ** (bits - 1), 2 ** bits + 1))
    Ns = list(Ns)
    M = prod(Ns)
    gens = [[M if i == j else 0 for j in range(n)] for i in range(n)]
    for side, N in zip(parts, Ns):
        Mi = M // N
        E = Mi * pow(Mi, -1, N)
        gens.append([E if i in side else 0 for i in range(n)])
    H = hermite_normal_form(Matrix(gens).T).T  # rows generate lattice
    basis = [list(H.row(i)) for i in range(H.rows) if any(H.row(i))]
    basis = [[int(x) for x in v] for v in basis]
    red = lll(basis)
    best = (0, None)
    primes = sorted(set(Ns) | set(SMALL))
    rng_c = range(-1, 2) if n > 6 else range(-2, 3)
    for coef in product(rng_c, repeat=n):
        if not any(coef):
            continue
        y = [0] + [sum(c * v[i] for c, v in zip(coef, red)) for i in range(n)]
        if len(set(y)) < B:
            continue
        lo = min(y)
        X = [t - lo for t in y]
        diam = max(X)
        if diam < 1000:
            continue
        Q = qmax_lb(X, primes)
        if Q <= 1:
            continue
        E = log(Q) / log(diam)
        if E > best[0]:
            best = (E, (Q, diam, sorted(X)))
    return best


overall = (0, None)
for t in range(TRIALS):
    E, info = trial()
    print("trial", t, E, flush=True)
    if E > overall[0]:
        overall = (E, info)
print("B", B, "best", overall)
