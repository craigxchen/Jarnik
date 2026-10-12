"""Common exact helpers for the fully flipped Hadamard phase lemma checks (round 4).

Gaussian integers are (a, b) pairs of Python ints.  Hadamard matrices are lists of lists of +-1
with FIRST COLUMN ALL ONES (normalised as in walsh-extraction.md); column 0 is the constant label.
"""
import math, random

# ---------------------------------------------------------------- primes and Gaussian arithmetic
def is_prime(n):
    """Deterministic Miller-Rabin for n < 3.3e24 (bases = first 12 primes)."""
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
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


def split_primes_from(start, count):
    out, p = [], max(5, start)
    while len(out) < count:
        if p % 4 == 1 and is_prime(p):
            out.append(p)
        p += 1
    return out


def gauss_prime(p):
    for b in range(1, math.isqrt(p) + 1):
        a2 = p - b * b
        a = math.isqrt(a2)
        if a * a == a2 and a > 0:
            return (max(a, b), min(a, b))
    raise ValueError(p)


def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gconj(x):
    return (x[0], -x[1])


def gnorm(x):
    return x[0] * x[0] + x[1] * x[1]


def gpow(x, e):
    assert e >= 0
    r, b = (1, 0), x
    while e:
        if e & 1:
            r = gmul(r, b)
        b = gmul(b, b)
        e >>= 1
    return r


def gsignpow(pi, e):
    """pi^(e) with the convention pi^(-k) := conj(pi)^k (a Gaussian INTEGER)."""
    return gpow(pi, e) if e >= 0 else gpow(gconj(pi), -e)


UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]


def gprod(lst):
    r = (1, 0)
    for t in lst:
        r = gmul(r, t)
    return r


def conj_primitive_part(exps, pis):
    """Given exponents e_j (ints) on distinct split Gaussian primes pi_j, return the conj-primitive
    Gaussian integer prod pi_j^(e_j) (no rational factor) -- exact."""
    return gprod([gsignpow(pi, e) for pi, e in zip(pis, exps) if e != 0])


# ---------------------------------------------------------------- Hadamard matrices
def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def paley_I(q):
    """Normalised (first column all ones) Paley type I Hadamard matrix of order q+1, q=3 mod 4 prime.
    Row 0 = infinity, rows 1..q = elements x=0..q-1 of F_q.  Column 0 constant, column 1+a = label a."""
    assert q % 4 == 3 and is_prime(q)
    M = q + 1
    H = [[1] * M]
    for x in range(q):
        row = [1]
        for a in range(q):
            row.append(-legendre(a - x, q) - (1 if a == x else 0))
        H.append(row)
    check_hadamard(H)
    return H


def sylvester(t):
    M = 2 ** t
    H = [[(-1) ** bin(x & a).count("1") for a in range(M)] for x in range(M)]
    check_hadamard(H)
    return H


def check_hadamard(H):
    M = len(H)
    for x in range(M):
        assert H[x][0] == 1
        for y in range(M):
            s = sum(H[x][k] * H[y][k] for k in range(M))
            assert s == (M if x == y else 0), (x, y, s)
    return True


# ---------------------------------------------------------------- profiles
class Profile:
    """Two-valued Hadamard-core profile with flips.

    columns: list of dicts {label a (1..M-1, or 0 for content-type), sigma (+-1), p, pi, F (set of flipped rows)}
    row x sign at column j: s_xj = sigma_j * H[x][a_j] * (-1)^[x in F_j].
    z_x = units[x] * prod_j pi_j^(s_xj)   (content g = 1).
    """

    def __init__(self, H, columns, units):
        self.H, self.cols, self.units = H, columns, units
        self.M = len(H)

    def sign(self, x, j):
        c = self.cols[j]
        s = c["sigma"] * self.H[x][c["label"]]
        return -s if x in c["F"] else s

    def point(self, x):
        z = UNITS[self.units[x] % 4]
        for j, c in enumerate(self.cols):
            z = gmul(z, gsignpow(c["pi"], self.sign(x, j)))
        return z

    def W(self):
        return sum(math.log(c["p"]) for c in self.cols)

    def flip_mass(self):
        return sum(len(c["F"]) * math.log(c["p"]) for c in self.cols)


def fully_flipped_profile(H, b, primes, assignment=None, rng=None, orient_random=True):
    """b copies per nonconstant label; row x flipped at one copy of label assignment[x];
    capacity: copies 0,1,... used in order (distinct physical columns).  primes: list of >= b(M-1) split primes."""
    rng = rng or random.Random(1)
    M = len(H)
    if assignment is None:
        # canonical: row x>=1 -> label x, row 0 -> label 1
        assignment = [1] + list(range(1, M))
    cols = []
    idx = 0
    colmap = {}
    for a in range(1, M):
        for c in range(b):
            p = primes[idx]
            idx += 1
            cols.append({"label": a, "sigma": rng.choice([1, -1]) if orient_random else 1,
                         "p": p, "pi": gauss_prime(p), "F": set()})
            colmap[(a, c)] = len(cols) - 1
    used = {}
    for x in range(M):
        a = assignment[x]
        c = used.get(a, 0)
        assert c < b, "capacity exceeded"
        used[a] = c + 1
        cols[colmap[(a, c)]]["F"].add(x)
    units = [rng.randrange(4) for _ in range(M)]
    return Profile(H, cols, units), assignment
