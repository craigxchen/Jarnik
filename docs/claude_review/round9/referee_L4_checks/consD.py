"""Referee's independent implementation of Construction D (L4.md §5) and of genuine labels.

Labels are stored at the top level Q = q^(A_q) as integers mod m_Q; lower levels are reductions.
"""
import math, itertools
from glib import *

S0 = (17, 31, 47)

def v2(n):
    c = 0
    while n % 2 == 0:
        n //= 2; c += 1
    return c

def lemma51_map(N, b):
    """Odd encoding iota_b: Z/(N/4) -> Z/N of Lemma 5.1 (N = 0 mod 4, N >= 16)."""
    n4 = N // 4
    targets = [t for t in range(N) if t % 2 == b and t not in (0, N // 2) and t < (-t) % N]
    seen, src = set(), []
    for x in range(n4):
        if x in seen:
            continue
        y = (-x) % n4
        seen.add(x); seen.add(y)
        src.append((x, y))
    assert len(src) <= len(targets), (N, b, len(src), len(targets))
    mp = {}
    for (x, y), t in zip(src, targets):
        mp[x] = t
        if y != x:
            mp[y] = (-t) % N
    return mp

def check_lemma51(N, b):
    mp = lemma51_map(N, b)
    n4 = N // 4
    for x in range(n4):
        t = mp[x]
        if t % 2 != b or t in (0, N // 2):
            return False
    for x in range(n4):
        for y in range(n4):
            for s in (1, -1):
                if mp[x] == (s * mp[y]) % N and x != (s * y) % n4:
                    return False
    return True

class System:
    """Holds tori and label vectors for a profile."""
    def __init__(self, primes, X):
        self.p = list(primes)
        self.r = len(primes)
        self.X = X
        self.qs = [q for q in sieve(X) if q > 2 and q not in self.p]
        self.T = {}
        for q in self.qs:
            A = int(math.floor(math.log(X) / math.log(q) + 1e-12))
            while q ** (A + 1) <= X:
                A += 1
            while q ** A > X:
                A -= 1
            self.T[q] = Torus(q, A)
        self.nu = {q: [1 if legendre(pj, q) == -1 else 0 for pj in self.p] for q in self.qs}

    def coherent(self, zetas):
        """labels log(zeta_j / conj zeta_j) at each q."""
        L = {}
        for q in self.qs:
            T = self.T[q]
            row = []
            for z in zetas:
                assert (z[0] ** 2 + z[1] ** 2) % q != 0, ("zeta not a unit", q)
                row.append(T.log(ratio(z, T.Q)))
            L[q] = row
        return L

    def construct_D(self, zetas):
        L = self.coherent(zetas)
        ell, case = {}, {}
        for q in self.qs:
            T = self.T[q]
            m, chi = T.m, T.chi
            qm = q - chi
            qa = q ** (T.A - 1)
            nu = self.nu[q]
            row = []
            if q % 8 in (3, 5):
                assert T.iota % 2 == 1
                for j, Lj in enumerate(L[q]):
                    w = (Lj - nu[j]) % 2
                    row.append((Lj + w * T.iota) % m)
                case[q] = 1
            elif q in S0:
                for j, Lj in enumerate(L[q]):
                    assert Lj % 2 == nu[j], ("S0 parity", q, j)
                    row.append(Lj)
                case[q] = 2
            else:
                b = v2(qm)
                Nq = qm // 4
                for j, Lj in enumerate(L[q]):
                    p1 = Lj % qm
                    pq = Lj % qa if qa > 1 else 0
                    if b == 3:
                        o = qm // 8
                        # pi_1 = (nu mod 8, p1 mod o) by CRT in Z/qm
                        t = crt2(nu[j], 8, p1 % o, o)
                        case[q] = 3
                    else:
                        mp = lemma51_map(Nq, nu[j])
                        x = p1 % (Nq // 4)
                        t = mp[x]  # lift: the representative in [0, Nq) of Z/qm
                        case[q] = 4
                    row.append(crt2(t, qm, pq, qa) if qa > 1 else t % m)
            ell[q] = row
        return ell, case, L

def crt2(a, m1, b, m2):
    if m2 == 1:
        return a % m1
    t = ((b - a) * pow(m1, -1, m2)) % m2
    return (a + m1 * t) % (m1 * m2)

def collision_profile(sysm, ell, v, qs=None):
    """returns dict s -> list of (q, a) with a = a_{q,s}(v) > 0."""
    out = {0: [], 1: [], 2: [], 3: []}
    for q in (qs if qs is not None else sysm.qs):
        T = sysm.T[q]
        S = sum(c * ell[q][j] for j, c in v.items())
        best = None
        for a in range(1, T.A + 1):
            ma = (q - T.chi) * q ** (a - 1)
            io = T.iota % ma
            Sa = S % ma
            hit = None
            for s in range(4):
                if Sa == (s * io) % ma:
                    hit = s
                    break
            if hit is None:
                break
            best = (hit, a)
        if best is not None:
            out[best[0]].append((q, best[1]))
    return out

def mod_from(lst):
    m = 1
    for q, a in lst:
        m *= q ** a
    return m

def Zv(zetas, v):
    W = (1, 0)
    for j, c in v.items():
        z = zetas[j] if c > 0 else gconj(zetas[j])
        W = gmul(W, gpow(z, abs(c)))
    return W

def monomials(r, k):
    """all v (dict) with ||v||_1 = k, up to global sign (first nonzero coefficient > 0)."""
    res = []
    for supp_size in range(1, k + 1):
        for supp in itertools.combinations(range(r), supp_size):
            for comp in compositions(k, supp_size):
                for signs in itertools.product((1, -1), repeat=supp_size - 1):
                    v = {supp[0]: comp[0]}
                    for jj, c, sg in zip(supp[1:], comp[1:], signs):
                        v[jj] = sg * c
                    res.append(v)
    return res

def compositions(k, n):
    if n == 1:
        if k >= 1:
            yield (k,)
        return
    for first in range(1, k - n + 2):
        for rest in compositions(k - first, n - 1):
            yield (first,) + rest
