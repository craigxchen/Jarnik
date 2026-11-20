"""Core exact routines for the two-class problem (own code, exact integers).

q2(y): the largest q such that at every exact prime power p^f || q the integers y occupy at most
two classes mod p^f.  Formula: q2(y) = lcm over all bipartitions sigma of {0..B-1} (including the
trivial one) of g_sigma = gcd of the differences inside each side.  [At a prime p the maximal f with
<=2 classes is max_sigma v_p(g_sigma).]  q2_bruteforce checks this definition directly.
"""
from math import gcd, log
from fractions import Fraction
import itertools, sympy

def split_gcds(y):
    B = len(y); out = {}
    for m in range(1 << (B - 1)):
        sides = [0] + [(m >> (i - 1)) & 1 for i in range(1, B)]
        g = 0; rep = [None, None]
        for i in range(B):
            s = sides[i]
            if rep[s] is None: rep[s] = y[i]
            else: g = gcd(g, abs(y[i] - rep[s]))
        out[m] = g
    return out

def q2(y):
    q = 1
    for g in split_gcds(y).values():
        if g > 1: q = q * g // gcd(q, g)
    return q

def q2_bruteforce(y):
    diam = max(y) - min(y); q = 1
    for p in sympy.primerange(2, diam + 1):
        f = 0
        while len({v % p ** (f + 1) for v in y}) <= 2: f += 1
        q *= p ** f
    return q

def ratio(y):
    y = sorted(y); L = y[-1] - y[0]
    return log(q2(y)) / log(L)

# ---------- lattice tools ----------
def hnf_rows(rows):
    """Row-style HNF basis of the Z-span of integer rows (simple, exact)."""
    A = [list(r) for r in rows]; n = len(A[0]); basis = []; r0 = 0
    for c in range(n):
        piv = [i for i in range(r0, len(A)) if A[i][c] != 0]
        if not piv: continue
        while True:
            piv = [i for i in range(r0, len(A)) if A[i][c] != 0]
            i0 = min(piv, key=lambda i: abs(A[i][c]))
            A[r0], A[i0] = A[i0], A[r0]
            done = True
            for i in range(r0 + 1, len(A)):
                if A[i][c]:
                    qt = A[i][c] // A[r0][c]
                    A[i] = [a - qt * b for a, b in zip(A[i], A[r0])]
                    if A[i][c]: done = False
            if done: break
        r0 += 1
    return [r for r in A[:r0]]

def lll(B, delta=Fraction(99, 100)):
    B = [list(b) for b in B]; n = len(B)
    def dot(u, v): return sum(a * b for a, b in zip(u, v))
    def gso():
        Bs = []; mu = [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = dot(B[i], Bs[j]) / dot(Bs[j], Bs[j])
                v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu
    Bs, mu = gso(); k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            qq = round(mu[k][j])
            if qq:
                B[k] = [a - qq * b for a, b in zip(B[k], B[j])]
                Bs, mu = gso()
        if dot(Bs[k], Bs[k]) >= (delta - mu[k][k - 1] ** 2) * dot(Bs[k - 1], Bs[k - 1]):
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu = gso(); k = max(k - 1, 1)
    return B

def crt_idem(m, Q):
    """e = 1 mod m, 0 mod Q/m."""
    r = Q // m
    return (r * pow(r, -1, m)) % Q

def split_lattice(B, splits):
    """splits: list of (S, m) with S a subset of {1..B-1} (side not containing 0), m pairwise coprime.
    Returns HNF basis of Lambda = {y in Z^{B-1}: y = t_sigma chi_S mod m_sigma}."""
    Q = 1
    for _, m in splits: Q *= m
    rows = []
    for S, m in splits:
        e = crt_idem(m, Q)
        rows.append([e if i in S else 0 for i in range(1, B)])
    for i in range(B - 1):
        rows.append([Q if j == i else 0 for j in range(B - 1)])
    return hnf_rows(rows), Q

if __name__ == "__main__":
    import random
    random.seed(1)
    for _ in range(300):
        B = random.randint(3, 7)
        y = sorted(random.sample(range(0, 60), B))
        assert q2(y) == q2_bruteforce(y), y
    # check exhaustive optima found by exh.c
    for y in [(0,104,405,860),(0,115,444,1410),(0,15,48,55,88,103)]:
        print(y, q2(y), q2_bruteforce(y), ratio(y))
    print("q2 formula agrees with brute force on 300 random tuples")
