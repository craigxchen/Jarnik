"""Gaussian triangle checks for X.5(1). Exact integer arithmetic only."""
from itertools import combinations, product
from math import prod, isqrt
from collections import Counter
import random

def heron2(A, B, C):
    # H(a,b,c) expressed in the squares A=a^2, B=b^2, C=c^2
    return 2*(A*B + B*C + C*A) - A*A - B*B - C*C

def gmul(x, y):
    return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])

def conj(x):
    return (x[0], -x[1])

def gsub(x, y):
    return (x[0]-y[0], x[1]-y[1])

def norm(x):
    return x[0]*x[0] + x[1]*x[1]

def gdivides(d, x):
    # does Gaussian d divide x ?
    n = norm(d)
    y = gmul(x, conj(d))
    return y[0] % n == 0 and y[1] % n == 0

def gdiv(x, d):
    n = norm(d)
    y = gmul(x, conj(d))
    assert y[0] % n == 0 and y[1] % n == 0
    return (y[0]//n, y[1]//n)

def split_gen(p):
    for a in range(1, isqrt(p)+1):
        b2 = p - a*a
        b = isqrt(b2)
        if b*b == b2:
            return (a, b)
    raise ValueError

def is_prime(n):
    if n < 2: return False
    i = 2
    while i*i <= n:
        if n % i == 0: return False
        i += 1
    return True

def check_config(primes, eta=(1, 0), verbose=False):
    gens = {p: split_gen(p) for p in primes}
    q = prod(primes)
    marks = []
    for bits in product((0, 1), repeat=len(primes)):
        z = eta
        for p, b in zip(primes, bits):
            g = gens[p] if b == 0 else conj(gens[p])
            z = gmul(z, g)
        assert norm(z) == q
        marks.append((bits, z))
    stats = Counter()
    sep_sq = {p: set() for p in primes}
    sep_leg = {p: set() for p in primes}
    for (bi, zi), (bj, zj) in combinations(marks, 2):
        d = gsub(zj, zi)
        h = prod(p for p, x, y in zip(primes, bi, bj) if x == y)
        H = (1, 0)
        for p, x, y in zip(primes, bi, bj):
            if x == y:
                H = gmul(H, gens[p] if x == 0 else conj(gens[p]))
        # z_j - z_i = 2 i eta H t with t a rational integer
        den = gmul(gmul((0, 2), eta), H)
        t = gdiv(d, den)
        assert t[1] == 0 and t[0] != 0
        assert norm(d) == 4*h*t[0]**2
        for p, x, y in zip(primes, bi, bj):
            if x != y:
                assert t[0] % p != 0
                sep_sq[p].add(norm(d) % p)
                sep_leg[p].add(pow(norm(d) % p, (p-1)//2, p))
            else:
                assert norm(d) % p == 0
    for (b1, z1), (b2, z2), (b3, z3) in combinations(marks, 3):
        A, B, C = norm(gsub(z1, z2)), norm(gsub(z2, z3)), norm(gsub(z1, z3))
        H = heron2(A, B, C)
        det = (z2[0]-z1[0])*(z3[1]-z1[1]) - (z2[1]-z1[1])*(z3[0]-z1[0])
        # 16 area^2 = 4 det^2
        assert H == 4*det*det
        assert H*q == A*B*C            # (abc)^2 / q
        g = prod(p for p, x, y, w in zip(primes, b1, b2, b3) if x == y == w)
        ts = []
        for (bi, zi), (bj, zj) in (((b1, z1), (b2, z2)), ((b2, z2), (b3, z3)), ((b1, z1), (b3, z3))):
            h = prod(p for p, x, y in zip(primes, bi, bj) if x == y)
            ts.append(norm(gsub(zj, zi)) // (4*h))
        T2 = prod(ts)   # (t12 t23 t13)^2
        assert H == 64*g*g*T2
        assert abs(det) == 4*g*isqrt(T2) and isqrt(T2)**2 == T2
        assert H > 0
        stats["triples"] += 1
        # the consequence at a nonconstant prime: one agreeing pair (A'), two separated (B', C')
        for k, p in enumerate(primes):
            bits = (b1[k], b2[k], b3[k])
            if len(set(bits)) == 1:
                continue
            if b1[k] == b2[k]:
                Aa, Bb, Cc = A, B, C
            elif b2[k] == b3[k]:
                Aa, Bb, Cc = B, A, C
            else:
                Aa, Bb, Cc = C, A, B
            assert Aa % p == 0 and Bb % p and Cc % p
            lhs = (-(Bb - Cc)**2 * (q//p)) % p
            rhs = ((Aa//p) * Bb * Cc) % p
            assert lhs == rhs
            stats["prime-triple congruences"] += 1
            stats["separated squares equal"] += (Bb - Cc) % p == 0
    return stats, {p: len(v) for p, v in sep_sq.items()}, {p: len(v) for p, v in sep_leg.items()}, len(marks)

if __name__ == "__main__":
    # the norm-65 example
    z = [(4, 7), (8, -1), (8, 1)]
    assert gmul((2, 1), (3, 2)) == (4, 7) and gmul((2, 1), (3, -2)) == (8, -1) and gmul((2, -1), (3, 2)) == (8, 1)
    ch = {(i, j): norm(gsub(z[i], z[j])) for i, j in combinations(range(3), 2)}
    print("norm-65 chord squares", ch, "mod 5", {k: v % 5 for k, v in ch.items()}, "mod 13", {k: v % 13 for k, v in ch.items()})
    A, B, C = ch[0, 1], ch[1, 2], ch[0, 2]
    det = (z[1][0]-z[0][0])*(z[2][1]-z[0][1]) - (z[1][1]-z[0][1])*(z[2][0]-z[0][0])
    print("H", heron2(A, B, C), "area", abs(det)/2, "(abc)^2/q", A*B*C/65)
    splits = [p for p in range(5, 400) if p % 4 == 1 and is_prime(p)]
    random.seed(1)
    tests = [[5, 13], [5, 13, 17], [5, 13, 17, 29], [13, 17, 29, 37, 41]]
    for _ in range(6):
        tests.append(sorted(random.sample(splits, random.choice([3, 4, 5, 6]))))
    tot = Counter()
    for primes in tests:
        for eta in ((1, 0), (0, 1)):
            st, sq, leg, n = check_config(primes, eta)
            tot += st
        print(primes, "marks", n, st["triples"], "triples; distinct separated chord^2 mod p:", sq, "distinct Legendre symbols:", leg)
    print("total", dict(tot))
