"""X.5(2a): squarefree 2<=q<=300, multisets of three positive lengths <= 2 sqrt q."""
from itertools import combinations_with_replacement
from math import isqrt, prod
from collections import Counter

def factor(n):
    f, p = [], 2
    while p*p <= n:
        e = 0
        while n % p == 0:
            n //= p; e += 1
        if e: f.append((p, e))
        p += 1
    if n > 1: f.append((n, 1))
    return f

def H(a, b, c):
    return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)

def admissible(primes, e):
    """signs free at each prime: all three = 0 mod p (constant cut), or exactly one = 0 and the two others = +-each other mod p."""
    g = 1
    for p in primes:
        z = [x % p == 0 for x in e]
        if all(z):
            g *= p; continue
        if sum(z) != 1:
            return None
        u, v = [x for x in e if x % p]
        if (u - v) % p and (u + v) % p:
            return None
    return g

cnt = Counter(); nonflat = []; capped = []
for q in range(2, 301):
    f = factor(q)
    if any(e > 1 for _, e in f):
        continue
    cnt["moduli"] += 1
    primes = [p for p, _ in f]
    top = isqrt(4*q)
    for e in combinations_with_replacement(range(1, top+1), 3):
        g = admissible(primes, e)
        if g is None: continue
        cnt["admissible"] += 1
        h = H(*e)
        assert h % (q*g)**2 == 0
        T = h // (q*g)**2
        Q = q // g
        assert T % 4 in (0, 3)
        if Q % 2 == 0: assert T % 4 == 0
        if h:
            cnt["nonzero H"] += 1
            nonflat.append((q, e, g, T))
            assert max(e)**2 >= q*g
        if max(e)**2 <= q:
            cnt["admissible, max^2<=q"] += 1
            capped.append((q, e))
            assert h == 0
        if max(e)**2 < q*g:
            assert h == 0
print(dict(cnt))
print("nonzero-H triples (q, lengths, g, H/(qg)^2):")
for x in nonflat: print("  ", x)
print("admissible with max^2<=q:", capped)
