"""X.5(2c): weighted PTE families; separated-square sets and the triple divisor formula."""
from itertools import combinations, product
from math import gcd, prod, isqrt

def rad_factors(n):
    n = abs(n); out = set(); p = 2
    while p*p <= n:
        while n % p == 0: out.add(p); n //= p
        p += 1
    if n > 1: out.add(n)
    return out

def H(a, b, c):
    return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)

def check(label, marks, primes):
    primes = sorted(primes)
    assert len(set(marks)) == len(marks)
    q = prod(primes)
    for p in primes:
        cls = {y % p for y in marks}
        assert len(cls) <= 2, (label, p, len(cls))
        sq = set(); 
        for i, j in combinations(range(len(marks)), 2):
            d = marks[j] - marks[i]
            if d % p: sq.add(d*d % p)
            else: assert marks[i] % p == marks[j] % p
        assert len(sq) <= 1
        if len(cls) == 2:
            r, s = cls; assert sq == {(r-s)**2 % p}
    n = 0
    for x, y, z in combinations(sorted(marks), 3):
        e = (y-x, z-y, z-x)
        g = prod(p for p in primes if all(t % p == 0 for t in e))
        assert g == gcd(q, gcd(*e))
        assert H(*e) == 0
        # the sign-free local test of X.1 at each prime
        for p in primes:
            zr = [t % p == 0 for t in e]
            assert all(zr) or (sum(zr) == 1 and len({t*t % p for t in e if t % p}) == 1)
        n += 1
    L = max(marks) - min(marks)
    import math
    print(f"{label}: {len(marks)} marks, {len(primes)} primes, {n} triples, log q/log L = {math.log(q)/math.log(L):.3f}")
    return n

def node_family(label, nodes, r, T):
    # nodes: (a, b) meaning node b/a with moment vector (a, b); modulus a T - b; coefficient a^2
    mods = [a*T - b for a, b in nodes]
    tot = prod(mods)
    marks = []
    for eps in product((0, 1), repeat=len(nodes)):
        if sum(e*a for e, (a, b) in zip(eps, nodes)): continue
        if r == 2 and sum(e*b for e, (a, b) in zip(eps, nodes)): continue
        marks.append(sum(e*a*a*(tot//k) for e, (a, b), k in zip(eps, nodes, mods)))
    primes = set().union(*(rad_factors(k) for k in mods))
    return check(label, marks, primes)

def maximal_modulus(marks):
    # Lemma 1.1 of ref8L_attack: the two-class modulus is lcm over bipartitions of internal gcds; take its radical
    primes = set()
    for flags in product((0, 1), repeat=len(marks)-1):
        bits = (0,) + flags; g = 0
        for i, j in combinations(range(len(marks)), 2):
            if bits[i] == bits[j]: g = gcd(g, marks[i]-marks[j])
        if g > 1: primes |= rad_factors(g)
    return primes

if __name__ == "__main__":
    tot = 0
    T = 1009
    for k in (4, 5, 8):
        w = [1]*((k+1)//2) + [-1]*(k//2)
        tot += node_family(f"{k}-node one-moment family", [(x, x*j) for j, x in enumerate(w)], 1, T)
    vec = [(-3, -3), (-3, -2), (0, -1), (1, 2), (1, 3), (2, 0), (2, 1)]
    tot += node_family("7-node two-moment family", [(a+7*b, b) for a, b in vec], 2, T)
    a, b = 1009, 431
    fams = {"parallelogram": [0, a, b, a+b],
            "5 projective marks": [0, a*(a+b), b*(a+2*b), a*(a+b)+b*(a+2*b), (a+b)*(a+2*b)],
            "6 projective marks": [0, (a-b)*(a-2*b), a*(a-2*b), -2*b*(a-b), -a*b, a*a-4*a*b+2*b*b]}
    nd = (1009, 1013, 1019)
    fams["cube"] = [sum((prod(nd)//v)*e for v, e in zip(nd, bits)) for bits in product((0, 1), repeat=3)]
    for lab, mk in fams.items():
        tot += check(lab + " (maximal modulus, radical)", mk, maximal_modulus(mk))
    print("total triples", tot)
