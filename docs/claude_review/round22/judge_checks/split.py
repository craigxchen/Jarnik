"""Judge check J4 (exact): the split analogue of A_5 on the split Losev-Manin space LM_k over Q.
(1) For a primitive integer k-tuple u (gcd 1, all u_j > 0), the toric -K height of [u] in G_m^k/G_m, with the
    root-polytope Weil functions lambda_v = max_{a,b} log|u_a/u_b|_v, equals (max u / min u) * lcm(u). Checked
    exactly on random tuples by computing the product over all places directly.
(2) Balanced split block tuples u_j = prod_{A ni j} n_A (A proper nonempty subsets of [k], n_A distinct primes):
    gcd(u) = 1, lcm = prod n_A, gcd(u_a,u_b) = prod_{A contains a,b} n_A. Exponent bookkeeping at n_A = e^w:
    U = e^{(2^{k-1}-1)w}, L = e^{(2^k-2)w}, pair gcd e^{(2^{k-2}-1)w}; Liouville |u_a-u_b| >= gcd gives
    dist >= L^{-1/alpha}, alpha = (2^k-2)/2^{k-2} = 4(1-2^{1-k}): the same balanced Liouville constant as the twisted
    form (15/4 at k=5). Verified on actual prime choices (exact integers).
"""
import random, itertools
from math import gcd, lcm, log
from fractions import Fraction
from sympy import factorint, prime

def toric_height(u):
    # archimedean factor
    H = Fraction(max(u), min(u))
    ps = set()
    for x in u: ps |= set(factorint(x).keys())
    for p in ps:
        v = [factorint(x).get(p, 0) for x in u]
        # max_{a,b} |u_a/u_b|_p = p^{max_b v_b - min_a v_a}
        H *= p ** (max(v) - min(v))
    return H

random.seed(5)
ok = 0
for t in range(300):
    k = random.randint(3, 6)
    while True:
        u = [random.randint(1, 10**6) for _ in range(k)]
        g = 0
        for x in u: g = gcd(g, x)
        if g == 1: break
    H = toric_height(u)
    L = 1
    for x in u: L = lcm(L, x)
    assert H == Fraction(max(u), min(u)) * L, (u, H, L)
    ok += 1
print("(1) toric -K height == (max/min)*lcm on", ok, "random primitive tuples: OK")

for k in (3, 4, 5, 6):
    subs = [A for r in range(1, k) for A in itertools.combinations(range(k), r)]
    n = {A: prime(100 + i) for i, A in enumerate(subs)}
    u = [1] * k
    for A in subs:
        for j in A: u[j] *= n[A]
    g = 0
    for x in u: g = gcd(g, x)
    L = 1
    for x in u: L = lcm(L, x)
    assert g == 1 and L == eval('*'.join(str(n[A]) for A in subs))
    for a, b in itertools.combinations(range(k), 2):
        pg = 1
        for A in subs:
            if a in A and b in A: pg *= n[A]
        assert gcd(u[a], u[b]) == pg
    cU = 2 ** (k - 1) - 1; cL = 2 ** k - 2; cg = 2 ** (k - 2) - 1
    assert sum(1 for A in subs if 0 in A) == cU and len(subs) == cL and sum(1 for A in subs if 0 in A and 1 in A) == cg
    alpha = Fraction(cL, cU - cg)
    print(f"(2) k={k}: #A ni j = {cU}, #A = {cL}, #A contains a,b = {cg}; balanced Liouville alpha = {alpha} = 4(1-2^(1-k)) = {4*(1-Fraction(1,2**(k-1)))}")
