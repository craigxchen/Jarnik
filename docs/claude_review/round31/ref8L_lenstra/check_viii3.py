#!/usr/bin/env python3
"""Independent exact check of the divisor matrices (VIII.5) of Part VIII.3.
For n with prime-power layers, q, r with gcd(r,q)=1 and eight distinct divisors
d_0..d_7 of n all = r mod q, check:
 * h_ij prime to q; t_ij=(d_j-d_i)/(q h_ij) nonzero integer
 * columns P_i=(a/g_i,(d_i-a)/(q g_i)) integral, primitive, P_0=(1,0)
 * det(P_i,P_j) = t_ij * prod_{T contains i,j} n_T, n_T = prod of p over layers (p,s)
   whose cut (slots differing from slot 0) is T
 * valuation formula b0+min(bi,bj)-min(b0,bi)-min(b0,bj) = #layers with both opposite
 * chain condition: at each p the T with p|n_T form two chains with disjoint unions
 * D_ij (layer mass) = log(d_i d_j/h_ij^2) exactly as integer identity
 * log identity  log|t| = D/2 - log q + log(|dj-di|/sqrt(di dj))  (squared, exact)
 * slopes Y_i/X_i = (d_i/a - 1)/q (common affine image)
"""
import random, itertools
from math import gcd, prod
from fractions import Fraction
from sympy import factorint, divisors, isprime, nextprime

def vp(x, p):
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v

def is_chain(sets):
    s = sorted(sets, key=len)
    return all(a <= b for a, b in zip(s, s[1:]))

def two_chains_disjoint(sets):
    # brute force: partition the family into two chains with disjoint unions
    sets = list(set(sets))
    m = len(sets)
    for mask in range(1 << m):
        A = [sets[i] for i in range(m) if mask >> i & 1]
        B = [sets[i] for i in range(m) if not mask >> i & 1]
        if is_chain(A) and is_chain(B):
            uA = frozenset().union(*A) if A else frozenset()
            uB = frozenset().union(*B) if B else frozenset()
            if not (uA & uB):
                return True
    return False

def check(n, q, ds):
    F = factorint(n)
    assert len(set(ds)) == 8 and all(n % d == 0 for d in ds)
    r = ds[0] % q
    assert gcd(r, q) == 1 and all(d % q == r for d in ds)
    a = ds[0]
    g = [gcd(a, d) for d in ds]
    X = [a // g[i] for i in range(8)]
    Y = []
    for i in range(8):
        num = ds[i] - a
        assert num % (q * g[i]) == 0
        Y.append(num // (q * g[i]))
    assert (X[0], Y[0]) == (1, 0)
    for i in range(8):
        assert gcd(X[i], Y[i]) == 1, "primitivity"
        # slope = (d_i/a - 1)/q
        if i:
            assert Fraction(Y[i], X[i]) == (Fraction(ds[i], a) - 1) / q
    # layers and cuts
    nT = {}
    layer_cut = {}
    for p, e in F.items():
        b = [vp(d, p) for d in ds]
        for s in range(1, e + 1):
            sign = [1 if b[i] >= s else 0 for i in range(8)]
            T = frozenset(i for i in range(1, 8) if sign[i] != sign[0])
            layer_cut[(p, s)] = T
            if T:
                nT[T] = nT.get(T, 1) * p
    # chain condition per prime
    for p, e in F.items():
        fam = [layer_cut[(p, s)] for s in range(1, e + 1) if layer_cut[(p, s)]]
        assert two_chains_disjoint(fam), ("chain", p, fam)
        # Astra's explicit chains
        b = [vp(d, p) for d in ds]
        low = [layer_cut[(p, s)] for s in range(1, b[0] + 1)]
        high = [layer_cut[(p, s)] for s in range(b[0] + 1, e + 1)]
        assert is_chain(low) and is_chain(high)
        assert all(T <= {i for i in range(1, 8) if b[i] < b[0]} for T in low)
        assert all(T <= {i for i in range(1, 8) if b[i] > b[0]} for T in high)
    for i, j in itertools.combinations(range(8), 2):
        h = gcd(ds[i], ds[j])
        assert gcd(h, q) == 1
        diff = ds[j] - ds[i]
        assert diff % (q * h) == 0
        t = diff // (q * h)
        assert t != 0
        det = X[i] * Y[j] - X[j] * Y[i]
        blk = prod(nT[T] for T in nT if i in T and j in T) if i > 0 else 1
        assert det == t * blk, (i, j)
        assert det * q * g[i] * g[j] == t * a * h * q   # (VIII.5) second form
        if i == 0:
            assert det == Y[j]
        # valuation formula
        for p, e in F.items():
            b0, bi, bj = vp(a, p), vp(ds[i], p), vp(ds[j], p)
            formula = b0 + min(bi, bj) - min(b0, bi) - min(b0, bj)
            cnt = sum(1 for s in range(1, e + 1)
                      if i in layer_cut[(p, s)] and j in layer_cut[(p, s)])
            if i == 0:
                assert formula == 0 and cnt == 0
            assert formula == cnt == vp(a * h // (g[i] * g[j]), p)
        # D_ij as exact exponent: prod p^{#differing layers} = d_i d_j / h^2
        Dexp = 1
        for p, e in F.items():
            Dexp *= p ** sum(1 for s in range(1, e + 1)
                             if (vp(ds[i], p) >= s) != (vp(ds[j], p) >= s))
        assert Dexp * h * h == ds[i] * ds[j]
        # log identity squared: t^2 = (e^D) / q^2 * diff^2/(di dj)
        assert Fraction(t * t) == Fraction(Dexp) / q**2 * Fraction(diff * diff, ds[i] * ds[j])
    return True

def random_instance(rng):
    # choose q, then n as product of prime powers; then find divisors in a class
    while True:
        q = rng.choice([3, 4, 5, 7, 8, 9, 11, 12, 13, 16, 24, 30, 60])
        ps = []
        cand = [p for p in range(2, 200) if isprime(p) and q % p]
        k = rng.randint(3, 6)
        F = {p: rng.randint(1, 4) for p in rng.sample(cand, k)}
        # adversarially add a prime dividing q? not allowed in divisors coprime to q;
        # but n may contain primes of q: divisors in invertible class avoid them
        if rng.random() < 0.3:
            pq = min(factorint(q))
            F[pq] = rng.randint(1, 3)
        n = prod(p**e for p, e in F.items())
        D = divisors(n)
        classes = {}
        for d in D:
            if gcd(d, q) == 1:
                classes.setdefault(d % q, []).append(d)
        big = [v for v in classes.values() if len(v) >= 8]
        if big:
            v = rng.choice(big)
            return n, q, rng.sample(v, 8)

def adversarial(rng):
    # deep prime-power layers, divisors with very different valuations at one prime
    q = 7
    n = 2**12 * 3**9 * 5**6 * 11**3
    D = [d for d in divisors(n) if d % 7 == 1]
    return n, q, rng.sample(D, 8)

rng = random.Random(20261119)
cnt = 0
for _ in range(800):
    n, q, ds = random_instance(rng); check(n, q, ds); cnt += 1
for _ in range(300):
    n, q, ds = adversarial(rng); check(n, q, ds); cnt += 1
# all orderings of one fixed set (anchor in every slot)
n = 2**6 * 3**5 * 5**3 * 13**2
q = 11
D = [d for d in divisors(n) if d % 11 == 4]
base = D[:8] if len(D) >= 8 else None
if base:
    for perm in itertools.islice(itertools.permutations(base), 0, 40320, 97):
        check(n, q, list(perm)); cnt += 1
print("all checks passed:", cnt, "eight-tuples")
