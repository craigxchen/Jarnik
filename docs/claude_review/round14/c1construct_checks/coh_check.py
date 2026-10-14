# Checks for Lemma L and Prop. G.
# (1) Lemma L (evidence for the size of Linnik's prime): for n = p * Norm(zeta) not a square / twice a square,
#     the least prime q = -1 mod 2^(E+1), q not dividing n, with (n/q) = -1.
# (2) Prop G(b) (exact): for xi = pi * sigma^2 (sigma a Gaussian prime of norm l != p), max(|D_1|,|D_3|)(xi)
#     exceeds sqrt(2p) * g, g = N-part with N = {p} (here only p itself can divide; we report the ratio),
#     i.e. an exactly coherent twisted column always violates the norm-1 per-unit height at some s in {1,3}
#     once its D_1, D_3 have no p_l-factors.  Also check D_1^2 + D_3^2 = 2|xi|^2 and oddness.
import sys, math, random
sys.path.insert(0, '.')
from lib14 import *
from sympy import primerange, factorint

random.seed(5)
print("(1) least q = -1 mod 2^(E+1), q prime, q not | n, (n/q) = -1, for n = p*Norm(zeta):")
worst = {}
for E in [2, 4, 6, 8]:
    mx = 0; mxn = None
    for trial in range(300):
        p = random.choice(split_primes(10**3, 10**5))
        l = random.choice(split_primes(5, 10**4))
        n = p * l
        mod = 2 ** (E + 1)
        q = mod - 1
        while True:
            if q > 2 and all(q % d for d in range(2, int(q**0.5) + 1)) and n % q and legendre_symbol(n % q, q) == -1:
                break
            q += mod
        if q > mx: mx, mxn = q, n
    print(f"   E={E}: max over 300 random n ~ 1e3..1e9 of the least such q: {mx} (n={mxn}); 2^(E+1)={2**(E+1)}")
from sympy import legendre_symbol

print("(2) twisted exact coherence xi = pi*sigma^2: norm-1 height at s=1,3")
viol = 0; tot = 0; minratio = 1e9
for p in split_primes(10**3, 2*10**4)[:300]:
    a, b = gauss_prime(p)
    for l in split_primes(5, 200):
        if l == p: continue
        sg = gauss_prime(l)
        for sig in [sg, gconj(sg)]:
            xi = gmul((a, b), gmul(sig, sig))
            D1, D3 = Ds(xi, 1), Ds(xi, 3)
            assert D1 % 2 and D3 % 2 and D1*D1 + D3*D3 == 2 * (xi[0]**2 + xi[1]**2)
            g = 1
            for x in (D1, D3):
                pass
            # N-part relative to the column's own prime set {p, l}: only p could divide (l cannot: primitivity)
            gp = max(abs(D1), abs(D3))
            while gp % p == 0: gp //= p
            ratio = gp / math.sqrt(2 * p)
            tot += 1
            if ratio > 1: viol += 1
            minratio = min(minratio, ratio)
print(f"   {tot} columns (p in [1e3,2e4], sigma of norm <= 200, both orientations): height violated at s=1 or 3 in {viol}; min ratio max|D_s|/sqrt(2p) = {minratio:.3f}")
