"""Referee checks of p3Dlower Theorem 1.1 and Lemma 5.1 (independent code, exact).
A. Thm 1.1(1): at q = +-3 mod 8, ell = ell' + iota*u permutes units by <v,u> (all levels a=1 here).
A'. NEW: at q = +-3 mod 8 a parity-respecting ell forces the PARITY of the collision unit: s = <v,nu> mod 2.
B. Thm 1.1(2): exact enumeration over all parity-respecting ell in (Z/m)^2, q = +-1 mod 8 small.
C. NEW counterexample to 'only as a veto': v = (m_q/8) e_j collides for EVERY parity-respecting ell when nu_j(q)=0,
   while the unconstrained ell_j = 1 avoids it.  Genuine confirmation at q = 7: every p = 1 mod 4, (p/7)=1, p<20000
   has 7 | D_s(pi) for some s."""
import random, itertools, math
from sympy import primerange, legendre_symbol, isprime
from tq import T, two_sq, D
random.seed(7); ok = True
# A, A'
nA = 0
for q in [q for q in primerange(3, 140) if q % 8 in (3, 5)]:
    t = T(q); m = t.m; io = t.iota
    for _ in range(300):
        r = random.randint(1, 6)
        ellp = [random.randrange(m) for _ in range(r)]; nu = [random.randrange(2) for _ in range(r)]
        u = [(ellp[j] - nu[j]) % 2 for j in range(r)]
        ell = [(ellp[j] + io * u[j]) % m for j in range(r)]
        assert all(ell[j] % 2 == nu[j] for j in range(r))
        v = [random.randint(-4, 4) for _ in range(r)]
        if not any(v): continue
        a = sum(x * y for x, y in zip(v, ell)) % m; b = sum(x * y for x, y in zip(v, ellp)) % m
        sa, sb = t.unit(a), t.unit(b)
        vu = sum(x * y for x, y in zip(v, u))
        ok &= ((sa is None) == (sb is None)) and (sa is None or sa == (sb + vu) % 4)
        if sa is not None:
            ok &= (sa % 2 == sum(x * y for x, y in zip(v, nu)) % 2)       # A'
        nA += 1
print("A/A' (Thm 1.1(1) + forced unit parity at q=+-3 mod 8):", ok, "cases", nA)
# B
okB = True
for q in [7, 17, 23, 31, 41]:
    t = T(q); m = t.m
    for nu in itertools.product((0, 1), repeat=2):
        P = [(a, b) for a in range(nu[0], m, 2) for b in range(nu[1], m, 2)]
        for v in itertools.product(range(-5, 6), repeat=2):
            if math.gcd(math.gcd(v[0], v[1]), m // 2) != 1: continue
            hits = sum(1 for (a, b) in P if t.unit(v[0] * a + v[1] * b) is not None)
            par = (v[0] * nu[0] + v[1] * nu[1]) % 2
            exp = 0 if par else 8 * len(P) // m
            okB &= (hits * m == 8 * len(P)) if par == 0 else (hits == 0)
print("B (Thm 1.1(2) exact probability 8/m or 0, r=2):", okB)
# C
okC = True
for q in [7, 17, 23, 31, 41, 47, 71, 73, 79]:
    t = T(q); m = t.m; c = m // 8
    forced = all(t.unit(c * l) is not None for l in range(0, m, 2))      # nu_j = 0: all even labels
    avoid = t.unit(c * 1) is None                                       # unconstrained label 1
    okC &= forced and avoid
    print(f"  q={q:3d} m={m:3d}: v=({c})e_j collides for all even labels: {forced}; label 1 avoids: {avoid}")
gen = [p for p in primerange(5, 20000) if p % 4 == 1 and legendre_symbol(p % 7, 7) == 1]
allhit = all(any(D(*two_sq(p), s) % 7 == 0 for s in range(4)) for p in gen)
print(f"C genuine q=7: {len(gen)} primes p=1 mod 4 with (p/7)=1 below 2e4; all have 7 | D_s(pi) for some s: {allhit}")
print("examples:", [(p, two_sq(p), [D(*two_sq(p), s) for s in range(4)]) for p in gen[:3]])
okC &= allhit
print("ALL OK" if ok and okB and okC else "FAIL")
