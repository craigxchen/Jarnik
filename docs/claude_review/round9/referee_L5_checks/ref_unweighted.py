#!/usr/bin/env python3
"""Referee check of L5.md section 1, Remark 1 ("the unweighted reading is false for every system").

(a) At q = 3, 5 every value of Z/m_q (m_q = 4) lies in T, so every monomial collides there: that is
    2 moduli, which an exception count O(||w||_1 + 1) with constant >= 2 absorbs.  Not a refutation.
(b) What does break the unweighted reading for GENUINE labels: a Gaussian prime pi = a + b i whose
    imaginary part b has many odd prime factors q.  Then pi = conj(pi) mod q, i.e. the single-label
    monomial (twin monomial (0, eps e_y) if pi sits on the flip column f(y)) collides with unit 1 at
    every such q.  The number of colliding moduli is unbounded (~ log p /(2 loglog p)), while the
    log q-weighted demand stays <= log b <= (1/2) log p, as the height bound requires.
Collisions are tested directly in (Z[i]/q)^*: pi * conj(pi)^(-1) == i^s.
"""
import random
random.seed(7)

def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if a % n == 0: continue
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def gmul(u, v, q): return ((u[0]*v[0] - u[1]*v[1]) % q, (u[0]*v[1] + u[1]*v[0]) % q)
def ginv(u, q):
    nrm = (u[0]*u[0] + u[1]*u[1]) % q
    t = pow(nrm, -1, q)
    return (u[0]*t % q, (-u[1])*t % q)
UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def unit_collisions(a, b, q):
    """s with (a+bi)/(a-bi) = i^s mod q (q odd prime, q not dividing a^2+b^2)."""
    z = gmul((a % q, b % q), ginv((a % q, (-b) % q), q), q)
    return [s for s, u in enumerate(UNITS) if (u[0] % q, u[1] % q) == z]

odd_primes = [q for q in range(3, 200) if is_prime(q)]

# (a): q = 3, 5: every element of T_q is a 4th root of unity, so every monomial value is a unit.
for q in (3, 5):
    ok = all(len(unit_collisions(a, b, q)) == 1 for a in range(q) for b in range(q)
             if (a*a + b*b) % q != 0)
    print('q=%d: every norm-one ratio is a unit (exactly one s):' % q, ok)

# (b): genuine Gaussian primes with many collisions at small moduli.
print('\n k  #moduli(unit 1)  log p   (1/2)log p  weighted demand sum_q log q   p')
import math
for k in range(1, 13):
    B = 1
    for q in odd_primes[1:k + 1]:    # skip 3 and use 5,7,11,... so q=3 is not even counted
        B *= q
    # find a with p = a^2 + B^2 prime, gcd(a, B) = 1
    a = 1
    while True:
        p = a*a + B*B
        if math.gcd(a, B) == 1 and is_prime(p):
            break
        a += 1
    coll = [q for q in odd_primes if p % q != 0 and 0 in unit_collisions(a, B, q)]
    wd = sum(math.log(q) for q in coll)
    print('%2d  %3d  %8.2f  %8.2f  %8.2f   %d' % (k, len(coll), math.log(p), math.log(p)/2, wd, p))
