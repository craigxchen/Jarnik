# Exact check for Lemma D.2 / GPA item 2: for every split prime p <= PMAX and pi = a+bi (a>b>0), find the
# Gaussian primes on the lines Im(alpha conj(pi)) = +1 and = -1 with Re(alpha conj(pi)) > 0 of smallest norm.
# alpha = x+iy:  Im(alpha conj pi) = a*y - b*x,  Re(alpha conj pi) = a*x + b*y.  Integer arithmetic only.
from sympy import isprime
import sys
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100000

def gauss_prime(x, y):
    n = x*x + y*y
    if x == 0 or y == 0:
        m = abs(x + y)
        return m % 4 == 3 and isprime(m)
    return isprime(n)

def egcd(u, v):
    if v == 0: return (u, 1, 0)
    g, s, t = egcd(v, u % v)
    return (g, t, s - (u // v) * t)

def rep(p):
    # a^2+b^2=p, a>b>0
    import math
    for b in range(1, int(math.isqrt(p // 2)) + 1):
        a2 = p - b*b; a = math.isqrt(a2)
        if a*a == a2: return (max(a, b), min(a, b))

worst = (0, None); worst_t = (0, None); count = 0; fails = []
for p in range(5, PMAX + 1, 4):
    if not isprime(p): continue
    a, b = rep(p); count += 1
    for m in (1, -1):
        # solve a*y - b*x = m : g = gcd(a,b) = 1
        g, s, t = egcd(a, -b)   # a*s + (-b)*t = g = +-1
        y0, x0 = s * m * g, t * m * g   # a*y0 - b*x0 = m
        assert a*y0 - b*x0 == m
        # general solution x = x0 + k*a, y = y0 + k*b ; Re = a*x+b*y = Re0 + k*p
        re0 = a*x0 + b*y0
        k = -(re0 // p) if re0 <= 0 else -(re0 // p)
        while a*(x0 + k*a) + b*(y0 + k*b) <= 0: k += 1
        found = None
        for step in range(0, 100000):
            x, y = x0 + (k + step)*a, y0 + (k + step)*b
            if gauss_prime(x, y):
                found = (x, y, step); break
        if found is None:
            fails.append((p, m)); continue
        x, y, step = found
        nrm = x*x + y*y
        assert a*y - b*x == m and a*x + b*y > 0
        if nrm / p > worst[0]: worst = (nrm / p, (p, m, x, y))
        if step > worst_t[0]: worst_t = (step, (p, m, x, y))
print("split primes p <= %d checked: %d (both lines each)" % (PMAX, count))
print("failures (no prime within 1e5 steps):", fails[:10], len(fails))
print("max Norm(varpi)/p over all p, lines:", worst)
print("max number of steps along the line:", worst_t)
