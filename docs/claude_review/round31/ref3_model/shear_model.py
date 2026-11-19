"""Uncomplemented pair-cut models with unimodular shear (Astra III.11, 'why this construction stops')."""
import random
from math import gcd, prod
from itertools import combinations
from fractions import Fraction
from sympy import factorint
from sympy.ntheory.modular import crt
from mpmath import mp, mpf, atan2, log
for m in (5, 6):
    s = [3,7,11,19,23,31][:m]; u = [x*x for x in s]
    pairs = list(combinations(range(m), 2))
    base = 2*5*13*17*29*31*37*47 if m == 5 else 2*5*13*17*29*37*47
    t0, md = crt(s, [pow(base, -1, p) for p in s]); t0, md = int(t0), int(md)
    rng = random.Random(7)
    while True:
        K = base*(t0+md*rng.randrange(1, 10**6))
        J = {e: K*K+u[e[0]]*u[e[1]] for e in pairs}
        if all(gcd(J[a],J[b]) == 1 for a,b in combinations(pairs,2)) and \
           all(gcd(J[e], abs(u[a]-u[b])) == 1 for e in pairs for a,b in pairs):
            break
    N = prod(J.values()); mp.dps = len(str(N))+60
    n = [prod(J[e] for e in pairs if i in e) for i in range(m)]
    V = [(u[i] + K**(m-3)*(u[i]**2-K*K), u[i]**2-K*K) for i in range(m)]
    assert all(gcd(*v) == 1 for v in V)
    ok = all((V[i][0]*V[j][1]-V[i][1]*V[j][0]) == J[i,j]*(u[j]-u[i]) for i,j in pairs)
    rel = max(abs(float(Fraction(V[i][0]**2+V[i][1]**2, n[i])-1)) for i in range(m))
    ang = [2*atan2(mpf(y), mpf(x)) for x,y in V]
    diam = max(ang)-min(ang)
    print("m=%d K=%d: primitive, all %d minors = J_ij*(u_j-u_i): %s; max rel norm err %.2e;"
          " diam*N^(1/4)=%s; log diam/log N = %s; diam*K^(2(m-2)) = %s" % (
          m, K, len(pairs), ok, rel, mp.nstr(diam*mpf(N)**0.25, 6), mp.nstr(log(diam)/log(mpf(N)), 8),
          mp.nstr(diam*mpf(K)**(2*(m-2)), 8)))
    sepfrac = Fraction(2*(m-2), len(pairs))
    print("   separated fraction of blocks per pair = %s (needs > 1/2)" % sepfrac)
for m in range(5, 12):
    print(m, Fraction(2*(m-2), m*(m-1)//2) > Fraction(1,2))
