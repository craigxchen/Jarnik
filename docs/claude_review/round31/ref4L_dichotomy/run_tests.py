import random, sympy
from math import prod, isqrt, gcd
from itertools import combinations
from dich import analyse, edges, CYC
random.seed(31)

def iv2(x, P, q):
    s = [sum(prod(c) for c in combinations(x, k)) for k in range(1, 6)]
    s1, s2, s3, s4, s5 = s
    assert s4 == 0 and isqrt(-s1*s5)**2 == -s1*s5 and (3*s5) % 2 == 0
    l = isqrt(-s1*s5); b0 = 3*s5//2
    b = {e: x[e[0]]*x[e[1]]*(s1-x[e[0]]-x[e[1]]) for e in edges}
    hh = [t**3-s1*t*t+s2*t for t in x]
    U = [P+q*hi for hi in hh]; B = b0*q*q; L = l*q
    V = [(P*U[i]+B*x[i], L*U[i]) for i in range(5)]
    J = {e: P+q*b[e] for e in edges}
    return V, J

def rad(n): return prod(sympy.factorint(abs(n)).keys())

def iv2_instance(x, q_extra=1, digits=40):
    s = [sum(prod(c) for c in combinations(x, k)) for k in range(1, 6)]
    l = isqrt(-s[0]*s[4]); b0 = 3*s[4]//2
    q = rad(l*b0*prod(x))*q_extra
    while True:
        P = sympy.nextprime(random.randrange(10**digits, 10**(digits+1)))
        V, J = iv2(x, P, q)
        if all(gcd(J[e], J[f]) == 1 for e, f in combinations(edges, 2)):
            return V, J, P, q

print("== Test 1: IV.2 family, Astra's quintuple, several P, q ==")
x = (-27, -270, 135, 120, 40)
for digits, qe in [(30, 1), (60, 1), (60, 7), (120, 11*13)]:
    V, J, P, q = iv2_instance(x, qe, digits)
    print(" P ~ 10^%d, q = %d" % (digits, q))
    analyse(V, J, set(), "IV.2")

print("== Test 2: Astra's balanced numerical parent (IV.3 certificate) ==")
from sympy import Rational
k, step = 3999318414, 924
shifts = (145, 79, -31, 109, 85, -47, 113, 19, -41, 25)
origin = k*k+1-329*step
J = {e: origin+step*bb for e, bb in zip(edges, shifts)}
R = (4680, 3861, 2860, 3168, 2808); d = (-76, -55, -6, -34, -4)
f = (26095, Rational(43223, 2), 16715, 16813, 16303)
T = [int(origin**2+a*step*origin+bq*step**2) for a, bq in zip(d, f)]
V = [(k*T[i]-R[i], T[i]) for i in range(5)]
analyse(V, J, CYC, "IV.3 cert")
