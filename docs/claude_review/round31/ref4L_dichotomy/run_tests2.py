import random, sympy
from math import prod, gcd
from itertools import combinations
from dich import analyse, edges, CYC
random.seed(7)
def gauss_prime(p):  # p = a^2 + b^2
    for a in range(1, int(p**0.5)+1):
        b2 = p-a*a; b = int(round(b2**0.5))
        if b*b == b2: return (a, b)
def sqrtm2_prime(p):  # p = a^2 + 2 b^2
    for b in range(1, int((p/2)**0.5)+1):
        a2 = p-2*b*b; a = int(round(a2**0.5))
        if a*a == a2: return (a, b)
def gmul(u, v, c=1): return (u[0]*v[0]-c*u[1]*v[1], u[0]*v[1]+u[1]*v[0])
def star(F, pr, c):
    S = [{e for e in edges if (i in e) != (e in F)} for i in range(5)]
    V = []
    for i in range(5):
        z = (1, 0)
        for e in sorted(S[i]): z = gmul(z, pr[e], c)
        V.append(z)
    return V
def primes(mod, res, lo, k):
    out = []; p = lo
    while len(out) < k:
        p = sympy.nextprime(p)
        if p % mod in res: out.append(p)
    return out
for F, name in [(set(), "uncompl"), (CYC, "balanced")]:
    ps = primes(8, {1}, 10**6, 10); random.shuffle(ps)
    J = dict(zip(edges, ps))
    print("== genuine Gaussian star,", name)
    pr = {e: gauss_prime(J[e]) for e in edges}
    pr = {e: (v if random.random() < .5 else (v[1], v[0])) for e, v in pr.items()}  # random choice of pi vs conj
    V = star(F, pr, 1); analyse(V, J, F, "genuine")
    print("== stretched genuine star (X, 3Y),", name)
    dd = next(d for d in range(3, 200, 2) if all(gcd(v[0], d) == 1 for v in V)); print('  stretch d =', dd); Vs = [(v[0], dd*v[1]) for v in V]
    if all(gcd(*v) == 1 for v in Vs): analyse(Vs, J, F, "stretched")
    print("== Z[sqrt-2] star,", name)
    pr2 = {e: sqrtm2_prime(J[e]) for e in edges}
    V2 = star(F, pr2, 2); analyse(V2, J, F, "sqrt-2")
    print("== perturbed IV-type: V_i + n_i u_i on genuine star,", name)
    S = [{e for e in edges if (i in e) != (e in F)} for i in range(5)]
    n = [prod(J[e] for e in S[i]) for i in range(5)]
    Vp = [(V[i][0]+n[i]*random.randint(-3,3), V[i][1]+n[i]*random.randint(-3,3)) for i in range(5)]
    analyse(Vp, J, F, "perturbed")
