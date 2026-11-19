# Index of the coefficient lattice {Q : Q(V_i) = 0 mod n_i} and directions mod p, small exact instances.
import random, sympy
from math import prod, gcd
from itertools import combinations
from sympy.matrices.normalforms import smith_normal_form
from sympy import Matrix, ZZ
from dich import edges, CYC
from run_tests2 import gauss_prime, gmul, star, primes
random.seed(3)
def index(V, n):
    rows = [[v[0]**2 for v in V], [v[0]*v[1] for v in V], [v[1]**2 for v in V]]
    rows += [[n[i] if j == i else 0 for j in range(5)] for i in range(5)]
    snf = smith_normal_form(Matrix(rows), domain=ZZ)
    d = prod(abs(snf[i, i]) for i in range(5))
    return prod(n)//d
def dirs_mod(V, p):
    s = set()
    for v in V:
        if v[1] % p == 0: s.add('inf')
        else: s.add(v[0]*pow(v[1], -1, p) % p)
    return s
for F, name in [(set(), 'uncompl'), (CYC, 'balanced')]:
    ps = primes(4, {1}, 50, 10); random.shuffle(ps); J = dict(zip(edges, ps)); N = prod(ps)
    S = [{e for e in edges if (i in e) != (e in F)} for i in range(5)]
    n = [prod(J[e] for e in S[i]) for i in range(5)]
    pr = {e: gauss_prime(J[e]) for e in edges}
    V = star(F, pr, 1)
    ok = [index(V, n) == N]
    # a non-genuine model with the same divisibilities: V_i + n_i u_i
    for trial in range(200):
        Vp = [(V[i][0]+n[i]*random.randint(-5,5), V[i][1]+n[i]*random.randint(-5,5)) for i in range(5)]
        if all(gcd(*v) == 1 for v in Vp): ok.append(index(Vp, n) == N)
    print(name, 'N =', N, 'index == N in all tests:', ok)
    # number of distinct directions mod each block prime for genuine star
    print('  #directions mod p (genuine):', [len(dirs_mod(V, J[e])) for e in edges])
