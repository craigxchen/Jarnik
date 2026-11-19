"""Exact one-list fibre test in lexicographic K5 edge order, over a prime field.

Run python3 p3_unbarred.py for deterministic tests.  To test a list, use
python3 p3_unbarred.py PRIME u01 u02 u03 u04 u12 u13 u14 u23 u24 u34.
For kernel dimension >1, the test enumerates the base field's projective
kernel directions. A dimension-one kernel has just one geometric direction.
"""
from itertools import combinations, product
from math import prod
from random import Random
from sys import argv

from flint import nmod_mat, nmod_poly

EDGES = list(combinations(range(5), 2))
STARS = [[k for k,e in enumerate(EDGES) if i in e] for i in range(5)]


def star_data(p, u):
    assert len(u) == len(set(z % p for z in u)) == 10
    factors = [nmod_poly([-z, 1],p) for z in u]
    polynomials = [prod((factors[k] for k in star), start=nmod_poly([1],p)) for star in STARS]
    matrix = nmod_mat([[int(P[j]) for j in range(5)] for P in polynomials], p)
    basis, size = matrix.nullspace()
    return factors, polynomials, [[int(basis[i,j]) for i in range(5)] for j in range(size)]


def directions(p, basis):
    for pivot in range(len(basis)):
        for rest in product(range(p), repeat=len(basis)-1-pivot):
            c = (0,)*pivot + (1,) + rest
            yield [sum(x * row[j] for x,row in zip(c,basis)) % p for j in range(5)]


def reconstruct(p, u, factors, polys, a):
    R = [nmod_poly([sum(a[l-1]*int(P[k+l]) for l in range(1,5-k)) % p
                        for k in range(4)],p) for P in polys]
    v = []
    for edge,(i,j) in enumerate(EDGES):
        N, remainder = divmod(polys[i]*R[j] - polys[j]*R[i], factors[edge])
        assert not remainder and N.degree() <= 1
        if N.degree() != 1:
            return ('pair-leading-coefficient', edge)
        v.append(-int(N[0]) * pow(int(N[1]),-1,p) % p)
    if len(set(u+v)) < 20:
        return ('collision', v)
    scalars = []
    for edge,(i,j) in enumerate(EDGES):
        z = v[edge]
        if not R[i](z):
            assert not R[j](z) and polys[i](z) and polys[j](z)
            return ('no-finite-scalar', edge)
        scalar = -int(polys[i](z)) * pow(int(R[i](z)), -1, p) % p
        assert (polys[j]+scalar*R[j])(z) == 0
        scalars.append(scalar)
    if len(set(scalars)) != 1:
        return ('inconsistent-scalars', v, scalars)
    for i,star in enumerate(STARS):
        assert polys[i] + scalars[0]*R[i] == prod(
            (nmod_poly([-v[k],1],p) for k in star), start=nmod_poly([1],p))
    for m in range(1,6):
        sums = [sum(pow(u[k],m,p)-pow(v[k],m,p) for k in star)%p for star in STARS]
        assert len(set(sums)) == 1
    return ('solution',v,scalars[0])


def check_lift():
    from sympy import Matrix, Poly, symbols
    x,z = symbols('x z')
    edges = [35,76,91,75,13,88,93,5,3,z]
    polynomials = [Poly(prod(x-edges[k] for k in star),x) for star in STARS]
    det = Poly(Matrix([[P.nth(k) for k in range(5)] for P in polynomials]).det(),z)
    assert det == Poly(23040*(690809615*z*z-85455679574*z+2146786577637),z)
    assert det.eval(92)%101 == 0 and det.diff().eval(92)%101 == 39


def test(p,u):
    u = [z%p for z in u]
    factors,polys,basis = star_data(p,u)
    print('p',p,'unbarred',u,'kernel',basis)
    for a in directions(p,basis):
        print('direction',a,'result',reconstruct(p,u,factors,polys,a))


if __name__ == '__main__':
    if len(argv)>1:
        test(int(argv[1]), list(map(int,argv[2:])))
    else:
        check_lift()
        for p in [101,103,107]:
            rng = Random(p)
            for attempt in range(10000):
                u = rng.sample(range(p),10)
                f,P,b = star_data(p,u)
                if len(b)==1:
                    result = reconstruct(p,u,f,P,b[0])
                    if result[0]=='inconsistent-scalars':
                        test(p,u)
                        break
            else:
                raise ValueError('failed to find deterministic witness')
