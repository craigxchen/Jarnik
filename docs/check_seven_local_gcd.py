"""Exact seven-point square-star and triple-content audit; no height claim."""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, isqrt, lcm, prod
from random import Random

from check_low_moment_plucker_independent import det2, primitive, primes_of, valuation


def check(form, vectors):
    aa, bb, cc = form
    determinant = aa * cc - bb * bb
    assert isqrt(determinant) ** 2 == determinant > 0
    qs = [aa*x*x + 2*bb*x*y + cc*y*y for x, y in vectors]
    stars = [prod(det2(v, w) for j, w in enumerate(vectors) if i != j)
             for i, v in enumerate(vectors)]
    assert all(stars)
    gale = [[F(qs[i] ** 2, stars[i]) * z for z in v]
            for i, v in enumerate(vectors)]
    edges = list(combinations(range(7), 2))
    pp = dict(zip(edges, primitive([det2(gale[i], gale[j]) for i, j in edges])))
    def bracket(i, j):
        return pp[i, j] if i < j else -pp[j, i]
    ds = [prod(bracket(i, j) for j in range(7) if i != j) for i in range(7)]
    assert all(d < 0 and isqrt(-d) ** 2 == -d for d in ds)
    hs = [isqrt(-d) for d in ds]
    ss = reduce(gcd, hs)
    cs = [h // ss for h in hs]
    assert reduce(gcd, cs) == 1
    assert prod(hs) == prod(abs(p) for p in pp.values())
    triples, roots = [], []
    for i, j, k in combinations(range(7), 3):
        x, y, z = cs[i]*pp[j,k]**2, cs[j]*pp[i,k]**2, cs[k]*pp[i,j]**2
        value = 2*(x*y+x*z+y*z)-x*x-y*y-z*z
        assert value > 0 and isqrt(value) ** 2 == value
        triples.append(abs(pp[i,j]*pp[i,k]*pp[j,k]))
        roots.append(isqrt(value))
    gg, rr = reduce(gcd, triples), reduce(gcd, roots)
    assert all(F(root, triple) == F(rr, gg) for root, triple in zip(roots, triples))
    original = [F(4*determinant*det2(vectors[i],vectors[j])**2, qs[i]*qs[j]) for i,j in edges]
    recovered = [F(rr**2*pp[i,j]**2, gg**2*cs[i]*cs[j]) for i,j in edges]
    assert original == recovered
    radius = lcm(*(chord.denominator for chord in original))
    formula = lcm(*(gg**2*cs[i]*cs[j] // gcd(gg**2*cs[i]*cs[j],rr**2*pp[i,j]**2) for i,j in edges))
    assert radius == formula
    assert all(p % 4 == 1 for p in primes_of(F(rr, gg).denominator))
    primes = {2}
    for value in list(pp.values()) + cs + [rr, gg, radius]:
        primes |= primes_of(value)
    for prime in primes:
        es = {edge: valuation(value,prime) for edge,value in pp.items()}
        av = [valuation(c,prime) for c in cs]
        gv, rv = valuation(gg,prime), valuation(rr,prime)
        assert gv == min(es[i,j]+es[i,k]+es[j,k] for i,j,k in combinations(range(7),3))
        assert valuation(radius,prime) == max(0,max(av[i]+av[j]-2*es[i,j] for i,j in edges)+2*gv-2*rv)
        assert all(sum(es[min(i,j),max(i,j)] for j in range(7) if i != j) % 2 == 0 for i in range(7))
    return gg, len(primes)


def main():
    rng = Random(719)
    node_sets = [[(i,1) for i in range(7)], [(1,0)]+[(i,1) for i in range(6)]]
    node_sets += [[(i,1) for i in rng.sample(range(-12,13),7)] for _ in range(12)]
    count = prime_count = 0
    for form in [(1,0,1),(1,1,2),(2,1,5)]:
        for nodes in node_sets:
            gg, pc = check(form,nodes)
            count += 1
            prime_count += pc
    assert check((1,0,1),[(i,1) for i in range(7)])[0] == 10
    print(f'Passed {count} exact seven-point examples and {prime_count} all-prime checks.')


if __name__ == '__main__':
    main()
