"""Exact tests of small integer Plücker coefficients on fixed Paley parents."""
from collections import Counter
from itertools import combinations, product
from math import gcd, isqrt, lcm, prod
from random import Random

from ii_paley import double_core


def split_primes(lower, count):
    primes = []
    candidate = lower + (1-lower) % 4
    while len(primes) < count:
        if all(candidate % p for p in range(2, isqrt(candidate) + 1)):
            primes.append(candidate)
        candidate += 4
    return primes


def norm2(v):
    return sum(x*x for x in v)


def modular_minimum(a, b, c):
    """Shortest squared norm for a*x + c*z = 0 (mod b), in Z^2."""
    r = (-c * pow(a, -1, b)) % b
    v, w = (b, 0), (r, 1)
    while True:
        if norm2(v) > norm2(w):
            v, w = w, v
        vv = norm2(v)
        vw = sum(x*y for x,y in zip(v,w))
        if 2*abs(vw) <= vv:
            assert abs(v[0]*w[1]-v[1]*w[0]) == b
            assert (a*v[0]+c*v[1]) % b == (a*w[0]+c*w[1]) % b == 0
            return vv
        multiple = (2*vw+vv)//(2*vv)
        w = tuple(y-multiple*x for x,y in zip(v,w))


def quartet_types(columns, indices):
    groups = [[], [], []]
    for index, column in enumerate(columns):
        bits = [column[i] for i in indices]
        if sum(bits) == 2:
            partner = next(j for j in range(1,4) if bits[j] == bits[0])
            groups[partner-1].append(index)
    return groups


def run(q, seed, lower=10**5):
    columns = [col[1:] for family in double_core(q) for col in family]
    primes = split_primes(lower, len(columns))
    Random(seed).shuffle(primes)
    parent = prod(primes)
    largest = max(prod(p for p,col in zip(primes,columns) if col[i] != col[j])
                  for i,j in combinations(range(q),2))
    denom = 16*isqrt(parent)
    upper = (largest+denom-1)//denom
    counts, survivors = Counter(), []
    bound = 2*upper**2
    for marks in combinations(range(q),4):
        groups = quartet_types(columns, marks)
        sizes = list(map(len,groups))
        assert sizes[0] == sizes[1] == sizes[2] >= 1
        counts[sizes[0]] += 1
        a,b,c = [prod(primes[j] for j in group) for group in groups]
        max_minimum = max(modular_minimum(a,b,c), modular_minimum(b,c,a),
                          modular_minimum(c,a,b))
        if max_minimum <= bound:
            survivors.append(marks)
    print('M ordinary marks, fixed permutation seed, parent digits:', q,seed,len(str(parent)))
    print('min and max prime, relaxed bound on products of two cofactors:',min(primes),max(primes),upper)
    print('two-two cut count histogram:',sorted(counts.items()))
    print('quartets not excluded by the three Euclidean lower bounds:',len(survivors),'of',sum(counts.values()))
    return len(survivors)


def positive_control():
    from iii_tri_balanced import RESIDUE, evaluate, words
    ds, ts = [], []
    for re_poly, imaginary in words():
        real = evaluate(re_poly, RESIDUE)
        den = gcd(real, imaginary)
        t = imaginary//den
        assert real > 0
        ds.append((real//den)**2+t*t)
        ts.append(t)
    pairings = ((0,5),(1,4),(2,3))
    pair_norms = [ds[i]*ds[j] for i,j in pairings]
    L = lcm(*pair_norms)
    A = [isqrt(L//n) for n in pair_norms]
    assert all(z*z*n==L for z,n in zip(A,pair_norms))
    x,y,z = (ts[i]*ts[j] for i,j in pairings)
    assert x*A[0]-y*A[1]+z*A[2] == 0
    assert modular_minimum(*A) <= x*x+z*z
    print("Genuine four-cycle: primitive matching products times cofactor products pass signed Pluecker.")


def bridge_control():
    def kind(v):
        if sum(v) != 2:
            return '0'
        return 'ABC'[next(j for j in range(1,4) if v[j]==v[0])-1]
    counts = Counter((kind(v[:4]),kind(v[:3]+v[4:])) for v in product((0,1),repeat=5))
    assert counts == {('0','0'):14, **{k:2 for s in 'ABC' for k in [('0',s),(s,'0'),(s,s)]}}
    print('Five-mark bridge: common modulus mass 2/32; reduced numerator and denominator mass 4/32 each.')


if __name__ == '__main__':
    positive_control()
    bridge_control()
    for q in (7,11,19):
        assert run(q,2026)==0
