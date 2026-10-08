"""Prime-power checks for the moving-degree direct Runge remainder."""

from fractions import Fraction as F
from math import gcd, isqrt, prod

from check_gaussian_reflection_replacement import gnorm
from check_mobius_reciprocal_stretch_grid import realize
from check_reciprocal_squareclass_runge import trunc_sqrt


def residue(value, modulus):
    value = F(value)
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def remainder(roots, Z):
    k = len(roots)//2
    S = trunc_sqrt(roots, k)
    q = 2**(2*k-1)
    value = q*q*(prod(Z-c for c in roots)
                 - sum(a*Z**j for j, a in enumerate(S))**2)
    assert value.denominator == 1
    return int(value)


def coefficient(units, k):
    if not units:
        return F(int(k == 0))
    return trunc_sqrt(units, k)[0]


def prime_after(bound):
    for p in range(bound+1, 10*bound+100):
        if p % 4 == 1 and all(p % d for d in range(2, isqrt(p)+1)):
            return p
    raise AssertionError('search bound')


def choose_units(p, k, count):
    squares = sorted({b*b % p for b in range(1, p)})
    fixed = squares[:count-1]
    for candidate in squares[count-1:]:
        values = fixed + [candidate]
        if residue(coefficient(values, k), p):
            roots = [next(b for b in range(1, p) if b*b % p == u)
                     for u in values]
            return roots, [b*b for b in roots]
    raise AssertionError('nonzero polynomial counting failed')


def sqrt_mod_power(value, p, precision):
    r = next(a for a in range(1, p) if (a*a-value) % p == 0)
    modulus = p
    for _ in range(1, precision):
        digit = ((value-r*r)//modulus)*pow(2*r, -1, p) % p
        r += digit*modulus
        modulus *= p
    assert (r*r-value) % modulus == 0
    return r


def general_congruences():
    cases = homogeneous = 0
    for k in range(1, 7):
        q = 2**(2*k-1)
        for p in (5, 13):
            for e in (1, 2, 3):
                modulus = p**e
                Z = 2*modulus
                for j in range(1, 2*k+1):
                    units = [u for u in range(1, 4*k*p) if u % p][:2*k-j]
                    roots = [modulus*t for t in range(1, j+1)] + units
                    value = remainder(roots, Z)
                    expected = -residue(q*coefficient(units, k), modulus)**2
                    assert (value-expected) % modulus == 0
                    cases += 1
                base = list(range(1, 2*k+1))
                assert remainder([modulus*c for c in base], modulus*(4*k+3)) \
                    == modulus**(2*k)*remainder(base, 4*k+3)
                homogeneous += 1
    return cases, homogeneous


def variable_degree_local_models():
    cases = 0
    for k in range(2, 11):
        p = prime_after(8*k)
        for j in range(2, 2*k, 2):
            bs, units = choose_units(p, k, 2*k-j)
            assert residue(coefficient(units, k), p) != 0
            for e in (1, 2, 3):
                precision = e+2
                modulus = p**precision
                Z = 2*p**e
                aa = [sqrt_mod_power(Z-b*b, p, precision) for b in bs]
                roots, labels = list(units), [1]*len(units)
                used = set()
                for t in range(1, p):
                    if (1+t*t) % p == 0:
                        continue
                    inv = pow(1+t*t, -1, modulus)
                    a, b = (1-t*t)*inv % modulus, 2*t*inv % modulus
                    if b*b % p in used:
                        continue
                    used.add(b*b % p)
                    assert (a*a+b*b-1) % modulus == 0
                    roots.append(Z*b*b % modulus)
                    labels.append(Z)
                    aa.append(a)
                    if len(used) == j:
                        break
                assert len(roots) == 2*k and len(set(roots)) == 2*k
                assert prod(labels) == Z**j
                assert all((Z-c-d*a*a) % modulus == 0
                           for c, d, a in zip(roots, labels, aa))
                assert (prod(Z-c for c in roots)
                        - (Z**(j//2)*prod(aa))**2) % modulus == 0
                assert remainder(roots, Z) % p != 0
                cases += 1
    return cases


def literal_fixture():
    rows = [(1, 0), (3, 2), (4, 7), (10, 11), (9, 32)]
    points, N = realize(rows)
    assert N == 1105
    assert points == [(32, 9), (4, 33), (-24, 23), (-12, 31), (-32, 9)]
    assert all(gcd(*point) == 1 for point in points)
    anchor = points[0]
    roots = [N-anchor[0]*z[0]-anchor[1]*z[1] for z in points[1:]]
    labels = [2*N//gnorm(h) for h in rows[1:]]
    assert roots == [680, 1666, 1210, 2048]
    assert labels == [170, 34, 10, 2] and prod(labels) == 340**2
    assert trunc_sqrt(roots, 2) == [1701512, -2802, 1]
    assert prod(2*N-c for c in roots) == 367200**2
    assert remainder(roots, 2*N) == -1264902967296
    assert remainder(roots, 2*N) % 5 == 4


def main():
    congruences, homogeneous = general_congruences()
    models = variable_degree_local_models()
    literal_fixture()
    print(f'PASS: {congruences} exact prime-power remainder congruences; '
          f'{homogeneous} full-occupancy scaling identities; {models} '
          'variable-degree local norm-square models; literal binary fixture.')


if __name__ == '__main__':
    main()
