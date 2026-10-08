"""Independent exact audit of positive square-determinant 4|4 cut lifts.

Uses only integer and Fraction arithmetic. Checks arbitrary-depth sample
lifts, the actual Gale formula, boundary reduction, and primitive scaling.
"""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, lcm, prod


P = 1009


def valuation(value):
    value = F(value)
    if not value:
        return float('inf')
    answer = 0
    numerator, denominator = abs(value.numerator), value.denominator
    while numerator % P == 0:
        answer += 1
        numerator //= P
    while denominator % P == 0:
        answer -= 1
        denominator //= P
    return answer


def root_mod_prime_power(target, depth):
    root = next(a for a in range(1, P) if (a*a-target) % P == 0)
    modulus = P
    for _ in range(1, depth):
        error = ((root*root-target) // modulus) % P
        root += ((-error * pow(2*root, -1, P)) % P) * modulus
        modulus *= P
        assert (root*root-target) % modulus == 0
    assert root > 0 and gcd(root, P) == 1
    return root


def determinant(rows):
    a, b, c = rows
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


def check(depth):
    accuracy = depth + 2
    aa, bb, cc = 2, 3, 7
    epsilon = P ** depth
    target = -bb*bb + 4*epsilon*aa*cc
    ss = root_mod_prime_power(target, depth+accuracy)
    c_tilde = F(ss*ss+bb*bb, 4*epsilon*aa)
    assert valuation(c_tilde-cc) >= accuracy
    assert aa > 0 and c_tilde > 0
    metric_det = aa*epsilon*c_tilde-F(bb*bb,4)
    assert metric_det == F(ss,2)**2 > 0
    assert valuation(metric_det) == 0

    xs, ys = list(range(1,5)), list(range(5,9))
    nodes = [F(x) for x in xs] + [F(epsilon,y) for y in ys]
    assert len(set(nodes)) == 8
    qs = [aa*r*r+bb*r+epsilon*c_tilde for r in nodes]
    assert all(q > 0 for q in qs)
    assert [valuation(q) for q in qs] == [0]*4+[depth]*4
    stars = [prod(r-other for j,other in enumerate(nodes) if i != j)
             for i,r in enumerate(nodes)]
    actual = [tuple(q*q/star*z for z in (r,r*r,F(epsilon)))
              for r,q,star in zip(nodes,qs,stars)]
    assert all(sum(row[col] for row in actual) == 0 for col in range(3))
    assert all(valuation(value) >= 0 for row in actual for value in row)

    boundary = []
    for x in xs:
        scale = F((aa*x+bb)**2, x*prod(other-x for other in xs if other != x))
        boundary.append((scale,scale*x,F(0)))
    for y in ys:
        scale = (F(cc)+F(bb,y))**2 / (y*prod(xs)*prod(F(1,other)-F(1,y) for other in ys if other != y))
        boundary.append((scale,F(0),scale*y))
    assert all(sum(row[col] for row in boundary) == 0 for col in range(3))
    # Both groups have the same minus sign relative to the displayed formula.
    assert all(valuation(actual[i][col]+boundary[i][col]) >= min(depth,accuracy)
               for i in range(8) for col in range(3))

    triples = list(combinations(range(8),3))
    minors = [determinant([actual[i] for i in triple]) for triple in triples]
    assert min(valuation(minor) for minor in minors) == 0
    assert min(valuation(determinant([actual[i] for i in triple]))
               for triple in combinations(range(7),3)) == 0
    common_denominator = lcm(*(minor.denominator for minor in minors))
    integer_minors = [int(common_denominator*minor) for minor in minors]
    content = reduce(gcd,integer_minors)
    primitive = [minor//content for minor in integer_minors]
    scalar = F(common_denominator,content)
    assert valuation(scalar) == 0
    assert reduce(gcd,primitive) == 1

    chord_denominators = []
    for i,j in combinations(range(8),2):
        chord = 4*metric_det*(nodes[i]-nodes[j])**2/(qs[i]*qs[j])
        assert 0 < chord <= 4
        assert -valuation(chord) == (depth if (i < 4) != (j < 4) else 0)
        chord_denominators.append(chord.denominator)
    conductor = lcm(*chord_denominators)
    assert valuation(conductor) == depth
    return len(triples)


def main():
    count = sum(check(depth) for depth in (1,2,4))
    print(f'PASS: positive rational lifts at depths 1, 2, 4; {count} minors, boundary convergence, exact chords, and primitive unit scaling.')


if __name__ == '__main__':
    main()
