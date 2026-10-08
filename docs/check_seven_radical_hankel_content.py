#!/usr/bin/env python3
"""Exact simultaneous-form contents and middle-conic determinant checks."""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations, product
from math import comb, gcd, isqrt, prod

from check_gale_square_quartic_characterization import determinant, primitive_pluckers
from check_seven_positive_radical_height import rational_gcd, check_saturated_basis_scaling
from check_seven_positive_radical_height import multiply as polynomial_multiply


INDICES = range(7)
PAIRS = list(combinations(INDICES, 2))


def sqrt_rational(value):
    numerator, denominator = isqrt(value.numerator), isqrt(value.denominator)
    assert numerator*numerator == value.numerator and denominator*denominator == value.denominator
    return F(numerator, denominator)


def binomial(value, degree):
    return prod((value-i)/F(i+1) for i in range(degree))


def check_content_and_coalescence():
    expected = [0, F(-35,2), -63, -140, F(-497,2), -385, -539]
    constants = []
    for k in range(1, 8):
        local = []
        for s in (1, 2, 3):
            local.append(min(a*(a+F(7,2)-s-k)
                             for a in range(max(0,k-7+s), min(k,s)+1)))
        assert sum(comb(7,s)*local[s-1] for s in (1,2,3)) == expected[k-1]
        matrix = [[binomial(F(7,2)-k, (6-a-b)//2)
                   if 6-a-b >= 0 and (6-a-b)%2 == 0 else F(0)
                   for b in range(k)] for a in range(k)]
        constants.append(determinant(matrix))
        assert constants[-1]
    assert constants[2] == F(-1,512)
    assert [sum(comb(7,s)*min(F(0),F(7-2*s-2*r,2)) for s in (1,2,3))
            for r in range(4)] == [0,F(-35,2),-63,F(-245,2)]
    print("PASS: exact content shifts and nonzero coalescent constants:", constants)


def configuration(form, vectors, check_all_branches=False):
    a, b, c = map(F, form)
    qs = [a*x*x+2*b*x*y+c*y*y for x,y in vectors]
    brackets = {(i,j): F(vectors[i][0]*vectors[j][1]-vectors[j][0]*vectors[i][1]) for i,j in PAIRS}
    def bracket(i,j):
        return brackets[i,j] if i<j else -brackets[j,i]
    stars = [prod(bracket(i,j) for j in INDICES if i!=j) for i in INDICES]
    tau = rational_gcd([qs[i]**5/stars[i]**2 for i in INDICES])
    cs = [int(qs[i]**5/stars[i]**2/tau) for i in INDICES]
    roots = [sqrt_rational(F(value)) for value in cs]
    signs = [1 if value>0 else -1 for value in stars]
    radius = sqrt_rational(a*c-b*b)
    all_brackets = abs(prod(brackets.values()))
    parameters = [F(y,x) for x,y in vectors]
    affine_values = [qs[i]/vectors[i][0]**2 for i in INDICES]
    for r in range(4):
        terms = [[signs[i]*roots[i]*F(vectors[i][0])**(2*r-j)*F(vectors[i][1])**j/qs[i]**r
                  for i in INDICES] for j in range(2*r+1)]
        entry_content = rational_gcd([value*value for row in terms for value in row if value])
        primitive_contents = [gcd(abs(x),abs(y)) for x,y in vectors]
        expected_content = rational_gcd([F(cs[i])*primitive_contents[i]**(4*r)/qs[i]**(2*r)
                                          for i in INDICES])
        assert entry_content == expected_content
        if r == 3:
            moments = [sum(row) for row in terms]
            for i in INDICES:
                polynomial = [F(1)]
                denominator = F(1)
                for j in INDICES:
                    if j != i:
                        polynomial = polynomial_multiply(polynomial,[-parameters[j],F(1)])
                        denominator *= parameters[i]-parameters[j]
                inverse_value = affine_values[i]**3*sum(a*b for a,b in zip(polynomial,moments))/denominator
                assert inverse_value == signs[i]*roots[i]
    contents, determinants, normalized = [], [], []
    for k in range(1,8):
        features = [[F(x)**(k-1-j)*F(y)**j for j in range(k)] for x,y in vectors]
        matrix = [[sum(signs[i]*roots[i]*features[i][j]*features[i][ell]/qs[i]**(k-1)
                       for i in INDICES) for ell in range(k)] for j in range(k)]
        actual = determinant(matrix)
        terms = []
        for subset in combinations(INDICES,k):
            within = prod(brackets[i,j] for i,j in combinations(subset,2))
            terms.append(prod(signs[i]*roots[i] for i in subset)*within**2/prod(qs[i]**(k-1) for i in subset))
        assert actual == sum(terms)
        content = rational_gcd([term*term for term in terms])
        normalized_value = actual/sqrt_rational(content)
        assert normalized_value.denominator == 1
        determinants.append(actual)
        contents.append(content)
        normalized.append(normalized_value)
    assert all(normalized[k-1] == -normalized[6-k] for k in range(1,7))
    assert normalized[6] == -1
    middle_content = contents[2]
    base = 8*radius**3/(tau**3*all_brackets*middle_content)
    assert base.denominator == 1 and base > 0
    gcd_terms = []
    for subset in combinations(INDICES,3):
        complement = [i for i in INDICES if i not in subset]
        gcd_terms.append(prod(qs[i] for i in subset)
                         *prod(brackets[i,j]**2 for i,j in combinations(subset,2))
                         *prod(brackets[i,j]**2 for i,j in combinations(complement,2)))
    big_gcd = rational_gcd(gcd_terms)
    assert middle_content == big_gcd/(tau**3*all_brackets**2)
    assert base == 8*radius**3*all_brackets/big_gcd
    pell_product = F(1)
    for i,j in PAIRS:
        x,y = vectors[i]; z,t = vectors[j]
        bilinear = a*x*z+b*(x*t+y*z)+c*y*t
        root_product = sqrt_rational(qs[i]*qs[j])
        factor = (bilinear+root_product)/(radius*abs(brackets[i,j]))
        conjugate = (bilinear-root_product)/(radius*abs(brackets[i,j]))
        assert factor > 0 and factor*conjugate == -1
        pell_product *= factor
    assert normalized[2]**2 == base/pell_product
    # Independently reconstruct the saturated Gale basis and transported metric.
    gale = [[qs[i]**2/stars[i]*z for z in vectors[i]] for i in INDICES]
    _, pp = primitive_pluckers(gale)
    def p(i,j):
        return pp[i,j] if i<j else -pp[j,i]
    hs = [isqrt(-prod(p(i,j) for j in INDICES if i!=j)) for i in INDICES]
    s = reduce(gcd,hs)
    assert [h//s for h in hs] == cs
    _, fitted, _, _ = check_saturated_basis_scaling(pp,cs,s)
    aa,bb,cc = fitted
    fitted_radius = sqrt_rational(aa*cc-bb*bb)
    fitted_gcd = rational_gcd([F((p(i,j)*p(i,k)*p(j,k))**4,(cs[i]*cs[j]*cs[k])**3)
                               for i,j,k in combinations(INDICES,3)])
    assert base == 8*fitted_radius**3*s**5/(prod(cs)*fitted_gcd)
    if check_all_branches:
        branch_product = F(1)
        for flips in product((-1,1),repeat=7):
            matrix = [[sum(flips[i]*signs[i]*roots[i]
                           *F(vectors[i][0])**(4-j-ell)*F(vectors[i][1])**(j+ell)/qs[i]**2
                           for i in INDICES) for ell in range(3)] for j in range(3)]
            value = determinant(matrix)/sqrt_rational(middle_content)
            assert value and value.denominator == 1
            branch_product *= value
        assert branch_product == base**64
    return int(base)


def check_laurent_constant():
    us = list(map(F,range(1,8)))
    ts = [(u-1/u)/2 for u in us]
    ys = [(u+1/u)/2 for u in us]
    pprime = [prod(ts[i]-ts[j] for j in INDICES if i!=j) for i in INDICES]
    matrix = [[sum(ys[i]*ts[i]**(a+b)/pprime[i] for i in INDICES) for b in range(3)] for a in range(3)]
    evaluation = [[1,t,t*t,t**3,y,t*y,t*t*y] for t,y in zip(ts,ys)]
    vt = prod(ts[j]-ts[i] for i,j in PAIRS)
    vu = prod(us[j]-us[i] for i,j in PAIRS)
    assert determinant(evaluation) == F(1,512)*vu/prod(u**3 for u in us)
    assert determinant(matrix) == -determinant(evaluation)/vt
    assert determinant(matrix) == -2**12*prod(u**3 for u in us)/prod(1+us[i]*us[j] for i,j in PAIRS)


def valuation(value, prime):
    exponent=0
    while value%prime==0:
        exponent+=1; value//=prime
    return exponent


def check_unbounded_base():
    for prime in (11,13):
        for exponent in (1,2):
            m=prime**exponent
            ts=[m*i for i in INDICES]
            qs=[1+t*t for t in ts]
            d={(i,j):ts[j]-ts[i] for i,j in PAIRS}
            D=prod(d.values())
            values=[]
            for subset in combinations(INDICES,3):
                other=[i for i in INDICES if i not in subset]
                values.append(prod(qs[i] for i in subset)
                              *prod(d[i,j]**2 for i,j in combinations(subset,2))
                              *prod(d[i,j]**2 for i,j in combinations(other,2)))
            G=reduce(gcd,values)
            base=F(8*D,G)
            assert base.denominator==1 and valuation(base.numerator,prime)==3*exponent


if __name__ == "__main__":
    check_content_and_coalescence()
    check_laurent_constant()
    vectors=[(2*u,u*u-1) for u in range(1,8)]
    bases=[configuration((1,0,1),vectors,True)]
    changed=[(x+2*y,y) for x,y in vectors]
    bases.append(configuration((1,-2,5),changed))
    bases.append(configuration((F(2,3),F(-4,3),F(10,3)),changed))
    assert len(set(bases))==1
    check_unbounded_base()
    print("PASS: exact Cauchy--Binet, all integral contents, complement duality, and 128 branch identities.")
    print("PASS: simultaneous linear-form denominator contents and exact seven-form Lagrange inverse.")
    print("PASS: affine constant, pair-factor identity, saturated-basis formula and metric covariance; base",bases[0])
    print("PASS: unbounded integer base has valuation 3e at the selected split and inert primes.")
