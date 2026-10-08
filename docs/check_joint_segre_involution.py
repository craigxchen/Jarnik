"""Exact certificate for the Igusa deck involution pulled to the Segre cubic."""

from fractions import Fraction
from itertools import permutations

from check_near_balanced_invariant_congruences import ZERO, add, mul, scale, variable


ONE = {ZERO: 1}


def sub(a, b):
    return add(a, scale(b, -1))


def product(*polynomials):
    result = ONE
    for polynomial in polynomials:
        result = mul(result, polynomial)
    return result


def total(polynomials):
    result = {}
    for polynomial in polynomials:
        result = add(result, polynomial)
    return result


def forms(m):
    a, b, c, d, e = m
    f = sub(c, b)
    h = add(b, c)
    u = add(h, scale(d, 2))
    v = add(h, scale(a, 2))
    g = add(add(v, scale(d, 2)), scale(e, 2))
    return f, g, h, u, v


def cubic(m):
    a, b, c, d, e = m
    return sub(product(b, c, e), product(a, d, total(m)))


def gradient(m):
    a, b, c, d, e = m
    l = total(m)
    return [scale(mul(d, add(l, a)), -1),
            sub(mul(c, e), mul(a, d)),
            sub(mul(b, e), mul(a, d)),
            scale(mul(a, add(l, d)), -1),
            sub(mul(b, c), mul(a, d))]


def transform(m):
    a, b, c, d, e = m
    f, g, h, u, v = forms(m)
    return [scale(product(a, u, f, g), -1),
            scale(product(b, u, v, g), -1),
            product(c, u, v, g),
            scale(product(d, v, f, g), -1),
            scale(product(e, f, f, f), -1)]


def remainder(polynomial, divisor):
    """One-divisor polynomial division in lexicographic order."""
    polynomial = dict(polynomial)
    result = {}
    leading_divisor = max(divisor)
    coefficient_divisor = divisor[leading_divisor]
    assert coefficient_divisor in (1, -1)
    while polynomial:
        leading = max(polynomial)
        coefficient = polynomial[leading]
        if all(a >= b for a, b in zip(leading, leading_divisor)):
            shift = tuple(a - b for a, b in zip(leading, leading_divisor))
            factor = coefficient // coefficient_divisor
            polynomial = sub(polynomial, mul({shift: factor}, divisor))
        else:
            result[leading] = coefficient
            del polynomial[leading]
    return result


def chart(a, b, c):
    f = c * (a + 1 - b) - a
    g = c * (a + b - 1) - a
    h = c * (b + 1 - a) + a - 2*b
    u = c * (b + 1 - a) - a
    v = c * (a + b - 1) - a * (2*b - 1)
    return (a*h/v, b*f/v, c*h/f), (f, g, h, u, v)


def normalized_coordinate(point, infinity, zero, one):
    def vector(x):
        return (Fraction(1), Fraction(0)) if x is None else (Fraction(x), Fraction(1))
    def determinant(x, y):
        return x[0]*y[1] - x[1]*y[0]
    p, inf, z, o = map(vector, (point, infinity, zero, one))
    numerator = determinant(p, z) * determinant(o, inf)
    denominator = determinant(p, inf) * determinant(o, z)
    return None if not denominator else numerator / denominator


def inverse_normalization(t, first, second, third):
    if t is None:
        return first
    numerator = t*(third-second)*first-(third-first)*second
    denominator = t*(third-second)-(third-first)
    return None if denominator == 0 else numerator/denominator


def main():
    m = [variable(i) for i in range(5)]
    phi = cubic(m)
    tm = transform(m)
    f, g, h, u, v = forms(m)
    assert not remainder(cubic(tm), phi)
    transformed_forms = forms(tm)
    expected_forms = [product(h,u,v,g), product(f,h,u,v), product(f,u,v,g),
                      product(f,h,v,g), product(f,h,u,g)]
    for actual, expected in zip(transformed_forms, expected_forms):
        assert not remainder(sub(actual, expected), phi)
    grad = gradient(m)
    deck = [product(grad[i], f, f) for i in range(4)] + [product(grad[4], g, g)]
    common = scale(product(h,u,v,g), -1)
    for actual, desired in zip(gradient(tm), deck):
        assert not remainder(sub(actual, mul(common, desired)), phi)
    p, gu, gv, q, gw = grad
    assert not remainder(sub(product(sub(gu,gv), g, gw), product(sub(mul(p,q),mul(gu,gv)), f)), phi)

    a, b, c = [variable(i) for i in range(3)]
    cf = sub(product(c, sub(add(a,ONE),b)), a)
    cg = sub(product(c, sub(add(a,b),ONE)), a)
    ch = sub(add(product(c,sub(add(b,ONE),a)),a),scale(b,2))
    cu = sub(product(c,sub(add(b,ONE),a)),a)
    cv = sub(product(c,sub(add(a,b),ONE)),product(a,sub(scale(b,2),ONE)))
    # Cleared five-factor transformation identities.
    assert sub(product(c,ch,add(sub(product(a,ch),product(b,cf)),cv)), product(a,ch,cf)) == product(ch,cu,cg)
    assert sub(product(c,ch,sub(add(product(a,ch),product(b,cf)),cv)), product(a,ch,cf)) == product(cf,ch,cu)
    assert add(product(c,ch,add(sub(product(b,cf),product(a,ch)),cv)), product(cf,sub(product(a,ch),scale(product(b,cf),2)))) == product(cf,cu,cg)
    assert sub(product(c,ch,add(sub(product(b,cf),product(a,ch)),cv)), product(a,ch,cf)) == product(cf,ch,cg)
    assert sub(product(c,ch,cv,sub(add(product(a,ch),product(b,cf)),cv)), product(a,ch,cf,sub(scale(product(b,cf),2),cv))) == product(cf,ch,cu,cg)
    # Three anchor differences and three moving-point differences.
    assert sub(product(a,ch),cv) == scale(product(sub(a,ONE),cf),-1)
    assert sub(product(b,cf),cv) == scale(product(sub(b,ONE),cu),-1)
    assert sub(product(c,ch),cf) == product(sub(c,ONE),cu)
    assert sub(product(a,ch),product(b,cf)) == product(sub(b,a),cg)
    assert sub(product(c,ch,cv),product(b,cf,cf)) == product(sub(c,b),cu,cg)
    assert sub(product(c,cv),product(a,cf)) == product(sub(c,a),cg)

    for sample in ((2,3,5),(2,4,7),(3,5,8),(-1,2,4),(2,3,4)):
        original = tuple(map(Fraction,sample))
        moved, fs = chart(*original)
        restored, _ = chart(*moved)
        assert restored == original
        x = original[0]*(original[2]-1)/(original[2]*(original[1]-1))
        xp = moved[0]*(moved[2]-1)/(moved[2]*(moved[1]-1))
        assert xp == -x
    original = (None,Fraction(0),Fraction(1),Fraction(2),Fraction(3),Fraction(5))
    moved, _ = chart(Fraction(2),Fraction(3),Fraction(5))
    assert moved == (Fraction(6,5),Fraction(-3,5),Fraction(-15))
    target = frozenset((None,Fraction(0),Fraction(1)) + moved)
    normalizations = 0
    for anchors in permutations(original,3):
        image = frozenset(normalized_coordinate(point,*anchors) for point in original)
        assert image != target
        normalizations += 1
    assert normalizations == 120
    positive_edges = ((3,1),(5,2),(0,4))
    negative_edges = ((0,3),(5,1),(4,2))
    scores = []
    for mask in range(64):
        def inside(edge):
            return all(mask >> vertex & 1 for vertex in edge)
        scores.append(sum(inside(edge) for edge in positive_edges)
                      - sum(inside(edge) for edge in negative_edges))
    assert sum(scores) == 0
    assert sum(max(0, score) for score in scores) == 12
    pole_example = [Fraction(0),Fraction(8,9),Fraction(10,11),
                    Fraction(40,43),Fraction(16,17),Fraction(1)]
    abc = tuple(normalized_coordinate(x,*pole_example[:3]) for x in pole_example[3:])
    assert abc == (Fraction(2),Fraction(5,2),Fraction(5))
    moved, fs = chart(*abc)
    assert fs[0] == Fraction(1,2)
    assert inverse_normalization(moved[2],*pole_example[:3]) is None
    reverse = pole_example[::-1]
    abc = tuple(normalized_coordinate(x,*reverse[:3]) for x in reverse[3:])
    assert abc == (Fraction(9,4),Fraction(3),Fraction(6))
    moved, fs = chart(*abc)
    assert fs[0] == Fraction(-3,4)
    lifted = [inverse_normalization(x,*reverse[:3]) for x in moved]
    assert lifted == [Fraction(250,269),Fraction(312,331),Fraction(160,161)]
    assert all(reverse[3] < x < reverse[0] for x in lifted)
    print("Passed: quartic Segre map, five transformed forms, deck-gradient composition,")
    print("six collision identities, generic involution, and X'=-X.")
    print("Example (2,3,5) -> (6/5,-3/5,-15); all120 original normalizations excluded.")
    print("Fixed-label pole example and reversed arc-preserving output passed.")
    print("No radius descent or preservation of the full endpoint profile is asserted.")


if __name__ == '__main__':
    main()
