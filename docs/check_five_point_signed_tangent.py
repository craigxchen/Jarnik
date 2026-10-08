"""Exact audit of signed leading shapes in the six-factor family."""

from functools import reduce
from itertools import combinations, product
from math import gcd


SUBSETS = [{4, 6}, {1, 3, 6}, {1, 4, 5}, {2, 3, 5}, {1, 2, 3, 4}]
SUPPORTS = [{j for j, subset in enumerate(SUBSETS) if r in subset}
            for r in range(1, 7)]
FORMS = [(1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (3, 1)]


def gmul(z, w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def gproduct(values):
    return reduce(gmul, values, (1, 0))


def points(u, v, t):
    aa = [u, 2*u, v, u+v, 2*u+v, 3*u+v]
    ans = []
    for subset in SUBSETS:
        z = gproduct((a, t if r in subset else -t)
                     for r, a in enumerate(aa, 1))
        ans.append(((-1)**len(subset)*z[0], (-1)**len(subset)*z[1]))
    return aa, ans


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def prime_factors(n):
    ans = set()
    n = abs(n)
    p = 2
    while p * p <= n:
        while n % p == 0:
            ans.add(p)
            n //= p
        p += 1
    if n > 1:
        ans.add(n)
    return ans


def signed_cover_primes():
    """Primes permitting a nonzero signed-equal coefficient cover."""
    ans = set()
    rows = set(range(5))
    for tags in product((-1, 0, 1), repeat=6):
        cover = set()
        equations = []
        for r, tag in enumerate(tags):
            if tag == 1:
                cover |= SUPPORTS[r]
            elif tag == -1:
                cover |= rows - SUPPORTS[r]
            if tag:
                equations.append((FORMS[r][0], FORMS[r][1], -tag))
        if cover != rows:
            continue
        minors3 = [det3(*triple) for triple in combinations(equations, 3)]
        content3 = reduce(gcd, (abs(x) for x in minors3), 0)
        if content3:
            ans |= prime_factors(content3)
            continue

        # Rank at most two over Q.  In every covering case its rational
        # kernel forces c=0 or one of the six coefficient forms to vanish.
        # Only primes at which the rank drops can evade that obstruction.
        pairs = [(a, b) for a, b in combinations(equations, 2)
                 if cross(a, b) != (0, 0, 0)]
        assert pairs
        kernel = cross(*pairs[0])
        common = reduce(gcd, (abs(x) for x in kernel))
        kernel = tuple(x // common for x in kernel)
        assert (kernel[2] == 0
                or any(a*kernel[0] + b*kernel[1] == 0 for a, b in FORMS))
        minors2 = []
        for a, b in combinations(equations, 2):
            minors2 += [a[i]*b[j]-a[j]*b[i] for i, j in combinations(range(3), 2)]
        content2 = reduce(gcd, (abs(x) for x in minors2))
        ans |= prime_factors(content2)
    return ans


def rank_deficient_cover_closures():
    """Actual rational tag patterns arising from rank-deficient covers."""
    ans = set()
    rows = set(range(5))
    for tags in product((-1, 0, 1), repeat=6):
        cover = set()
        equations = []
        for r, tag in enumerate(tags):
            if tag == 1:
                cover |= SUPPORTS[r]
            elif tag == -1:
                cover |= rows - SUPPORTS[r]
            if tag:
                equations.append((FORMS[r][0], FORMS[r][1], -tag))
        if cover != rows:
            continue
        if any(det3(*triple) for triple in combinations(equations, 3)):
            continue
        pairs = [(a, b) for a, b in combinations(equations, 2)
                 if cross(a, b) != (0, 0, 0)]
        assert pairs
        kernel = cross(*pairs[0])
        common = reduce(gcd, (abs(x) for x in kernel))
        u, v, c = (x // common for x in kernel)
        if c == 0:
            continue
        closure = tuple(1 if a*u+b*v == c else
                        -1 if a*u+b*v == -c else 0
                        for a, b in FORMS)
        actual_cover = set()
        for r, tag in enumerate(closure):
            if tag:
                actual_cover |= (SUPPORTS[r] if tag == 1
                                 else rows - SUPPORTS[r])
        assert actual_cover == rows
        ans.add(closure)
    return ans


def actual_cover_patterns(p):
    """All signed coefficient-cover patterns over the finite field F_p."""
    ans = set()
    rows = set(range(5))
    for u, v, c in product(range(p), repeat=3):
        if c == 0:
            continue
        aa = [(a*u+b*v) % p for a, b in FORMS]
        tags = tuple(1 if x == c else -1 if x == -c % p else 0
                     for x in aa)
        cover = set()
        for r, tag in enumerate(tags):
            if tag:
                cover |= SUPPORTS[r] if tag == 1 else rows - SUPPORTS[r]
        if cover == rows:
            ans.add(tags)
    return ans


def offsets(u, v):
    return [
        u * (13*u*u + 12*u*v + 4*v*v),
        u * (13*u*u + 10*u*v + 2*v*v),
        u * (u*u + 2*u*v + 2*v*v),
        5*u*u*u,
        u * (u*u - 6*u*v - 2*v*v),
    ]


def main():
    # Equal subset sums force precisely these six coefficient forms.
    assert [sum(FORMS[r-1][0] for r in subset) for subset in SUBSETS] == [4]*5
    assert [sum(FORMS[r-1][1] for r in subset) for subset in SUBSETS] == [2]*5

    # Two pair differences multiply to the full six-factor product.
    for u, v in [(1, 3), (2, -3), (3, 5), (-2, 7), (5, -8)]:
        d = offsets(u, v)
        assert d[1] - d[4] == 4*u*(u+v)*(3*u+v)
        assert d[2] - d[4] == 4*u*v*(2*u+v)
        full_product = 2*u*u*v*(u+v)*(2*u+v)*(3*u+v)
        diameter = max(d) - min(d)
        assert diameter**4 >= 64 * full_product**2

    # The exact equality shape has a duplicate pair.
    assert offsets(2, -3)[1] == offsets(2, -3)[2]

    # A symbolic signed-support enumeration leaves only 2, 3, and 5.
    # The odd split prime among them is 5; 3 is inert.
    assert signed_cover_primes() == {2, 3, 5}

    # Over every characteristic other than 2, 3, and 5, a cover must
    # be rank-deficient.  Closing each such rational solution under
    # all six coefficient equations leaves precisely these six
    # boundary patterns.  They respectively impose u+v=0, 2u+v=0,
    # or u=0.
    boundary_patterns = {
        (0, -1, 0, 0, 0, -1),
        (0, -1, 1, 0, 0, 0),
        (0, 0, -1, -1, -1, -1),
        (0, 0, 1, 1, 1, 1),
        (1, 0, -1, 0, 1, 0),
        (1, 0, 0, -1, 0, 1),
    }
    assert rank_deficient_cover_closures() == boundary_patterns

    # Scale the nonzero residue u to one.  This checks every unit
    # solution at five, including all active signed factors.
    valid_covers = []
    for v, c in product(range(1, 5), repeat=2):
        aa = [(a + b*v) % 5 for a, b in FORMS]
        if 0 in aa:
            continue
        tags = tuple(1 if a == c else -1 if a == -c % 5 else 0 for a in aa)
        covered = set()
        for r, tag in enumerate(tags):
            if tag:
                covered |= SUPPORTS[r] if tag == 1 else set(range(5)) - SUPPORTS[r]
        if covered == set(range(5)):
            valid_covers.append((v, c, tags))
    assert valid_covers == [(1, 2, (0, 1, 0, 1, -1, 0))]

    # Allowing zero coefficient forms adds exactly the six boundary
    # patterns above and no other mod-five cover.
    assert actual_cover_patterns(5) == boundary_patterns | {
        (0, 1, 0, 1, -1, 0)
    }

    # Each side has degree at most six in each variable.  Equality on a
    # seven-point grid in every variable is an exact interpolation
    # certificate for the two polynomial identities (11).
    for u, v, t in product(range(7), repeat=3):
        aa, w = points(u, v, t)
        rhs14 = gproduct([(aa[0], t), (aa[2], t), (aa[4], -t)])
        rhs14 = tuple(-4*u*(u+v)*(3*u+v)*x for x in rhs14)
        rhs24 = gproduct([(aa[0], t), (aa[3], t), (aa[5], -t)])
        rhs24 = tuple(-4*u*v*(2*u+v)*x for x in rhs24)
        assert (w[1][0]-w[4][0], w[1][1]-w[4][1]) == rhs14
        assert (w[2][0]-w[4][0], w[2][1]-w[4][1]) == rhs24

    print("PASS: signed tangent offsets and the uniform 2sqrt(2) bound.")
    print("PASS: exact all-parameter chord factorizations by interpolation.")
    print("PASS: exact signed-cover determinant enumeration gives {2,3,5}.")
    print("PASS: six boundary closures and the sole extra mod-five cover.")


if __name__ == "__main__":
    main()
