"""Exact certificates for cubic norm descent and cyclic-field nonsplitting.

Run: python3 docs/check_cubic_norm_descent_mod17.py
Only the standard library is used. The finite residue loop has 16 parameters.
"""

from fractions import Fraction as F
from itertools import permutations


ZERO = (F(0), F(0))
ONE = (F(1), F(0))
I = (F(0), F(1))


def zadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def zneg(a):
    return (-a[0], -a[1])


def zmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def zpadd(a, b):
    out = [ZERO] * max(len(a), len(b))
    for i in range(len(out)):
        out[i] = zadd(a[i] if i < len(a) else ZERO, b[i] if i < len(b) else ZERO)
    return out


def zpmul(a, b):
    out = [ZERO] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = zadd(out[i + j], zmul(x, y))
    return out


def det3(matrix):
    out = [ZERO]
    for perm in permutations(range(3)):
        inversions = sum(perm[i] > perm[j] for i in range(3) for j in range(i + 1, 3))
        term = [ONE]
        for row, col in enumerate(perm):
            term = zpmul(term, matrix[row][col])
        if inversions % 2:
            term = [zneg(x) for x in term]
        out = zpadd(out, term)
    return out


def kmul(a, b, field):
    """Multiply 1,alpha,alpha^2 coordinates over Q(i)."""
    raw = [ZERO] * 5
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            raw[i + j] = zadd(raw[i + j], zmul(x, y))
    if field == "pure":  # alpha^3=2, alpha^4=2alpha
        return (zadd(raw[0], zmul((F(2), F(0)), raw[3])),
                zadd(raw[1], zmul((F(2), F(0)), raw[4])), raw[2])
    assert field == "cyclic"  # alpha^3=3alpha-1, alpha^4=3alpha^2-alpha
    return (zadd(raw[0], zneg(raw[3])),
            zadd(zadd(raw[1], zmul((F(3), F(0)), raw[3])), zneg(raw[4])),
            zadd(raw[2], zmul((F(3), F(0)), raw[4])))


def kadd(a, b):
    return tuple(zadd(x, y) for x, y in zip(a, b))


def kscale(a, scalar):
    return tuple(zmul(x, scalar) for x in a)


def kneginv(a):
    return tuple(zneg(x) for x in a)


def field_norm_polynomial(poly, field):
    """Determinant of multiplication by a K(i)[t] element."""
    matrix = [[[] for _ in range(3)] for _ in range(3)]
    basis = [(ONE, ZERO, ZERO), (ZERO, ONE, ZERO), (ZERO, ZERO, ONE)]
    for col, e in enumerate(basis):
        products = [kmul(coef, e, field) for coef in poly]
        for row in range(3):
            matrix[row][col] = [x[row] for x in products]
    return det3(matrix)


def exact_fixture(field, delta, r):
    """Check monicity, the fixed root, and the exact imaginary norm."""
    if field == "pure":
        w = (ZERO, ONE, (F(1, 2), F(0)))
        v = (ZERO, ONE, (F(-1, 2), F(0)))
    else:
        w = ((F(2), F(0)), ZERO, (F(-1), F(0)))
        v = ((F(-2), F(0)), (F(2), F(0)), ONE)
    e = (ONE, ZERO, ZERO)
    u = 2 * r / (1 + delta * r * r)
    s = (1 - delta * r * r) / (1 + delta * r * r)
    h = delta * u * u
    k = s / (4 * h)
    # a,b,c in ascending polynomial order.
    a = [k, F(-3, 2), 2 * k, F(-1, 2), k]
    b = [k * (1 - 2 * h), F(-1, 2), 2 * k * (1 - h), F(-1, 2), k]
    c = [-(h + 1) * u / (2 * h), s * u / (2 * h),
         -(h + 1) * u / (2 * h), s * u / (2 * h), F(0)]
    U = [kadd(kadd(kscale(e, (a[j], F(0))), kscale(w, (b[j], F(0)))),
              kscale(v, (c[j], F(0)))) for j in range(5)]
    U[0] = kadd(U[0], kscale(e, I))
    # Y=[k(1+w)]^-1 is found by exact 3-by-3 rational linear algebra.
    lead = U[4]
    assert lead == kscale(kadd(e, w), (k, F(0)))
    # In these explicit fields, invert by the adjugate polynomial identity.
    # Compute inverse as a degree-two element by Gaussian elimination.
    columns = [kmul(lead, basis, field) for basis in
               ((ONE, ZERO, ZERO), (ZERO, ONE, ZERO), (ZERO, ZERO, ONE))]
    mat = [[columns[col][row][0] for col in range(3)] + [F(row == 0)]
           for row in range(3)]
    for col in range(3):
        pivot = next(row for row in range(col, 3) if mat[row][col])
        mat[col], mat[pivot] = mat[pivot], mat[col]
        scale = mat[col][col]
        mat[col] = [x / scale for x in mat[col]]
        for row in range(3):
            if row != col:
                factor = mat[row][col]
                mat[row] = [x - factor * y for x, y in zip(mat[row], mat[col])]
    Y = tuple((mat[row][3], F(0)) for row in range(3))
    P = [kmul(Y, coef, field) for coef in U]
    assert P[4] == e
    assert all(all(component[1] == 0 for component in coef) for coef in P[1:])
    assert tuple(component[1] for component in P[0]) == tuple(component[0] for component in Y)
    at_i = (ZERO, ZERO, ZERO)
    power = ONE
    for coef in P:
        at_i = kadd(at_i, kscale(coef, power))
        power = zmul(power, I)
    assert at_i == (ZERO, ZERO, ZERO)
    norm = field_norm_polynomial(P, field)
    norm_y = field_norm_polynomial([Y], field)[0][0]
    assert norm[-1] == ONE
    assert all(coef[1] == 0 for coef in norm[1:])
    assert norm[0][1] == -4 * norm_y
    return norm_y, len(norm) - 1


def cubic_mod17_basis(w, v, r, i=4):
    """Unnormalized residual cubic at a chosen basis and local i."""
    p = 17
    delta = 3
    den = (1 + delta * r * r) % p
    s_num = (1 - delta * r * r) % p
    assert den and s_num and r
    u = 2 * r * pow(den, -1, p) % p
    s = s_num * pow(den, -1, p) % p
    h = delta * u * u % p
    k = s * pow(4 * h, -1, p) % p
    half = pow(2, -1, p)
    ih = pow(2 * h, -1, p)
    return (
        k * (1 + w) % p,
        (-(1 + w) * half + i * k * (1 + w) + v * s * u * ih) % p,
        (k * (1 + w * (1 - 2 * h)) - i * (1 + w) * half
         + v * u * ih * (-h - 1 + i * s)) % p,
        (-1 + i * (k * (1 + w * (1 - 2 * h)) - v * u * (h + 1) * ih)) % p,
    )


def cubic_mod17(alpha, r):
    """Unnormalized residual cubic at alpha and i=4, r in F_17^*."""
    return cubic_mod17_basis((2 - alpha * alpha) % 17,
                             (alpha * alpha + 2 * alpha - 2) % 17, r)


def discriminant_mod17(coef):
    a, b, c, d = coef
    return (b * b * c * c - 4 * a * c**3 - 4 * b**3 * d
            - 27 * a * a * d * d + 18 * a * b * c * d) % 17


def root_count_mod17(coef):
    a, b, c, d = coef
    return sum((a * t**3 + b * t * t + c * t + d) % 17 == 0
               for t in range(17))


def reduction_splits(coef):
    """Factor a degree-at-most-three reduction, retaining multiplicities."""
    poly = list(coef)
    while poly and poly[0] == 0:
        poly.pop(0)
    assert poly
    for root in range(17):
        while len(poly) > 1:
            quotient = [poly[0]]
            for a in poly[1:-1]:
                quotient.append((a + root * quotient[-1]) % 17)
            if (poly[-1] + root * quotient[-1]) % 17:
                break
            poly = quotient
    return len(poly) == 1


def rotated_unit_certificate():
    p = 17
    old = ((4, 10), (3, 6), (10, 1))  # alpha=7,13,14; i=4
    rotations = [(x, y) for x in range(p) for y in range(p)
                 if (x * x + 3 * y * y) % p == 1]
    assert len(rotations) == 18
    bad_rotations = []
    for x, y in rotations:
        bases = [((x * w + y * v) % p, (-3 * y * w + x * v) % p)
                 for w, v in old]
        bad = [j for j, (w, _) in enumerate(bases) if (1 + w) % p == 0]
        if not bad:
            for r in range(1, p):
                assert any(not reduction_splits(cubic_mod17_basis(w, v, r))
                           for w, v in bases)
            continue
        bad_rotations.append((x, y))
        assert sorted(bases) == [(2, 0), (16, 3), (16, 14)]
        regular = next(j for j in range(3) if j not in bad)
        candidates = [r for r in range(1, p)
                      if reduction_splits(cubic_mod17_basis(*bases[regular], r))]
        assert candidates == [2, 5, 12, 15]
        # A split primitive polynomial over Z_17 reduces to linear or
        # constant factors. These degree-two reductions are irreducible.
        bad_basis = {basis for j, basis in enumerate(bases) if j in bad}
        assert bad_basis == {(16, 3), (16, 14)}
        witnesses = {2: ((16, 3), 4, (0, 5, 0, 4)),
                     5: ((16, 14), 13, (0, 2, 16, 5)),
                     12: ((16, 3), 13, (0, 2, 16, 5)),
                     15: ((16, 14), 4, (0, 5, 0, 4))}
        for r, (basis, local_i, expected) in witnesses.items():
            reduced = cubic_mod17_basis(*basis, r, local_i)
            assert reduced == expected
            assert not reduction_splits(reduced)
    assert bad_rotations == [(2, 13), (5, 3), (10, 1)]


def local_split_197_example():
    """A clean six-place split reduction, not a global splitting claim."""
    prime = 197
    rotate_p, rotate_q = F(-971, 973), F(36, 973)
    assert rotate_p * rotate_p + 3 * rotate_q * rotate_q == 1
    def reduce_fraction(x):
        return x.numerator * pow(x.denominator, -1, prime) % prime
    rot_p, rot_q = reduce_fraction(rotate_p), reduce_fraction(rotate_q)
    assert (rot_p, rot_q) == (163, 194)
    r = 102
    den, snum = (1 + 3 * r * r) % prime, (1 - 3 * r * r) % prime
    assert (den, snum) == (87, 112)
    u = 2 * r * pow(den, -1, prime) % prime
    s = snum * pow(den, -1, prime) % prime
    h = 3 * u * u % prime
    k = s * pow(4 * h, -1, prime) % prime
    ih, half = pow(2 * h, -1, prime), pow(2, -1, prime)
    expected = {
        (34, 14): ((179, 54, 42, 105), (8, 72, 120)),
        (34, 183): ((179, 164, 142, 90), (70, 177, 178)),
        (169, 14): ((56, 24, 59, 65), (3, 126, 152)),
        (169, 183): ((56, 32, 63, 130), (158, 165, 183)),
        (191, 14): ((68, 124, 181, 129), (24, 73, 75)),
        (191, 183): ((68, 190, 119, 66), (10, 63, 127)),
    }
    bases = []
    for alpha in (34, 169, 191):
        assert (alpha**3 - 3 * alpha + 1) % prime == 0
        assert (3 * alpha * alpha - 3) % prime
        w, v = (2 - alpha * alpha) % prime, (alpha * alpha + 2 * alpha - 2) % prime
        w, v = (rot_p * w + rot_q * v) % prime, (-3 * rot_q * w + rot_p * v) % prime
        bases.append((w, v))
        assert (1 + w) % prime
        for local_i in (14, 183):
            assert (local_i * local_i + 1) % prime == 0
            a = k * (1 + w) % prime
            b = (-(1 + w) * half + local_i * k * (1 + w)
                 + v * s * u * ih) % prime
            c = (k * (1 + w * (1 - 2 * h)) - local_i * (1 + w) * half
                 + v * u * ih * (-h - 1 + local_i * s)) % prime
            d = (-1 + local_i * (k * (1 + w * (1 - 2 * h))
                                 - v * u * (h + 1) * ih)) % prime
            coeff, roots = expected[(alpha, local_i)]
            assert (a, b, c, d) == coeff
            assert tuple(t for t in range(prime)
                         if (a * t**3 + b * t * t + c * t + d) % prime == 0) == roots
    assert bases == [(110, 74), (179, 192), (105, 128)]


def splitting_certificate():
    # alpha=7 and 13 are simple roots of x^3-3x+1 modulo 17.
    assert all((a**3 - 3 * a + 1) % 17 == 0 for a in (7, 13))
    assert all((3 * a * a - 3) % 17 for a in (7, 13))
    assert 4 * 4 % 17 == 16
    exceptions = []
    for r in range(1, 17):
        coef = cubic_mod17(7, r)
        disc = discriminant_mod17(coef)
        if not (disc and root_count_mod17(coef) < 3):
            exceptions.append(r)
    assert exceptions == [2, 3, 16]
    for r in exceptions:
        coef = cubic_mod17(13, r)
        disc = discriminant_mod17(coef)
        assert disc and pow(disc, 8, 17) == 16
        assert root_count_mod17(coef) == 1
    # At alpha=7, i=4: w=4, v=10. The tangent quadratics for r->0
    # and r->infinity have these nonsquare discriminants. This alone
    # determines their squareclasses, but the following formal series
    # independently checks the leading coefficient of the FULL monic
    # cubic discriminant. x is r at zero and 1/r at infinity.
    w, v, delta, i = 4, 10, 3, 4
    monic_scale = pow((1 + w) ** 2, -1, 17)
    d0 = 32 * i * (v * v - 4 * delta * (1 + w)) * monic_scale % 17
    dinf = (-32 * i * (v * v - 4 * delta * (1 + w))
            * pow(delta * delta, -1, 17) * monic_scale) % 17
    assert (d0, dinf) == (11, 12)
    assert pow(d0, 8, 17) == pow(dinf, 8, 17) == 16

    zero = (0, 0, 0)
    one = (1, 0, 0)

    def add(x, y):
        return tuple((a + b) % 17 for a, b in zip(x, y))

    def scale(x, c):
        return tuple(c * a % 17 for a in x)

    def sub(x, y):
        return add(x, scale(y, -1))

    def mul(x, y):
        return tuple(sum(x[j] * y[k - j] for j in range(k + 1)) % 17
                     for k in range(3))

    def full_discriminant_series(at_infinity):
        # Exact expansions modulo x^3 of u, v0, h=delta*u^2,
        # A=2h/v0, B=2u/v0. Every omitted term is divisible by x^3.
        if not at_infinity:
            u = (0, 2, 0)
            s = (1, 0, -2 * delta)
            h = (0, 0, 4 * delta)
            A = (0, 0, 8 * delta)
            B = (0, 4, 0)
        else:
            di = pow(delta, -1, 17)
            u = (0, 2 * di, 0)
            s = (-1, 0, 2 * di)
            h = (0, 0, 4 * di)
            A = (0, 0, -8 * di)
            B = (0, -4 * di, 0)
        lead = 1 + w
        # Qhat=(U+i)/(k(t-i)) has leading coefficient lead.
        # Its coefficients below are divided by lead to make it monic.
        b = add(sub(scale(one, i * lead), scale(A, lead)), scale(u, 2 * v))
        c = add(sub(sub(scale(one, lead), scale(h, 2 * w)),
                    scale(A, i * lead)),
                scale(mul(B, add(sub(scale(s, i), one), scale(h, -1))), v))
        d = add(scale(A, -2),
                scale(sub(sub(scale(one, lead), scale(h, 2 * w)),
                          scale(mul(B, add(one, h)), v)), i))
        inv_lead = pow(lead, -1, 17)
        b, c, d = (scale(x, inv_lead) for x in (b, c, d))
        b2, c2, d2 = mul(b, b), mul(c, c), mul(d, d)
        disc = add(sub(sub(mul(b2, c2), scale(mul(c2, c), 4)),
                       scale(mul(mul(b2, b), d), 4)),
                   add(scale(d2, -27), scale(mul(mul(b, c), d), 18)))
        return disc

    assert full_discriminant_series(False) == (0, 0, 7)
    assert full_discriminant_series(True) == (0, 0, 3)


if __name__ == "__main__":
    for field, delta in (("pure", -1), ("cyclic", 3)):
        norm_y, degree = exact_fixture(field, delta, F(1, 2))
        print(field, "r=1/2: deg norm", degree, "Norm(Y)", norm_y)
    splitting_certificate()
    rotated_unit_certificate()
    local_split_197_example()
    print("cyclic cubic: mod-17 nonsplitting certified for every rational r!=0")
    print("all 18 rotated bases: mod-17 nonsplitting certified for 17-adic unit r")
    print("p=197: one rotated rational parameter has six split local reductions")
