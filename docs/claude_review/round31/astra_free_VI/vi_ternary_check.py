"""Exact fixed-field F_3 norm construction in Part VI (python-flint 0.9).

The auxiliary root labels are pairs in F_27. Polynomial arithmetic for the
actual marks and their common norm uses only F_9[t]. No random choices.
"""
from itertools import combinations, product

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx


def chi(a, q):
    if a == 0:
        return 0
    return 1 if a ** ((q - 1) // 2) == 1 else -1


def orientation(forms):
    base = fmpz_mod_poly_ctx(3)
    # x^3 - x + 1 has no root in F_3.
    field = fq_default_ctx(modulus=base([1, -1, 0, 1]), var="a")
    a = field.gen()
    elements = [b0 + b1*a + b2*a**2
                for b0, b1, b2 in product(range(3), repeat=3)]
    labels = []
    for x, y in product(elements, repeat=2):
        values = [u*x + v*y for u, v in forms]
        signs = [chi(z, 27) for z in values]
        assert [chi(z**9, 27) for z in values] == signs
        assert [chi(-z, 27) for z in values] == [-s for s in signs]
        labels.append(signs)
    possibilities = []
    for flip in product((1, -1), repeat=len(forms)):
        common = sum(min(1 + f*s for f, s in zip(flip, row))
                     for row in labels)
        possibilities.append((common, flip))
    return max(possibilities)


def raw_residues(flip, q, L):
    fp3 = fmpz_mod_poly_ctx(3)
    modulus = [1, -1, 0, 0, 0, 1]
    assert fp3(modulus).is_irreducible()
    f9 = fq_default_ctx(modulus=fp3([1, 0, 1]), var="i")
    R = fq_default_poly_ctx(f9)
    t = R.gen()
    pi = R(modulus)
    assert pi.is_irreducible()
    X = (t**L) % pi
    forms = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1)]
    monomials = [R(1), t, t**2]
    residues = []
    for form, sign in zip(forms, flip):
        pivot = form.index(1)
        probe = f9.gen()*sign*monomials[pivot]
        basis = [f9.gen()*(u-form[j]*monomials[pivot])
                 for j, u in enumerate(monomials) if j != pivot]
        slope = (basis[0]**(q-1)) % pi

        def line(w):
            return (w**q - slope*w) % pi

        slope2 = (line(basis[1])**(q-1)) % pi

        def plane(w):
            u = line(w)
            return (u**q - slope2*u) % pi

        h, c = plane(X), plane(probe)
        raw = (h*(h**((q-1)//2)-c**((q-1)//2))**2) % pi
        residues.append(raw)
    assert all(x != y for x, y in combinations(residues, 2))
    return residues


def five_mark_factor_check():
    # Five F_3-linear forms in three variables over F_27. Record each exponent
    # vector with its exact multiplicity; this is a compact splitting-field
    # factorization certificate, distinct from the expanded four-mark test.
    from collections import Counter
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(3)([1, -1, 0, 1]), var="a")
    a = field.gen()
    els = [b0 + b1*a + b2*a**2 for b0, b1, b2 in product(range(3), repeat=3)]
    add = [[els.index(x+y) for y in els] for x in els]
    square_frob = [els.index(x**9) for x in els]
    neg = [els.index(-x) for x in els]
    chars = [chi(x, 27) for x in els]
    assert all(chars[j] == chars[square_frob[j]] == -chars[neg[j]]
               for j in range(27))
    table = Counter()
    for x, y, z in product(range(27), repeat=3):
        values = (x, y, z, add[x][y], add[x][z])
        table[tuple(chars[j] for j in values)] += 1
    common, flip = max(
        (sum(freq*min(1+f*c for f, c in zip(flips, row))
             for row, freq in table.items()), flips)
        for flips in product((1, -1), repeat=5)
    )
    assert all(table[row] == table[tuple(-c for c in row)] for row in table)
    n, h = 27**3, 27**2
    assert common > h
    for j in range(5):
        assert sum(freq*(1+flip[j]*row[j]) for row, freq in table.items()) == n
        for k in range(j):
            assert any(flip[j]*row[j] != flip[k]*row[k] for row in table)
    # If j is the X-degree in a formal reduced point, its t-coefficient
    # has degree <= 2*(n-common-j), by expanding its linear factors.
    # Choose the least L giving a strict gap by this bound and verify
    # distinctness after evaluation in the finite field F_9[t]/pi.
    L = (2*(n-h))//(common-h) + 1
    residues = raw_residues(flip, 27, L)
    D = (n-common)*L
    gap_upper = ((n+h)//2 - common)*L + (n-h)
    assert 2*gap_upper < D
    print("m=5 formal certificate: labels, h, common degree:", n, h, common)
    print("orientations, evaluation degree:", flip, L)
    print("certified point degree, parent degree, difference upper bound:", D, 2*D, gap_upper)
    print("raw residues modulo t^5-t+1:", tuple(str(x) for x in residues))


def main():
    q = 27
    forms = ((1, 0), (0, 1), (1, 1), (1, -1))
    common, flip = orientation(forms)
    assert common > q
    # For a linear form (a,b), its kernel line has basis vector (-b,a).
    # The finite-field subspace polynomial of this line is X^q - w^(q-1)X.
    fp3 = fmpz_mod_poly_ctx(3)
    f9 = fq_default_ctx(modulus=fp3([1, 0, 1]), var="i")
    ring = fq_default_poly_ctx(f9)
    i = f9.gen()
    t = ring.gen()
    power = max(2, ((q - 1)*q + common - q - 1) // (common - q))
    X = t**power
    raw = []
    for f, (a, b) in zip(flip, forms):
        w = i * (-b + a*t)
        # v0 has oriented form value 1, since f^2 = 1.
        v0 = ring(i*f) if a else i*f*t
        slope = w**(q - 1)
        c = v0**q - slope*v0
        h = X**q - slope*X
        raw.append(h * (h**((q-1)//2) - c**((q-1)//2))**2)
    common_poly = raw[0]
    for point in raw[1:]:
        common_poly = common_poly.gcd(point)
    assert common_poly.degree() >= common*power
    points = [p.exact_division(common_poly) for p in raw]
    degrees = {point.degree() for point in points}
    assert len(degrees) == 1
    degree = degrees.pop()
    norms = [point * ring([c**3 for c in point.coeffs()]) for point in points]
    assert all(n == norms[0] for n in norms)
    assert norms[0].degree() == 2*degree
    assert all(c**3 == c for c in norms[0].coeffs())
    differences = [(x-y).degree() for x, y in combinations(points, 2)]
    assert min(differences) >= 0 and 2*max(differences) <= degree
    print("m, q, k, root labels, common formal degree:", len(points), q, 2, q*q, common)
    print("orientations:", flip)
    print("evaluation degree, gcd degree:", power, common_poly.degree())
    print("point degree, parent degree:", degree, norms[0].degree())
    print("pair-difference degrees:", differences)


if __name__ == "__main__":
    main()
    five_mark_factor_check()
