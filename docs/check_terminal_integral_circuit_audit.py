"""Exact terminal circuit/reconstruction checks for m=6,12.

This is independent of the general Specht argument.  It expands products of
disjoint triangle determinants over Z, specializes every terminal cut, checks
the two differential identities and reconstruction, and measures coefficient
norm conversion to noncrossing matching bases on four and eight rows.
"""

from fractions import Fraction
from itertools import combinations
from math import prod


def add(a, b, scale=1):
    out = dict(a)
    for mon, c in b.items():
        out[mon] = out.get(mon, 0) + scale * c
        if not out[mon]:
            del out[mon]
    return out


def multiply(a, b):
    out = {}
    for u, cu in a.items():
        for v, cv in b.items():
            w = tuple(x + y for x, y in zip(u, v))
            out[w] = out.get(w, 0) + cu * cv
    return {u: c for u, c in out.items() if c}


def variable(nvars, i):
    e = [0] * nvars
    e[i] = 1
    return {tuple(e): 1}


def determinant(i, j, m):
    x = [variable(2*m, a) for a in range(m)]
    y = [variable(2*m, m+a) for a in range(m)]
    return add(multiply(x[i], y[j]), multiply(x[j], y[i]), -1)


def triangle(a, b, c, m):
    return multiply(multiply(determinant(a, b, m), determinant(b, c, m)),
                    determinant(c, a, m))


def graph_product(m, triangles):
    out = {(0,) * (2*m): 1}
    for tri in triangles:
        out = multiply(out, triangle(*tri, m))
    return out


def restrict_cut(qpoly, m, inside):
    """Set inside rows=(0,1), keep outside=(P_i,Y_i), then divide P_out."""
    inside = set(inside)
    out = {}
    for mon, c in qpoly.items():
        if any(mon[i] for i in inside):
            continue
        if any(mon[m+i] == 0 for i in inside):
            continue
        new = list(mon)
        for i in inside:
            new[m+i] = 0
        for i in range(m):
            if i not in inside:
                assert new[i] >= 1, (inside, mon)
                new[i] -= 1
        key = tuple(new)
        out[key] = out.get(key, 0) + c
    return {u: c for u, c in out.items() if c}


def multiply_variable(poly, i, m):
    v = variable(2*m, i)
    return multiply(poly, v)


def terminal_audit(m, k, triangles):
    q = graph_product(m, triangles)
    s = 2*k
    coeffs = {}
    for inside in combinations(range(m), s):
        r = restrict_cut(q, m, inside)
        if r:
            coeffs[frozenset(inside)] = r

    # For R(x)=sum r_S x_S, check D_P R=D_Y R=0 coefficientwise.
    for weights_offset in (0, m):
        derived = {}
        for inside, r in coeffs.items():
            for i in inside:
                support = inside - {i}
                derived.setdefault(support, {})
                term = multiply_variable(r, weights_offset+i, m)
                derived[support] = add(derived[support], term)
        assert all(not poly for poly in derived.values()), (m, k, weights_offset)

    # Reconstruct Q(P,Y)=sum_S r_S prod_(i outside S)P_i prod_(i inside S)Y_i^2.
    reconstructed = {}
    for inside, r in coeffs.items():
        term = dict(r)
        for i in range(m):
            if i in inside:
                term = multiply(term, multiply(variable(2*m, m+i),
                                              variable(2*m, m+i)))
            else:
                term = multiply(term, variable(2*m, i))
        reconstructed = add(reconstructed, term)
    assert reconstructed == q, (m, k)

    # Restriction divides each surviving original monomial by one P outside
    # and deletes its Y^2 inside.  For these terminal graph products this is
    # an injective reassignment of monomials, hence exact aggregate l1 norm.
    original_l1 = sum(abs(c) for c in q.values())
    circuit_l1 = sum(abs(c) for r in coeffs.values() for c in r.values())
    assert original_l1 == circuit_l1, (m, k, original_l1, circuit_l1)
    expected_nonzero_cuts = 3**(2*k)
    assert len(coeffs) == expected_nonzero_cuts
    assert original_l1 == 6**(2*k)
    assert all(sum(abs(c) for c in r.values()) == 2**(2*k)
               for r in coeffs.values())
    return original_l1, len(coeffs), circuit_l1


def all_matchings(vertices):
    if not vertices:
        yield ()
        return
    i = vertices[0]
    for j in vertices[1:]:
        rest = tuple(x for x in vertices if x not in (i, j))
        for tail in all_matchings(rest):
            yield ((i, j),) + tail


def noncrossing(match):
    edges = [tuple(sorted(e)) for e in match]
    return not any(a < c < b < d or c < a < d < b
                   for (a, b), (c, d) in combinations(edges, 2))


def matching_poly(match, n):
    x = [variable(2*n, i) for i in range(n)]
    y = [variable(2*n, n+i) for i in range(n)]
    out = {(0,) * (2*n): 1}
    for i, j in match:
        out = multiply(out, add(multiply(x[i], y[j]), multiply(x[j], y[i]), -1))
    return out


def matching_basis_bound(n):
    basis = [m for m in all_matchings(tuple(range(n))) if noncrossing(m)]
    allm = list(all_matchings(tuple(range(n))))
    polys = [matching_poly(m, n) for m in basis]
    mons = sorted(set().union(*(p.keys() for p in polys)))
    def matrix_rank(mat):
        a = [[Fraction(x) for x in line] for line in mat]
        rank = 0
        for col in range(len(basis)):
            pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
            if pivot is None:
                continue
            a[rank], a[pivot] = a[pivot], a[rank]
            z = a[rank][col]
            a[rank] = [x/z for x in a[rank]]
            for i in range(rank+1, len(a)):
                if a[i][col]:
                    z = a[i][col]
                    a[i] = [x-z*y for x, y in zip(a[i], a[rank])]
            rank += 1
        return rank

    # Select a unimodular set of monomial coordinates. Its inverse gives a
    # certified l1 conversion bound from raw polynomial coefficients.
    pivot_rows = []
    for mon in mons:
        candidate = pivot_rows + [[p.get(mon, 0) for p in polys]]
        if matrix_rank(candidate) > len(pivot_rows):
            pivot_rows = candidate
        if len(pivot_rows) == len(basis):
            break
    assert len(pivot_rows) == len(basis)
    d = len(basis)
    aug = [[Fraction(pivot_rows[i][j]) for j in range(d)]
           + [Fraction(int(i == j)) for j in range(d)] for i in range(d)]
    det = Fraction(1)
    for col in range(d):
        pivot = next(i for i in range(col, d) if aug[i][col])
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]
            det = -det
        value = aug[col][col]
        det *= value
        aug[col] = [x/value for x in aug[col]]
        for i in range(d):
            if i != col and aug[i][col]:
                z = aug[i][col]
                aug[i] = [x-z*y for x, y in zip(aug[i], aug[col])]
    assert abs(det) == 1
    inv = [line[d:] for line in aug]
    raw_to_basis_l1 = max(sum(abs(inv[i][j]) for i in range(d))
                          for j in range(d))
    # Each matching expands in this basis by exact rational elimination.
    max_l1, max_coeff = 0, 0
    for match in allm:
        target = matching_poly(match, n)
        A = [[Fraction(p.get(mon, 0)) for p in polys]
             + [Fraction(target.get(mon, 0))] for mon in mons]
        row = 0
        for col in range(len(basis)):
            pivot = next((i for i in range(row, len(A)) if A[i][col]), None)
            assert pivot is not None
            A[row], A[pivot] = A[pivot], A[row]
            z = A[row][col]
            A[row] = [x/z for x in A[row]]
            for i in range(len(A)):
                if i != row and A[i][col]:
                    z = A[i][col]
                    A[i] = [x-z*y for x, y in zip(A[i], A[row])]
            row += 1
        assert all(not any(line[:len(basis)]) and not line[-1]
                   for line in A[row:])
        coords = [A[j][-1] for j in range(len(basis))]
        assert all(c.denominator == 1 for c in coords)
        rebuilt = {}
        for c, p in zip(coords, polys):
            if c:
                rebuilt = add(rebuilt, {u: int(v*c) for u, v in p.items()})
        assert rebuilt == target
        l1 = sum(abs(c) for c in coords)
        max_l1 = max(max_l1, l1)
        max_coeff = max(max_coeff, max(abs(c) for c in coords))
    return len(basis), len(allm), max_l1, max_coeff, raw_to_basis_l1


def gaussian_mul(z, w):
    return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])


def gaussian_divides(z, w):
    """Whether Gaussian integer z divides w, by exact norm division."""
    norm = z[0]*z[0]+z[1]*z[1]
    assert norm > 0
    a = w[0]*z[0]+w[1]*z[1]
    b = w[1]*z[0]-w[0]*z[1]
    return a % norm == 0 and b % norm == 0


def gaussian_power(z, e):
    out = (1, 0)
    for _ in range(e):
        out = gaussian_mul(out, z)
    return out


def cancellation_valuation_fixtures():
    # pi=2+i has norm p=5. For H=pi^e and kappa=pi^a, an
    # ordinary integer r satisfies H | kappa*r iff p^(e-min(e,a)) | r.
    pi = (2, 1)
    p = 5
    for e in range(1, 6):
        H = gaussian_power(pi, e)
        for a in range(0, e+2):
            kappa = gaussian_power(pi, a)
            need = max(0, e-a)
            modulus = p**need
            for r in (0, modulus, 2*modulus, modulus//p if need else 1):
                product_value = (kappa[0]*r, kappa[1]*r)
                assert gaussian_divides(H, product_value) == (r % modulus == 0)
    # A conjugate factor in kappa does not cancel the pi-oriented conductor.
    pib = (2, -1)
    for e in range(1, 5):
        H = gaussian_power(pi, e)
        for a in range(e+1):
            for b in range(4):
                kappa = gaussian_mul(gaussian_power(pi, a), gaussian_power(pib, b))
                modulus = p**max(0, e-a)
                r = modulus
                value = (kappa[0]*r, kappa[1]*r)
                assert gaussian_divides(H, value)


def main():
    cases = (
        (6, 1, ((0, 1, 2), (3, 4, 5))),
        (12, 2, ((0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11))),
    )
    for m, k, triangles in cases:
        print("m,k,original_l1,nonzero_cuts,circuit_l1 =", (m, k, *terminal_audit(m, k, triangles)))
    for n in (4, 8):
        print("outside_rows,basis_size,matching_count,max_matching_l1,max_matching_coeff,raw_to_basis_l1 =",
              (n, *matching_basis_bound(n)))
    cancellation_valuation_fixtures()
    print("PASS: exact split-prime cancellation rule after denominator clearing.")
    print("PASS: exact terminal circuit, derivative, reconstruction and matching-basis checks.")


if __name__ == "__main__":
    main()
