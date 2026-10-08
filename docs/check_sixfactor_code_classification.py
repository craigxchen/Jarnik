"""Exact finite classification of five equal-sum six-column bit codes.

One row is normalized to zero by independent column complements.  The four
remaining rows have weight at least three and pairwise Hamming distance at
least three whenever they can be realized by six nonzero coefficients with
distinct absolute values.  We enumerate the resulting 4-cliques, retain
exactly those whose rational kernel is not contained in any forbidden
hyperplane ``a_i=0`` or ``a_i +/- a_j=0``, and quotient by row permutations,
column permutations, and independent column complements.

The exhaustive classification and exact polynomial sign analysis cover
arbitrary real coefficients with distinct nonzero absolute values in this
six-column problem. They do not classify arbitrary circle configurations.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd, isqrt, lcm, prod


N = 6
ROW_PERMS = tuple(permutations(range(5)))


def hamming(a, b):
    return bin(a ^ b).count("1")


def kernel_basis(rows):
    """Return exact rank and a Fraction basis of ker(rows over Q)."""
    matrix = [
        [Fraction((mask >> column) & 1) for column in range(N)]
        for mask in rows
    ]
    rank = 0
    pivots = []
    for column in range(N):
        pivot = next(
            (i for i in range(rank, len(matrix)) if matrix[i][column]), None
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = matrix[rank][column]
        matrix[rank] = [value / scale for value in matrix[rank]]
        for i in range(len(matrix)):
            if i != rank and matrix[i][column]:
                scale = matrix[i][column]
                matrix[i] = [
                    x - scale * y for x, y in zip(matrix[i], matrix[rank])
                ]
        pivots.append(column)
        rank += 1

    free = [column for column in range(N) if column not in pivots]
    basis = []
    for free_column in free:
        vector = [Fraction(0)] * N
        vector[free_column] = 1
        for row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -matrix[row][free_column]
        basis.append(tuple(vector))
    return rank, tuple(basis)


def forbidden_forms():
    forms = []
    for i in range(N):
        form = [0] * N
        form[i] = 1
        forms.append(tuple(form))
    for i, j in combinations(range(N), 2):
        for sign in (1, -1):
            form = [0] * N
            form[i], form[j] = 1, sign
            forms.append(tuple(form))
    return tuple(forms)


FORBIDDEN = forbidden_forms()


def admissible(rows):
    """Test that ker(rows) is not contained in a forbidden hyperplane."""
    rank, basis = kernel_basis(rows)
    for form in FORBIDDEN:
        if all(
            sum(coefficient * coordinate for coefficient, coordinate in zip(form, vector))
            == 0
            for vector in basis
        ):
            return False, rank
    return True, rank


def canonical_code(rows):
    """Canonical columns modulo row/column permutations and column flips."""
    candidates = []
    for row_permutation in ROW_PERMS:
        columns = []
        for column in range(N):
            column_mask = sum(
                ((rows[source_row] >> column) & 1) << output_row
                for output_row, source_row in enumerate(row_permutation)
            )
            # Independent complement of this column identifies mask and its
            # 5-bit complement. Sorting then quotients column permutations.
            columns.append(min(column_mask, column_mask ^ 31))
        candidates.append(tuple(sorted(columns)))
    return min(candidates)


def old_rows():
    subsets = ((4, 6), (1, 3, 6), (1, 4, 5), (2, 3, 5), (1, 2, 3, 4))
    return tuple(sum(1 << (i - 1) for i in subset) for subset in subsets)


def multiply(p, q):
    result = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            result[i+j] += x*y
    return result


def old_orbit_order():
    """Exact sign intervals, including every possible tie and extra row."""
    rows = old_rows()
    coefficients = ((1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (3, 1))
    cubes = [multiply(multiply(a, a), a) for a in coefficients]
    cubic_sums = [tuple(sum(cubes[j][k] for j in range(N) if row >> j & 1)
                        for k in range(4)) for row in rows]
    expected = ((28,30,12,2), (28,27,9,2), (10,15,9,2),
                (16,12,6,2), (10,3,3,2))
    assert tuple(cubic_sums) == expected
    for row in rows:
        assert tuple(sum(coefficients[j][k] for j in range(N) if row >> j & 1)
                     for k in (0,1)) == (4,2)
    # The two displayed coefficient vectors span the full old-code kernel.
    rank, _ = kernel_basis(tuple(row ^ rows[0] for row in rows))
    assert rank == 4
    assert coefficients[0] == (1,0) and coefficients[2] == (0,1)

    roots = set()
    for i, j in combinations(range(5), 2):
        a, b, c = (cubic_sums[i][k]-cubic_sums[j][k] for k in range(3))
        assert cubic_sums[i][3] == cubic_sums[j][3]
        if c:
            disc = b*b-4*a*c
            assert disc.denominator == 1 and disc >= 0
            square = isqrt(disc.numerator)
            assert square*square == disc
            roots.update((Fraction(-b+square,2*c), Fraction(-b-square,2*c)))
        else:
            assert b
            roots.add(Fraction(-a,b))
    assert roots == {Fraction(x) for x in (-4,-3,-2,-1,0,1)} | {Fraction(-3,2)}

    def valid_ratio(r):
        values = [a+b*r for a,b in coefficients]
        return all(values) and len(set(map(abs,values))) == N

    assert not any(valid_ratio(r) for r in roots)
    roots = sorted(roots)
    samples = [roots[0]-1] + [(a+b)/2 for a,b in zip(roots,roots[1:])] + [roots[-1]+1]
    orders = []
    for r in samples:
        values = [sum(a*r**k for k,a in enumerate(p)) for p in cubic_sums]
        assert len(set(values)) == 5
        order = tuple(sorted(range(5),key=values.__getitem__))
        orders.append(order)
        for j in range(N):
            bits = tuple((rows[i] >> j) & 1 for i in order)
            assert bits not in ((0,1,1,1,0),(1,0,0,0,1))
    assert orders == [(4,3,1,2,0),(4,1,3,2,0),(1,4,3,0,2),(1,2,0,3,4),
                      (2,1,0,3,4),(2,4,3,0,1),(4,2,3,1,0),(4,3,2,1,0)]

    exceptional = set()
    persistent = set()
    for row in range(64):
        a,b = (sum(coefficients[j][k] for j in range(N) if row >> j & 1)
               for k in (0,1))
        if b == 2:
            if a == 4:
                persistent.add(row)
        else:
            r = Fraction(4-a,b-2)
            exceptional.add(r)
            assert not valid_ratio(r)
    assert persistent == set(rows)
    assert exceptional == {Fraction(x) for x in (-5,-4,-3,-2,-1,0,1,2)} | {
        Fraction(-5,2),Fraction(-3,2),Fraction(-1,2)}
    print("PASS: all cubic-order intervals avoid the outer-pair cut; no admissible sixth row")
    return orders


def collision_bound_checks(orders):
    rows = old_rows()
    representatives = {tuple(sorted((order[0],order[-1]))): order for order in orders}
    count = 0
    for order in representatives.values():
        for depths in product(range(-2,3), repeat=N):
            values = [sum(depths[j] for j in range(N) if rows[i] >> j & 1)
                      for i in order]
            outer_lo,outer_hi = sorted((values[0],values[4]))
            inner_lo,inner_hi = min(values[1:4]),max(values[1:4])
            gap = max(0,outer_lo-inner_hi,inner_lo-outer_hi)
            assert gap <= sum(map(abs,depths))-max(map(abs,depths))
            count += 1
    print(f"PASS: {count} signed-allocation collision bounds")


def actual_family_checks():
    from functools import reduce
    from check_ordered_affine_circuit_conductor_index import gmul, gconj, ggcd, gnorm, gsub
    from check_ordered_pair_triple_gap_denominator import check_tuple

    checked = 0
    for u,v in ((2,-11),(4,-13),(4,-9),(4,-7),(4,-5),(4,-3),(4,1),(4,7)):
        a = (u,2*u,v,u+v,2*u+v,3*u+v)
        assert all(a) and len(set(map(abs,a))) == N
        cubes = [sum(a[j]**3 for j in range(N) if row >> j & 1) for row in old_rows()]
        expected_order = sorted(range(5),key=cubes.__getitem__)
        resultant = prod(abs(a[j]**2-a[k]**2) for j,k in combinations(range(N),2))
        for t in (1000,1001,1002):
            points = []
            for row in old_rows():
                z = ((-1)**bin(row).count('1'),0)
                for j in range(N):
                    factor = (a[j],t)
                    z = gmul(z,factor if row >> j & 1 else gconj(factor))
                points.append(z)
            assert all(x<0 for x,y in points)
            actual_order = sorted(range(5),key=lambda i:Fraction(points[i][1],points[i][0]))
            assert actual_order == expected_order
            k,_ = check_tuple([points[i] for i in actual_order])
            common_norm = gnorm(reduce(ggcd,points))
            primitive_norm = gnorm(points[0])//common_norm
            chord2 = max(gnorm(gsub(points[i],points[4])) for i in (1,2))//common_norm
            assert 100*chord2**2 >= primitive_norm
            norms = [t*t+x*x for x in a]
            collision = prod(norms)//lcm(*norms)
            pair_gcds = prod(gcd(norms[j],norms[h]) for j,h in combinations(range(N),2))
            assert collision%k == 0 and pair_gcds%collision == 0 and resultant%pair_gcds == 0
            if t == 1000:
                for flip in range(64):
                    flipped_a = [(-x if flip >> j & 1 else x) for j,x in enumerate(a)]
                    for row,point in zip(old_rows(),points):
                        new_row = row ^ flip
                        z = ((-1)**bin(new_row).count('1'),0)
                        for j in range(N):
                            factor = (flipped_a[j],t)
                            z = gmul(z,factor if new_row >> j & 1 else gconj(factor))
                        assert z == point
            checked += 1
    print(f"PASS: {checked} actual fixed-coefficient families, arc order and denominator divisibility")
    print("PASS: exact column-flip transport and primitive small-arc exclusion fixtures")


def main():
    candidates = tuple(
        mask for mask in range(1, 1 << N) if bin(mask).count("1") >= 3
    )
    all_four_sets = compatible_four_sets = retained = 0
    ranks = Counter()
    classes = defaultdict(int)

    for four_rows in combinations(candidates, 4):
        all_four_sets += 1
        if any(hamming(a, b) < 3 for a, b in combinations(four_rows, 2)):
            continue
        compatible_four_sets += 1
        ok, rank = admissible((0,) + four_rows)
        if not ok:
            continue
        retained += 1
        ranks[rank] += 1
        classes[canonical_code((0,) + four_rows)] += 1

    old_class = canonical_code(old_rows())
    assert len(candidates) == 42
    assert all_four_sets == 111930
    assert compatible_four_sets == 3150
    assert retained == 1800
    assert len(classes) == 1
    assert ranks == Counter({4: 1800})
    assert old_class in classes
    assert classes[old_class] == 1800
    orders = old_orbit_order()
    collision_bound_checks(orders)
    actual_family_checks()

    print("PASS: 42 candidates, 111930 four-sets, 3150 compatible cliques")
    print("PASS: 1800 admissible codes, all rank 4 / nullity 2")
    print("PASS: one symmetry orbit; old six-factor code is its representative")
    print("canonical column masks:", old_class)


if __name__ == "__main__":
    main()
