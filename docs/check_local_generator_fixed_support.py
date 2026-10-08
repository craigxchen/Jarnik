"""Exact checks for fixed-support determinant heights and EN dimensions.

The finite checks audit the combinatorial content and the proposed
Eagon--Northcott resolution dimensions. They do not replace the general
proofs in local_generator_support_family_diagnostic.md.
"""

from fractions import Fraction
from math import comb, factorial, gcd, isqrt
from pathlib import Path
import importlib.util
from itertools import combinations


HERE = Path(__file__).resolve().parent
FUSION_PATH = HERE / "check_higher_rank_terminal_fusion_grid.py"
SPEC = importlib.util.spec_from_file_location("fusion_grid", FUSION_PATH)
FUSION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FUSION)


def hook_dimension(partition):
    """Dimension of the S_N irreducible indexed by a partition."""
    rows = tuple(partition)
    size = sum(rows)
    conjugate = [sum(row > col for row in rows)
                 for col in range(max(rows, default=0))]
    hooks = []
    for i, row in enumerate(rows):
        for col in range(row):
            hooks.append(row - col + conjugate[col] - i - 1)
    numerator = factorial(size)
    denominator = 1
    for hook in hooks:
        denominator *= hook
    assert numerator % denominator == 0
    return numerator // denominator


def en_a_values(n, d):
    """Return a_j from the proposed resolution formula."""
    e = d - 2
    values = []
    for j in range(2 * e + 1):
        total = 0
        for a in range(e + 1):
            b = j - a
            if 0 <= b <= a:
                shape = (e,) * (n - 2) + (e - b, e - a)
                assert all(shape[i] >= shape[i + 1]
                           for i in range(len(shape) - 1))
                total += (a - b + 1) * hook_dimension(shape)
        values.append(total)
    return tuple(values)


def en_dimensions(n, d):
    m = n * d
    a_values = en_a_values(n, d)
    h_alt = sum((-1) ** j * comb(m, 2 * n + j) * a
                for j, a in enumerate(a_values)
                if 0 <= 2 * n + j <= m)
    delta_alt = sum((-1) ** j * comb(m - 1, 2 * n + j - 1) * a
                    for j, a in enumerate(a_values)
                    if 0 <= 2 * n + j - 1 <= m - 1)
    return h_alt, delta_alt, a_values


def matching_internal_minimum(size, inside):
    """Enumerate perfect matchings and return minimum internal edges."""
    labels = tuple(range(size))
    selected = set(range(inside))
    minima = []

    def visit(remaining, internal):
        if not remaining:
            minima.append(internal)
            return
        first = remaining[0]
        for pos in range(1, len(remaining)):
            second = remaining[pos]
            rest = remaining[1:pos] + remaining[pos + 1:]
            visit(rest, internal + int(first in selected and second in selected))

    visit(labels, 0)
    return min(minima)


def check_occupancy_and_gcd():
    for n in range(2, 5):
        size = 2 * n
        for inside in range(size + 1):
            assert matching_internal_minimum(size, inside) == max(0, inside - n)
        combinatorial = sum(comb(size, t) * max(0, t - n)
                            for t in range(size + 1))
        # For any fixed 2n-support I among m labels, choose the t labels of
        # T inside I and any subset of the remaining labels outside I.
        for m in range(size, size + 5):
            support = tuple(range(size))
            brute = 0
            for mask in range(1 << m):
                T = {i for i in range(m) if mask & (1 << i)}
                brute += max(0, len(T.intersection(support)) - n)
            formula = (1 << (m - size)) * combinatorial
            assert brute == formula
    # Concrete full-profile (3,4) accounting.
    n, d, m = 3, 4, 12
    core_gcd_exp = (1 << (m - 2 * n)) * sum(
        comb(2 * n, t) * (t - n) for t in range(n + 1, 2 * n + 1))
    raw_matching_exp = n * (1 << (m - 2))
    primitive_exp = raw_matching_exp - core_gcd_exp
    B32 = n * (1 << (2 * n - 2)) - sum(
        comb(2 * n, t) * (t - n) for t in range(n + 1, 2 * n + 1))
    assert (core_gcd_exp, raw_matching_exp, primitive_exp, B32) == (1920, 3072, 1152, 18)
    supports = comb(12, 6)
    complement_rank = FUSION.source_dimension(3, 2)
    assert (supports, complement_rank, supports * complement_rank) == (924, 5, 4620)
    assert FUSION.case(3, 4)["h"] == 341
    assert FUSION.case(3, 4)["threshold"] == 90


def determinant(matrix):
    """Small exact determinant by fraction-free elimination."""
    a = [list(map(int, row)) for row in matrix]
    n = len(a)
    sign = 1
    previous = 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
        for i in range(k + 1, n):
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def gram(columns):
    return [[sum(x * y for x, y in zip(left, right))
             for right in columns] for left in columns]


def tensor(left, right):
    return [x * y for x in left for y in right]


def check_tensor_saturation_and_gram():
    # Primitive v and a non-orthogonal saturated rank-two lattice W.
    v = (2, 3, 5)
    assert gcd(gcd(*v[:2]), v[2]) == 1
    W = ((1, 0), (0, 1), (1, 1))  # coordinate rows, basis columns
    w_columns = tuple(tuple(row[j] for row in W) for j in range(2))
    gram_w = gram(w_columns)
    assert gram_w == [[2, 1], [1, 2]]
    det_w = determinant(gram_w)
    norm_v2 = sum(x * x for x in v)
    product_columns = [tensor(v, w) for w in w_columns]
    gram_product = gram(product_columns)
    assert gram_product == [[norm_v2 * x for x in row] for row in gram_w]
    assert determinant(gram_product) == norm_v2 ** 2 * det_w
    assert norm_v2 == 38 and det_w == 3

    # The matrix defining W is primitive/saturated: one 2x2 minor is 1.
    minors = [determinant((W[i], W[j])) for i, j in combinations(range(3), 2)]
    assert gcd(gcd(*minors[:2]), minors[2]) == 1

    # Bézout recovery: if v tensor q is integral then q is integral.
    # Here 2*(-1)+3*1=1, so the same combination of tensor slices
    # recovers every coordinate of q.
    for denominator in range(1, 8):
        for x in range(-4, 5):
            q = Fraction(x, denominator)
            if all((coord * q).denominator == 1 for coord in v):
                recovered = -((v[0] * q).numerator // (v[0] * q).denominator) + \
                    ((v[1] * q).numerator // (v[1] * q).denominator)
                assert recovered == q and q.denominator == 1


def check_en_against_fusion():
    tested = 0
    for n in range(2, 9):
        for d in range(2, 13):
            h_en, delta_en, _ = en_dimensions(n, d)
            fusion = FUSION.case(n, d)
            assert h_en == fusion["h"], ("h mismatch", n, d, h_en, fusion["h"])
            assert delta_en == fusion["delta"], (
                "delta mismatch", n, d, delta_en, fusion["delta"])
            tested += 1
    assert tested == 77


def main():
    check_occupancy_and_gcd()
    check_tensor_saturation_and_gram()
    check_en_against_fusion()
    print("PASS: matching occupancy minima and cut-content exponent sums.")
    print("PASS: (3,4) exact values: raw=3072, gcd=1920, primitive=1152; "
          "supports=924, fixed-support rank=5, total generators=4620, "
          "kernel rank=341, threshold=90.")
    print("PASS: primitive tensor saturation recovery and exact Gram scaling.")
    print("PASS: EN resolution h and delta match fusion in all 77 cases "
          "(2<=n<=8, 2<=d<=12).")


if __name__ == "__main__":
    main()
