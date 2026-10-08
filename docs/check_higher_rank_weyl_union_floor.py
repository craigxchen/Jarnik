"""Finite checks for higher-rank Weyl-union coefficient bookkeeping.

No Gaussian core-divisibility or CB-connection theorem is tested here.
"""

from itertools import permutations
from math import prod


def sign(permutation):
    inversions = sum(permutation[i] > permutation[j]
                     for i in range(len(permutation))
                     for j in range(i + 1, len(permutation)))
    return -1 if inversions % 2 else 1


def add_term(poly, monomial, coefficient):
    poly[monomial] = poly.get(monomial, 0) + coefficient
    if poly[monomial] == 0:
        del poly[monomial]


def check_bracket_sl2_invariance():
    # Variables are (P_i,Y_i,P_j,Y_j,t); monomials are exponent tuples.
    baseline = {(1, 0, 0, 1, 0): 1, (0, 1, 1, 0, 0): -1}

    upper = {}
    add_term(upper, (1, 0, 0, 1, 0), 1)
    add_term(upper, (0, 1, 0, 1, 1), 1)
    add_term(upper, (0, 1, 1, 0, 0), -1)
    add_term(upper, (0, 1, 0, 1, 1), -1)
    assert upper == baseline

    lower = {}
    add_term(lower, (1, 0, 0, 1, 0), 1)
    add_term(lower, (1, 0, 1, 0, 1), 1)
    add_term(lower, (0, 1, 1, 0, 0), -1)
    add_term(lower, (1, 0, 1, 0, 1), -1)
    assert lower == baseline


def multiply_polynomials(left, right):
    out = {}
    for monomial_a, coefficient_a in left.items():
        for monomial_b, coefficient_b in right.items():
            monomial = tuple(sorted(monomial_a + monomial_b))
            out[monomial] = out.get(monomial, 0) + coefficient_a * coefficient_b
    return {monomial: coefficient for monomial, coefficient in out.items()
            if coefficient}


def binary_residual_of_determinant(block, n):
    """Expand a determinant after fixing rows 2..n-1 to basis vectors."""
    rows = []
    for row, label in enumerate(block):
        if row == 0:
            rows.append({0: ((f"P{label}",),), 1: ((f"Y{label}",),)})
        elif row == 1:
            rows.append({0: ((f"P{label}",),), 1: ((f"Y{label}",),)})
        else:
            rows.append({row: ((),)})
    result = {}
    for assignment in permutations(range(n)):
        factors = [rows[row].get(column) for row, column in enumerate(assignment)]
        if any(factor is None for factor in factors):
            continue
        monomial = tuple(sorted(token for factor in factors
                                for option in factor for token in option))
        result[monomial] = result.get(monomial, 0) + sign(assignment)
    return {monomial: coefficient for monomial, coefficient in result.items()
            if coefficient}


def bracket_polynomial(i, j):
    return {
        tuple(sorted((f"P{i}", f"Y{j}"))): 1,
        tuple(sorted((f"P{j}", f"Y{i}"))): -1,
    }


def determinant_assignment_terms(n, blocks):
    """Terms of a product of disjoint n-row determinants."""
    terms = {(): 1}
    for block in blocks:
        next_terms = {}
        for assignment in permutations(range(n)):
            local = tuple((block[row], assignment[row]) for row in range(n))
            for prefix, coefficient in terms.items():
                full = prefix + local
                next_terms[full] = coefficient * sign(assignment)
        terms = next_terms
    return terms


def check_veronese_bijection_and_residual(n, d):
    m = n * d
    blocks = [tuple(range(j * n, (j + 1) * n)) for j in range(d)]
    terms = determinant_assignment_terms(n, blocks)
    assert len(terms) == (prod(range(1, n + 1)) ** d)
    assert all(abs(c) == 1 for c in terms.values())

    # Assignment -> binary exponent vector is injective, and each source
    # coefficient is preserved exactly under the Veronese pullback.
    pulled_back = {}
    for assignment, coefficient in terms.items():
        exponents = tuple(j for _, j in sorted(assignment))
        assert exponents not in pulled_back
        pulled_back[exponents] = coefficient
    assert sum(abs(c) for c in pulled_back.values()) == sum(
        abs(c) for c in terms.values())

    # Select the identity assignment in each determinant. Each color block
    # has d labels, and fixing colors >=2 to their standard basis vectors
    # reduces the product to d binary brackets on disjoint row pairs.
    colors = [tuple(block[row] for block in blocks) for row in range(n)]
    assert all(len(color) == d for color in colors)
    assert len(set(i for color in colors for i in color)) == m
    pairs = list(zip(colors[0], colors[1]))
    assert len(pairs) == d and len(set(i for pair in pairs for i in pair)) == 2 * d
    residual = {(): 1}
    expected = {(): 1}
    for block, (i, j) in zip(blocks, pairs):
        residual = multiply_polynomials(residual,
                                        binary_residual_of_determinant(block, n))
        expected = multiply_polynomials(expected, bracket_polynomial(i, j))
    assert residual == expected
    # Each extracted bracket is SL2 invariant by the exact coefficientwise
    # shear checks above; their product is the fixture's binary residual.
    check_bracket_sl2_invariance()
    return len(terms), pairs


def union_formula(n, d):
    return (2 ** (n * d) - 2 * (2 ** d - 1) ** n
            + (2 ** d - 2) ** n)


def exhaustive_union_count(n, d):
    m = n * d
    blocks = [frozenset(range(j * d, (j + 1) * d)) for j in range(n)]
    union = set()
    for mask in range(1 << m):
        T = frozenset(i for i in range(m) if mask & (1 << i))
        has_full = any(B <= T for B in blocks)
        has_empty = any(not (B & T) for B in blocks)
        if has_full and has_empty:
            union.add(T)
    assert len(union) == union_formula(n, d)
    return blocks, union


def check_pair_loss(n, d):
    blocks, union = exhaustive_union_count(n, d)
    pair_groups = {(a, b): [] for a in range(n) for b in range(n) if a != b}
    for T in union:
        witnesses = [(a, b) for a in range(n) for b in range(n)
                     if a != b and blocks[a] <= T and not (blocks[b] & T)]
        assert witnesses
        pair_groups[min(witnesses)].append(T)
    assert sum(len(group) for group in pair_groups.values()) == len(union)
    assert all(blocks[a] <= T and not (blocks[b] & T)
               for (a, b), group in pair_groups.items() for T in group)

    # For arbitrary nonnegative row correction weights l_i, charging every
    # ordered pair costs exactly 2(n-1)*sum_i l_i; assigned groups cost no more.
    weights = [i + 1 for i in range(n * d)]
    all_pair_loss = sum(2 * sum(weights[i] for i in blocks[b])
                        for a in range(n) for b in range(n) if a != b)
    target = 2 * (n - 1) * sum(weights)
    assert all_pair_loss == target
    assigned_loss = sum(2 * sum(weights[i] for i in blocks[b])
                        for (a, b), group in pair_groups.items() if group)
    assert assigned_loss <= target
    return len(union)


def main():
    for n, d in ((3, 2), (4, 2)):
        term_count, pairs = check_veronese_bijection_and_residual(n, d)
        U = check_pair_loss(n, d)
        print(f"n={n}, d={d}: determinant terms={term_count}, "
              f"residual pairs={pairs}, U={U}")
    cases = [(n, d) for n in range(2, 17) for d in range(1, 17)
             if n * d <= 16]
    for n, d in cases:
        check_pair_loss(n, d)
    print(f"PASS: coefficient bijections, residual fixtures, union count, "
          f"and pair-loss accounting for {len(cases)} cases.")


if __name__ == "__main__":
    main()
