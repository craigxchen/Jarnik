"""Exact, standard-library certificate for the six-row invariant example.

This checks polynomial restrictions only. It does not assert a numerical
zero relation at any lattice configuration.
"""

from itertools import combinations


N = 6
ZERO = (0,) * N
PARTITIONS = [
    ((0,) + pair, tuple(j for j in range(N) if j not in (0,) + pair))
    for pair in combinations(range(1, N), 2)
]


def add(a, b):
    result = dict(a)
    for monomial, coefficient in b.items():
        result[monomial] = result.get(monomial, 0) + coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def scale(a, coefficient):
    return {monomial: coefficient * value for monomial, value in a.items() if coefficient * value}


def mul(a, b):
    result = {}
    for monomial_a, coefficient_a in a.items():
        for monomial_b, coefficient_b in b.items():
            monomial = tuple(x + y for x, y in zip(monomial_a, monomial_b))
            result[monomial] = result.get(monomial, 0) + coefficient_a * coefficient_b
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def variable(j):
    monomial = tuple(int(i == j) for i in range(N))
    return {monomial: 1}


def bracket(i, j, inside):
    """det(v_i,v_j), with v_i=e1 inside and v_j=(z_j,1) outside."""
    if i in inside and j in inside:
        return {}
    if i in inside:
        return {ZERO: 1}
    if j in inside:
        return {ZERO: -1}
    return add(variable(i), scale(variable(j), -1))


def restriction(inside):
    result = {}
    for index, (first, second) in enumerate(PARTITIONS):
        term = {ZERO: 1}
        for triple in (first, second):
            a, b, c = triple
            for i, j in ((a, b), (b, c), (c, a)):
                term = mul(term, bracket(i, j, inside))
        result = add(result, scale(term, 2**index))
    return result


def quadratic_coefficient(polynomial, i, j):
    monomial = tuple(int(k == i) + int(k == j) for k in range(N))
    return polynomial.get(monomial, 0)


def main():
    assert len(PARTITIONS) == 10
    for inside in combinations(range(N), 3):
        assert not restriction(set(inside))

    rows = []
    for inside in combinations(range(N), 2):
        polynomial = restriction(set(inside))
        outside = [j for j in range(N) if j not in inside]
        a, b, c, d = outside
        alpha = quadratic_coefficient(polynomial, a, c)
        beta = quadratic_coefficient(polynomial, a, b)
        first = mul(add(variable(a), scale(variable(b), -1)),
                    add(variable(c), scale(variable(d), -1)))
        second = mul(add(variable(a), scale(variable(c), -1)),
                     add(variable(b), scale(variable(d), -1)))
        assert polynomial == add(scale(first, alpha), scale(second, beta))
        assert alpha * beta * (alpha + beta) != 0
        for i in outside:
            assert quadratic_coefficient(polynomial, i, i) == 0
            assert sum(quadratic_coefficient(polynomial, i, j) for j in outside if j != i) == 0
        rows.append((tuple(i + 1 for i in inside), alpha, beta, alpha + beta))

    assert 64 * sum(2**j for j in range(10)) == 65472
    print("20 balanced restrictions vanish; all 15 smaller restrictions have rank three.")
    print("S, alpha, beta, alpha+beta (one-based row labels):")
    for row in rows:
        print(row)
    print("No numerical zero relation or lattice construction is asserted.")


if __name__ == "__main__":
    main()
