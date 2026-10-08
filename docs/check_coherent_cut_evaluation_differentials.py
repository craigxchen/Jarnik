"""Exact finite checks for coherent cut differential identities.

The general proof uses the Specht filtration. These checks evaluate actual
triangle/pair graph invariants, differentiate their coherent cut arrays,
check the terminal integral identities, and verify the height comparison.
No short relation or full-profile counterexample is asserted.
"""

from itertools import combinations, product
from math import comb, factorial
from fractions import Fraction
from random import Random


def derivative(polynomial, weights):
    result = {}
    for support, coefficient in polynomial.items():
        for i in support:
            smaller = support - {i}
            result[smaller] = result.get(smaller, 0) + coefficient * weights[i]
    return {support: coefficient for support, coefficient in result.items()
            if coefficient}


def evaluate(polynomial, values):
    total = 0
    for support, coefficient in polynomial.items():
        for i in support:
            coefficient *= values[i]
        total += coefficient
    return total


def graph_blocks(m, k, ordering):
    triangles = [tuple(ordering[3*i:3*i+3]) for i in range(2*k)]
    pairs = [tuple(ordering[i:i+2]) for i in range(6*k, m, 2)]
    return triangles, pairs


def coherent_array(triangles, pairs, b):
    blocks = triangles + pairs
    result = {}
    for choices in product(*blocks):
        coefficient = 1
        for (a, c, d), inside in zip(triangles, choices):
            coefficient *= {
                a: b[d] - b[c],
                c: -(b[d] - b[a]),
                d: b[c] - b[a],
            }[inside]
        result[frozenset(choices)] = coefficient
    return result


def bracket(rows, i, j):
    return rows[i][0]*rows[j][1] - rows[j][0]*rows[i][1]


def terminal_arrays(triangles, rows):
    r = {}
    e = {}
    for choices in product(*triangles):
        coefficient = 1
        for (a, b, c), inside in zip(triangles, choices):
            coefficient *= {
                a: bracket(rows, b, c),
                b: -bracket(rows, a, c),
                c: bracket(rows, a, b),
            }[inside]
        support = frozenset(choices)
        r[support] = coefficient
        for i, row in enumerate(rows):
            if i not in support:
                coefficient *= row[0]
        e[support] = coefficient
    return r, e


def graph_identity_checks():
    rng = Random(438)
    count = 0
    for m in (6, 8, 10, 12, 14, 18):
        b = [i*i + 3*i + 1 for i in range(m)]
        for k in range(1, m//6 + 1):
            p = m//2 - 3*k
            for _ in range(8):
                ordering = list(range(m))
                rng.shuffle(ordering)
                triangles, pairs = graph_blocks(m, k, ordering)
                polynomial = coherent_array(triangles, pairs, b)
                assert all(len(support) == m//2-k for support in polynomial)
                for j in range(p+2):
                    differentiated = polynomial
                    for _ in range(j):
                        differentiated = derivative(differentiated, b)
                    for _ in range(p+1-j):
                        differentiated = derivative(differentiated, [1]*m)
                    assert not differentiated, (m, k, j)
                count += 1
    print(f"PASS: all mixed differential identities on {count} graph fixtures.")


def terminal_identity_checks():
    for k in (1, 2, 3):
        m = 6*k
        rows = [(i+1, i*i+2*i+3) for i in range(m)]
        triangles, pairs = graph_blocks(m, k, list(range(m)))
        assert not pairs
        r, e = terminal_arrays(triangles, rows)
        assert not derivative(r, [row[0] for row in rows])
        assert not derivative(r, [row[1] for row in rows])
        assert not derivative(e, [row[0]**2 for row in rows])
        assert not derivative(e, [row[0]*row[1] for row in rows])
        actual = 1
        for a, b, c in triangles:
            actual *= (bracket(rows, a, b)*bracket(rows, a, c)
                       * bracket(rows, b, c))
        assert evaluate(e, [row[1]**2 for row in rows]) == actual
    print("PASS: terminal denominator-free derivatives and evaluation identity.")


def matrix_identity_checks():
    for m, k in ((6, 1), (8, 1), (10, 1), (12, 2)):
        b = [i*i+3*i+1 for i in range(m)]
        p = m//2-3*k
        triangles, pairs = graph_blocks(m, k, list(range(m)))
        polynomial = coherent_array(triangles, pairs, b)
        for U in combinations(range(m), 2*k-1):
            U = frozenset(U)
            for j in range(p+2):
                value = 0
                for S, coefficient in polynomial.items():
                    if not U <= S:
                        continue
                    elementary = 0
                    for selected in combinations(S-U, j):
                        term = 1
                        for i in selected:
                            term *= b[i]
                        elementary += term
                    value += coefficient*elementary
                assert value == 0
    print("PASS: explicit elementary-symmetric constraint matrix.")


def height_checks():
    for m in range(6, 402, 2):
        q = m//2
        weights = [comb(m, j) for j in range(m+1)]
        for k in range(1, m//6+1):
            weighted_max = sum(weights[j]*max(k, abs(j-q))
                               for j in range(m+1))
            A = m*2**(m-2)-weighted_max
            base = m*2**(m-2)-q*comb(m, q)
            discount = sum(weights[j]*max(0, k-abs(j-q))
                           for j in range(m+1))
            assert A == base-discount
            assert 24*A >= 2**m*(2*m-9)
            N = comb(m, q-k)
            assert A > N
            if m >= 8:
                assert A > 2*N
    print("PASS: all-depth height identities and deficit bounds through m=400.")


def circuit_height_checks():
    expected = {1: 18, 2: 2600, 3: 266560}
    for k in range(1, 13):
        tail = sum(comb(4*k+1, a)*max(0, a-2*k-1)
                   for a in range(4*k+2))
        assert 2*tail == (4*k+1)*comb(4*k, 2*k)-2**(4*k)
        J = 2**(2*k-1)*tail
        H = 2*k*2**(6*k-2)-J
        assert H == 2**(2*k-2)*((2*k+1)*2**(4*k)
                                -(4*k+1)*comb(4*k, 2*k))
        if k in expected:
            assert H == expected[k]
    print("PASS: universal circuit-content sums and primitive-height upper bounds.")


def coherent_kernel_covolume_checks():
    for k in range(2, 13):
        m = 6*k
        r = 2*factorial(m)//(factorial(2*k+2)*factorial(2*k+1)*factorial(2*k))
        rho = comb(m, 2*k)-comb(m, 2*k-1)
        H = k*2**(2*k)*(2**(4*k-1)-comb(4*k, 2*k))
        assert 0 < rho < r
        ratio = Fraction(rho*H, (r-rho)*2**(2*k+1))
        assert ratio > 1
    print("PASS: basic coherent-kernel covolume comparison for k=2,...,12.")


if __name__ == "__main__":
    graph_identity_checks()
    terminal_identity_checks()
    matrix_identity_checks()
    height_checks()
    circuit_height_checks()
    coherent_kernel_covolume_checks()
