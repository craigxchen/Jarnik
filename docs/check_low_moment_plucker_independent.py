"""Independent exact audit of the seven/eight-point low-moment kernel.

Uses a direct rational nullspace, homogeneous quartic evaluation minors,
and the weighted Gale matrix as three independent presentations. Formal
fair-cut enumeration is kept separate from the actual circle examples.
"""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, isqrt, lcm, prod

from check_distance_laplacian_denominator import determinant


def det2(v, w):
    return v[0] * w[1] - v[1] * w[0]


def primitive(values):
    denominator = lcm(*(x.denominator for x in values))
    raw = [int(x * denominator) for x in values]
    content = reduce(gcd, raw)
    return [x // content for x in raw]


def primes_of(value):
    value = abs(value)
    primes = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            primes.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        primes.add(value)
    return primes


def valuation(value, prime):
    value = abs(value)
    assert value
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def nullspace(matrix):
    a = [[F(value) for value in row] for row in matrix]
    height, width = len(a), len(a[0])
    pivots = []
    row = 0
    for col in range(width):
        pivot = next((i for i in range(row, height) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        divisor = a[row][col]
        a[row] = [value / divisor for value in a[row]]
        for i in range(height):
            if i != row:
                scale = a[i][col]
                a[i] = [x - scale * y for x, y in zip(a[i], a[row])]
        pivots.append(col)
        row += 1
        if row == height:
            break
    assert row == 5
    basis = []
    for free in (col for col in range(width) if col not in pivots):
        vector = [F(0)] * width
        vector[free] = F(1)
        for i, pivot in enumerate(pivots):
            vector[pivot] = -a[i][free]
        basis.append(vector)
    return basis


def rational_sqrt(value):
    assert value > 0
    a, b = isqrt(value.numerator), isqrt(value.denominator)
    assert a * a == value.numerator and b * b == value.denominator
    return F(a, b)


def check_seven_inverse(kernel, original_chords):
    """Recover the circle metric from just the rational rank-two kernel."""
    rows = list(zip(*kernel))
    brackets = {(i, j): det2(rows[i], rows[j]) for i, j in combinations(range(7), 2)}
    stars = [prod(det2(rows[i], rows[j]) for j in range(7) if j != i) for i in range(7)]
    values = [rational_sqrt(star / stars[0]) for star in stars]
    primitive_brackets = dict(zip(brackets, primitive(list(brackets.values()))))
    def oriented(i, j):
        return primitive_brackets[i, j] if i < j else -primitive_brackets[j, i]
    integer_stars = [prod(oriented(i, j) for j in range(7) if i != j) for i in range(7)]
    common_content = reduce(gcd, integer_stars)
    metric_values = [isqrt(abs(star) // common_content) for star in integer_stars]
    assert all(value * value == abs(star) // common_content for value, star in zip(metric_values, integer_stars))
    assert reduce(gcd, metric_values) == 1
    assert values == [F(value, metric_values[0]) for value in metric_values]
    coefficients = [[x * x, 2 * x * y, y * y] for x, y in rows[:3]]
    denominator = determinant(coefficients)
    assert denominator
    metric = []
    for col in range(3):
        replaced = [[values[i] if j == col else coefficients[i][j] for j in range(3)] for i in range(3)]
        metric.append(determinant(replaced) / denominator)
    aa, bb, cc = metric
    metric_det = aa * cc - bb * bb
    assert aa > 0 and metric_det > 0
    rational_sqrt(metric_det)
    assert all(aa * x * x + 2 * bb * x * y + cc * y * y == value for (x, y), value in zip(rows, values))
    reconstructed = {(i, j): 4 * metric_det * bracket ** 2 / (values[i] * values[j]) for (i, j), bracket in brackets.items()}
    assert reconstructed == original_chords
    # The quartic is stars[0] times the square of the recovered quadratic.
    assert all(stars[0] * value ** 2 == star for value, star in zip(values, stars))


def check_case(form, vectors):
    n = len(vectors)
    rank = n - 5
    aa, bb, cc = form
    radius_form = isqrt(aa * cc - bb * bb)
    assert radius_form ** 2 == aa * cc - bb * bb > 0
    qs = [aa * x * x + 2 * bb * x * y + cc * y * y for x, y in vectors]
    delta = {(i, j): det2(vectors[i], vectors[j]) for i, j in combinations(range(n), 2)}
    assert all(delta.values())
    evaluation = [[F(x ** (4 - power) * y ** power, q ** 2) for (x, y), q in zip(vectors, qs)] for power in range(5)]
    subsets = list(combinations(range(n), 5))
    raw = [F(prod(delta[i, j] for i, j in combinations(subset, 2)), prod(qs[i] ** 2 for i in subset)) for subset in subsets]
    assert raw == [determinant([[row[i] for i in subset] for row in evaluation]) for subset in subsets]
    p_eval = primitive(raw)
    lookup = {tuple(sorted(set(range(n)) - set(subset))): abs(p) for subset, p in zip(subsets, p_eval)}
    kernel = nullspace(evaluation)
    small_subsets = list(combinations(range(n), rank))
    raw_kernel = [determinant([[row[i] for i in subset] for row in kernel]) for subset in small_subsets]
    p_kernel = primitive(raw_kernel)
    assert [abs(p) for p in p_kernel] == [lookup[subset] for subset in small_subsets]

    phases = [(F((aa * x + bb * y) ** 2 - radius_form ** 2 * y ** 2, aa * q), F(2 * radius_form * y * (aa * x + bb * y), aa * q)) for (x, y), q in zip(vectors, qs)]
    assert all(x * x + y * y == 1 for x, y in phases)
    circle_columns = [(F(1), x, y, x * x - y * y, 2 * x * y) for x, y in phases]
    circle_rows = [[column[j] for column in circle_columns] for j in range(5)]
    circle_kernel = nullspace(circle_rows)
    circle_pluckers = primitive([determinant([[row[i] for i in subset] for row in circle_kernel]) for subset in small_subsets])
    assert [abs(p) for p in circle_pluckers] == [abs(p) for p in p_kernel]

    stars = [prod(det2(vectors[i], vectors[j]) for j in range(n) if i != j) for i in range(n)]
    gale = [[F(qs[i] ** 2, stars[i]) * x ** (rank - 1 - power) * y ** power for power in range(rank)] for i, (x, y) in enumerate(vectors)]
    assert all(sum(evaluation[a][i] * gale[i][b] for i in range(n)) == 0 for a in range(5) for b in range(rank))
    raw_gale = [determinant([gale[i] for i in subset]) for subset in small_subsets]
    p_gale = primitive(raw_gale)
    assert [abs(p) for p in p_gale] == [abs(p) for p in p_kernel]
    if n == 8:
        central = primitive([F(q ** 3, star) for q, star in zip(qs, stars)])
        t = [[F(central[i], qs[i]) * x ** (2 - power) * y ** power for power in range(3)] for i, (x, y) in enumerate(vectors)]
        assert [abs(p) for p in primitive([determinant([t[i] for i in subset]) for subset in small_subsets])] == [abs(p) for p in p_kernel]

    primes = {2}
    for value in qs + list(delta.values()) + [radius_form]:
        primes |= primes_of(value)
    for prime in sorted(primes):
        rho = valuation(radius_form, prime)
        lam = {(i, j): valuation(qs[i], prime) + valuation(qs[j], prime) - 2 * rho - 2 * valuation(d, prime) for (i, j), d in delta.items()}
        energies = [sum(lam[i, j] for i, j in combinations(subset, 2)) for subset in subsets]
        maximum = max(energies)
        for p, energy in zip(p_eval, energies):
            assert maximum % 2 == energy % 2
            assert 2 * valuation(p, prime) == maximum - energy
        actual_chord_energies = [energy - (20 if prime == 2 else 0) for energy in energies]
        assert [max(actual_chord_energies) - e for e in actual_chord_energies] == [maximum - e for e in energies]

    chords = {(i, j): F(4 * radius_form ** 2 * d ** 2, qs[i] * qs[j]) for (i, j), d in delta.items()}
    if n == 7:
        check_seven_inverse(kernel, chords)
    n_intrinsic = lcm(*(chord.denominator for chord in chords.values()))
    h_squared = sum(p * p for p in p_kernel)
    return n_intrinsic, h_squared, len(primes)


def fair_coefficients(n):
    all_fives = list(combinations(range(n), 5))
    maximum_total = 0
    fixed_total = 0
    for mask in range(1, 2 ** (n - 1)):
        cut = {i + 1 for i in range(n - 1) if mask >> i & 1}
        energies = [len(cut.intersection(subset)) * (5 - len(cut.intersection(subset))) for subset in all_fives]
        maximum_total += max(energies)
        fixed_total += energies[0]
    assert maximum_total == 6 * (2 ** (n - 1) - 1) - 2 * n
    assert fixed_total == 10 * 2 ** (n - 2)
    assert (maximum_total - fixed_total) // 2 == 2 ** (n - 2) - n - 3
    return maximum_total, fixed_total, (maximum_total - fixed_total) // 2


def main():
    count = 0
    prime_checks = 0
    for n in (7, 8):
        vector_sets = [[(i, 1) for i in range(n)], [(1, 0)] + [(i, 1) for i in range(n - 1)], [(1, 2), (2, 3), (3, 5), (5, 7), (7, 11), (11, 13), (13, 17), (17, 19)][:n]]
        for form in ((1, 0, 1), (1, 1, 2), (2, 1, 5)):
            for vectors in vector_sets:
                _, _, checked = check_case(form, vectors)
                count += 1
                prime_checks += checked
    n, h2, _ = check_case((1, 0, 1), [(i, 1) for i in range(8)])
    assert (n, h2) == (1022125, 1053718749231052896)
    assert h2 % n == 879896
    print(f"Toy t=0..7: N={n}; H^2={h2}; H^2 mod N={h2 % n}.")
    for size in range(6, 11):
        print(f"Formal fair n={size}: max-energy sum, fixed-five energy sum, log-height coefficient = {fair_coefficients(size)}")
    print(f"Passed {count} actual seven/eight-point examples and {prime_checks} complete prime checks, including 2.")


if __name__ == '__main__':
    main()
