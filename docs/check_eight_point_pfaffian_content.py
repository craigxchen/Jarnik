"""Exact supplementary checks for eight_point_pfaffian_content.md."""
from collections import Counter
from functools import reduce
from fractions import Fraction
from itertools import combinations
from math import factorial, gcd


def prod(values):
    return reduce(lambda a, b: a * b, values, 1)


def rational_rank(rows):
    basis = {}
    for original in rows:
        row = list(map(Fraction, original))
        for pivot in sorted(basis):
            if row[pivot]:
                scalar = row[pivot]
                row = [a - scalar * b for a, b in zip(row, basis[pivot])]
        pivot = next((j for j, entry in enumerate(row) if entry), None)
        if pivot is not None:
            scalar = row[pivot]
            basis[pivot] = [a / scalar for a in row]
    return len(basis)


def matchings(vertices):
    if not vertices:
        yield (), 1
        return
    a = vertices[0]
    for pos in range(1, len(vertices)):
        b = vertices[pos]
        rest = vertices[1:pos] + vertices[pos + 1:]
        for edges, sign in matchings(rest):
            yield ((a, b),) + edges, sign * (-1) ** (pos - 1)


EDGES = tuple(combinations(range(8), 2))
MATCHINGS = tuple(matchings(tuple(range(8))))
assert len(MATCHINGS) == 105


def mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def div_exact(a, b):
    num = mul(a, (b[0], -b[1]))
    norm = b[0] * b[0] + b[1] * b[1]
    assert num[0] % norm == 0 and num[1] % norm == 0
    return (num[0] // norm, num[1] // norm)


def gauss_gcd(a, b):
    while b != (0, 0):
        num = mul(a, (b[0], -b[1]))
        norm = b[0] * b[0] + b[1] * b[1]
        q = ((num[0] + norm // 2) // norm, (num[1] + norm // 2) // norm)
        qb = mul(q, b)
        a, b = b, (a[0] - qb[0], a[1] - qb[1])
    return a


def circle_from_rows(rows):
    gauss_product = reduce(mul, rows, (1, 0))
    conjugate_product = (gauss_product[0], -gauss_product[1])
    return [mul(h, div_exact(conjugate_product, (h[0], -h[1]))) for h in rows]


def gaussian_v(a, pi):
    power = 0
    norm = pi[0] * pi[0] + pi[1] * pi[1]
    while True:
        num = mul(a, (pi[0], -pi[1]))
        if num[0] % norm or num[1] % norm:
            return power
        a = (num[0] // norm, num[1] // norm)
        power += 1


def pair_data(rows):
    result = {}
    for i, j in EDGES:
        a, b = mul(rows[j], (rows[i][0], -rows[i][1]))
        common = gcd(abs(a), abs(b))
        assert a and b
        result[i, j] = (a // common, b // common)
    return result


def v(n, p):
    n = abs(n)
    assert n
    result = 0
    while n % p == 0:
        result += 1
        n //= p
    return result


def check_pfaffian(data):
    lhs = 0
    common_lcm = 1
    for matching, sign in MATCHINGS:
        chosen = set(matching)
        lhs += sign * prod(data[e][1] if e in chosen else data[e][0] for e in EDGES)
        value = prod(abs(data[e][0]) for e in matching)
        common_lcm = common_lcm // gcd(common_lcm, value) * value
    rhs = prod(data[e][1] for e in EDGES)
    assert lhs == rhs
    real_product = prod(abs(data[e][0]) for e in EDGES)
    assert real_product % common_lcm == 0
    content = real_product // common_lcm
    assert rhs % content == 0
    return content


def check_occupancies(xs, data, p):
    assert all((x * x + 1) % p for x in xs)
    real_v = {edge: v(data[edge][0], p) for edge in EDGES}
    imag_v = {edge: v(data[edge][1], p) for edge in EDGES}
    defect = sum(real_v.values()) - max(sum(real_v[e] for e in m) for m, _ in MATCHINGS)
    expected = 0
    max_depth = max(tuple(real_v.values()) + tuple(imag_v.values()))
    for depth in range(1, max_depth + 1):
        modulus = p ** (depth + (p == 2))
        classes = Counter()
        for x in xs:
            inv = pow(x * x + 1, -1, modulus)
            classes[((x * x - 1) * inv % modulus, 2 * x * inv % modulus)] += 1
        seen = set()
        for c in classes:
            if c in seen:
                continue
            opposite = ((-c[0]) % modulus, (-c[1]) % modulus)
            seen.update((c, opposite))
            diff = abs(classes[c] - classes[opposite])
            expected += diff * (diff - 1) // 2
    assert sum(imag_v.values()) - defect == expected


def check_conductor_occupancies(xs, data, pi):
    p = pi[0] * pi[0] + pi[1] * pi[1]
    rows = [(x, 1) for x in xs]
    circle = circle_from_rows(rows)
    e = v(prod(x * x + 1 for x in xs), p)
    local = []
    for z in circle:
        a = gaussian_v(z, pi)
        u = z
        for _ in range(a):
            u = div_exact(u, pi)
        for _ in range(e - a):
            u = div_exact(u, (pi[0], -pi[1]))
        assert (u[0] * u[0] + u[1] * u[1]) % p
        local.append((a, u))
    real_v = {edge: v(data[edge][0], p) for edge in EDGES}
    imag_v = {edge: v(data[edge][1], p) for edge in EDGES}
    defect = sum(real_v.values()) - max(sum(real_v[x] for x in m) for m, _ in MATCHINGS)
    expected = 0
    for depth in range(1, max(tuple(real_v.values()) + tuple(imag_v.values())) + 1):
        modulus = p ** depth
        classes = Counter((a, u[0] % modulus, u[1] % modulus) for a, u in local)
        seen = set()
        for c in classes:
            if c in seen:
                continue
            opposite = (c[0], -c[1] % modulus, -c[2] % modulus)
            seen.update((c, opposite))
            diff = abs(classes[c] - classes[opposite])
            expected += diff * (diff - 1) // 2
    assert sum(imag_v.values()) - defect == expected
    return e > 0


def check_plus_content(xs, content):
    circle = circle_from_rows([(x, 1) for x in xs])
    pair_cores = {edge: gauss_gcd(circle[edge[0]], circle[edge[1]]) for edge in EDGES}
    pair_sums = {}
    for i, j in EDGES:
        numerator = (circle[i][0] + circle[j][0], circle[i][1] + circle[j][1])
        assert numerator[0] % 2 == numerator[1] % 2 == 0
        pair_sums[i, j] = (numerator[0] // 2, numerator[1] // 2)
    core_products = [reduce(mul, (pair_cores[e] for e in m), (1, 0)) for m, _ in MATCHINGS]
    sum_products = [reduce(mul, (pair_sums[e] for e in m), (1, 0)) for m, _ in MATCHINGS]
    gamma0 = reduce(gauss_gcd, core_products)
    gamma_plus = reduce(gauss_gcd, sum_products)
    delta_gauss = div_exact(gamma_plus, gamma0)
    assert delta_gauss[0] == 0 or delta_gauss[1] == 0
    delta = abs(delta_gauss[0]) + abs(delta_gauss[1])
    assert content % delta == 0
    primitive = [div_exact(value, gamma_plus) for value in sum_products]
    assert all(a == 0 or b == 0 for a, b in primitive)
    assert reduce(gcd, (abs(a) + abs(b) for a, b in primitive)) == 1
    assert len(set(primitive)) > 1


def main():
    cuts = [set(c) for c in combinations(range(8), 4) if 0 in c]
    balanced = [[int(all((i in cut) != (j in cut) for i, j in matching))
                 for matching, _ in MATCHINGS] for cut in cuts]
    assert rational_rank(balanced) == 35
    for i, cut in enumerate(cuts):
        for j, other in enumerate(cuts):
            intersection = len(cut & other)
            assert sum(a * b for a, b in zip(balanced[i], balanced[j])) == factorial(intersection) * factorial(4 - intersection)
    quadratic = [[int(edge not in matching) for matching, _ in MATCHINGS] for edge in EDGES]
    cubic = [[1 - sum(i in triple and j in triple for i, j in matching)
              for matching, _ in MATCHINGS] for triple in map(set, combinations(range(8), 3))]
    assert rational_rank(quadratic) == 21
    assert rational_rank([[1] * 105] + quadratic + cubic) == 21
    windows = 0
    local_checks = 0
    conductor_checks = 0
    for step in (2, 4, 6, 10):
        for start in (2, 8, 26):
            xs = [start + step * j for j in range(8)]
            data = pair_data([(x, 1) for x in xs])
            content = check_pfaffian(data)
            check_plus_content(xs, content)
            windows += 1
            for p in (2, 3, 7, 11, 13, 17, 19, 23, 29, 31):
                if all((x * x + 1) % p for x in xs):
                    check_occupancies(xs, data, p)
                    local_checks += 1
            for pi in ((2, 1), (3, 2), (4, 1), (5, 2), (6, 1)):
                conductor_checks += check_conductor_occupancies(xs, data, pi)

    xs = (26, 2, 4, 10, 52, 28, 30, 36)
    base_data = pair_data([(x, 1) for x in xs])
    cross_forms = {}
    for i, j in EDGES:
        if i < 4 <= j:
            coefficients = base_data[i, j]
            if coefficients[0] < 0:
                coefficients = tuple(-a for a in coefficients)
            cross_forms.setdefault(coefficients, []).append((i, j))
    for group in cross_forms.values():
        assert len(set(sum((list(edge) for edge in group), []))) == 2 * len(group)
    determinants = [abs(a[0] * b[1] - a[1] * b[0]) for a, b in combinations(cross_forms, 2)]
    assert all(determinants)
    fixed_real = prod(abs(base_data[i, j][0]) for i, j in EDGES if (i < 4) == (j < 4))
    resultant_bound = fixed_real * prod(determinants) ** 16
    pi_power = (1, 0)
    contents = []
    for exponent in range(1, 21):
        pi_power = mul(pi_power, (3, 2))
        assert gcd(abs(pi_power[0]), abs(pi_power[1])) == 1
        rows = [mul(pi_power, (x, 1)) if j < 4 else (x, 1) for j, x in enumerate(xs)]
        data = pair_data(rows)
        assert all(a % 13 and b % 13 for a, b in data.values())
        for i, j in EDGES:
            if i < 4 <= j:
                c, d = base_data[i, j]
                assert data[i, j][0] == c * pi_power[0] + d * pi_power[1]
        contents.append(check_pfaffian(data))
        assert contents[-1] <= resultant_bound
    print("PASS:", windows, "eight-row Pfaffian/LCM and pair-sum gcd windows;", local_checks,
          "outside-prime occupancy checks;", conductor_checks,
          "conductor occupancy checks; 20 balanced-prime exponents.")
    print("Observed B_X range in the balanced-prime family:", min(contents), max(contents))
    print("PASS: balanced Gram rank35, Taylor ranks1/21/35, and fixed-resultant bound.")


if __name__ == "__main__":
    main()
