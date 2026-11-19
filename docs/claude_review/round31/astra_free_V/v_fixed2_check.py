"""Exact F4[t] check of the fixed-characteristic hyperplane construction.

A pair of bit masks (a,b) represents a(t)+alpha*b(t), alpha^2=alpha+1.
No third-party packages are needed.  Run with python3 v_fixed2_check.py.
"""


def carryless(x, y):
    if x.bit_count() < y.bit_count():
        x, y = y, x
    result = 0
    while y:
        lowest = y & -y
        result ^= x << (lowest.bit_length() - 1)
        y ^= lowest
    return result


def product(x, y):
    a, b = x
    c, d = y
    ac = carryless(a, c)
    bd = carryless(b, d)
    cross = carryless(a ^ b, c ^ d)
    return ac ^ bd, ac ^ cross


def multiply_factors(indices, power):
    result = (1, 0)
    for v in indices:
        result = product(result, (power ^ v, 1))
    return result


def degree(x):
    return (x[0] | x[1]).bit_length() - 1


def check(d):
    q = 1 << d
    marks = list(range(1, q, 2))
    kernels = [
        {v for v in range(q) if not ((v & mark).bit_count() & 1)}
        for mark in marks
    ]
    # At each pair, the coefficient is the product over 1+(H intersect H').
    period_degrees = [
        sum((v ^ 1).bit_length() - 1 for v in (h & k))
        for i, h in enumerate(kernels) for k in kernels[:i]
    ]
    exponent = max(2 * max(period_degrees), d)
    power = 1 << exponent
    points = [multiply_factors(sorted(h - {0}), power) for h in kernels]
    parent = multiply_factors(range(2, q), power)
    assert parent[1] == 0
    point_degree = (q // 2 - 1) * exponent
    assert degree(parent) == 2 * point_degree
    assert len(set(points)) == len(marks)
    for pt in points:
        assert degree(pt) == point_degree
        a, b = pt
        assert product(pt, (a ^ b, b)) == parent
    degrees = [
        degree((x[0] ^ y[0], x[1] ^ y[1]))
        for i, x in enumerate(points) for y in points[:i]
    ]
    assert all(deg == (q // 4 - 1) * exponent + period
               for deg, period in zip(degrees, period_degrees))
    assert 2 * max(degrees) <= point_degree
    print(d, q, len(points), exponent, point_degree, max(degrees),
          len(degrees), sep="\t")


def remainder(x, divisor):
    divisor_degree = degree(divisor)
    lead_a = (divisor[0] >> divisor_degree) & 1
    lead_b = (divisor[1] >> divisor_degree) & 1
    inverse_lead = (lead_a ^ lead_b, lead_b)  # F4 conjugate
    while degree(x) >= divisor_degree:
        top = degree(x)
        shift = top - divisor_degree
        lead = ((x[0] >> top) & 1, (x[1] >> top) & 1)
        coefficient = product(lead, inverse_lead)
        multiple = product(coefficient, divisor)
        x = (x[0] ^ (multiple[0] << shift), x[1] ^ (multiple[1] << shift))
    return x


def gcd(x, y):
    while y != (0, 0):
        x, y = y, remainder(x, y)
    return x  # Only its degree will be used.


def check_residual_gcd():
    # Translation by an element of V permutes the sides and their intersection.
    # These four representatives cover all T of degree < 6, including zero.
    d = 4
    q = 1 << d
    forms = range(1, q, 2)
    full = set(range(q))
    peak = 0
    cases = 0
    for evaluation in (0, 16, 32, 48):
        factor_degrees = [degree((evaluation ^ v, 1)) for v in range(q)]
        sides = []
        for form in forms:
            kernel = {v for v in full if not ((v & form).bit_count() & 1)}
            sides.append([(side, multiply_factors(sorted(side), evaluation))
                          for side in (kernel, full - kernel)])

        def visit(index, count, common_indices, common_gcd):
            nonlocal cases, peak
            if index == len(sides):
                if count >= 2:
                    cases += 1
                    residual_degree = degree(common_gcd) - sum(
                        factor_degrees[v] for v in common_indices
                    )
                    assert 0 <= residual_degree <= q*d//4
                    peak = max(peak, residual_degree)
                return
            visit(index + 1, count, common_indices, common_gcd)
            for subset, raw in sides[index]:
                visit(index + 1, count + 1, common_indices & subset,
                      gcd(raw, common_gcd))

        visit(0, 0, full, (0, 0))
    print(f"d=4 residual gcd: {cases} cases up to translation, max degree {peak}")


if __name__ == "__main__":
    print("d\tq\tmarks\tM\tdeg_point\tmax_deg_difference\tpairs")
    for dimension in range(3, 7):
        check(dimension)
    check_residual_gcd()
