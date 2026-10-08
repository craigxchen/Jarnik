"""Exact whole-block ray and reflection-orbit checks for the Pell family."""

from itertools import product


# (a,b,c,d) represents a+b*sqrt(5)+i(c+d*sqrt(5)).
def real_mul(x, y):
    a, b = x
    c, d = y
    return a * c + 5 * b * d, a * d + b * c


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def negate(x):
    return tuple(-a for a in x)


def multiply(x, y):
    xr, xi = x[:2], x[2:]
    yr, yi = y[:2], y[2:]
    real = add(real_mul(xr, yr), negate(real_mul(xi, yi)))
    imag = add(real_mul(xr, yi), real_mul(xi, yr))
    return real + imag


def conjugate(x):
    return x[0], x[1], -x[2], -x[3]


def unit_times(unit, x):
    if unit == 1:
        return x
    if unit == -1:
        return negate(x)
    if unit == 1j:
        return -x[2], -x[3], x[0], x[1]
    assert unit == -1j
    return x[2], x[3], -x[0], -x[1]


def block_product(blocks, bits):
    answer = (1, 0, 0, 0)
    for block, bit in zip(blocks, bits):
        answer = multiply(answer, conjugate(block) if bit else block)
    return answer


def main():
    # Exact leading coefficients of P,A/V,B/V,C/V.
    blocks = [(-1, 0, 2, 0), (1, 1, 2, 0),
              (2, 0, -1, 1), (8, 2, 1, 3)]
    base = block_product(blocks, (0, 0, 0, 0))
    units = (1, -1, 1j, -1j)
    aligned = []
    for bits in product((0, 1), repeat=4):
        value = block_product(blocks, bits)
        matches = [unit for unit in units if unit_times(unit, value) == base]
        if matches:
            aligned.append(bits)
            assert len(matches) == 1
    expected = {(0, 0, 0, 0), (1, 0, 1, 1),
                (1, 1, 0, 1), (1, 1, 1, 0)}
    assert set(aligned) == expected

    # Binary integral points in the full affine reflection orbit.
    integral_patterns = set()
    for t_a, t_b, t_c in product((0, 1), repeat=3):
        numerators = (t_b + t_c - t_a,
                      t_a + t_c - t_b,
                      t_a + t_b - t_c)
        if any(x % 2 for x in numerators):
            continue
        c1, c2, c3 = (x // 2 for x in numerators)
        t_p = c1 + c2 + c3
        if t_p in (0, 1):
            integral_patterns.add((t_p, t_a, t_b, t_c))
    assert integral_patterns == expected
    print("PASS: exactly four Gaussian-unit-aligned whole-block rays")
    print("PASS: the full rational reflection orbit has exactly the original four integral points")


if __name__ == "__main__":
    main()
