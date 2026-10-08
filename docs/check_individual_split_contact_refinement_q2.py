"""Exact arithmetic certificates supporting the restricted contact lemma."""


def add(z, w):
    return (z[0] + w[0], z[1] + w[1])


def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def divides(z, w):
    norm = z[0] ** 2 + z[1] ** 2
    if norm == 0:
        return False
    numerator = mul(w, (z[0], -z[1]))
    return numerator[0] % norm == numerator[1] % norm == 0


def sharp_polynomial(z):
    return add(add(mul(z, z), z), (1, 0))


def main():
    assert sharp_polynomial((0, 1)) == (0, 1)
    assert sharp_polynomial((-1, -1)) == (0, 1)

    # Every nonzero Gaussian divisor of 2i has norm dividing 4, so
    # both coordinates have absolute value at most 2. This is an
    # exhaustive finite certificate, not a bounded polynomial search.
    divisors = {
        (x, y)
        for x in range(-2, 3)
        for y in range(-2, 3)
        if divides((x, y), (0, 2))
    }
    assert all(4 % (x * x + y * y) == 0 for x, y in divisors)
    assert {(x, y) for x, y in divisors if abs(y) == 2} == {(0, -2), (0, 2)}
    assert {y for y in (-2, -1, 1, 2) if divides((0, 2 * y), (0, 2))} == {-1, 1}
    print("PASS: the two-point bound is sharp for F=t^2+t+1.")
    print("PASS: exhaustive Gaussian-divisor certificate for the same-sign exclusion.")
    print("The general result and degree-eight residual-fiber statement are proved in the note.")


if __name__ == "__main__":
    main()
