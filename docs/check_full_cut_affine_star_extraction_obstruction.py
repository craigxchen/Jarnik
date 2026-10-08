"""Exact checks for the full-cut affine star extraction obstruction."""

from itertools import combinations, product


def regular(v):
    k = len(v) // 2
    return sorted(v[:k]) == sorted(v[k:])


def star(v):
    k = len(v) // 2
    return len(set(v[:k])) == len(set(v[k:])) == 1


def allowed(v):
    return regular(v) or star(v)


def thresholds_allowed(v):
    k = len(v) // 2
    for level in sorted(set(v))[1:]:
        a = sum(x >= level for x in v[:k])
        b = sum(x >= level for x in v[k:])
        if a != b and (a, b) not in ((k, 0), (0, k)):
            return False
    return True


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def power(z, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, z)
    return out


def prod(xs):
    out = 1
    for x in xs:
        out *= x
    return out


def main():
    threshold_checks = shift_checks = 0
    for k in (2, 3, 4):
        for v in product(range(-1, 2), repeat=2 * k):
            assert allowed(v) == thresholds_allowed(v)
            threshold_checks += 1
            if regular(v) and not star(v):
                for alpha, beta in ((-1, 0), (0, 1), (-1, 1)):
                    s = (alpha,) * k + (beta,) * k
                    for sign in (-1, 1):
                        shifted = tuple(x + sign * y for x, y in zip(s, v))
                        assert not allowed(shifted)
                        shift_checks += 1

    rows = [(-2, 1, 2, 0), (-1, 2, 0, 0), (0, 0, 1, 0),
            (-1, 0, 2, 0), (0, 1, 0, 0), (-2, 2, 1, 0)]
    assert len(set(rows)) == 6 and all(sum(row) == 1 for row in rows)
    for mask in range(16):
        values = tuple(sum(row[j] for j in range(4) if mask >> j & 1)
                       for row in rows)
        assert regular(values)
    cuts = [mask for mask in range(1, 15) if mask & 1]
    primes = [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1), (5, 4), (7, 2)]
    norms = [a * a + b * b for a, b in primes]
    assert norms == [5, 13, 17, 29, 37, 41, 53]
    cols, widths = [], []
    for mask in cuts:
        values = [sum(row[j] for j in range(4) if mask >> j & 1)
                  for row in rows]
        low, high = min(values), max(values)
        assert regular(values)
        cols.append([value - low for value in values])
        widths.append(high - low)
    norm = prod(p ** e for p, e in zip(norms, widths))
    points = []
    for i in range(6):
        z = (1, 0)
        for gaussian, col, e in zip(primes, cols, widths):
            z = mul(z, power(gaussian, col[i]))
            z = mul(z, power((gaussian[0], -gaussian[1]), e - col[i]))
        points.append(z)
        assert z[0] ** 2 + z[1] ** 2 == norm
    assert len(set(points)) == 6
    assert all(min(col) == 0 and max(col) == e
               for col, e in zip(cols, widths))  # Complete Gaussian gcd is one.
    cross_norms = []
    for i in range(3):
        for j in range(3, 6):
            cross_norms.append(prod(p ** abs(col[i] - col[j])
                                    for p, col in zip(norms, cols)))
    assert prod(cross_norms) ** 2 <= norm ** 9
    assert any(p ** 2 <= norm for p in cross_norms)
    print(f"PASS: {threshold_checks} threshold profiles, {shift_checks} shifted "
          "multisets; all 16 subset sums and six distinct primitive Gaussian "
          "outputs with exact cross-pair norm budget.")


if __name__ == "__main__":
    main()
