"""Exact six-pair Gaussian-norm cross-ratio countermodel."""

from fractions import Fraction


def product_quadratics(a, b):
    """Ascending coefficients of (t^2+a)(t^2+b)."""
    return a * b, 0, a + b, 0, 1


def main():
    b = {'12': 1, '34': 3, '13': 6, '24': 7, '14': 5, '23': 15}
    assert len({abs(x) for x in b.values()}) == 6
    assert all(x != 0 for x in b.values())
    n = {key: value * value for key, value in b.items()}
    terms = [product_quadratics(n['12'], n['34']),
             product_quadratics(n['13'], n['24']),
             product_quadratics(n['14'], n['23'])]
    assert all(11 * x - 16 * y + 5 * z == 0 for x, y, z in zip(*terms))
    a = [Fraction(1), Fraction(38, 27), Fraction(3, 2), Fraction(2)]
    mu = lambda i, j: a[j - 1] - a[i - 1]
    assert [mu(1, 2) * mu(3, 4), mu(1, 3) * mu(2, 4),
            mu(1, 4) * mu(2, 3)] == [Fraction(11, 54),
                                     Fraction(16, 54), Fraction(5, 54)]
    print('Six distinct Gaussian-linear norm blocks satisfy the pair cross ratio.')


if __name__ == '__main__':
    main()
