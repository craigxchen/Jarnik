"""Exact specialization checks for three profiles of different block degrees.

The general theorem is proved in the accompanying note. These checks
exercise its rational base point, primitive row scalars, whole-progression
prime certificate, and actual circle construction on the known three-row
family composed with t, t^2, and t^3.
"""

from math import gcd, prod
from itertools import combinations

from check_balanced_three_row_gaussian_polynomial_family import (
    FACTORS, INCIDENT, conjugate, determinant, exact_divide, factor,
    gaussian_product, mul, norm,
)


def resultant_monic_quadratics(a, b):
    """Norm(w+a) against Norm(w+b), for Gaussian integer constants."""
    linear_a, constant_a = 2 * a[0], norm(a)
    linear_b, constant_b = 2 * b[0], norm(b)
    delta_linear = linear_b - linear_a
    delta_constant = constant_b - constant_a
    return (constant_a * delta_linear**2
            - linear_a * delta_linear * delta_constant + delta_constant**2)


def primitive(z):
    return gcd(*z) == 1 and norm(z) % 2 == 1


def check_degree(e):
    # H_T((1/2+u)^e), scaled by 2^e, equals A_T+v(u), where
    # v(u)=(1+2u)^e-1.  All A_T have odd real and even imaginary part.
    scale = 2**e
    base = {key: (1 + scale * a, scale * b)
            for key, (a, b) in FACTORS.items()}
    constants = {key: norm(z) for key, z in base.items()}
    contents, corrections = [], []
    for row in INCIDENT:
        z = gaussian_product(base[key] for key in row)
        content = gcd(*z)
        correction = tuple(x // content for x in z)
        assert primitive(correction)
        contents.append(content)
        corrections.append(correction)

    certificates = [2]
    certificates += [abs(z[1]) for z in base.values()]
    certificates += list(constants.values())
    certificates += [norm(z) for z in corrections]
    certificates += [abs(resultant_monic_quadratics(base[a], base[b]))
                     for a, b in combinations(base, 2)]
    assert all(certificates)
    bad = set().union(*(set(factor(n)) for n in certificates))
    constant_factors = [factor(n) for n in constants.values()]
    modulus = prod(p**(1 + max(f.get(p, 0) for f in constant_factors))
                   for p in bad)
    assert all(modulus % c == 0 for c in constants.values())
    # If u=Mh then p^(1+max v_p(c_T)) divides v(u), proving
    # v_p(Norm(A_T+v(u)))=v_p(c_T) for the entire progression.
    assert all(modulus % (p**(1 + max(f.get(p, 0)
                                    for f in constant_factors))) == 0
               for p in bad)

    residuals = {(i, j): determinant(corrections[i], corrections[j])
                 for i, j in combinations(range(3), 2)}
    assert all(residuals.values())
    correction_product = gaussian_product(conjugate(z) for z in corrections)
    circle_corrections = [exact_divide(mul(correction_product, z), conjugate(z))
                          for z in corrections]

    for h in range(1, 7):
        value = (1 + 2 * modulus * h)**e - 1
        blocks = {key: exact_divide((z[0] + value, z[1]), z)
                  for key, z in base.items()}
        norms = {key: norm(z) for key, z in blocks.items()}
        assert all(primitive(z) for z in blocks.values())
        assert all(n > 1 for n in norms.values())
        assert all(gcd(a, b) == 1 for a, b in combinations(norms.values(), 2))
        assert all(n % p for n in norms.values() for p in bad)
        rows = []
        for i, incident in enumerate(INCIDENT):
            row = mul(corrections[i], gaussian_product(blocks[key] for key in incident))
            raw = gaussian_product((base[key][0] + value, base[key][1])
                                   for key in incident)
            assert raw == tuple(contents[i] * x for x in row)
            assert primitive(row)
            assert row[1] == corrections[i][1] != 0
            rows.append(row)
        for (i, j), residual in residuals.items():
            shared = set(INCIDENT[i]) & set(INCIDENT[j])
            assert determinant(rows[i], rows[j]) == residual * prod(norms[k] for k in shared)

        b = gaussian_product(blocks.values())
        anchor = mul(conjugate(b), correction_product)
        points = [anchor]
        for i, incident in enumerate(INCIDENT):
            oriented = gaussian_product(z if key in incident else conjugate(z)
                                        for key, z in blocks.items())
            point = mul(oriented, circle_corrections[i])
            assert norm(point) == norm(anchor)
            assert mul(point, conjugate(rows[i])) == mul(anchor, rows[i])
            points.append(point)
        assert len(set(points)) == 4
    return len(bad)


def main():
    sizes = [check_degree(e) for e in (1, 2, 3)]
    print('PASS: three block degrees, whole-progression certificates, '
          '18 primitive profiles and equal-radius circle realizations; '
          f'bad-prime set sizes {sizes}.')


if __name__ == '__main__':
    main()
