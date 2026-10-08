"""Exact modular certificates for irreducibility and J(Q)[2]=0.

UV_CACHE_DIR=/tmp/uv-cache uv run --offline --with sympy python \
    docs/check_six_point_branch_two_torsion.py

The finite-field factors are independently checked by the Frobenius/gcd
irreducibility criterion. No Galois group identification is assumed.
"""
import sympy as sp
from check_six_point_isotropic_circle_cover import probe


def modular_power(base, exponent, modulus):
    result = sp.Poly(1, base.gens[0], modulus=base.get_modulus())
    while exponent:
        if exponent & 1:
            result = (result * base).rem(modulus)
        base = (base * base).rem(modulus)
        exponent //= 2
    return result


def check_irreducible(f):
    p, n = f.get_modulus(), f.degree()
    x = sp.Poly(f.gens[0], f.gens[0], modulus=p)
    assert (modular_power(x, p**n, f) - x).rem(f).is_zero
    for q in sp.factorint(n):
        assert sp.gcd(f, modular_power(x, p**(n // q), f) - x).degree() == 0


def certificate(f, p, expected):
    fp = sp.Poly(f, f.gens[0], modulus=p)
    assert fp.degree() == 12
    assert sp.gcd(fp, fp.diff()).degree() == 0
    unit, factors = fp.factor_list()
    product = sp.Poly(unit, f.gens[0], modulus=p)
    degrees = []
    for factor, exponent in factors:
        assert exponent == 1
        check_irreducible(factor)
        product *= factor
        degrees.append(factor.degree())
    assert product == fp
    assert sorted(degrees) == list(expected)
    return sorted(degrees)


def subset_sums(degrees):
    sums = {0}
    for d in degrees:
        sums |= {s + d for s in sums}
    return sums


def fixed_even_partition_classes(degrees):
    """Enumerate invariants of one cycle permutation on even subsets / complement."""
    permutation, start = [], 0
    for length in degrees:
        permutation.extend(list(range(start + 1, start + length)) + [start])
        start += length
    assert start == 12
    full = (1 << 12) - 1
    fixed = set()
    for mask in range(1 << 12):
        if mask.bit_count() % 2:
            continue
        moved = sum(1 << permutation[i] for i in range(12) if mask & (1 << i))
        if moved in (mask, full ^ mask):
            fixed.add(min(mask, full ^ mask))
    return fixed


def main():
    x = sp.Symbol('t')
    cases = [
        ((0, 1, 2, 3, 4, 5), [(59, (1, 11)), (73, (5, 7))]),
        ((0, 1, 2, 4, 7, 11), [(43, (3, 9)), (29, (4, 4, 4))]),
    ]
    for params, certificates in cases:
        weights, _, _, _, coefficients = probe(params)
        f = sp.Poly.from_dict(
            {(i,): sp.Rational(c.numerator, c.denominator)
             for i, c in enumerate(coefficients)}, x, domain=sp.QQ
        ).clear_denoms()[1].primitive()[1]
        if params == (0, 1, 2, 3, 4, 5):
            # The S12 witness's sign field is not Q(i): -disc(f) is
            # a nonsquare at the good prime 59. This also rules out
            # zero branch translation kernel on a pointed two-cover.
            discriminant_mod_59 = int(sp.Poly(f, x, modulus=59).discriminant()) % 59
            assert discriminant_mod_59 == 22
            assert pow(discriminant_mod_59, 29, 59) == 1
            assert pow((-discriminant_mod_59) % 59, 29, 59) == 58
            print('PASS: first witness disc=22 mod 59; -disc is nonsquare.')
        common_degrees = set(range(13))
        for p, expected in certificates:
            degrees = certificate(f, p, expected)
            common_degrees &= subset_sums(degrees)
            if p == certificates[0][0]:
                assert fixed_even_partition_classes(degrees) == {0}
            print(f'PASS: p={p}, squarefree factor degrees={degrees}')
        assert common_degrees == {0, 12}
        print(f'PASS: weights={weights}: irreducible over Q by modular subset sums; '
              'J(Q)[2]=0 by a single Frobenius permutation.')
    # Regression: irreducibility alone must not imply zero invariants.
    # A transitive 12-cycle fixes a nonzero unordered alternating 6+6 partition.
    assert len(fixed_even_partition_classes([12])) > 1
    print('PASS: transitive 12-cycle retains an unordered 6+6 class; '
          'the irreducibility-only shortcut is rejected.')


if __name__ == '__main__':
    main()
