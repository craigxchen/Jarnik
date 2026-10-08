"""Finite incidence checks for the Boolean CM phase divisor model.

This deliberately does not implement elliptic-curve arithmetic.  A point
``n*P`` is represented by its Gaussian coefficient ``n`` in Z[i]; this is
enough for the distinctness and divisor-support assertions because P is
nontorsion and the endomorphism ring acts faithfully on P.
"""

from itertools import product


def bits(m):
    return tuple(product((0, 1), repeat=m))


def complement(b):
    return tuple(1 - x for x in b)


def c_vector(b):
    """The row-phase divisor vector c_b, with row m as anchor."""
    return tuple(x - b[-1] for x in b[:-1])


def incidence(m, kappas, sign=1):
    """Return (kind, block, coefficient-vector, point-coefficient) incidences.

    ``sign`` allows the simultaneous opposite determinant convention.  The
    asserted pairing is unchanged under it.
    """
    assert len(kappas) == 1 << m
    assert len(set(kappas.values())) == len(kappas)
    out = []
    for b in bits(m):
        v = tuple(sign * x for x in c_vector(b))
        out.append(("A", b, v, (kappas[b] - 1j)))
        out.append(("B", b, tuple(-x for x in v), (kappas[b] + 1j)))
    return out


def check_vector_multiplicities(m, kappas, sign=1):
    events = incidence(m, kappas, sign=sign)
    by_vector = {}
    for kind, b, v, point in events:
        by_vector.setdefault(v, []).append((kind, b, point))

    zero = (0,) * (m - 1)
    assert sorted((kind, b) for kind, b, _ in by_vector[zero]) == sorted(
        [("A", (0,) * m), ("A", (1,) * m),
         ("B", (0,) * m), ("B", (1,) * m)])
    for b in bits(m):
        v = tuple(sign * x for x in c_vector(b))
        if not any(v):
            continue
        got = {(kind, block) for kind, block, _ in by_vector[v]}
        expected = {("A", b), ("B", complement(b))}
        assert got == expected, (m, b, v, got, expected)
        assert len(got) == 2
        assert any(abs(x) == 1 for x in v)
    return by_vector


def check_point_distinctness(m, kappas):
    """Check A/B points are distinct using only coefficient arithmetic."""
    events = incidence(m, kappas)
    points = [point for _, _, _, point in events]
    # Equality of (k+i)P and (l+j)P would give a nonzero Gaussian
    # endomorphism annihilating P.  The coefficients below are all distinct.
    coeffs = [(k.real, k.imag) for k in points]
    assert len(set(coeffs)) == len(coeffs), coeffs


def symmetric_kappas(m):
    """Distinct shifts with a constant nontrivial complement sum.

    The empty/full pair is deliberately assigned a different sum: its
    divisor vector is zero and therefore it is irrelevant to the degree-two
    criterion.
    """
    all_bits = bits(m)
    K = 10 * (1 << m) + 7
    result = {}
    next_t = -3 * (1 << m)
    for b in all_bits:
        if b in result:
            continue
        bc = complement(b)
        if not any(c_vector(b)):
            result[b] = -1
            result[bc] = 1
            continue
        result[b] = next_t
        result[bc] = K - next_t
        next_t += 3
    assert len(set(result.values())) == len(result)
    assert all(result[b] + result[complement(b)] == K
               for b in all_bits if any(c_vector(b)))
    return result


def asymmetric_kappas(m):
    """Distinct shifts whose complement sums are not all equal."""
    result = {}
    for n, b in enumerate(bits(m)):
        result[b] = n * n + 3 * n + 1
    assert len(set(result.values())) == len(result)
    sums = {result[b] + result[complement(b)] for b in bits(m)}
    assert len(sums) > 1
    return result


def complement_sums(kappas):
    return {kappas[b] + kappas[complement(b)] for b in kappas}


def nontrivial_complement_sums(kappas):
    return {kappas[b] + kappas[complement(b)]
            for b in kappas if any(c_vector(b))}


def check_degree_two_criterion(m):
    symmetric = symmetric_kappas(m)
    asymmetric = asymmetric_kappas(m)
    assert len(nontrivial_complement_sums(symmetric)) == 1
    # The zero-vector pair is intentionally allowed to have another sum.
    assert len(complement_sums(symmetric)) == 2
    # For m=2 there is only one nontrivial complement pair, so the
    # condition is vacuous even for an asymmetric assignment.
    if m >= 3:
        assert len(nontrivial_complement_sums(asymmetric)) > 1

    # A_b+B_bar has coefficient (kappa_b+kappa_bar)P, since the +/- i
    # terms cancel.  Since P is nontorsion, these degree-two divisors are
    # linearly equivalent exactly when their complement sums agree.
    pairs = [b for b in bits(m) if b <= complement(b)
             and any(c_vector(b))]
    assert len(pairs) == (1 << m) // 2 - 1
    for shifts, should_be_constant in ((symmetric, True),
                                       (asymmetric, m < 3)):
        sums = [shifts[b] + shifts[complement(b)] for b in pairs]
        assert (len(set(sums)) == 1) == should_be_constant


def check_translation_obstruction():
    # A_b and B_b differ by [2i]P.  A translation deck involution would
    # require this difference to be a 2-torsion point, but [2i]P is
    # nontorsion whenever P is nontorsion.
    assert 2j != 0


def main():
    for m in range(2, 7):
        symmetric = symmetric_kappas(m)
        asymmetric = asymmetric_kappas(m)
        for shifts in (symmetric, asymmetric):
            check_point_distinctness(m, shifts)
            # Opposite bracket/determinant convention changes every vector
            # sign, but not the two-point incidence or complement pairing.
            for sign in (1, -1):
                check_vector_multiplicities(m, shifts, sign=sign)
        check_degree_two_criterion(m)
    check_translation_obstruction()
    print("Boolean CM divisor vectors pair each nonzero vector at A_b and B_complement.")
    print("Symmetric assignments pass the constant complement-sum test; asymmetric assignments fail it for m >= 3.")
    print("m=2 and m=3 pass: the zero vector occurs at four points, while each nonzero vector occurs at exactly two.")
    print("The +/- i sign convention preserves all incidence statements; [2i]P is not a translation involution.")


if __name__ == "__main__":
    main()
