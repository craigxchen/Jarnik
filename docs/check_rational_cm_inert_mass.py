"""Exact finite checks of the rational CM Boolean denominator decomposition.

Uses only integer and Fraction arithmetic.  It does not numerically test or
claim to certify the compactness/diagonal subsequence theorem.
"""

from fractions import Fraction
from math import prod

from check_elliptic_split_inert_norm_route import (
    C, P, add, conj, denominator_generator, divide, exact, ggcd,
    multiple, norm, times,
)


def integral(z):
    return all(Fraction(a).denominator == 1 for a in z)


def trim(z, primes):
    for pi in primes:
        while integral(divide(z, pi)):
            z = exact(z, pi)
    return z


def subsets(mask):
    s = mask
    while True:
        yield s
        if s == 0:
            return
        s = (s - 1) & mask


def conjugate_mask(mask, pairs):
    return sum(((mask >> (j ^ 1)) & 1) << j for j in range(2 * pairs))


def cm_multiple(alpha, r):
    ip = (times(C(-1), r[0]), times(C(0, 1), r[1]))
    return add(multiple(alpha[0], r), multiple(alpha[1], ip))


def check_cells(split_generators, t):
    primes = [pi for generator in split_generators
              for pi in (generator, conj(generator))]
    m = prod(int(norm(pi)) for pi in split_generators)
    masks = range(1 << len(primes))
    full = len(masks) - 1
    r = multiple(t, P)
    denominators = {}
    cells = {}
    for mask in masks:
        alpha = C(1)
        for j, pi in enumerate(primes):
            if mask >> j & 1:
                alpha = times(alpha, pi)
        alpha = tuple(int(a) for a in alpha)
        point = cm_multiple(alpha, r)
        denominators[mask] = trim(denominator_generator(point), [(1, 1)] + primes)
        numerator, denominator = C(1), C(1)
        for s in subsets(mask):
            if (bin(mask).count("1") - bin(s).count("1")) % 2:
                denominator = times(denominator, denominators[s])
            else:
                numerator = times(numerator, denominators[s])
        cells[mask] = exact(numerator, denominator)

    for mask in masks:
        rebuilt = C(1)
        for s in subsets(mask):
            rebuilt = times(rebuilt, cells[s])
        assert norm(exact(rebuilt, denominators[mask])) == 1
        other = conjugate_mask(mask, len(split_generators))
        assert norm(exact(conj(cells[mask]), cells[other])) == 1
        for previous in range(mask):
            assert norm(ggcd(cells[mask], cells[previous])) == 1
        if other != mask:
            # An inert prime would divide the conjugate gcd.  This exact
            # gcd test proves its absence without integer factorization.
            assert norm(ggcd(cells[mask], conj(cells[mask]))) == 1

    rational_full = trim(denominator_generator(multiple(m * t)), [(1, 1)] + primes)
    assert norm(exact(denominators[full], rational_full)) == 1

    weights = {
        mask: prod(int(norm(pi)) - 1 for j, pi in enumerate(primes) if mask >> j & 1)
        for mask in masks
    }
    stable_weight = sum(w for mask, w in weights.items()
                        if conjugate_mask(mask, len(split_generators)) == mask)
    assert sum(weights.values()) == m * m
    assert stable_weight == prod(1 + (int(norm(pi)) - 1) ** 2 for pi in split_generators)
    rho = prod(Fraction(q*q - 2*q + 2, q*q)
               for q in (int(norm(pi)) for pi in split_generators))
    assert Fraction(stable_weight, m*m) == rho
    largest_weight = prod(int(norm(pi)) - 1 for pi in primes)
    assert max(weights.values()) == largest_weight
    assert Fraction(largest_weight, m*m) <= rho
    # Whole nonstable conjugate cells are mutually coprime rational atoms.
    # Their norms multiply to the retained rational denominator, exactly.
    chosen = [mask for mask in masks
              if mask < conjugate_mask(mask, len(split_generators))]
    retained = C(1)
    for mask in masks:
        if mask != conjugate_mask(mask, len(split_generators)):
            retained = times(retained, cells[mask])
    retained_integer = prod(int(norm(cells[mask])) for mask in chosen)
    assert norm(exact(retained, C(retained_integer))) == 1
    for mask in chosen:
        for other in chosen:
            assert norm(ggcd(cells[mask], conj(cells[other]))) == 1
    print(f"m={m}, t={t}: {len(masks)} integral coprime cells; reconstruction, conjugation, rho={rho} pass.")


def main():
    for t in range(1, 7):
        check_cells([(2, 1)], t)
    for t in (1, 2):
        check_cells([(3, 2)], t)
    check_cells([(2, 1), (3, 2)], 1)
    print("All exact finite checks pass; the subsequence theorem uses the written proof.")


if __name__ == "__main__":
    main()
