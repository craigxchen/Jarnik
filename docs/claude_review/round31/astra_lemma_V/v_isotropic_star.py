"""Exact integer specialization of the four-block isotropic star in Part V."""

from fractions import Fraction
from itertools import combinations
from math import gcd, prod

from sympy import isprime
from sympy.ntheory.modular import crt


def main():
    roots = (1, -9, -7, 2, 4)
    marks = [r * r for r in roots]
    edges = list(combinations(range(len(marks)), 2))
    scale = 6_350_400
    parameter = 3_175_200_023
    t = Fraction(parameter, scale)

    # The last four roots have elementary symmetric functions -10,-25,250,504.
    d_real = t**4 + 25 * t**2 + 504
    d_imag = 10 * t**3 + 250 * t
    columns = []
    for u in marks:
        multiplier = Fraction(u - 1, u) * (u + t**2)
        real = scale**4 * (d_real + multiplier * t / 10)
        imag = scale**4 * (d_imag + multiplier)
        assert real.denominator == imag.denominator == 1
        column = (int(real), int(imag))
        assert gcd(*column) == 1
        columns.append(column)

    blocks = [parameter**2 + (scale * roots[i] * roots[j])**2 for i, j in edges]
    # All primes lie below 2^64, the deterministic range for SymPy isprime.
    factors = [
        (10085161720052160529,),
        (37, 1297, 210128805298861),
        (14321, 33221957, 21190957),
        (17, 61, 1249, 4091873, 1902421),
        (3821, 36793, 617509, 117977),
        (29, 353, 20752289, 47518853),
        (13, 113, 6898679189888741),
        (10089799391770560529,),
        (21728981, 465438853709),
        (10084476151189440529,),
    ]
    all_primes = [p for entry in factors for p in entry]
    assert len(set(all_primes)) == len(all_primes)
    assert all(p < 2**64 and p % 4 == 1 and isprime(p) for p in all_primes)
    assert all(prod(entry) == block for entry, block in zip(factors, blocks))

    norms = [x**2 + y**2 for x, y in columns]
    for (i, j), block in zip(edges, blocks):
        x0, y0 = columns[i]
        x1, y1 = columns[j]
        det = x0 * y1 - y0 * x1
        cofactor = 504 * scale**6 * (Fraction(1, marks[i]) - Fraction(1, marks[j]))
        assert cofactor.denominator == 1 and cofactor
        assert det == cofactor * block
        for endpoint in (i, j):
            assert (norms[endpoint] % block == 0) == (i == 0)

    parent = prod(blocks)
    print(f"{len(all_primes)} distinct split primes, {len(str(parent))}-digit parent")
    print("Ten signed determinants verified. Isotropy holds on exactly four blocks.")

    # Send the marked line at each p to the line of (K,L*r_i*r_j).
    lower_residues = []
    for (i, _), entry in zip(edges, factors):
        lower_residues.extend(int(columns[i][1] % p == 0) for p in entry)
    lower = int(crt(all_primes, lower_residues)[0])
    tilted = [(x, y + lower * x) for x, y in columns]
    upper_residues = []
    for (i, j), entry in zip(edges, factors):
        x, y = tilted[i]
        for p in entry:
            target = parameter * pow(scale * roots[i] * roots[j], -1, p)
            upper_residues.append((target - x * pow(y, -1, p)) % p)
    upper = int(crt(all_primes, upper_residues)[0]) + parent**4
    moved = [(x + upper * y, y) for x, y in tilted]
    marked = [set(edge) for edge in edges]

    def verify_full(points, sides):
        values = [x*x + y*y for x, y in points]
        assert all(gcd(x, y) == 1 for x, y in points)
        for block, side in zip(blocks, sides):
            assert all(values[i] % block == 0 for i in side)
        for h, (i, j) in enumerate(edges):
            shared = prod(block for block, side in zip(blocks, sides) if {i, j} <= side)
            x, y = points[i]
            xx, yy = points[j]
            minor = x * yy - y * xx
            signed_cofactor = (columns[i][0]*columns[j][1] - columns[i][1]*columns[j][0]) // blocks[h]
            assert minor == signed_cofactor * shared
            assert (x * xx + y * yy) % shared == 0
            assert 4 * minor**2 * parent**6 < values[i] * values[j]

    verify_full(moved, marked)
    for h, (i, j) in enumerate(edges):
        if j-i != 1 and (i, j) != (0, 4):
            continue
        imag = scale * roots[i] * roots[j]
        adjusted = []
        for k, (x, y) in enumerate(moved):
            re, im = parameter * x + imag * y, parameter * y - imag * x
            divisor = blocks[h] if k in (i, j) else 1
            assert re % divisor == im % divisor == 0
            adjusted.append((re // divisor, im // divisor))
        moved = adjusted
        marked[h] = set(range(5)) - marked[h]
    verify_full(moved, marked)
    assert all(sum(i in side for side in marked) == 5 for i in range(5))
    print("CRT and shears give width < N^(-3) with all 20 endpoint incidences.")
    print("Gaussian reorientation verifies all 25 balanced incidences and all signed cofactors.")


if __name__ == "__main__":
    main()
