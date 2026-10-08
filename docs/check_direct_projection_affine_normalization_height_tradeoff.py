"""Exact affine gcd, Gaussian residue, and denominator-clearing fixtures.

No synthetic large endpoint family or finite test of the large-q diameter
argument is asserted; that part of the note is a Gram/sine proof.
"""

from itertools import combinations
from math import gcd, isqrt, lcm, log, pi, prod, sqrt

from check_direct_projection_monic_height_bridge import points_on_circle, primitive
from check_gaussian_reflection_replacement import ggcd, gnorm
from check_integer_cotangent_offset_height_cut_residue_dictionary import (
    factor, gaussian_prime, gaussian_v, vp,
)


def halfangle(z0: tuple[int, int], z: tuple[int, int], n: int) -> tuple[int, int, int, int]:
    dot = z0[0] * z[0] + z0[1] * z[1]
    cross = z0[0] * z[1] - z0[1] * z[0]
    if dot == -n:
        a, b = 0, 1
    else:
        common = gcd(n + dot, cross)
        a, b = (n + dot) // common, cross // common
    assert gcd(a, b) == 1 and b != 0
    f = a * a + b * b
    assert 2 * n % f == 0
    d = 2 * n // f
    c = n - dot
    assert c == d * b * b and 2 * n - c == d * a * a
    return a, b, d, c


def aff_data(z0: tuple[int, int], selected: tuple[tuple[int, int], ...], n: int,
             prime_check: bool = True) -> tuple[int, int, int, int]:
    assert len(selected) >= 2 and len({z0, *selected}) == len(selected) + 1
    h = [halfangle(z0, z, n) for z in selected]
    aa, bb, dd, cc = zip(*h)
    assert len(set(cc)) == len(cc)
    gaff = gcd(*(c - cc[0] for c in cc))
    v = 2 * n - cc[0]
    geff = gcd(gaff, v)
    assert geff == gcd(*(2 * n - c for c in cc))
    dcommon, acommon = gcd(*dd), gcd(*(abs(a) for a in aa))
    assert geff == dcommon * acommon * acommon

    ps = (set().union(*(factor(value) for value in
                       (gaff, geff, n, 2*n, *dd, *bb, *cc)))
          if prime_check else set())
    for p in ps:
        s = [vp(c, p) for c in cc]
        assert s == [vp(d, p) + 2 * vp(b, p)
                     for d, b in zip(dd, bb)]
        valuations = []
        for i in range(1, len(cc)):
            if s[i] != s[0]:
                predicted = min(s[i], s[0])
            else:
                ui = cc[i] // p**s[i]
                ur = cc[0] // p**s[0]
                predicted = s[0] + vp(ui - ur, p)
            valuations.append(predicted)
            assert predicted == vp(cc[i] - cc[0], p)
        assert min(valuations) == vp(gaff, p)
        assert vp(geff, p) == min(vp(d, p) for d in dd) + 2 * min(vp(abs(a), p) if a else 10**9 for a in aa)

    for i, j in combinations(range(len(h)), 2):
        ai, bi = aa[i], bb[i]
        aj, bj = aa[j], bb[j]
        normg = gnorm(ggcd((ai, bi), (aj, bj)))
        re_num = aj * ai + bj * bi
        im_num = bj * ai - aj * bi
        assert re_num % normg == im_num % normg == 0
        re, im = re_num // normg, im_num // normg
        assert im != 0 and im % acommon == 0
        pair_g = gnorm(ggcd(selected[i], selected[j]))
        assert n % pair_g == 0
        pair_norm = n // pair_g
        assert re*re + im*im in (pair_norm, 2*pair_norm)

    bden = gaff // geff
    anum = v // geff
    assert gcd(anum, bden) == 1
    for size in (2, 4):
        for idx in combinations(range(len(h)), size):
            label_product = prod(dd[i] for i in idx)
            root = isqrt(label_product)
            if root * root != label_product:
                continue
            original = prod(2 * n - cc[i] for i in idx)
            assert original == (root * prod(aa[i] for i in idx)) ** 2
            integer_cleared = prod(anum - bden * ((cc[i] - cc[0]) // gaff)
                                   for i in idx)
            assert integer_cleared == original // geff**size
            assert isqrt(integer_cleared)**2 == integer_cleared
    return gaff, geff, dcommon, acommon


def gamma_exact(points: tuple[tuple[int, int], ...], n: int) -> float:
    """Common-threshold weight for all these physical rows."""
    out = 0.0
    for p, exponent in factor(n).items():
        assert p % 4 == 1
        pi = gaussian_prime(p)
        es = [gaussian_v(z, pi) for z in points]
        assert all(0 <= e <= exponent for e in es)
        out += (exponent - max(es) + min(es)) * log(p)
    return out


def pell_points(power: int) -> tuple[tuple[int, int], ...]:
    u, v = 1, 0
    for _ in range(power):
        u, v = 9*u + 20*v, 4*u + 9*v
    kx = 40*u*v*v + 32*u*u*v
    ky = 60*u*v*v + 4*u*u*v
    return (
        (kx + u - 17*v, ky + 8*u + 6*v),
        (kx + u + 15*v, ky - 8*u - 10*v),
        (kx - u - 15*v, ky + 8*u + 10*v),
        (kx + 7*u - 15*v, ky + 4*u - 10*v),
    )


def main() -> None:
    z0 = (-3, -2)
    selected = ((3, -2), (2, -3), (-3, 2))
    assert aff_data(z0, selected, 13) == (5, 1, 1, 1)
    assert [halfangle(z0, z, 13)[3] for z in selected] == [18, 13, 8]

    cases = 0
    for n in (5, 13, 25, 65, 85, 125, 325, 625):
        points = points_on_circle(n)
        for anchor in points:
            for rows in combinations([z for z in points if z != anchor], 3):
                if not primitive((anchor, *rows)):
                    continue
                deficits = [n - anchor[0]*z[0] - anchor[1]*z[1] for z in rows]
                if len(set(deficits)) != 3:
                    continue
                gaff, geff, dcommon, acommon = aff_data(anchor, rows, n)
                gamma = gamma_exact((anchor, *rows), n)
                gamma_i = gamma_exact(rows, n)
                label_lcm = lcm(*(halfangle(anchor, z, n)[2] for z in rows))
                assert dcommon <= 2 * gcd(*anchor) * (2.718281828459045**gamma) + 1e-8
                assert log(label_lcm) + 1e-8 >= log(n) - log(gcd(*anchor)) - gamma_i
                cases += 1
    assert cases > 1000

    short = pell_points(11)
    n = gnorm(short[0])
    assert primitive(short) and all(gnorm(z) == n for z in short)
    max_chord = max(sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2)
                    for x, y in combinations(short, 2))
    c_upper = (pi/2) * max_chord / n**0.25
    assert c_upper < sqrt(2)
    for anchor in short:
        rows = tuple(z for z in short if z != anchor)
        deficits = [n - anchor[0]*z[0] - anchor[1]*z[1] for z in rows]
        if len(set(deficits)) != 3:
            continue
        _, geff, dcommon, acommon = aff_data(anchor, rows, n,
                                           prime_check=False)
        label_lcm = lcm(*(halfangle(anchor, z, n)[2] for z in rows))
        assert log(label_lcm) + 1e-8 >= log(n)*(1-1/len(rows)) - log(gcd(*anchor))
        assert log(acommon) <= log(n) / (4*(len(rows)-1)) + 1e-8
        assert log(geff) <= log(2*gcd(*anchor)) + log(n)/len(short) + log(n)/(2*(len(rows)-1)) + 1e-8
    print(f"PASS: {cases} exact affine/source subsets, rational denominator fixture, and short Pell tuple")


if __name__ == "__main__":
    main()
