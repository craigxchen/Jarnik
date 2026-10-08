"""Exact checks of common scaling for pair-midpoint reflection unions."""

from functools import reduce
from itertools import combinations

from check_complement_reflection_invariant_relations import (
    ONE, conjugate, exact_div, gaussian_gcd, mul, norm, power,
    primitive_pair, product,
)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def lcm(a, b):
    return exact_div(mul(a, b), gaussian_gcd(a, b))


def main():
    prime_elements = []
    used_norms = set()
    for a in range(1, 35):
        for b in range(1, a):
            nn = a*a + b*b
            if prime(nn) and nn not in used_norms:
                prime_elements.append((a, b))
                used_norms.add(nn)
    assert len(prime_elements) >= 20
    cases = pairs = 0
    for m in (3, 4):
        for variant in range(8):
            blocks = {
                mask: power(prime_elements[mask-1], 1 + (mask + variant) % 3)
                for mask in range(1, 1 << m)
            }
            core = product(blocks.values())
            aa = [ONE] + [
                product(value for mask, value in blocks.items() if mask >> i & 1)
                for i in range(m)
            ]
            kk = [ONE]
            for i in range(m):
                inactive = next(mask for mask in blocks if not (mask >> i & 1))
                correction = power(conjugate(prime_elements[inactive-1]), 1 + variant % 2)
                if variant % 2:
                    correction = mul(correction, prime_elements[(1 << m)-1+i])
                kk.append(correction)
            hh = [mul(a, k) for a, k in zip(aa, kk)]
            assert all(norm(gaussian_gcd(h, conjugate(h))) == 1 for h in hh)
            ll = reduce(lcm, hh, ONE)
            jj = exact_div(ll, core)
            zz = [exact_div(mul(conjugate(ll), h), conjugate(h)) for h in hh]
            assert all(norm(z) == norm(ll) for z in zz)
            assert norm(reduce(gaussian_gcd, zz)) == 1
            assert len(set(zz)) == len(zz)

            def inside(mask, row):
                return row > 0 and bool(mask >> (row-1) & 1)

            for i, j in combinations(range(m+1), 2):
                sep = product(v for mask, v in blocks.items()
                              if inside(mask, i) != inside(mask, j))
                agree = exact_div(core, sep)
                numerator = mul(zz[i], zz[j])
                denominator = (norm(ll), 0)
                common = gaussian_gcd(numerator, denominator)
                rho_num = exact_div(numerator, common)
                rho_den = exact_div(denominator, common)
                assert norm(rho_num) == norm(rho_den)
                assert norm(gaussian_gcd(rho_num, rho_den)) == 1
                union = ([mul(rho_den, z) for z in zz]
                         + [mul(rho_num, conjugate(z)) for z in zz])
                assert norm(reduce(gaussian_gcd, union)) == 1
                assert len({norm(z) for z in union}) == 1
                assert mul(rho_num, conjugate(zz[i])) == mul(rho_den, zz[j])
                assert (norm(rho_den) * norm(kk[i]) * norm(kk[j]) * norm(jj)
                        >= norm(agree))
                vv = primitive_pair(hh[i], hh[j])
                assert vv[1] != 0
                assert norm(vv) <= norm(kk[i]) * norm(kk[j]) * norm(sep)
                chord = (zz[i][0]-zz[j][0], zz[i][1]-zz[j][1])
                assert norm(chord) * norm(vv) == 4 * norm(ll) * vv[1]**2

                # Fourth powers remove every square root in (8).
                budget = norm(sep) * norm(kk[i]) * norm(kk[j])
                assert (norm(rho_den) * norm(chord)**2 * budget**3
                        >= 16 * norm(ll) * norm(core)**2)
                pairs += 1
            cases += 1
    print(f'PASS: {cases} exact independent-block reconstructions and {pairs} midpoint reflections.')
    print('PASS: primitive integral unions, reduced denominator budgets, and endpoint-constant bound.')


if __name__ == '__main__':
    main()
