"""Exact audits of Möbius content transfer and the compressed radius bound."""

from itertools import combinations, product
from math import gcd, lcm, prod

from check_gaussian_reflection_replacement import (
    gconj, gmul, gnorm, ggcd, gdivexact,
)
from check_gaussian_near_real_divisor_reduction import pair_factor


def least_norm(rows):
    out = 1
    for x, y in combinations(rows, 2):
        numerator = gmul(x, gconj(y))
        denominator = gmul(gconj(x), y)
        common = ggcd(numerator, denominator)
        out = lcm(out, gnorm(gdivexact(denominator, common)))
    return out


def residues(rows):
    out = []
    for x, y in combinations(rows, 2):
        determinant = abs(x[0]*y[1]-x[1]*y[0])
        common = gnorm(ggcd(x, y))
        assert determinant and determinant % common == 0
        out.append(determinant // common)
    return out


def valuation_audit():
    cases = 0
    for m in range(3, 6):
        for ts in product(range(-2, 3), repeat=m):
            for k in range(3):
                beta = sum(max(0, min(abs(u), abs(v))-k)
                           for u, v in combinations(ts, 2) if u*v > 0)
                loss = sum(map(abs, ts)) - (max(ts)-min(ts))
                assert loss <= (m-2)*k+beta+2*abs(ts[0])
                cases += 1
    return cases


def fixtures():
    sets = [
        [(1, 0), (3, 2), (5, 2)],
        [(1, 0), (3, 2), (5, 2), (7, 2), (9, 2)],
        [(1, 0), (2, 1), (1, 2), (3, 2), (2, 3), (5, 2)],
    ]
    # Actual four-point C<2 Pell cluster; the common source point is anchor 0.
    u, v = 9, 4
    blocks = [(-1, 2), (u+v, 2*v), (2*v, u-v), (2*u+8*v, 3*u+v)]
    exponents = [(1, 1, 1, 1), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    pell = [(1, 0)] + [pair_factor(blocks, row, exponents[0])
                       for row in exponents[1:]]
    sets.append(pell)
    return sets


def matrix_audit():
    count = conformal = 0
    source_sets = fixtures()
    for entries in product(range(-3, 4), repeat=4):
        a, b, c, d = entries
        delta = a*d-b*c
        if not delta or gcd(*entries) != 1:
            continue
        A, B, D = a*a+c*c, a*b+c*d, b*b+d*d
        trace = A+D
        K = (A-D)**2+4*B*B
        assert K == trace*trace-4*delta*delta
        for source_index, rows in enumerate(source_sets):
            assert rows[0] == (1, 0)
            assert all(gcd(*h) == 1 and gnorm(h) % 2 for h in rows)
            raw = [(a*x+b*y, c*x+d*y) for x, y in rows]
            contents = [gcd(*h) for h in raw]
            assert all(abs(delta) % r == 0 for r in contents)
            images = [(x//r, y//r) for (x, y), r in zip(raw, contents)]
            old_norm, new_norm = least_norm(rows), least_norm(images)
            if K == 0:
                # All rows, including the transformed anchor, must be retained.
                assert new_norm == old_norm
                conformal += 1
                continue
            old_b = residues(rows)
            for h, raw_h, image_h in zip(rows, raw, images):
                assert K % gcd(gnorm(h), gnorm(raw_h)) == 0
                assert trace*gnorm(image_h) >= gnorm(h)
            Q = gcd(K, abs(delta)*old_norm)  # Source-capped transfer.
            for (h, k), (x, y), residue in zip(combinations(rows, 2), combinations(images, 2), old_b):
                pair_cap = gcd(K, abs(delta)*gnorm(ggcd(h, k)))
                assert pair_cap*residue % gnorm(ggcd(x, y)) == 0
                assert Q*residue % gnorm(ggcd(x, y)) == 0
            m = len(rows)
            residual_product = prod(old_b)
            F0 = gnorm(images[0])
            compressed_denominator = 2*F0*F0*Q**(m-2)*residual_product
            assert new_norm*compressed_denominator >= prod(map(gnorm, images))
            anchored_denominator = 2*trace**m*Q**(m-2)*residual_product
            assert new_norm*anchored_denominator >= prod(map(gnorm, rows))
            # Joint accounting at determinant primes coprime to K, including 2.
            good_determinant = abs(delta)
            while gcd(good_determinant, K) > 1:
                good_determinant //= gcd(good_determinant, K)
            joint_numerator = good_determinant**(2*m-3)*prod(map(gnorm, rows))
            assert new_norm*anchored_denominator*trace**2 >= joint_numerator
            if abs(delta) == 1:
                assert old_norm % Q == 0
                exceptions = sum(8*old_norm*gnorm(z) < trace*gnorm(h)
                                 for z, h in zip(images, rows))
                assert exceptions <= 1
            all_pairs_denominator = 2*trace**m*K**(m*(m-1)//2)*residual_product
            assert new_norm*all_pairs_denominator >= joint_numerator
            if source_index == 3:
                # C<=2 consequence; remove the old residue product using Plotkin.
                assert residual_product**8 <= old_norm**m
                assert all(gnorm(h)**2 >= old_norm for h in rows[1:])
                endpoint_denominator = 2*trace**m*Q**(m-2)
                assert (new_norm*endpoint_denominator)**8 >= old_norm**(3*m-4)
                assert (new_norm*endpoint_denominator*trace**2)**8 >= good_determinant**(8*(2*m-3))*old_norm**(3*m-4)
                if abs(delta) == 1:
                    # The large-trace product lower bound at C=2, squared.
                    assert prod(map(gnorm, images))**2*(8*old_norm)**(2*(m-1)) >= trace**(2*(m-1))*old_norm**(m-2)
            count += 1
    return count, conformal


def joint_local_sharp_audit():
    def valuation(value, prime):
        assert value
        value = abs(value)
        answer = 0
        while value % prime == 0:
            value //= prime
            answer += 1
        return answer

    count = 0
    for ell in range(1, 9):
        for prime in (2, 17):
            q = prime**ell
            if prime == 2:
                matrix = (q+1,0,0,q)
                rows = [(1,0),(q,1),(q,3)]
            else:
                matrix = (1,q,1,0)
                rows = [(1,0),(q,2),(2*q,1)]
            a,b,c,d = matrix
            delta = a*d-b*c
            trace = a*a+b*b+c*c+d*d
            K = trace*trace-4*delta*delta
            assert K % prime and valuation(delta,prime) == ell
            assert all(gcd(*z) == 1 and gnorm(z) % 2 for z in rows)
            raw = [(a*x+b*y,c*x+d*y) for x,y in rows]
            rs = [gcd(*z) for z in raw]
            images = [(x//r,y//r) for (x,y),r in zip(raw,rs)]
            cs = [valuation(r,prime) for r in rs]
            assert cs == [0,ell,ell]
            beta = [valuation(v,prime) for v in residues(rows)]
            if prime == 2:
                parities = [gnorm(z) % 2 == 0 for z in images]
                assert sum(parities) == 2
                assert 2*sum(cs)+sum(parities) == 3*ell+sum(beta)+1
                for (i,j), val in zip(combinations(range(3),2),beta):
                    assert val >= min(cs[i],cs[j])+int(parities[i] and parities[j])
            else:
                assert all(gnorm(z) % prime for z in images)
                assert 2*sum(cs) == 3*ell+sum(beta)
            count += 1
    return count


def main():
    valuations = valuation_audit()
    matrices, conformal = matrix_audit()
    sharp_cases = joint_local_sharp_audit()
    print(f"PASS: {valuations} compressed signed-valuation cases.")
    print(f"PASS: {matrices} nonconformal matrix/configuration cases, including ordinary contents and powers of two.")
    print(f"PASS: {conformal} conformal cases preserve the full configuration's least norm.")
    print("PASS: direct Gaussian all-pair radius reconstruction and the actual Pell endpoint consequence.")
    print("PASS: joint content gain from every determinant prime coprime to K, including two.")
    print("PASS: source-capped pair and compressed bounds, with unimodular expansion and its Pell endpoint consequence.")
    print(f"PASS: {sharp_cases} sharp local examples with determinant valuations up to eight, including two.")


if __name__ == '__main__':
    main()
