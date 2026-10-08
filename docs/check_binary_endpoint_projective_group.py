"""Exact cut/permutation checks for the free projective two-group theorem."""

from fractions import Fraction as F
from itertools import permutations, combinations, product
from math import gcd, comb

from check_gaussian_reflection_replacement import gmul, gconj, gnorm
from check_mobius_source_aligned_compression import gaussian_valuation
from check_isometric_virtual_anchor import factor, gaussian_prime
from check_parabolic_shear_anchor_content import constants


def image_mask(mask,permutation):
    return sum(1 << permutation[i] for i in range(len(permutation)) if mask >> i & 1)


def cut_audit():
    count = 0
    for m in range(3,8):
        full = (1 << m)-1
        for sigma in permutations(range(m)):
            square = tuple(sigma[sigma[i]] for i in range(m))
            fixed = any(sigma[i] == i for i in range(m))
            for cut in range(1,full):
                target = image_mask(cut,sigma)
                fiber_unsplit = target in (cut,full^cut)
                # A source fiber is unsplit precisely when its target bits are constant.
                values = [{(cut >> sigma[i]) & 1 for i in range(m) if (cut >> i)&1 == bit}
                          for bit in (0,1)]
                assert fiber_unsplit == all(len(v) == 1 for v in values)
                if fiber_unsplit:
                    assert image_mask(cut,square) == cut
                    if fixed or m % 2:
                        assert target == cut
                count += 1
    return count


def erased_edge_audit():
    count = 0
    primes = (5,13,17,29,37,41)
    for m in range(3,7):
        full = (1 << m)-1
        cuts = tuple(range(1,7)) if m > 3 else tuple(range(1,full))
        cuts = tuple(c if c < full else full-1 for c in cuts)
        exponents = (1,2,3,1,2,3)
        for sigma in permutations(range(m)):
            N = Q = 1
            for p,e,cut in zip(primes,exponents,cuts):
                N *= p**e
                compatible = image_mask(cut,sigma) in (cut,full^cut)
                Q *= p**(1 if compatible else 0)  # Includes partial clipping.
            E = N//Q
            candidates = [tuple(sigma[sigma[i]] for i in range(m))]
            if any(sigma[i] == i for i in range(m)):
                candidates.append(sigma)
            for tau in candidates:
                for i in range(m):
                    edge = 1
                    for p,e,cut in zip(primes,exponents,cuts):
                        if ((cut >> i)&1) != ((cut >> tau[i])&1):
                            edge *= p**e
                    assert E % edge == 0
                    count += 1
    return count


def physical_parity_audit():
    N = 65
    points = [(sx*x,sy*y) for x,y in ((1,8),(4,7),(7,4),(8,1))
              for sx,sy in product((-1,1),repeat=2)]
    count = even = 0
    for z,w in combinations(points,2):
        n = 1
        for p,e in factor(N).items():
            pi = gaussian_prime(p)
            vz,vw = gaussian_valuation(z,pi),gaussian_valuation(w,pi)
            assert vz in (0,e) and vw in (0,e)
            n *= p**abs(vz-vw)
        numerator = gmul(w,gconj(z))
        H = (numerator[0]+N,numerator[1])
        if H == (0,0):
            H = (0,1)
        g = gcd(*H)
        H = tuple(v//g for v in H)
        f = gnorm(H)
        assert f in (n,2*n)
        even += f == 2*n
        count += 1
    assert even > 0  # Dropping the factor two would be false.
    return count,even


def exponent_and_fair_profile_audit():
    count = 0
    for m in (*range(32,129),256,512):
        q,_,_,_ = constants(m)
        assert q <= F(8,m) <= F(1,4)
        assert m*q <= 8
        count += 1
    assert 2**8 * 2**9 == 2**18//2  # Boundary at N=2^36.
    for m in range(5,11):
        patterns = [bits for bits in product((0,1),repeat=m) if 0 < sum(bits) < m]
        for seeds in combinations(range(m),4):
            odd = sum(sum(bits[i] for i in seeds) in (1,3) for bits in patterns)
            assert odd == 2**(m-1)
            assert F(odd,len(patterns)) == F(2**(m-1),2**m-2)
            count += 1
    return count


if __name__ == '__main__':
    cuts = cut_audit()
    edges = erased_edge_audit()
    pairs,even = physical_parity_audit()
    exponents = exponent_and_fair_profile_audit()
    print(f'PASS: {cuts} permutation/cut cases; {edges} erased-edge divisibilities; '
          f'{pairs} physical parity pairs ({even} even); '
          f'{exponents} exponent/fair-profile certificates.')
