"""Residues of separated Gaussian differences modulo pi_p and conj(pi_p)."""
import sys
sys.path.insert(0, '.')
from gauss import split_gen, gmul, conj, gsub, norm
from itertools import combinations, product
from math import prod

def res_mod_pi(z, pi, p):
    # Z[i]/pi = F_p via i -> r where r = -a/b mod p for pi = a+bi (b != 0 mod p): a + b r = 0
    a, b = pi
    r = (-a * pow(b, -1, p)) % p
    return (z[0] + z[1]*r) % p

for primes in ([5, 13, 17, 29], [13, 17, 29, 37, 41], [53, 113, 257, 293, 317, 337]):
    gens = {p: split_gen(p) for p in primes}
    q = prod(primes)
    marks = []
    for bits in product((0, 1), repeat=len(primes)):
        z = (1, 0)
        for p, b in zip(primes, bits):
            z = gmul(z, gens[p] if b == 0 else conj(gens[p]))
        marks.append((bits, z))
    out = {}
    for k, p in enumerate(primes):
        pi = gens[p]
        diffpi, diffpibar, t2 = set(), set(), set()
        for (bi, zi), (bj, zj) in combinations(marks, 2):
            if bi[k] == bj[k]:
                continue
            if bi[k] == 1:
                (bi, zi), (bj, zj) = (bj, zj), (bi, zi)
            # now pi | z_i, conj(pi) | z_j
            assert res_mod_pi(zi, pi, p) == 0 and res_mod_pi(zj, conj(pi), p) == 0
            d = gsub(zi, zj)
            diffpi.add(res_mod_pi(d, pi, p))
            diffpibar.add(res_mod_pi(d, conj(pi), p))
            h = prod(pp for pp, x, y in zip(primes, bi, bj) if x == y)
            t2.add((norm(d)//(4*h)) % p)
        out[p] = (len(diffpi), len(diffpibar), len(t2))
    print(primes, "per prime: #residues of z_i-z_j mod pi, mod conj(pi), #t^2 mod p (separated pairs):", out)
