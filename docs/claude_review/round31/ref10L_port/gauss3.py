"""Rank-one (vertex-label) form of separated chord^2 modulo a split prime."""
import sys
sys.path.insert(0, '.')
from gauss import split_gen, gmul, conj, gsub, norm
from gauss2 import res_mod_pi
from itertools import combinations, product
from math import prod
import random
random.seed(2)
cnt = cyc = 0
for primes in ([5, 13, 17, 29], [13, 17, 29, 37, 41], [53, 113, 257, 293, 317, 337], [5, 13, 17, 29, 37, 41, 53]):
    gens = {p: split_gen(p) for p in primes}
    for eta in ((1, 0), (0, 1), (-1, 0), (0, -1)):
        marks = []
        for bits in product((0, 1), repeat=len(primes)):
            z = eta
            for p, b in zip(primes, bits):
                z = gmul(z, gens[p] if b == 0 else conj(gens[p]))
            marks.append((bits, z))
        for k, p in enumerate(primes):
            pi = gens[p]
            P0 = [z for b, z in marks if b[k] == 0]   # pi | z
            P1 = [z for b, z in marks if b[k] == 1]   # conj(pi) | z
            for zi in P0:
                for zj in P1:
                    lhs = norm(gsub(zi, zj)) % p
                    rhs = (-res_mod_pi(zi, conj(pi), p) * res_mod_pi(zj, pi, p)) % p
                    assert lhs == rhs and lhs != 0
                    cnt += 1
            for _ in range(200):
                a, b = random.sample(P0, 2); c, d = random.sample(P1, 2)
                assert (norm(gsub(a, c))*norm(gsub(b, d)) - norm(gsub(a, d))*norm(gsub(b, c))) % p == 0
                cyc += 1
print("rank-one identity checked on", cnt, "separated pairs;", cyc, "4-cycle cross-ratio checks")
