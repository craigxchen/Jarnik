"""Theorem C (direction pinning for even M at the CC scale): exact checks.

(1) Layer bookkeeping on random actual clusters: prod_i z_i equals a unit times
    prod_{layers} p^{min(h,M-h)} times a Gaussian integer of norm
    prod_{layers} p^{|M-2h|}  (h = number of points on the high side).
(2) Four-point Pell family (L=24): the primitive part of prod z_i (after removing
    its rational content) has bounded norm along the family, and 4*arg(z_0)
    converges to its argument, arg(-3+4i).
"""
import math, random
from itertools import combinations
from cot_lib import (gmul, gnorm, gexact, ggcd_all, gdivides, factor, split_prime_pi,
                     circle_points, angle, primitive_tuple, all_edge_lcm, is_clique)
from check_dictionary import SPLIT


def layer_check(seed=5, trials=600):
    random.seed(seed)
    cnt = 0
    for _ in range(trials):
        r = random.randint(1, 4)
        pp = [(p, random.randint(1, 3)) for p in random.sample(SPLIT, r)]
        N = 1
        for p, e in pp:
            N *= p ** e
        if N > 10 ** 12:
            continue
        pts = sorted(circle_points(pp), key=angle)
        n = len(pts)
        M = random.randint(2, min(8, n // 4))
        s = random.randrange(n)
        cl = [pts[(s + t) % n] for t in range(M)]
        g = ggcd_all(cl)
        cl = [gexact(z, g) for z in cl]
        Np = gnorm(cl[0])
        prod = (1, 0)
        for z in cl:
            prod = gmul(prod, z)
        real_part = 1
        rem_norm = 1
        for p, e in factor(Np).items():
            pi = split_prime_pi(p)
            al = []
            for z in cl:
                t = 0
                while gdivides(pi, z):
                    z = gexact(z, pi); t += 1
                al.append(t)
            for tau in range(1, e + 1):
                h = sum(1 for a in al if a >= tau)
                real_part *= p ** min(h, M - h)
                rem_norm *= p ** abs(M - 2 * h)
        assert prod[0] % real_part == 0 and prod[1] % real_part == 0
        rest = (prod[0] // real_part, prod[1] // real_part)
        assert gnorm(rest) == rem_norm
        cnt += 1
    print('layer bookkeeping identity verified on %d clusters' % cnt)


def pell(n):
    U, V = 1, 0
    for _ in range(n):
        U, V = 9 * U + 20 * V, 4 * U + 9 * V
    return U, V


def pell_check():
    target = math.atan2(4, -3)
    for n in (1, 11, 21, 31, 41, 51, 61):
        U, V = pell(n)
        L = 24
        X = sorted([90 * U * V + 30 * V * V + 3, 96 * U * V, 160 * V * V + 128 * U * V + 16])
        assert is_clique(X, L)
        rows = primitive_tuple(X, L)
        N = gnorm(rows[0])
        prod = (1, 0)
        for z in rows:
            prod = gmul(prod, z)
        c = math.gcd(prod[0], prod[1])
        phi = (prod[0] // c, prod[1] // c)
        # 4*arg z_0 - arg(phi), computed exactly enough via logs of big ints
        z = rows[0]
        def arg(w):
            # robust for huge ints
            s = max(abs(w[0]), abs(w[1]))
            k = max(0, s.bit_length() - 900)
            return math.atan2(w[1] >> k if w[1] >= 0 else -((-w[1]) >> k),
                              w[0] >> k if w[0] >= 0 else -((-w[0]) >> k))
        a4 = (4 * arg(z)) % (2 * math.pi)
        d = (a4 - arg(phi) + math.pi) % (2 * math.pi) - math.pi
        print('n=%2d digits(N)=%4d  primitive part of prod z_i = %s (norm %d)  '
              '4arg z0 - arg = %.3e  arg-target=%.2e' %
              (n, len(str(N)), phi, gnorm(phi), d, (arg(phi) - target)))


if __name__ == '__main__':
    layer_check()
    pell_check()
