"""Referee check of the Part IV trace construction over F_q(t), q = 2^d.

Everything is recomputed from the definitions, with gcds by the Euclidean algorithm
(not from the predicted common marked set).
Usage: python3 trace_check.py d1 d2 ...
"""
import sys
from fractions import Fraction
from itertools import combinations
from gf2d import GF


def run(d):
    F = GF(d)
    q = F.q
    k = q // 2 - 1
    nz = list(range(1, q))
    # parent N = prod_{x != 0} (t - x); check = t^(q-1) + 1
    N = F.from_roots(nz)
    assert N == [1] + [0] * (q - 2) + [1]
    data = {}
    for a in nz:
        # L_a(t) = a^(-q/2) * sum_j (a t)^(2^j)
        La = [0] * (q // 2 + 1)
        c0 = F.pw(F.inv(a), q // 2)
        for j in range(d):
            e = 1 << j
            La[e] ^= F.mul(c0, F.pw(a, e))
        La = F.strip(La)
        assert La[-1] == 1 and F.deg(La) == q // 2 and La[0] == 0  # monic, t | L_a
        # f_a = (t^(q/2) - L_a)/t
        num = F.add([0] * (q // 2) + [1], La)
        assert num[0] == 0
        fa = F.strip(num[1:])
        assert F.deg(fa) == q // 4 - 1
        marked = [x for x in nz if F.tr(F.mul(a, x)) == 0]
        assert len(marked) == q // 2 - 1
        na = F.from_roots(marked)
        # n_a = L_a / t  (roots of L_a are exactly the trace hyperplane)
        assert na == F.strip(La[1:])
        # n_a | N
        assert F.pdivmod(N, na)[1] == []
        # f_a(x) = x^(q/2-1) at marked x, and at UNmarked x it is different
        for x in nz:
            v = F.peval(fa, x)
            if x in marked:
                assert v == F.pw(x, k)
            else:
                assert v != F.pw(x, k)
        # identity f_a = t^(q/2-1) - n_a
        assert fa == F.add([0] * k + [1], na)
        # norm-type identities
        lead = F.pw(F.inv(a), q // 2)
        tfa2p1 = F.add(F.pmul([0, 1], F.pmul(fa, fa)), [1])
        assert tfa2p1 == F.scal(lead, na)          # t f_a^2 + 1 = a^(-q/2) n_a
        co, rem = F.pdivmod(N, na)
        assert rem == [] and co == F.add(La, [lead])  # N/n_a = L_a + a^(-q/2)
        data[a] = dict(f=fa, n=na, co=co, marked=set(marked))
    L = q - 1
    npairs = 0
    min_slack = None
    for a, b in combinations(nz, 2):
        A, B = data[a], data[b]
        D = F.add(A['f'], B['f'])                  # minor [P_a,P_b] = f_a - f_b (columns (f,1))
        g = F.pgcd(A['n'], B['n'])
        assert F.deg(g) == q // 4 - 1
        assert g == F.from_roots(sorted(A['marked'] & B['marked']))
        assert F.deg(D) == q // 4 - 1
        quo, rem = F.pdivmod(D, g)
        assert rem == [] and F.deg(quo) == 0       # nonzero constant cofactor
        # degree budget: deg D <= (alpha_a+alpha_b)/2 - L/4
        budget = Fraction(F.deg(A['n']) + F.deg(B['n']), 2) - Fraction(L, 4)
        slack = budget - F.deg(D)
        assert slack == Fraction(1, 4)
        # pair separation: deg d_ab = |marked_a xor marked_b| = q/2 > L/2
        sep = len(A['marked'] ^ B['marked'])
        assert sep == q // 2 and 2 * sep > L
        # hyperbola (split-form) cluster: both coordinates of Z_a=(n_a, N/n_a) are close
        assert F.deg(F.add(A['n'], B['n'])) == q // 4 - 1
        assert F.deg(F.add(A['co'], B['co'])) == q // 4
        # unimodular change: f_a - f_b = n_a - n_b
        assert D == F.add(A['n'], B['n'])
        npairs += 1
    return q, npairs


if __name__ == '__main__':
    for d in map(int, sys.argv[1:]):
        q, np_ = run(d)
        print(f"q={q}: m={q-1} marks, {np_} pairs: all checks pass "
              f"(deg D = deg g = {q//4-1}, constant cofactor, slack 1/4, separation {q//2} > {(q-1)/2}, "
              f"W = (q-1) log q, hyperbola degrees {q//4-1},{q//4} vs size {q//2})")
