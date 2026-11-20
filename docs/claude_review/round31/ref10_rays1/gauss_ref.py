"""Referee's check of the local X.2 data with actual Gaussian residues modulo l^a.

For an exterior odd prime l and a level a with l^(a-1) >= U, build norm-correct local factors
r_j mod l^a for 2 blocks of g = 2B physical labels, with the block products of r_j/conj(r_j)
equal to those of given "old" data, by the positional digit construction at the l-component
and plain block sums (with the norm parity) at the other components of d_a = (l - chi(l)) l^(a-1).
Then test, with Gaussian arithmetic only (no logarithms), every word of mass <= H on a label
subset, every n <= NMAX and every unit i^s: the exact collision level a' of
(prod_j (r_j/conj r_j)^{v_j})^n with i^s satisfies d_{a'} | 4 n gcd(Lambda, d_a), the l-part of
gcd(Lambda, d_a) is < U for words with nonzero class sum, and l^{a'} <= 6 n G_eff.
"""
import random, itertools
from math import gcd
from sympy import factorint, isprime, legendre_symbol
from sympy.ntheory import sqrt_mod

from digit_ref import build_block, all_class_vectors   # referee's own module


def gmul(x, y, Q):
    return ((x[0] * y[0] - x[1] * y[1]) % Q, (x[0] * y[1] + x[1] * y[0]) % Q)


def gconj(x, Q):
    return (x[0], (-x[1]) % Q)


def gnorm(x, Q):
    return (x[0] * x[0] + x[1] * x[1]) % Q


def ginv(x, Q):
    ni = pow(gnorm(x, Q), -1, Q)
    c = gconj(x, Q)
    return (c[0] * ni % Q, c[1] * ni % Q)


def gpow(x, n, Q):
    if n < 0:
        x, n = ginv(x, Q), -n
    r = (1, 0)
    while n:
        if n & 1:
            r = gmul(r, x, Q)
        x = gmul(x, x, Q)
        n >>= 1
    return r


def ratio(x, Q):
    return gmul(x, ginv(gconj(x, Q), Q), Q)


def vl(n, l):
    k = 0
    while n % l == 0 and n:
        n //= l
        k += 1
    return k


def coll_level(y, s, l, a):
    unit = [(1, 0), (0, 1), (-1, 0), (0, -1)][s]
    d0 = y[0] - unit[0]
    d1 = y[1] - unit[1]
    k = 0
    while k < a and d0 % l ** (k + 1) == 0 and d1 % l ** (k + 1) == 0:
        k += 1
    return k


def run(l, H, b, rng, nmax=12):
    B = 64 * H * H
    g = 2 * B
    rows = [[int(i == k) for k in range(b)] for i in range(b)]
    U = 100 * g * B ** len(rows)
    kl = 0
    while l ** kl < U:
        kl += 1
    a = kl + 1
    Q = l ** a
    chi = 1 if l % 4 == 1 else -1
    d = (l - chi) * l ** (a - 1)
    fd = factorint(d)
    # generator t0 = z0/conj z0 of T_Q
    while True:
        z0 = (rng.randrange(Q), rng.randrange(Q))
        if gnorm(z0, Q) % l == 0:
            continue
        t0 = ratio(z0, Q)
        if gpow(t0, d, Q) != (1, 0):
            raise AssertionError("order")
        if all(gpow(t0, d // p, Q) != (1, 0) for p in fd):
            break
    Nz0 = gnorm(z0, Q)
    assert legendre_symbol(Nz0 % l, l) == -1      # generator <-> nonsquare norm
    # physical primes
    primes = []
    p = 5
    while len(primes) < 2 * g:
        if isprime(p) and p % 4 == 1 and p != l:
            primes.append(p)
        p += 4
    nu = [int(legendre_symbol(p % l, l) == -1) for p in primes]
    classes = [rng.randrange(b) for _ in range(2 * g)]
    old = [2 * rng.randrange(d // 2) + nu[j] for j in range(2 * g)]
    L = [0] * (2 * g)
    for blk in range(2):
        idx = list(range(blk * g, (blk + 1) * g))
        S = sum(old[j] for j in idx) % d
        comps = []
        for r, e in fd.items():
            mod = r ** e
            Sr = S % mod
            par = [nu[j] for j in idx]
            if r == l:
                Lr, D, _ = build_block(H, rows, [classes[j] for j in idx], par, r, e, Sr)
            else:
                Lr = [rng.randrange(mod) for _ in idx]
                if r == 2:
                    Lr = [x - (x - pj) % 2 for x, pj in zip(Lr, par)]
                Lr[0] += (Sr - sum(Lr)) % mod if r != 2 else ((Sr - sum(Lr)) % mod)
                assert (sum(Lr) - Sr) % mod == 0
                if r == 2:
                    assert all((x - pj) % 2 == 0 for x, pj in zip(Lr, par))
            comps.append((mod, Lr))
        for pos, j in enumerate(idx):
            x, M = 0, 1
            for mod, Lr in comps:
                # CRT
                t = ((Lr[pos] - x) * pow(M, -1, mod)) % mod
                x += M * t
                M *= mod
            L[j] = x % d
    D = l ** kl
    # build r_j and check norms, ratios, parities, block products
    t = []
    for j, pj in enumerate(primes):
        assert L[j] % 2 == nu[j]
        r0 = gpow(z0, L[j], Q)
        q = pj * pow(gnorm(r0, Q), -1, Q) % Q
        sq = sqrt_mod(q, Q)
        assert sq is not None
        rj = (sq * r0[0] % Q, sq * r0[1] % Q)
        assert gnorm(rj, Q) == pj % Q
        tj = ratio(rj, Q)
        assert tj == gpow(t0, L[j], Q)
        # the other parity is impossible: p_j / N(z0)^(L+1) is a nonsquare mod l
        q2 = pj * pow(gnorm(gpow(z0, L[j] + 1, Q), Q), -1, Q) % l
        assert legendre_symbol(q2, l) == -1
        t.append(tj)
    for blk in range(2):
        idx = range(blk * g, (blk + 1) * g)
        newp, oldp = (1, 0), (1, 0)
        for j in idx:
            newp = gmul(newp, t[j], Q)
            oldp = gmul(oldp, gpow(t0, old[j], Q), Q)
        assert newp == oldp
    # word tests
    labs = sorted(set(list(range(0, 2 * g, 32)) + [1, 2 * g - 1]))
    nt = nu_ = 0
    maxlev = 0
    for h in range(1, H + 1):
        for combo in itertools.combinations_with_replacement(labs, h):
            for signs in itertools.product((1, -1), repeat=h):
                v = {}
                for j, sg in zip(combo, signs):
                    v[j] = v.get(j, 0) + sg
                v = {j: c for j, c in v.items() if c}
                if not v:
                    continue
                w = [0] * b
                for j, c in v.items():
                    w[classes[j]] += c
                lam = sum(c * L[j] for j, c in v.items())
                gg = gcd(lam, d)
                x = (1, 0)
                for j, c in v.items():
                    x = gmul(x, gpow(t[j], c, Q), Q)
                tested = any(w)
                if tested:
                    assert vl(gg, l) <= kl - 1
                    geff = gg
                y = (1, 0)
                for n in range(1, nmax + 1):
                    y = gmul(y, x, Q)
                    for s in range(4):
                        ap = coll_level(y, s, l, a)
                        if ap == 0:
                            continue
                        dap = (l - chi) * l ** (ap - 1)
                        assert (4 * n * gcd(lam, dap)) % dap == 0
                        if tested:
                            assert l ** ap <= 6 * n * geff
                            maxlev = max(maxlev, ap)
                if tested:
                    nt += 1
                else:
                    nu_ += 1
    return dict(l=l, a=a, kl=kl, Q_digits=len(str(Q)), labels=2 * g, tested=nt, untested=nu_,
                max_level_tested=maxlev)


if __name__ == '__main__':
    rng = random.Random(7)
    for l in (3, 5, 7, 13):
        print(run(l, 2, 3, rng))
    print("ALL GAUSSIAN LOCAL CHECKS PASSED")
