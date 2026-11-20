"""Referee check of the square-block encoding (VIII.3) of Part VIII.2.

Independent code (does not import Astra's script).  Exact integer arithmetic.
For random and adversarial moduli q and prescribed residue pairs at every
p^e || q, take subsets X of Omega(q), reduce (x - x_*)/g with
g = gcd(q, differences), Q = q/g, and check for many ordered distinct
8-tuples of the reduced grid Y (anchor = slot 0):
  * local structure: at every p^f || Q, Y occupies exactly two classes mod p^f
    and their difference is a unit mod p;
  * columns P_i = (Q/d_i, (y_i-y_0)/d_i), d_i = gcd(Q, y_i-y_0): primitive,
    P_0 = (1,0);
  * blocks n_T = prod p^{2f} over p^f || Q whose side opposite the anchor is T
    (T nonempty subset of {1..7});  pairwise coprime squares, chain condition;
  * det(P_i,P_j) = t_ij * prod_{T contains i,j} n_T with
    t_ij = (y_j-y_i)/gcd(Q,y_j-y_i) a nonzero integer, t_0j = Y_j;
  * the valuation of Q h_ij/(d_i d_j) at p is 2f iff both i,j lie in T_p,
    0 otherwise (anchor pairs included);
  * D_ij (mass where signs differ, atom mass 2 f log p) equals 2 log(Q/h_ij),
    checked as an integer identity  prod_{p separates i,j} p^{2f} = (Q/h_ij)^2;
  * cofactor identity  log|t_ij| = D_ij/2 - W/4 + log|y_j-y_i| - (1/2)log Q,
    checked as the integer identity |t_ij|^2 * Q * h_ij^2 ... (see code).
"""
import random
from math import gcd, prod
from itertools import combinations
from sympy import factorint

rng = random.Random(20261119)

PRIMES = [2, 3, 5, 7, 11, 13, 17]


def crt(residues, moduli):
    M = prod(moduli)
    x = 0
    for r, m in zip(residues, moduli):
        Mi = M // m
        x += r * Mi * pow(Mi, -1, m)
    return x % M


def random_instance():
    k = rng.randint(1, 5)
    ps = rng.sample(PRIMES, k)
    if rng.random() < 0.5 and 2 not in ps:
        ps[0] = 2
    mods, choices, modes = [], [], []
    for p in ps:
        e = rng.randint(1, 6 if p == 2 else 4)
        pe = p ** e
        mode = rng.choice(["single", "random", "coalesce", "coalesce", "unitdiff"])
        r1 = rng.randrange(pe)
        if mode == "single":
            rs = [r1]
        elif mode == "random":
            rs = sorted({r1, rng.randrange(pe)})
        elif mode == "coalesce" and e >= 2:
            s = rng.randint(1, e - 1)
            u = rng.randrange(1, p ** (e - s))
            while u % p == 0:
                u = rng.randrange(1, p ** (e - s))
            rs = [r1, (r1 + u * p ** s) % pe]
        else:
            u = rng.randrange(1, pe)
            while u % p == 0:
                u = rng.randrange(1, pe)
            rs = [r1, (r1 + u) % pe]
        mods.append(pe)
        choices.append(rs)
        modes.append(mode)
    q = prod(mods)
    # Omega(q) via CRT, then lift by random multiples of q (marks are integers,
    # not necessarily in [0,q)); keep it finite.
    omega = [0]
    cur = [[]]
    for rs in choices:
        cur = [c + [r] for c in cur for r in rs]
    omega = sorted({crt(c, mods) for c in cur})
    pts = set()
    for r in omega:
        for t in range(rng.randint(1, 4)):
            pts.add(r + q * rng.randint(-3, 3))
    pts = sorted(pts)
    return q, mods, choices, modes, pts


def check_tuple(Q, fac, ys, stats):
    y0 = ys[0]
    d = [gcd(Q, y - y0) for y in ys]
    assert d[0] == Q
    cols = [(Q // di, (y - y0) // di) for y, di in zip(ys, d)]
    assert cols[0] == (1, 0)
    for X, Yc in cols:
        assert gcd(X, Yc) == 1
    # slopes = affine image of marks
    for (X, Yc), y in zip(cols, ys):
        assert Yc * Q == (y - y0) * X
    blocks = {}
    side = {}
    for p, f in fac.items():
        pf = p ** f
        T = frozenset(i for i in range(1, 8) if (ys[i] - y0) % pf)
        side[p] = T
        if T:
            blocks[T] = blocks.get(T, 1) * pf * pf
    # pairwise coprime squares and chain condition (one set per prime)
    bl = list(blocks.items())
    for (T1, n1), (T2, n2) in combinations(bl, 2):
        assert gcd(n1, n2) == 1
    for p in fac:
        sets = [T for T, n in bl if n % p == 0]
        assert len(sets) <= 1
    for T, n in bl:
        r = int(round(n ** 0.5)) if n < 10 ** 30 else None
        fn = factorint(n)
        assert all(v % 2 == 0 for v in fn.values())
    if len(set(side.values()) - {frozenset()}) < len(fac):
        stats["shared_side"] += 1
    if any(len(T) == 7 for T in side.values()):
        stats["top_block"] += 1
    if any(len(T) == 1 for T in side.values()):
        stats["singleton"] += 1
    if any(len(T) == 0 for T in side.values()):
        stats["empty_side"] += 1
    for i, j in combinations(range(8), 2):
        (a, b), (c, e) = cols[i], cols[j]
        det = a * e - b * c
        forced = prod(n for T, n in bl if i in T and j in T)
        h = gcd(Q, ys[j] - ys[i])
        diff = ys[j] - ys[i]
        assert diff % h == 0
        t = diff // h
        assert t != 0
        assert det == t * forced, (i, j)
        if i == 0:
            assert t == cols[j][1] and forced == 1
        # Astra's displayed identity
        assert (Q * h) % (d[i] * d[j]) == 0
        lastf = Q * h // (d[i] * d[j])
        assert det == t * lastf
        for p, f in fac.items():
            v = 0
            z = lastf
            while z % p == 0:
                z //= p
                v += 1
            both = (i in side[p]) and (j in side[p])
            assert v == (2 * f if both else 0)
        # D_ij identity: atoms where i,j have different signs
        sep = prod(p ** (2 * f) for p, f in fac.items()
                   if ((i in side[p]) or i == 0 and False) != (j in side[p]))
        # sign of slot i at prime p: 1 if opposite anchor
        si = {p: (i in side[p]) for p in fac}
        sj = {p: (j in side[p]) for p in fac}
        sep = prod(p ** (2 * f) for p, f in fac.items() if si[p] != sj[p])
        assert sep * h * h == Q * Q
        # cofactor identity: |t|^2 * Q = (Q/h)^2 * (|diff|^2 / Q) i.e.
        # 2log|t| = D_ij - W/2 + 2log|diff| - log Q  <=> t^2 * h^2 = diff^2
        assert t * t * h * h == diff * diff
        stats["pairs"] += 1
    stats["tuples"] += 1


def main():
    stats = dict(instances=0, tuples=0, pairs=0, shared_side=0, top_block=0,
                 singleton=0, empty_side=0, oneclass_prime=0, coalesce_primes=0,
                 two_adic=0, Q1=0)
    for it in range(4000):
        q, mods, choices, modes, pts = random_instance()
        # random subset X of the marks
        if len(pts) < 8:
            continue
        m = rng.randint(8, min(len(pts), 40))
        X = rng.sample(pts, m)
        xs = X[0]
        g = q
        for x in X:
            g = gcd(g, x - xs)
        Q = q // g
        Y = [(x - xs) // g for x in X]
        assert len(set(Y)) == len(Y)
        fac = factorint(Q) if Q > 1 else {}
        if Q == 1:
            stats["Q1"] += 1
        # local structure
        for p, f in fac.items():
            pf = p ** f
            cls = sorted({y % pf for y in Y})
            assert len(cls) == 2, (q, mods, choices, p, f, cls)
            assert (cls[1] - cls[0]) % p != 0
            if p == 2:
                stats["two_adic"] += 1
        for pe, rs, mode in zip(mods, choices, modes):
            p = min(factorint(pe))
            if p not in fac:
                stats["oneclass_prime"] += 1
            if mode == "coalesce":
                stats["coalesce_primes"] += 1
        stats["instances"] += 1
        for _ in range(6):
            tup = rng.sample(Y, 8)  # ordered, distinct, random anchor
            check_tuple(Q, fac, tup, stats)
        # adversarial tuples: all non-anchor marks in anchor class at one prime,
        # or all opposite, when available
        for p, f in fac.items():
            pf = p ** f
            for want_same in (True, False):
                a = rng.choice(Y)
                pool = [y for y in Y if y != a and (((y - a) % pf == 0) == want_same)]
                if len(pool) >= 7:
                    tup = [a] + rng.sample(pool, 7)
                    check_tuple(Q, fac, tup, stats)
    print(stats)


if __name__ == "__main__":
    main()
