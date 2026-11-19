"""(b) The symmetry exclusion: integer prime assignments to the trace (Walsh) incidence table.

1. For random assignments of distinct split primes (p = 1 mod 4) to the cuts x in F_q^*, compute
   n_a, g_ab, d_ab = n_a n_b / g_ab^2 as integers and check d_ab depends only on a+b.
2. Translation by u on S_u = F_q^* \\ {u}: unoriented cuts preserved, d preserved, exponent array
   is the original one with the cuts at {x : Tr(ux)=1} complemented.
3. Pair separation for unequal weights: min_c sum_{Tr(cx)=1} log p_x vs W/2 (random weights;
   and the averaging identity sum_c = (q/2) W).
4. Group theory: elementary abelian 2-subgroups of PGL_2 have order <= 4. Symbolic centraliser of
   diag(1,-1) in PGL_2 over a field of char 0 (sympy), and brute force in PGL_2(F_p), p = 3,5,7,11,13.
"""
import math
import random
from itertools import combinations, product
import sympy as sp
from gf2d import GF

random.seed(31)


def split_primes(k, lo):
    out, p = [], lo
    while len(out) < k:
        p = sp.nextprime(p)
        if p % 4 == 1:
            out.append(p)
    return out


def part1_2(d, trials=3):
    F = GF(d)
    q = F.q
    nz = list(range(1, q))
    for trial in range(trials):
        lo = random.choice([10, 1000, 10 ** 6])
        pr = split_primes(4 * q, lo)
        random.shuffle(pr)
        px = dict(zip(nz, pr[:q - 1]))
        marked = {a: {x for x in nz if F.tr(F.mul(a, x)) == 0} for a in nz}
        n = {a: math.prod(px[x] for x in marked[a]) for a in nz}
        Nn = math.prod(px.values())
        dd = {}
        for a, b in combinations(nz, 2):
            g = math.gcd(n[a], n[b])
            assert g == math.prod(px[x] for x in marked[a] & marked[b])
            d_ab = n[a] * n[b] // g ** 2
            key = a ^ b
            if key in dd:
                assert dd[key] == d_ab
            else:
                dd[key] = d_ab
            assert d_ab == math.prod(px[x] for x in nz if F.tr(F.mul(key, x)) == 1)
        # translation by u
        for u in nz:
            S = [a for a in nz if a != u]
            for x in nz:
                cut = frozenset(a for a in S if a in marked[x] if False)  # placeholder
            for x in nz:
                side = frozenset(a for a in S if x in marked[a])
                side_t = frozenset(a for a in S if x in marked[a ^ u])
                comp = frozenset(S) - side
                if F.tr(F.mul(u, x)) == 0:
                    assert side_t == side
                else:
                    assert side_t == comp
            for a, b in combinations(S, 2):
                assert dd[(a ^ u) ^ (b ^ u)] == dd[a ^ b]
        # separation with these (unequal) weights
        W = sum(math.log(p) for p in px.values())
        sepmin = min(sum(math.log(px[x]) for x in nz if F.tr(F.mul(c, x)) == 1) for c in nz)
        tot = sum(sum(math.log(px[x]) for x in nz if F.tr(F.mul(c, x)) == 1) for c in nz)
        assert abs(tot - (q / 2) * W) < 1e-6 * W
        print(f"  q={q} trial {trial} (primes > {lo}): d_ab depends only on a+b; translations preserve "
              f"unoriented cuts and d; min separation/W = {sepmin / W:.4f} (> 1/2: {sepmin > W / 2})")


def adversarial_sep(d):
    """Can unequal weights break pair separation? Put large weight on a hyperplane complement."""
    F = GF(d)
    q = F.q
    nz = list(range(1, q))
    c0 = 1
    w = {x: (1.0 if F.tr(F.mul(c0, x)) == 1 else 3.0) for x in nz}
    W = sum(w.values())
    s = sum(w[x] for x in nz if F.tr(F.mul(c0, x)) == 1)
    print(f"  q={q}: weights 3 on ker Tr(c0 .), 1 elsewhere: separation of pairs with a+b=c0 is {s / W:.3f} W < W/2,"
          f" so pair separation is a genuine extra hypothesis")


def centraliser_symbolic():
    a, b, c, dd, lam = sp.symbols('a b c d lam')
    M = sp.Matrix([[a, b], [c, dd]])
    D = sp.diag(1, -1)
    sols = sp.solve(list(M * D - lam * D * M), [a, b, c, dd], dict=True)
    print("  symbolic: M D = lam D M  solutions:", sols, "; det(M D - lam D M) system forces lam^2=1 for invertible M")
    # lam = 1 -> diagonal; lam = -1 -> antidiagonal
    for l in (1, -1):
        s = sp.solve(list(M * D - l * D * M), [a, b, c, dd], dict=True)
        print(f"   lam={l}: {s}")
    # involutions among diagonal diag(x,1): (diag(x,1))^2 scalar iff x^2 = 1
    x = sp.symbols('x')
    print("   diag(x,1)^2 scalar iff", sp.solve(x ** 2 - 1, x), "; antidiag [[0,y],[1,0]]^2 = y I always")


def pgl2_brute(p):
    els = []
    seen = set()
    for a, b, c, d in product(range(p), repeat=4):
        if (a * d - b * c) % p == 0:
            continue
        # normalise projectively: first nonzero entry = 1
        v = (a, b, c, d)
        f = next(e for e in v if e)
        inv = pow(f, p - 2, p)
        v = tuple(e * inv % p for e in v)
        if v not in seen:
            seen.add(v)
            els.append(v)
    assert len(els) == p * (p * p - 1)

    def mul(X, Y):
        a, b, c, d = X; e, f, g, h = Y
        v = ((a * e + b * g) % p, (a * f + b * h) % p, (c * e + d * g) % p, (c * f + d * h) % p)
        fz = next(t for t in v if t)
        inv = pow(fz, p - 2, p)
        return tuple(t * inv % p for t in v)
    I = (1, 0, 0, 1)
    inv2 = [X for X in els if X != I and mul(X, X) == I]
    # max set of pairwise commuting involutions
    best = 1
    for X, Y in combinations(inv2, 2):
        if mul(X, Y) == mul(Y, X):
            best = max(best, 2)
            for Z in inv2:
                if Z not in (X, Y) and mul(X, Z) == mul(Z, X) and mul(Y, Z) == mul(Z, Y):
                    if Z != mul(X, Y):
                        best = max(best, 4)  # would generate order >= 8
                        return p, len(els), 'FOUND order >= 8 elementary abelian 2-group'
    return p, len(els), f"max elementary abelian 2-subgroup order = {2 ** best if best < 2 else 4}"


if __name__ == '__main__':
    print("Part 1-2: integer assignments")
    for d in (3, 4, 5, 6):
        part1_2(d)
    for d in (3, 5):
        adversarial_sep(d)
    print("Part 4: group theory")
    centraliser_symbolic()
    for p in (3, 5, 7, 11, 13):
        print("  PGL_2(F_%d):" % p, pgl2_brute(p))
