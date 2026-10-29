#!/usr/bin/env python3
"""Independent check of stab.md Theorem 3.3 (and construct.md Theorem B).

For finite S in Z^n:
  MW_n(S) = 2^{1-n} sum_{P subset [n]} w_P(S),  w_P = max x(P) - min x(P)
  A_t(S)  = measures c on S with sum_s c_s s^alpha = 0 for all |alpha| < t
  dim A_t(S) = |S| - rank(moment matrix of degree < t)
Claim (Thm 3.3): dim A_t(S) <= 1 for every integer t >= MW_n(S).
Claim (Thm B):   reg(S) <= MW_n(S), where reg(S) = max{t : A_t(S) != 0}.

Certification: rank mod p <= rank over Q, so dim_p >= dim_Q.  dim_p <= 1 certifies dim_Q <= 1.
If dim_p >= 2 we recompute exactly over Q (fmpq_mat rank).
reg_p >= reg_Q, so reg_p <= MW certifies Theorem B.
Written from scratch for the kill report (does not import the author's scripts).
"""
import itertools, random, sys, time
from fractions import Fraction
import flint

P1 = (1 << 61) - 1
P2 = 2305843009213693921  # another prime < 2^61? checked below

def is_prime(n):
    return flint.fmpz(n).is_prime()

if not is_prime(P2):
    P2 = flint.fmpz(1 << 60).next_prime() if hasattr(flint.fmpz, 'next_prime') else 1152921504606846883
assert is_prime(P1) and is_prime(P2), (P1, P2)

def monomials(n, deg_lt):
    out = []
    for d in range(deg_lt):
        for a in itertools.combinations_with_replacement(range(n), d):
            e = [0] * n
            for i in a:
                e[i] += 1
            out.append(tuple(e))
    return out

def moment_rank_mod(S, t, p):
    if t <= 0:
        return 0
    n = len(S[0])
    mons = monomials(n, t)
    rows = []
    for al in mons:
        row = []
        for s in S:
            v = 1
            for si, ai in zip(s, al):
                v = (v * pow(si % p, ai, p)) % p
            row.append(v)
        rows.append(row)
    M = flint.nmod_mat(rows, p)
    return M.rank()

def moment_rank_Q(S, t):
    if t <= 0:
        return 0
    n = len(S[0])
    mons = monomials(n, t)
    rows = []
    for al in mons:
        row = []
        for s in S:
            v = 1
            for si, ai in zip(s, al):
                v *= si ** ai
            row.append(v)
        rows.append(row)
    M = flint.fmpz_mat(rows)
    return M.rank()

def dimA(S, t, exact_if_ge=2):
    """certified upper bound on dim_Q A_t(S): min over two primes; exact if >= exact_if_ge"""
    d = len(S) - max(moment_rank_mod(S, t, P1), moment_rank_mod(S, t, P2))
    if d >= exact_if_ge:
        d = len(S) - moment_rank_Q(S, t)
    return d

def MW(S):
    n = len(S[0])
    tot = 0
    for r in range(n + 1):
        for P in itertools.combinations(range(n), r):
            vals = [sum(s[i] for i in P) for s in S]
            tot += max(vals) - min(vals)
    return Fraction(tot, 2 ** (n - 1))

def reg_upper(S):
    """certified upper bound for reg_Q(S): least d with rank_p(deg<=d) = |S| (rank_p<=rank_Q)."""
    d = 0
    while True:
        r = max(moment_rank_mod(S, d + 1, P1), moment_rank_mod(S, d + 1, P2))
        if r == len(S):
            return d
        d += 1

def check_set(S, stats):
    S = sorted(set(S))
    if len(S) == 0:
        return
    K = MW(S)
    stats['sets'] += 1
    # Theorem B (certified)
    r = reg_upper(S)
    if r > K:
        # recheck exactly
        n_ok = len(S) - moment_rank_Q(S, int(K) + 1) if K == int(K) else None
        stats['B_fail'].append((S, K, r))
    # Theorem 3.3 at t = ceil(K)
    t = -(-K.numerator // K.denominator)  # ceil
    d = dimA(S, t)
    if d >= 2:
        stats['T33_fail'].append((S, K, t, d))
    if d == 1:
        stats['tight'] += 1
        if K.denominator != 1:
            stats['nonint_tight'].append((S, K))

def boxes_exhaustive(dims, stats, label):
    pts = list(itertools.product(*[range(a + 1) for a in dims]))
    N = len(pts)
    t0 = time.time()
    for mask in range(1, 1 << N):
        S = [pts[i] for i in range(N) if mask >> i & 1]
        check_set(S, stats)
    print(f"{label}: box {dims} all {2**N-1} subsets, {time.time()-t0:.1f}s", flush=True)

def random_sets(n, box, size_rng, count, stats, label, seed):
    rnd = random.Random(seed)
    t0 = time.time()
    for _ in range(count):
        k = rnd.randint(*size_rng)
        S = set()
        while len(S) < k:
            S.add(tuple(rnd.randint(0, box) for _ in range(n)))
        check_set(list(S), stats)
    print(f"{label}: n={n} box [0,{box}] {count} random sets size {size_rng}, {time.time()-t0:.1f}s", flush=True)

def minkowski(A, B):
    return list({tuple(a + b for a, b in zip(x, y)) for x in A for y in B})

def random_zonotope_unions(n, count, stats, label, seed):
    """Adversarial: tight sets (zonotopes of root segments [0,e_i],[e_i,e_j]) and unions / perturbations
    of two zonotopes with the same MW; these are where a second tight measure would live."""
    rnd = random.Random(seed)
    segs = [((0,) * n, tuple(1 if k == i else 0 for k in range(n))) for i in range(n)]
    segs += [(tuple(1 if k == i else 0 for k in range(n)), tuple(1 if k == j else 0 for k in range(n)))
             for i in range(n) for j in range(i + 1, n)]
    t0 = time.time()
    for _ in range(count):
        def zono():
            Z = [(0,) * n]
            for _ in range(rnd.randint(1, 4)):
                a, b = rnd.choice(segs)
                Z = minkowski(Z, [a, b])
            return Z
        Z1 = zono()
        mode = rnd.randint(0, 3)
        if mode == 0:
            S = Z1
        elif mode == 1:
            Z2 = zono()
            S = list(set(Z1) | set(Z2))
        elif mode == 2:  # add one extra point near
            S = list(set(Z1) | {tuple(c + rnd.randint(-1, 1) for c in rnd.choice(Z1))})
        else:  # remove a point
            S = list(Z1)
            if len(S) > 1:
                S.pop(rnd.randrange(len(S)))
        check_set(S, stats)
    print(f"{label}: n={n} {count} zonotope-derived sets, {time.time()-t0:.1f}s", flush=True)

def cut_closure_tests(n, box, count, stats, label, seed):
    """Cut closure T of random S: all x in box-hull with min_S x(P) <= x(P) <= max_S x(P) for all P.
    Same widths as S, maximal set with them; the natural place for a 2-dim tight kernel."""
    rnd = random.Random(seed)
    Ps = [P for r in range(1, n + 1) for P in itertools.combinations(range(n), r)]
    t0 = time.time()
    for _ in range(count):
        k = rnd.randint(2, 5)
        S = list({tuple(rnd.randint(0, box) for _ in range(n)) for _ in range(k)})
        lo = {P: min(sum(s[i] for i in P) for s in S) for P in Ps}
        hi = {P: max(sum(s[i] for i in P) for s in S) for P in Ps}
        T = [x for x in itertools.product(range(box + 1), repeat=n)
             if all(lo[P] <= sum(x[i] for i in P) <= hi[P] for P in Ps)]
        check_set(T, stats)
    print(f"{label}: n={n} {count} cut closures in [0,{box}]^{n}, {time.time()-t0:.1f}s", flush=True)

if __name__ == '__main__':
    stats = dict(sets=0, tight=0, B_fail=[], T33_fail=[], nonint_tight=[])
    boxes_exhaustive((4, 1), stats, 'E1')            # 10 points
    boxes_exhaustive((1, 1, 1), stats, 'E2')         # 8 points (also in author's list; cheap)
    boxes_exhaustive((2, 1, 1), stats, 'E3')         # 12 points (author used [0,2]x{0,1}^2 too)
    boxes_exhaustive((5, 1), stats, 'E4')            # 12 points
    boxes_exhaustive((1, 1, 1, 1), stats, 'E5')      # 16 points: 65535 subsets
    print('after exhaustive:', stats['sets'], 'sets', stats['tight'], 'tight',
          len(stats['B_fail']), 'B fails', len(stats['T33_fail']), 'T3.3 fails', flush=True)
    random_sets(5, 2, (2, 14), 1500, stats, 'R5', 1)
    random_sets(6, 1, (2, 14), 800, stats, 'R6', 2)
    random_sets(3, 4, (2, 16), 2000, stats, 'R3', 3)
    random_zonotope_unions(3, 1500, stats, 'Z3', 4)
    random_zonotope_unions(4, 1000, stats, 'Z4', 5)
    random_zonotope_unions(5, 300, stats, 'Z5', 6)
    cut_closure_tests(3, 3, 400, stats, 'C3', 7)
    cut_closure_tests(4, 2, 200, stats, 'C4', 8)
    print('TOTAL sets', stats['sets'], 'tight(dim A_ceilMW = 1)', stats['tight'])
    print('Theorem B failures', len(stats['B_fail']), stats['B_fail'][:3])
    print('Theorem 3.3 failures', len(stats['T33_fail']), stats['T33_fail'][:3])
    print('tight with non-integer MW (impossible if Thm B):', len(stats['nonint_tight']))
