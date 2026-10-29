"""Theorem R (rigidity of mean-width extremal measures): exact checks.

For finite S in Z^n: MW_n(S) = 2^{1-n} sum_{P subset [n]} w_P(S);  A_t(S) = measures on S killing Poly_{<t};
dim A_t(S) = |S| - rank(Poly_{<t}|_S)  (exact integer rank, python-flint).
Claims checked:
  (R)  dim A_t(S) <= 1 for every integer t >= MW_n(S)   [Theorem R]
  (C)  critical-index lemma: slicing by x_1 and ordering slices by decreasing g_c = MW_{n-1}(S_c),
       #{i : g_(i) + i - 1 = MW_n(S)} <= 1
  (B)  Theorem B: reg(S) <= MW_n(S) (sanity)
Exhaustive families: all subsets of small boxes; random subsets in Z^3, Z^4; zonotope supports in the k=5 slice.
"""
import sys, random, itertools
from fractions import Fraction
from itertools import combinations, product
import flint

def monomials(n, d):
    out = []
    for deg in range(d+1):
        for c in itertools.combinations_with_replacement(range(n), deg):
            e = [0]*n
            for i in c:
                e[i] += 1
            out.append(tuple(e))
    return out

def rank_poly(S, n, d):
    """rank of Poly_{<=d} restricted to S"""
    if d < 0:
        return 0
    mons = monomials(n, d)
    M = flint.fmpz_mat(len(S), len(mons), [ (lambda s, m: eval_mon(s, m))(s, m) for s in S for m in mons])
    return M.rank()

def eval_mon(s, m):
    r = 1
    for si, mi in zip(s, m):
        r *= si**mi
    return r

def reg(S, n):
    S = list(S)
    if len(S) <= 1:
        return 0
    d = 0
    while rank_poly(S, n, d) < len(S):
        d += 1
    return d

def dimA(S, n, t):
    """dim of measures of order >= t on S"""
    S = list(S)
    return len(S) - rank_poly(S, n, t-1)

def MW(S, n):
    S = list(S)
    if n == 0:
        return Fraction(0)
    tot = 0
    for k in range(n+1):
        for P in combinations(range(n), k):
            vals = [sum(s[i] for i in P) for s in S]
            tot += max(vals) - min(vals)
    return Fraction(tot, 2**(n-1))

def slices(S):
    d = {}
    for s in S:
        d.setdefault(s[0], []).append(s[1:])
    return d

def critical_count(S, n):
    K = MW(S, n)
    sl = slices(S)
    g = sorted((MW(v, n-1) for v in sl.values()), reverse=True)
    return sum(1 for i, gi in enumerate(g, start=1) if gi + i - 1 == K)

def check(S, n, stats):
    K = MW(S, n)
    r = reg(S, n)
    stats['B_viol'] += (r > K)
    import math
    t = math.ceil(K)
    dA = dimA(S, n, t)
    stats['R_viol'] += (dA > 1)
    stats['tight'] += (r == K)
    if n >= 1:
        cc = critical_count(S, n)
        stats['C_viol'] += (cc > 1)
    stats['tot'] += 1
    if r == K:
        stats['tight_dim1'] += (dA == 1)

def new_stats():
    return dict(tot=0, B_viol=0, R_viol=0, C_viol=0, tight=0, tight_dim1=0)

if __name__ == '__main__':
    out = sys.stdout
    # exhaustive small boxes
    for n, box in [(1, [range(6)]), (2, [range(3), range(3)]), (2, [range(4), range(3)]), (3, [range(2)]*3),
                   (3, [range(3), range(2), range(2)])]:
        pts = list(product(*box))
        st = new_stats()
        for mask in range(1, 2**len(pts)):
            S = [pts[i] for i in range(len(pts)) if mask >> i & 1]
            check(S, n, st)
        out.write('exhaustive n=%d box=%s: %s\n' % (n, [len(b) for b in box], st)); out.flush()
    # random subsets in Z^3, Z^4
    rng = random.Random(1)
    for n, side, size, trials in [(3, 4, (3, 14), 4000), (4, 3, (3, 16), 3000)]:
        st = new_stats()
        for _ in range(trials):
            m = rng.randint(*size)
            S = list({tuple(rng.randrange(side) for _ in range(n)) for _ in range(m)})
            check(S, n, st)
        out.write('random n=%d side=%d: %s\n' % (n, side, st)); out.flush()
    # k=5 slice: zonotope supports of root-binomial products and random subsets of them (coordinates 1..4 after dropping 5th)
    roots = [(a, b) for a in range(5) for b in range(a+1, 5)]
    def vec(a, b):
        v = [0]*5; v[a] += 1; v[b] -= 1
        return tuple(v[:4])
    st = new_stats(); st2 = new_stats()
    for _ in range(1500):
        k = rng.randint(1, 4)
        gens = [vec(*rng.choice(roots)) for _ in range(k)]
        Z = {tuple([0]*4)}
        for gv in gens:
            Z = Z | {tuple(z[i] + gv[i] for i in range(4)) for z in Z}
        Z = sorted(Z)
        # full zonotope lattice points = sumset; also fill convex hull lattice points? use sumset (vertex sums)
        check(Z, 4, st)
        sub = [z for z in Z if rng.random() < 0.8] or Z[:1]
        check(sub, 4, st2)
    out.write('k=5 zonotope sumsets: %s\n' % st)
    out.write('k=5 random subsets of zonotope sumsets: %s\n' % st2)
