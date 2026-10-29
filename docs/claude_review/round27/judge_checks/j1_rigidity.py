# Judge check J1 (own code): stab Theorem 3.3  -- dim A_t(S) <= 1 for every integer t >= MW_n(S).
# MW_n(S) = 2^{1-n} sum_{P subset [n]} (max_S x(P) - min_S x(P)),  x(P) = sum_{i in P} x_i.
# A_t(S) = measures c on S with sum_x c_x x^alpha = 0 for all |alpha| < t.  Exact ranks via python-flint fmpq_mat.
import itertools, random, sys, math
from fractions import Fraction
import flint
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)

def MW(S, n):
    tot = 0
    for r in range(n+1):
        for P in itertools.combinations(range(n), r):
            vals = [sum(x[i] for i in P) for x in S]
            tot += max(vals) - min(vals)
    return Fraction(tot, 2**(n-1))

def monos(n, deg):
    out = []
    for d in range(deg+1):
        for c in itertools.combinations_with_replacement(range(n), d):
            a = [0]*n
            for i in c: a[i] += 1
            out.append(tuple(a))
    return out

def dimA(S, n, t):
    if t <= 0: return len(S)
    ms = monos(n, t-1)
    rows = []
    for a in ms:
        rows.append([math.prod(x[i]**a[i] for i in range(n)) for x in S])
    M = flint.fmpq_mat(len(rows), len(S), [v for r in rows for v in r])
    return len(S) - M.rank()

def binom_product_support(n, gens):
    # support of prod (Y^v - 1): Minkowski sum of {0, v}; cancellation is ignored (we take the full sumset, a superset)
    pts = {tuple([0]*n)}
    for v in gens:
        pts = pts | {tuple(p[i]+v[i] for i in range(n)) for p in pts}
    return sorted(pts)

def roots(n):
    R = []
    for i in range(n):
        e = [0]*n; e[i] = 1; R.append(tuple(e))
    for i in range(n):
        for j in range(i+1, n):
            e = [0]*n; e[i] = 1; e[j] = -1; R.append(tuple(e))
    return R

stats = dict(sets=0, tight=0, viol=0, nonzero_at_ceil=0)
def test(S, n):
    S = sorted(set(S))
    mw = MW(S, n)
    t0 = math.ceil(mw)
    d = dimA(S, n, t0)
    stats['sets'] += 1
    if d >= 1: stats['nonzero_at_ceil'] += 1
    if d >= 1 and mw == t0: stats['tight'] += 1
    if d > 1:
        stats['viol'] += 1
        print('VIOLATION', n, S, mw, d, flush=True)
    # Theorem B sanity: no nonzero measure of order > MW
    if dimA(S, n, t0+1) != 0 and t0+1 > mw:
        print('THM B VIOLATION', n, S, mw, flush=True); stats['viol'] += 1

for n in (1, 2, 3, 4):
    R = roots(n)
    N = {1: 300, 2: 600, 3: 400, 4: 120}[n]
    for it in range(N):
        k = random.randint(1, {1:6, 2:5, 3:4, 4:3}[n])
        gens = [random.choice(R) for _ in range(k)]
        S = binom_product_support(n, gens)
        mode = random.random()
        if mode < 0.4:
            pass
        elif mode < 0.8:
            # add 1-2 random nearby points
            for _ in range(random.randint(1, 2)):
                S.append(tuple(random.randint(-1, 3) for _ in range(n)))
        else:
            # remove one point and add one
            if len(S) > 2: S.pop(random.randrange(len(S)))
            S.append(tuple(random.randint(-1, 2) for _ in range(n)))
        test(S, n)
    # purely random small sets
    for it in range(N//2):
        m = random.randint(2, 9)
        S = [tuple(random.randint(0, 3) for _ in range(n)) for _ in range(m)]
        test(S, n)
    print('n =', n, stats, flush=True)
print('FINAL', stats)
