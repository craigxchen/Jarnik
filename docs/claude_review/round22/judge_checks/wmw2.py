"""Judge check J3b (evidence only): weighted mean-width inequality (W) of toric.md open point 1,
    reg(S) <= (1/6) sum_{|A|=2} width_S(1_A),
tested on lattice points S of lattice polytopes in M = {sum x = 0} subset Z^5 (where extremal sections live).
Modes: 'poly' random polytopes, 'zono' random zonotopes, 'cut' half-space pieces of m R_5, 'rootzono' zonotopes of
roots e_a - e_b (where (W) is an equality for products of differences).
reg computed by rank mod p (flint nmod_mat): an upper bound for the rational reg. A violation must be re-certified
over Q (lower bound needs an explicit rational kernel vector); no violation = evidence only.
Also checks the refereed Thm B bound (sum over all 32 subsets / 16) as a sanity test of the code.
"""
import itertools, random, sys
from fractions import Fraction
import numpy as np
from scipy.spatial import Delaunay
import flint

P = 2147483629

def monos(d, n=4):
    return [e for e in itertools.product(range(d + 1), repeat=n) if sum(e) <= d]

def reg(S):
    pts = [s[:4] for s in S]
    n = len(pts); d = 0
    while True:
        cols = monos(d)
        if len(cols) >= n:
            M = flint.nmod_mat(n, len(cols), [ (p[0]**e[0] * p[1]**e[1] * p[2]**e[2] * p[3]**e[3]) % P for p in pts for e in cols], P)
            if M.rank() == n: return d
        d += 1

subsets2 = list(itertools.combinations(range(5), 2))
allA = [A for r in range(0, 6) for A in itertools.combinations(range(5), r)]

def width(S, A):
    vals = [sum(s[a] for a in A) for s in S]
    return max(vals) - min(vals)

def lattice_points(verts):
    V = np.array([v[:4] for v in verts], dtype=float)
    if np.linalg.matrix_rank(V[1:] - V[0]) < 4: return None
    tri = Delaunay(V)
    lo = V.min(axis=0).astype(int); hi = V.max(axis=0).astype(int)
    pts = []
    for x in itertools.product(*[range(lo[i], hi[i] + 1) for i in range(4)]):
        if tri.find_simplex(np.array(x, dtype=float), tol=1e-9) >= 0:
            pts.append(tuple(list(x) + [-sum(x)]))
    return pts

def rand_vec(c):
    while True:
        x = [random.randint(-c, c) for _ in range(4)]
        if abs(sum(x)) <= c and any(x): return tuple(x + [-sum(x)])

roots = [tuple((1 if j == a else -1 if j == b else 0) for j in range(5)) for a in range(5) for b in range(5) if a != b]

def gen(mode):
    if mode == 'poly':
        c = random.choice([1, 2, 2, 3]); k = random.randint(5, 9)
        return lattice_points([rand_vec(c) for _ in range(k)])
    if mode in ('zono', 'rootzono'):
        j = random.randint(4, 7)
        vs = [random.choice(roots) if mode == 'rootzono' else rand_vec(random.choice([1, 1, 2])) for _ in range(j)]
        verts = set()
        for mask in itertools.product([0, 1], repeat=j):
            verts.add(tuple(sum(m * v[i] for m, v in zip(mask, vs)) for i in range(5)))
        return lattice_points(list(verts))
    if mode == 'cut':
        m = random.choice([2, 3])
        pts = [p for p in itertools.product(range(-m, m + 1), repeat=4)]
        S = [tuple(list(p) + [-sum(p)]) for p in pts if sum(abs(t) for t in p) + abs(sum(p)) <= 2 * m]
        a = rand_vec(2); t = random.randint(-m, m)
        S2 = [s for s in S if sum(x * y for x, y in zip(a, s)) <= t]
        return S2 if len(S2) >= 5 else None

def main(mode, seed, trials):
    random.seed(seed)
    worst = Fraction(0); worstS = None; viol = 0; bviol = 0; tight = 0; done = 0
    while done < trials:
        S = gen(mode)
        if not S or len(S) < 5 or len(S) > 400: continue
        done += 1
        r = reg(S)
        b = Fraction(sum(width(S, A) for A in subsets2), 6)
        tb = Fraction(sum(width(S, A) for A in allA), 16)
        if r > tb: bviol += 1; print("THM B VIOLATION (code bug?)", r, tb, len(S))
        if r > b: viol += 1; print("W VIOLATION", r, b, len(S), sorted(S)[:5])
        if r == b: tight += 1
        if Fraction(r) / b > worst: worst = Fraction(r) / b; worstS = (r, b, len(S))
        sys.stdout.flush()
    print(f"mode {mode} seed {seed}: {trials} sets, W-violations {viol}, ThmB-violations {bviol}, tight {tight}, max reg/bound {worst} ({float(worst):.4f}) at (reg, bound, |S|) = {worstS}")

if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
