"""c4: N(u) = #{boundary divisors D_T of Mbar_{0,8} with ord_{D_T}(u) != 0} for units u = prod_{i<j<=7}(x_i-x_j)^{e_ij}
(Kapranov chart, anchor 8), sum e = 0; ord_{D_T}(u) = e(T) := sum_{i<j in T} e_ij for T subset [7], 2<=|T|<=6.
Lemma (multilinear Schwartz-Zippel): q(x) = sum e_ij x_i x_j is a nonzero multilinear quadric => >= 2^{7-2} = 32
nonzeros on {0,1}^7, all at 2<=|T|<=6 => N(u) >= 32.  Checks (exact integers):
 (a) cross-ratio (x1-x2)/(x1-x3): N = 32 (bound sharp);
 (b) exhaustive over all e with entries in {-1,0,1} supported on <= 4 edges with sum 0, and 200000 random e: min N >= 32;
 (c) evidence: random search for min N over 'genuine' units (every one of the 8 points, incl. anchor, involved),
     where the anchor 8 is involved iff some vertex degree deg_i = sum_j e_ij != 0.
"""
import itertools, random
edges = [(i, j) for i in range(7) for j in range(i+1, 7)]
Ts = [frozenset(c) for r in range(2, 7) for c in itertools.combinations(range(7), r)]
assert len(Ts) == 119
Tmask = []
for T in Ts: Tmask.append([k for k, (i, j) in enumerate(edges) if i in T and j in T])
def N(e): return sum(1 for m in Tmask if sum(e[k] for k in m) != 0)
def involved(e):
    pts = set()
    for k, (i, j) in enumerate(edges):
        if e[k]: pts.add(i); pts.add(j)
    deg = [sum(e[k] for k, (i, j) in enumerate(edges) if v in (i, j)) for v in range(7)]
    if any(deg): pts.add(7)
    return len(pts)
e = [0]*21; e[edges.index((0, 1))] = 1; e[edges.index((0, 2))] = -1
print("(a) cross-ratio (x1-x2)/(x1-x3): N =", N(e), " points involved:", involved(e))
mn = 10**9; cnt = 0
for r in range(2, 5):
    for supp in itertools.combinations(range(21), r):
        for signs in itertools.product((1, -1), repeat=r):
            if sum(signs) != 0: continue
            e = [0]*21
            for k, s in zip(supp, signs): e[k] = s
            mn = min(mn, N(e)); cnt += 1
print("(b) exhaustive small-support e (<=4 edges, +-1, sum 0):", cnt, "vectors, min N =", mn)
random.seed(1); mn = 10**9
for _ in range(200000):
    e = [random.choice((-2, -1, 0, 0, 0, 1, 2)) for _ in range(21)]
    s = sum(e)
    if s: e[random.randrange(21)] -= s
    if any(e): mn = min(mn, N(e))
print("(b) random e: min N =", mn)
best = {}
random.seed(2)
for _ in range(300000):
    r = random.randint(2, 8)
    e = [0]*21
    for k in random.sample(range(21), r): e[k] = random.choice((-1, 1, -2, 2))
    s = sum(e)
    if s: e[random.randrange(21)] -= s
    if not any(e): continue
    k = involved(e); n_ = N(e)
    if k not in best or n_ < best[k][0]: best[k] = (n_, tuple(e))
for k in sorted(best):
    print(f"(c) evidence: min N found for units involving exactly {k} of the 8 points: {best[k][0]}  "
          f"(method bound gain <= Lambda_k.beta - N = {(16 if k==4 else 32-2**(8-k)) - best[k][0]})")
