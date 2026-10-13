"""zeta(A) = max sum_i L_i over lattice zonotopes k0 + sum_i [0,L_i] v_i (v_i nonzero primitive integer
vectors, pairwise non-parallel, L_i >= 1, any number of generators) whose lattice points k0 + sum t_i v_i
(0 <= t_i <= L_i) all lie in A = Q cap Z^d, Q = {|k_x| <= T, |k(S)| <= L}.
Exact: Q is convex and given by the constraint functionals lam with bound b_lam, so Z subset Q iff for every lam:
  lam(k0) + sum_i L_i min(0, lam(v_i)) >= -b  and  lam(k0) + sum_i L_i max(0, lam(v_i)) <= b.
Then every point k0 + sum t_i v_i lies in Q (convexity) and is a lattice point, hence in A.
Search: DFS over generators in a fixed order with width pruning sum_i L_i |lam(v_i)| <= 2 b_lam, then
existence of k0 checked by enumerating k0 in the box.  Non-primitive generators are subsumed (v = g v',
L -> gL).  Parallel generators are subsumed (merge).  So the search is exhaustive."""
import itertools, math, sys
sys.setrecursionlimit(1000000)
from core import build_A

def functionals(n, T, prof, L):
    d = n - 1
    lams = []
    for x in range(d):
        lams.append((tuple(1 if y == x else 0 for y in range(d)), T))
    for S in prof:
        lam = tuple((1 if x in S else 0) - (1 if (n - 1) in S else 0) for x in range(d))
        lams.append((lam, L))
    return lams

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))

def prim_vectors(d, R):
    out = []
    for v in itertools.product(range(-R, R + 1), repeat=d):
        if not any(v):
            continue
        if math.gcd(*[abs(x) for x in v]) != 1:
            continue
        if next(x for x in v if x) < 0:
            continue
        out.append(v)
    return out

def zeta_search(n, T, prof, L, target, d=None, lams=None):
    """return a zonotope with sum L_i >= target inside A, or None (exhaustive).
    General form: pass d, lams = [(lam, b)] (must include the coordinate functionals with bound T)."""
    if lams is None:
        d = n - 1
        lams = functionals(n, T, prof, L)
    vecs = prim_vectors(d, 2 * T)
    # per-vector cost on each functional
    cost = {v: [abs(dot(l, v)) for l, b in lams] for v in vecs}
    caps = [2 * b for l, b in lams]
    vecs = [v for v in vecs if all(c <= cap for c, cap in zip(cost[v], caps))]
    # order: cheap vectors first
    vecs.sort(key=lambda v: sum(cost[v]))
    nv = len(vecs)
    box = list(itertools.product(range(-T, T + 1), repeat=d))
    found = [None]

    def feasible(gens):
        # exists k0 with all constraints
        for k0 in box:
            ok = True
            for (l, b) in lams:
                base = dot(l, k0)
                lo = base + sum(Li * min(0, dot(l, v)) for v, Li in gens)
                hi = base + sum(Li * max(0, dot(l, v)) for v, Li in gens)
                if lo < -b or hi > b:
                    ok = False
                    break
            if ok:
                return k0
        return None

    def dfs(i, gens, used, total):
        if found[0] is not None:
            return
        if total >= target:
            k0 = feasible(gens)
            if k0 is not None:
                found[0] = (k0, list(gens))
            return
        if i >= nv:
            return
        # optimistic bound (valid): every generator has a nonzero coordinate, so
        # sum_i L_i <= sum_x (coordinate budget left); the first d functionals are the coordinates.
        if total + sum(caps[x] - used[x] for x in range(d)) < target:
            return
        v = vecs[i]
        c = cost[v]
        maxL = min([(cap - u) // cc for cc, cap, u in zip(c, caps, used) if cc > 0] + [target - total])
        for Lv in range(maxL, 0, -1):
            nu = [u + Lv * cc for u, cc in zip(used, c)]
            gens.append((v, Lv))
            dfs(i + 1, gens, nu, total + Lv)
            gens.pop()
            if found[0] is not None:
                return
        dfs(i + 1, gens, used, total)

    dfs(0, [], [0] * len(lams), 0)
    return found[0]

def zeta(n, T, prof, L, upper, d=None, lams=None):
    """largest s <= upper with a zonotope of total s (monotone: sub-zonotopes)."""
    for s in range(upper, 0, -1):
        z = zeta_search(n, T, prof, L, s, d=d, lams=lams)
        if z is not None:
            return s, z
    return 0, None

if __name__ == "__main__":
    import ast
    n = int(sys.argv[1]); T = int(sys.argv[2]); L = int(sys.argv[3]); prof = ast.literal_eval(sys.argv[4])
    up = int(sys.argv[5])
    print(zeta(n, T, prof, L, up))
