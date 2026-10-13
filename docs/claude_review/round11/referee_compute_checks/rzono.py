# Referee's zonotope witness search and generic slab-box utilities (independent of the author's zono.py).
import itertools, math
import flint

def functionals(n, profile, T, L, base):
    """constraint functionals on k' (coords x != base, in increasing order) with bounds.
    profile: list of frozensets S (subsets of range(n))."""
    others = [x for x in range(n) if x != base]
    d = len(others)
    cons = []
    for i in range(d):
        lam = [0] * d; lam[i] = 1
        cons.append((tuple(lam), T))
    for S in profile:
        lam = [(1 if x in S else 0) - (1 if base in S else 0) for x in others]
        cons.append((tuple(lam), L))
    return others, cons

def points(d, cons, T):
    out = []
    for kp in itertools.product(range(-T, T + 1), repeat=d):
        if all(abs(sum(a * b for a, b in zip(lam, kp))) <= b0 for lam, b0 in cons):
            out.append(kp)
    return out

def rank_poly(A, m, T, p=None):
    d = len(A[0])
    mons = [al for al in itertools.product(range(min(2 * T, m - 1) + 1), repeat=d) if sum(al) < m]
    rows = [[math.prod(k[i] ** al[i] for i in range(d)) for k in A] for al in mons]
    if p:
        return flint.nmod_mat(len(rows), len(A), [v % p for r in rows for v in r], p).rank()
    return flint.fmpz_mat(rows).rank()

def rho(A, T, exact=True):
    """largest m with dim V(A,m) > 0, i.e. rank Poly_{<m}|_A < |A|."""
    m = 1
    while True:
        r = rank_poly(A, m + 1, T, None if exact else 2**61 - 1)
        if r == len(A):
            return m
        m += 1

def prim_gens(d, E):
    gens = []
    for v in itertools.product(range(-E, E + 1), repeat=d):
        if not any(v): continue
        first = next(t for t in v if t)
        if first < 0: continue
        g = 0
        for t in v: g = math.gcd(g, abs(t))
        if g != 1: continue
        gens.append(v)
    gens.sort(key=lambda v: (sum(map(abs, v)), v))
    return gens

def find_zonotope(cons, d, target, E=2, node_limit=2_000_000):
    """search a lattice zonotope k0 + sum [0,L_i] v_i inside the polytope {|lam(k)| <= b} with sum L_i >= target.
    Returns (k0, [(v, L)]) or None (None = not found within limits; NOT a proof of absence)."""
    gens = prim_gens(d, E)
    lv = [[sum(a * b for a, b in zip(lam, v)) for lam, _ in cons] for v in gens]
    budget0 = [2 * b for _, b in cons]
    box_idx = list(range(d))   # the first d constraints are the coordinate box
    nodes = [0]
    chosen = []
    def k0_exists(chosen):
        # need integer k0 with lo_j <= lam_j(k0) <= hi_j for all j
        lo = []; hi = []
        for j, (lam, b) in enumerate(cons):
            pos = sum(Lc * max(0, lv[g][j]) for g, Lc in chosen)
            neg = sum(Lc * min(0, lv[g][j]) for g, Lc in chosen)
            lo.append(-b - neg); hi.append(b - pos)
        # box constraints give k0 ranges directly
        rng = [range(lo[i], hi[i] + 1) for i in range(d)]
        for k0 in itertools.product(*rng):
            if all(lo[j] <= sum(a * c for a, c in zip(cons[j][0], k0)) <= hi[j] for j in range(len(cons))):
                return k0
        return None
    def dfs(start, budget, total):
        nodes[0] += 1
        if nodes[0] > node_limit:
            return None
        if total >= target:
            k0 = k0_exists(chosen)
            if k0 is not None:
                return (k0, [(gens[g], Lc) for g, Lc in chosen])
            return None
        # upper bound on what can still be added: box budgets
        if total + sum(budget[i] for i in box_idx) < target:
            return None
        for g in range(start, len(gens)):
            costs = lv[g]
            Lmax = min((budget[j] // abs(c) for j, c in enumerate(costs) if c), default=0)
            for Lc in range(min(Lmax, target - total), 0, -1):
                nb = [budget[j] - Lc * abs(costs[j]) for j in range(len(costs))]
                chosen.append((g, Lc))
                r = dfs(g + 1, nb, total + Lc)
                chosen.pop()
                if r is not None:
                    return r
        return None
    return dfs(0, budget0, 0)
