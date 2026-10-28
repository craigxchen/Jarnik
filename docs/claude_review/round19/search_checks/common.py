# Shared combinatorics for the round-19 abc-gain search (unit pencils lam*A/C on Mbar_{0,8}).
# Conventions copied from round18/referee_theory_checks/ref_gain.py: labels 0..6 = x_1..x_7, edges (i,j), i<j,
# factor x_i - x_j; visible T: pairs (21) and Ts = subsets of size 3..6 (98).
import itertools, json, os
edges = [(i, j) for i in range(7) for j in range(i + 1, 7)]
Ts = [frozenset(c) for r in range(3, 7) for c in itertools.combinations(range(7), r)]
assert len(edges) == 21 and len(Ts) == 98

def canon(ms):
    """canonical representative of an edge multiset under S_7 (brute force relabelling, as ref_gain.py)."""
    verts = sorted(set(v for e in ms for v in e)); k = len(verts); best = None
    for perm in itertools.permutations(range(k)):
        mp = {verts[i]: perm[i] for i in range(k)}
        t = tuple(sorted(tuple(sorted((mp[a], mp[b]))) for (a, b) in ms))
        if best is None or t < best: best = t
    return best

def orbit_reps(D, cache_dir=None):
    fn = None
    if cache_dir:
        fn = os.path.join(cache_dir, 'reps_D%d.json' % D)
        if os.path.exists(fn):
            return [tuple(tuple(e) for e in r) for r in json.load(open(fn))]
    reps = sorted(set(canon(ms) for ms in itertools.combinations_with_replacement(edges, D)))
    if fn: json.dump(reps, open(fn, 'w'))
    return reps

def nu_unit(ms, T):
    return sum(1 for (a, b) in ms if a in T and b in T)

def N0(A, C):
    """#{T visible : a_T != c_T}; pairs: e in A or C (A, C disjoint)."""
    n = sum(1 for T in Ts if nu_unit(A, T) != nu_unit(C, T))
    n += len(set(A) | set(C))
    return n

def red_init(ms, T):
    """initial form of prod_{e in ms}(x_i - x_j) along L_T, |T|>=3, i0 = min T, local coords y_i = (x_i-x_i0)/eps.
    Returned as (sign, sorted multiset of normalised factors). Factor labels: ('y',i) for y_i (y_{i0} = 0, written
    ('y',i0)), ('x',i) for x_i.  x_i - x_j -> inside T: y_i - y_j ; i in T only: x_{i0} - x_j ; j in T only: x_i - x_{i0};
    neither: x_i - x_j.  A factor u - v is normalised to (min,max) with a sign."""
    i0 = min(T); sg = 1; out = []
    for (i, j) in ms:
        if i in T and j in T: u, v = ('y', i), ('y', j)
        elif i in T: u, v = ('x', i0), ('x', j)
        elif j in T: u, v = ('x', i), ('x', i0)
        else: u, v = ('x', i), ('x', j)
        if u > v: u, v = v, u; sg = -sg
        out.append((u, v))
    return sg, tuple(sorted(out))
