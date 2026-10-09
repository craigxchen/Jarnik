"""Combinatorics of REAL balanced curves (exact).

A real rational curve of class beta in Mbar_{0,n} meets each boundary divisor D_S once, transversally,
at a real point.  Its real locus (a circle) therefore traces a closed walk in the dual graph of the
associahedral tiling of Mbar_{0,n}(R): vertices = dihedral orders of [n], edges = reversal of a cyclic
interval S (2 <= |S| <= n-2), labelled by the split {S, S^c}.  Each split must be used exactly once.
This script searches for such closed walks (Euler-type tours) by randomized backtracking.
"""
import itertools, random, sys

def canon(seq):
    n = len(seq)
    best = None
    for s in (list(seq), list(reversed(seq))):
        k = s.index(0)
        r = tuple(s[k:] + s[:k])
        if best is None or r < best:
            best = r
    return best

def moves(order):
    n = len(order)
    out = []
    seen = set()
    for start in range(n):
        for length in range(2, n-1):
            idx = [(start+k) % n for k in range(length)]
            S = frozenset(order[i] for i in idx)
            key = S if 0 not in S else frozenset(range(n)) - S
            if key in seen:
                continue
            seen.add(key)
            new = list(order)
            vals = [order[i] for i in idx][::-1]
            for i, v in zip(idx, vals):
                new[i] = v
            out.append((key, canon(new)))
    return out

def search(n, seed, max_nodes=10**6):
    rng = random.Random(seed)
    allsplits = set()
    for r in range(2, n-1):
        for S in itertools.combinations(range(n), r):
            S = frozenset(S)
            key = S if 0 not in S else frozenset(range(n)) - S
            allsplits.add(key)
    N = len(allsplits)
    start = canon(list(range(n)))
    mcache = {}
    def M(o):
        if o not in mcache:
            mcache[o] = moves(o)
        return mcache[o]
    path = []
    used = set()
    nodes = [0]
    def dfs(o):
        nodes[0] += 1
        if nodes[0] > max_nodes:
            return None
        if len(used) == N:
            return o == start
        opts = [(S, o2) for (S, o2) in M(o) if S not in used]
        rng.shuffle(opts)
        for S, o2 in opts:
            used.add(S); path.append((S, o2))
            r = dfs(o2)
            if r:
                return True
            used.discard(S); path.pop()
            if r is None:
                return None
        return False
    r = dfs(start)
    return r, path, nodes[0], N

if __name__ == '__main__':
    n = int(sys.argv[1]); tries = int(sys.argv[2])
    for seed in range(tries):
        r, path, nn, N = search(n, seed)
        print('n', n, 'seed', seed, 'result', r, 'nodes', nn, 'splits', N, flush=True)
        if r:
            print([tuple(sorted(S)) for S, _ in path])
            break
