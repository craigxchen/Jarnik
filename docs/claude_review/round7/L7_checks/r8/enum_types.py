# Enumerate 3-uniform hypergraphs on [7] (0..6) with 12 distinct edges, every pair in <= 2 edges,
# up to S_7 isomorphism. Output canonical forms to types12.txt
import itertools, sys
V = range(7)
triples = list(itertools.combinations(V,3))
pairs = list(itertools.combinations(V,2))
perms = list(itertools.permutations(V))
def canon(E):
    best = None
    for p in perms:
        F = tuple(sorted(tuple(sorted((p[a],p[b],p[c]))) for (a,b,c) in E))
        if best is None or F < best: best = F
    return best
# backtracking with pair counts; canonical augmentation is costly; use simple orderly search with
# set of canonical forms at each level (levels up to 12; prune by pair multiplicity)
level = {()}
for k in range(12):
    nxt = set()
    for E in level:
        cnt = {}
        for e in E:
            for pr in itertools.combinations(e,2): cnt[pr] = cnt.get(pr,0)+1
        last = E[-1] if E else None
        for t in triples:
            if t in E: continue
            if any(cnt.get(pr,0) >= 2 for pr in itertools.combinations(t,2)): continue
            F = tuple(sorted(E + (t,)))
            nxt.add(canon(F))
    level = nxt
    print(k+1, len(level), flush=True)
with open(sys.argv[1] if len(sys.argv)>1 else 'types12.txt','w') as f:
    for E in sorted(level): f.write(repr(E)+'\n')
