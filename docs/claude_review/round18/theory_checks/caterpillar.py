"""Exact check of Lemma 1.1 (caterpillar lemma), k = 5..8.
For every tree on leaves {0..k-1} (enumerated as subsets of internal edges of all
binary trees), compare:
 (a) for anchor a: the visible sides (sides not containing a, size 2..k-2) form a union of
     two chains (under inclusion) with disjoint unions;
 (b) caterpillar: one can choose one side per split so that the chosen sides form a chain.
Claim: (a) holds for every anchor a  <=>  (a) holds for anchor 0  <=>  (b).
"""
import itertools, sys

def binary_trees(leaves):
    # returns list of trees, each as frozenset of splits (frozenset side not containing leaves[0])
    leaves = list(leaves)
    if len(leaves) == 3:
        return [frozenset()]
    out = []
    first = leaves[:-1]; new = leaves[-1]
    for t in binary_trees(first):
        # edges of t: internal splits + pendant edges; insert new leaf on any edge
        full = set(t) | {frozenset([l]) for l in first if l != first[0]} | {frozenset(first[1:])}
        # represent each edge by side not containing first[0]; pendant edge of first[0] is side first[1:]
        for e in full:
            # subdividing edge e, attaching new leaf on the side e
            ns = set()
            for s in t:
                if e <= s and s != e:
                    ns.add(s | {new})
                elif s == e:
                    ns.add(s | {new}); ns.add(s)
                else:
                    ns.add(s)
            if e not in t:  # pendant or root edge
                if len(e) == 1:
                    ns.add(e | {new})
                # root-pendant edge (side = all but first[0]): new split {first[0], new} vs rest -> side = rest
            # nontrivial splits only, sides not containing first[0]
            allv = set(first) | {new}
            fixed = set()
            for s in ns:
                s = frozenset(s)
                if first[0] in s: s = frozenset(allv - s)
                if 2 <= len(s) <= len(allv) - 2: fixed.add(s)
            if e == frozenset(first[1:]):
                fixed.add(frozenset(allv - {first[0], new}))
            out.append(frozenset(fixed))
    return list(set(out))

def chains_ok(splits, anchor, n):
    allv = frozenset(range(n))
    vis = [s if anchor not in s else allv - s for s in splits]
    vis = [s for s in vis if 2 <= len(s) <= n - 2]
    # try to 2-color into two chains with disjoint unions
    vis.sort(key=len)
    m = len(vis)
    for mask in range(1 << m):
        A = [vis[i] for i in range(m) if mask >> i & 1]
        B = [vis[i] for i in range(m) if not mask >> i & 1]
        ok = True
        for C in (A, B):
            for i in range(len(C) - 1):
                if not C[i] <= C[i + 1]: ok = False; break
            if not ok: break
        if not ok: continue
        uA = frozenset().union(*A) if A else frozenset()
        uB = frozenset().union(*B) if B else frozenset()
        if uA & uB: continue
        return True
    return False

def caterpillar(splits, n):
    allv = frozenset(range(n))
    sp = list(splits)
    for choice in itertools.product([0, 1], repeat=len(sp)):
        S = sorted([s if c == 0 else allv - s for s, c in zip(sp, choice)], key=len)
        if all(S[i] <= S[i + 1] for i in range(len(S) - 1)):
            return True
    return False

for n in range(5, 9):
    trees = binary_trees(list(range(n)))
    allsys = set()
    for t in trees:
        t = list(t)
        for r in range(len(t) + 1):
            for sub in itertools.combinations(t, r):
                allsys.add(frozenset(sub))
    bad = 0; ncat = 0
    for S in allsys:
        cat = caterpillar(S, n)
        a0 = chains_ok(S, 0, n)
        alla = all(chains_ok(S, a, n) for a in range(n))
        if not (cat == a0 == alla): bad += 1
        ncat += cat
    print(f"k={n}: binary trees {len(trees)}, distinct trees {len(allsys)}, caterpillars {ncat}, mismatches {bad}")
    sys.stdout.flush()
