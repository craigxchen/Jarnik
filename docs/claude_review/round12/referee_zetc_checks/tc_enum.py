import itertools
def blocks(maxd):
    out = []; seen = set()
    for dd in range(2, maxd+1):
        # columns as subsets of range(dd), each size>=2, distinct, union = range(dd); dedupe under S_dd x S_3
        subsets = [frozenset(c) for r in range(2, dd+1) for c in itertools.combinations(range(dd), r)]
        for trip in itertools.combinations(subsets, 3):
            if frozenset().union(*trip) != frozenset(range(dd)): continue
            for Ls in itertools.product((1, 2), repeat=3):
                if any(L >= len(c) for L, c in zip(Ls, trip)): continue
                # canonical: region-size vector under column permutations (coordinate perms are absorbed)
                best = None
                for pi in itertools.permutations(range(3)):
                    cs = [trip[pi[j]] for j in range(3)]; ls = tuple(Ls[pi[j]] for j in range(3))
                    sig = tuple(sum(1 for x in range(dd) if tuple(int(x in c) for c in cs) == T)
                                for T in itertools.product((0, 1), repeat=3))
                    key = (dd, sig, ls)
                    best = key if best is None or key < best else best
                if best in seen: continue
                seen.add(best); out.append((dd, [tuple(sorted(c)) for c in trip], Ls))
    return out
if __name__ == "__main__":
    b = blocks(6)
    from collections import Counter
    print(len(b), Counter(x[0] for x in b))
