from zonocheck import *
import itertools
regions = [T for T in itertools.product((0,1), repeat=3) if any(T)]   # membership in columns 1..3
seen = set(); cnt = gaps = 0
for dd in range(2, 7):
    for sizes in itertools.product(range(dd+1), repeat=7):
        if sum(sizes) != dd: continue
        # build columns
        coords = []; 
        for T, sz in zip(regions, sizes): coords += [T]*sz
        cols = [tuple(x for x in range(dd) if coords[x][j]) for j in range(3)]
        if any(len(c) < 2 for c in cols) or len(set(cols)) < 3: continue
        for Ls in itertools.product((1, 2), repeat=3):
            if any(L >= len(c) for L, c in zip(Ls, cols)): continue
            # canonical form under permutations of columns
            key = min(tuple(sorted((tuple(sizes[regions.index(tuple(T[pi[j]] for j in range(3)))] for T in regions)), )) + (tuple(Ls[pi[j]] for j in range(3)),)
                      for pi in itertools.permutations(range(3)))
            if key in seen: continue
            seen.add(key)
            cuts = list(zip(cols, Ls))
            A = box(dd, cuts)
            r = rho_p(A, dd); z, edges, k0 = zono(dd, cuts, A)
            cnt += 1
            if z != r:
                gaps += 1; print("GAP", sizes, Ls, r, z, flush=True)
    print(f"d={dd} blocks so far {cnt} gaps {gaps}", flush=True)
