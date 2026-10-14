# 3-cut blocks: is every T=1 block with three cuts zonotopal?  rho_p (upper bound on rho_Q) vs verified zonotope.
import itertools, sys, time
from zono import zeta_lower
from tclib import box_points, hilbert, cube_monomials
pats = [p for p in itertools.product((0,1), repeat=3) if any(p)]   # 7 membership patterns
dmax = int(sys.argv[1]); Lmax = int(sys.argv[2])
seen = set(); tot = gap = 0
for d in range(3, dmax+1):
    mons = cube_monomials(d)
    for comp in itertools.product(range(d+1), repeat=7):
        if sum(comp) != d: continue
        # canonical under permutations of the 3 columns
        def key(perm):
            return tuple(comp[pats.index(tuple(p[perm[j]] for j in range(3)))] for p in pats)
        can = min(key(perm) for perm in itertools.permutations(range(3)))
        if can in seen: continue
        seen.add(can)
        coords = []; x = 0
        for p, sz in zip(pats, comp):
            coords += [(p, x + i) for i in range(sz)]; x += sz
        cols = [[xx for p, xx in coords if p[j]] for j in range(3)]
        if any(len(c) < 2 for c in cols): continue
        if len({tuple(c) for c in cols}) < 3: continue
        for Ls in itertools.product(range(1, Lmax+1), repeat=3):
            if any(L >= len(c) for L, c in zip(Ls, cols)): continue
            pts = box_points(d, list(zip(cols, Ls)))
            H = hilbert(pts, mons=mons); rho = max(H)
            z, edges = zeta_lower(pts, cols, list(Ls), d)
            tot += 1
            if z < rho:
                gap += 1
                print('GAP d=%d comp=%s Ls=%s |A|=%d rho_p=%d zeta>=%d' % (d, comp, Ls, len(pts), rho, z), flush=True)
    print('d=%d done, blocks so far %d, gaps %d' % (d, tot, gap), flush=True)
