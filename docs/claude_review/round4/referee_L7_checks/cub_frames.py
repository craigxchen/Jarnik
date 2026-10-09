import sys, itertools
sys.path.insert(0, '../L7')
import classes7 as C
planes = list(itertools.combinations(range(1, 7), 3))
for choice in itertools.combinations(planes, 13):
    vec = C.kapranov(3, {}, {}, {p: 1 for p in choice})
    if min(vec.values()) >= 0: break
cub = vec
lab = list(range(1, 8))
def key(S):
    S = frozenset(S); return S if 7 not in S else frozenset(lab) - S
def crdeg(a, b, c, d):
    rest = [x for x in lab if x not in (a, b, c, d)]
    tot = 0
    for r in range(len(rest) + 1):
        for U in itertools.combinations(rest, r):
            S = set((a, b)) | set(U)
            if 2 <= len(S) <= 5: tot += cub[key(S)]
    return tot
print('nonzero D-values:')
for S in C.SPL:
    if cub[S]: print('  ', sorted(S), cub[S])
best = None
for anchors in itertools.combinations(lab, 3):
    movers = [x for x in lab if x not in anchors]
    a, b, c = anchors
    degs = [crdeg(a, b, c, m) for m in movers]   # cross-ratio (a,b;c,m)-type degree = deg of x_m in frame a,b,c
    tot = sum(degs)
    if best is None or tot < best[0]: best = (tot, anchors, movers, degs)
    print(anchors, movers, degs)
print('best frame', best)
