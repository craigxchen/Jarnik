# Search Gaussian-integer root configurations: 8 distinct roots (no conjugate pairs), with >=3 splits into two 4-sets
# whose monic quartics have constant imaginary part (e1,e2,e3 real <=> Im p1,p2,p3 = 0). Gives integer 8-point templates, c=0, D=8.
import itertools, sys
from collections import defaultdict
B=int(sys.argv[1])
pts=[complex(x,y) for x in range(-B,B+1) for y in range(-B,B+1) if y!=0]
def key(rs):
    return tuple(int(round(sum(r**m for r in rs).imag)) for m in (1,2,3))
idx={p:i for i,p in enumerate(pts)}
pairs=defaultdict(list)
for i,j in itertools.combinations(range(len(pts)),2):
    a,b=pts[i],pts[j]
    if a==b.conjugate(): continue
    pairs[key((a,b))].append((i,j))
good=set()
for k,L in pairs.items():
    nk=tuple(-x for x in k)
    if nk not in pairs: continue
    for (i,j) in L:
        for (u,v) in pairs[nk]:
            s=frozenset((i,j,u,v))
            if len(s)<4: continue
            ps=[pts[t] for t in s]
            if any(ps[a]==ps[b].conjugate() for a in range(4) for b in range(4)): continue
            good.add(s)
print("good 4-sets:",len(good),flush=True)
# normalize translation: we will just count 8-sets
splits=defaultdict(set)
goodl=list(good)
byel=defaultdict(list)
for g in goodl:
    for e in g: byel[e].append(g)
for g in goodl:
    # complement candidates: disjoint good sets h with no conj pairs across
    gp=[pts[t] for t in g]
    for h in goodl:
        if h<=g or g&h: continue
        if min(h)<min(g): continue
        hp=[pts[t] for t in h]
        if any(a==b.conjugate() for a in gp for b in hp): continue
        splits[g|h].add(frozenset((g,h)))
best=sorted(splits.items(),key=lambda kv:-len(kv[1]))
from collections import Counter
print("8-sets by #splits:",Counter(len(v) for v in splits.values()))
for S,sp in best[:5]:
    print(len(sp),sorted([pts[t] for t in S],key=lambda z:(z.real,z.imag)))
