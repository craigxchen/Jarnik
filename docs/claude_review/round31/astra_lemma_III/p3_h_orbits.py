"""Enumerate geometric obstructions for all distinct rational h in small fields."""
from itertools import combinations
from p3_h_plane import calc,elim
import sys
p=int(sys.argv[1])
assert p>5
sets={tuple(h) for h in combinations(range(p),5) if sum(h)%p==0}
while sets:
    h=min(sets)
    sets.difference_update(tuple(sorted(x*s%p for x in h)) for s in range(1,p))
    hs=list(h);hs[-1]-=sum(hs)
    out=[]
    try:
        for chart in range(3):
            ctx,x,ds,eqs,D,Delta=calc(p,H=hs,chart=chart)
            out.append((elim(list(ds.values()),'v').total_degree(),elim(eqs,'v').total_degree()))
    except Exception as e:out.append(type(e).__name__)
    print('p',p,'centered orbit',h,'rank-drop and solutions projection degree',out,flush=True)
