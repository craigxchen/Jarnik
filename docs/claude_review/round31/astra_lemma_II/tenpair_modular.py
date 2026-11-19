"""Enumerate good reductions of the degree-one K5 moment system for inert p<=7."""
from collections import defaultdict
from itertools import combinations, product
import sys

def power_moments(p, x, y):
    a,b=1,0
    out=[]
    for k in range(5):
        a,b=(a*x-b*y)%p,(a*y+b*x)%p
        out.append(b)
    return tuple(out)

def count_good(p, allow_disjoint_conjugates=False):
    data=[]
    for x in range(p):
        for y in range(1,(p+1)//2):
            data.extend([power_moments(p,x,y),power_moments(p,x,-y)])
    groups=defaultdict(list)
    for uns in combinations(range(len(data)//2),4):
        for signs in product(range(2),repeat=4):
            roots=[2*u+s for u,s in zip(uns,signs)]
            sig=tuple(sum(v)%p for v in zip(*(data[r] for r in roots)))
            val=sum(1<<r for r in roots)
            flip=sum(1<<(r^1) for r in roots)
            groups[sig].append((val,flip))
    def dfs(group, part, used, twice, forbidden):
        if len(part)==5:
            return 1
        count=0
        for n,(v,f) in enumerate(group):
            if v & twice or (not allow_disjoint_conjugates and v & forbidden):
                continue
            if any((v&w).bit_count()!=1 for w in part):
                continue
            count+=dfs(group[n+1:],part+[v],used|v,twice|v&used,forbidden|f)
        return count
    sizes=[len(g) for g in groups.values()]
    total=sum(dfs(g,[],0,0,0) for g in groups.values())
    print("p",p,"quartets",sum(sizes),"signatures",len(sizes),"largest class",max(sizes,default=0),"allow disjoint conjugate collisions",allow_disjoint_conjugates,"unordered systems",total,flush=True)

if __name__=='__main__':
    count_good(int(sys.argv[1]) if len(sys.argv)>1 else 7, "--allow-disjoint-conjugates" in sys.argv)
