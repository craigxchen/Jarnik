import itertools, flint
from collections import Counter
t=4
V=[a for a in itertools.product((0,1),repeat=t) if any(a)]
dot=lambda a,z: sum(x*y for x,y in zip(a,z))%2
M={a:[dot(a,z) for z in V] for a in V}
def rank2(F):
    basis=[]
    for a in F:
        x=int(''.join(map(str,a)),2)
        for bb in basis: x=min(x,x^bb)
        if x: basis.append(x)
    return len(basis)
planes=set()
for B in itertools.combinations(V,3):
    if rank2(B)==3:
        planes.add(frozenset(tuple(sum(ci*b[i] for ci,b in zip(c,B))%2 for i in range(t)) for c in itertools.product((0,1),repeat=3) if any(c)))
sizes=Counter(); ex={}
for r in range(1,16):
    for F in itertools.combinations(V,r):
        if rank2(F)!=t: continue
        A=flint.fmpz_mat([M[a] for a in F]); S=A.snf()
        inv=[abs(int(S[i,i])) for i in range(len(F))]
        idx=1
        for x in inv: idx*=x
        if idx!=2**(len(F)-t) and not any(P<=set(F) for P in planes):
            sizes[r]+=1
            if r not in ex: ex[r]=(F,inv)
print("non-binary block-free families by size:",dict(sizes))
for r,(F,inv) in sorted(ex.items()):
    # circuits of F
    print("size",r,"F=",[''.join(map(str,a)) for a in F],"SNF",inv,"index",eval('*'.join(map(str,inv))),"vs 2^(|F|-t)=",2**(r-t))
