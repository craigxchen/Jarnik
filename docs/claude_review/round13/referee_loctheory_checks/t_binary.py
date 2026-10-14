# exhaustive binarity census: Lambda_F = Lambda_{C(F)} iff index = 2^{|F|-t}; compare with "contains complete 3-block"
import itertools, flint, sys
from collections import Counter
for t in (2,3,4):
    V=[a for a in itertools.product((0,1),repeat=t) if any(a)]
    dot=lambda a,z: sum(x*y for x,y in zip(a,z))%2
    M={a:[dot(a,z) for z in V] for a in V}
    def rank2(F):
        rows=[int(''.join(map(str,a)),2) for a in F]; r=0
        basis=[]
        for x in rows:
            for bb in basis: x=min(x,x^bb)
            if x: basis.append(x); r+=1
        return r
    planes=[]  # complete 3-blocks: 7 nonzero vectors of a 3-dim subspace
    for B in itertools.combinations(V,3):
        if rank2(B)==3:
            S=set()
            for c in itertools.product((0,1),repeat=3):
                if any(c): S.add(tuple(sum(ci*b[i] for ci,b in zip(c,B))%2 for i in range(t)))
            planes.append(frozenset(S))
    planes=set(planes)
    stats=Counter(); bad=[]
    for r in range(1,len(V)+1):
        for F in itertools.combinations(V,r):
            if rank2(F)!=t: continue
            A=flint.fmpz_mat([M[a] for a in F])
            S=A.snf(); inv=[abs(int(S[i,i])) for i in range(len(F))]
            idx=1
            for x in inv: idx*=x
            binary=(idx==2**(len(F)-t))
            hasblock=any(P<=set(F) for P in planes)
            stats[(binary,hasblock)]+=1
            if binary==hasblock: bad.append(F)
    print("t",t,"counts (binary,has complete 3-block):",dict(stats),"counterexamples to 'binary iff no block':",len(bad)); sys.stdout.flush()
