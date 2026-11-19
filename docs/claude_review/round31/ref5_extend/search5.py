# dim-5 V=F_2[t]_{<5}: all transversal products over W = (V minus pair {0,1}) minus at most `maxdrop` further pairs.
import sys
from f2lib import *
V=list(range(32)); P=[v for v in V if v%2==0 and v!=0]  # reps of pairs {v,v+1}, excluding {0,1}
maxdrop=int(sys.argv[1]); best={}
Ts=[1<<L for L in range(5,70)]+[x<<5 for x in range(1,64)]
for T in Ts:
    groups={}
    def dfs(k,mask,drops,prod):
        if k==len(P):
            D=deg4(prod); K=D//2
            key=(mask,prod[0]>>(K+1),prod[1]>>(K+1)); groups[key]=groups.get(key,0)+1; return
        if drops<maxdrop: dfs(k+1,mask,drops+1,prod)
        w=P[k]; m=mask|1<<k
        dfs(k+1,m,drops,m4(prod,(T^w,1))); dfs(k+1,m,drops,m4(prod,(T^w^1,1)))
    dfs(0,0,0,(1,0))
    for key,m in groups.items():
        if m<3: continue
        D=sum(max(deg(T^P[k]),0) for k in range(len(P)) if key[0]>>k&1)
        if m not in best or D<best[m][0]: best[m]=(D,bin(T),bin(key[0]))
for m in sorted(best): print(m,best[m])
