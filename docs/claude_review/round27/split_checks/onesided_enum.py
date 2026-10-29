# Exact enumeration of one-sided near-tight profiles: choose one side B of each of the 15 cuts with X_a = #{B in F: a in B} equal.
import itertools
n=5
cuts=[]
for r in (1,2):
    for c in itertools.combinations(range(n),r):
        B=frozenset(c); cuts.append((B,frozenset(range(n))-B))
sols=[]
for mask in range(1<<15):
    F=[cuts[i][(mask>>i)&1] for i in range(15)]
    X=[sum(1 for B in F if a in B) for a in range(n)]
    if len(set(X))==1: sols.append((F,X[0]))
print("one-sided near-tight profiles:",len(sols)," X_c values:",sorted(set(x for _,x in sols)))
# code vertices
codes=[]
for r in range(1,6):
    for R in itertools.combinations(range(n),r):
        if r%2==1:
            F=frozenset(B for B,Bc in cuts for B in (B,Bc) if len(B&set(R))%2==0)
            codes.append((R,F))
S=set(frozenset(F) for F,_ in sols)
print("code families F_R (|R| odd) that are near-tight:",[len(R) for R,F in codes if F in S])
