import itertools
from collections import Counter
def best(G, add):
    c=Counter()
    for S in itertools.combinations(G,4):
        s=S[0]
        for x in S[1:]: s=add(s,x)
        c[s]+=1
    return c.most_common(3)
Z8=list(range(8)); print("Z/8", best(Z8, lambda a,b:(a+b)%8))
Z24=[(i,j) for i in range(2) for j in range(4)]; print("Z2xZ4", best(Z24, lambda a,b:((a[0]+b[0])%2,(a[1]+b[1])%4)))
Z222=[(i,j,k) for i in range(2) for j in range(2) for k in range(2)]; print("Z2^3 (char 2 only)", best(Z222, lambda a,b:tuple((x+y)%2 for x,y in zip(a,b))))
# best 8-subsets of a torsion-free group Z^2 inside a box (any finite set embeds in an elliptic curve group)
bestv=0
box=[(i,j) for i in range(4) for j in range(3)]
for T in itertools.combinations(box,8):
    c=Counter()
    for S in itertools.combinations(T,4):
        c[(sum(p[0] for p in S),sum(p[1] for p in S))]+=1
    v=max(c.values())
    if v>bestv: bestv=v; print("Z^2 box", v, T)
# subsets of Z/n x Z/n of size 8, n=4 (in elliptic curve), quick random
import random
for n in (4,6,8):
    G=[(i,j) for i in range(n) for j in range(n)]
    bv=0
    for _ in range(20000):
        T=random.sample(G,8)
        c=Counter()
        for S in itertools.combinations(T,4):
            c[(sum(p[0] for p in S)%n,sum(p[1] for p in S)%n)]+=1
        v=max(c.values()); 
        if v>bv: bv=v
    print("random 8-subsets of (Z/%d)^2:"%n,bv)
