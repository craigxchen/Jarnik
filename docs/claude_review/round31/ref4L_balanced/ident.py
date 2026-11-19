from sympy import *
from itertools import combinations
P,q=symbols('P q')
edges=list(combinations(range(1,6),2))
b=dict(zip(edges,[145,79,-31,109,85,-47,113,19,-41,25]))
print(b)
cyc={(1,2),(2,3),(3,4),(4,5),(1,5)}
S={i:set(e for e in edges if (i in e)!=(e in cyc)) for i in range(1,6)}
for i in S: print(i,sorted(S[i]),len(S[i]))
J={e:P+q*b[e] for e in edges}
n={i:expand(prod([J[e] for e in S[i]])) for i in S}
g={(i,j):prod([J[e] for e in S[i]&S[j]]) for (i,j) in edges}
print({k:len(S[k[0]]&S[k[1]]) for k in edges})
print({k:len(S[k[0]]^S[k[1]]) for k in edges})
R={1:4680,2:3861,3:2860,4:3168,5:2808}
d={1:-76,2:-55,3:-6,4:-34,5:-4}
f={1:26095,2:Rational(43223,2),3:16715,4:16813,5:16303}
T={i:P**2+d[i]*q*P+f[i]*q**2 for i in S}
ps=P+329*q
for (i,j) in edges:
    print((i,j), expand(T[i]*R[j]-T[j]*R[i]-(R[j]-R[i])*g[(i,j)]))
for i in S:
    e=Poly(expand(ps*T[i]**2-n[i]),P,q)
    print(i,{m:c for m,c in zip(e.monoms(),e.coeffs())})
# shift set disjointness
print(sorted(b.values()), 329 in b.values(), 0 in b.values())
