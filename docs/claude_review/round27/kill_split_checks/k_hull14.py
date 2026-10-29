# Affine dimension of the hull of all one-sided near-tight vertices except the x_[5] pair (classes 2^10 4^5 and 1^5 3^10)
import itertools, sympy as sp
full=frozenset(range(5))
Bs=[frozenset(s) for r in range(1,5) for s in itertools.combinations(range(5),r)]
idx={B:i for i,B in enumerate(Bs)}; cuts=[B for B in Bs if 0 in B]
verts=[]; excluded=0
for choice in range(1<<15):
    F=[cuts[i] if (choice>>i)&1 else full-cuts[i] for i in range(15)]
    X=[sum(1 for B in F if a in B) for a in range(5)]
    if len(set(X))!=1: continue
    sig=sorted(len(B) for B in F)
    if sig in ([2]*10+[4]*5,[1]*5+[3]*10): excluded+=1; continue
    v=[0]*30
    for B in F: v[idx[B]]=1
    verts.append(v)
print("vertices used:",len(verts),"excluded:",excluded)
D=sp.Matrix([[v[i]-verts[0][i] for i in range(30)] for v in verts[1:]])
print("affine dim of hull:",D.rank())
avg=[sp.Rational(sum(v[i] for v in verts),len(verts)) for i in range(30)]
print("average is 1/2:",all(a==sp.Rational(1,2) for a in avg))
