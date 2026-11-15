# Prop 2.5 check: Theorem-B fake (D0=8,c=0) and its double (D=16,c=0) violate the saturated-lattice Lemma 1.2(1).
import itertools, numpy as np, sympy as sp
k=3; U=list(itertools.product([0,1],repeat=k)); D0=8
chi=lambda a,u: (-1)**(sum(x*y for x,y in zip(a,u))%2)
B=[(1,0,0),(0,1,0),(0,0,1)]; A=[a for a in U if a!=(0,0,0) and a not in B]
W=sp.Matrix([[chi(b,u) for u in U] for b in B])
code01=[[0]*8,[1]*8]+[[ (1+s*chi(a,u))//2 for u in U] for a in A for s in (1,-1)]
C=[[1-2*c for c in w] for w in code01]
assert all(min(sum(x!=y for x,y in zip(p,q)) for q in C if q!=p)>=4 for p in C)
assert len({tuple(W*sp.Matrix(s)) for s in C})==1, 'not one fibre'
print('Theorem-B fake: points',len(C),'one fibre of W: yes; L=',len(B))
h=sp.Matrix([1 if u in [(0,0,0),(1,1,1)] else 0 for u in U])
print('1_H in ker W:',list(W*h),' weight',sum(h),' vs L+1=',len(B)+1)
diffs=sp.Matrix([[ (p[j]-C[0][j])//2 for j in range(8)] for p in C[1:]])
print('1_H in span of differences:',diffs.rank()==diffs.col_join(h.T).rank())
# doubled fake, rows (x, x+e), e perp span(C-C), L=D0-1=7 (c=0 at D=16)
rng=np.random.default_rng(3)
perp=diffs.nullspace(); E=np.array([[float(x) for x in v] for v in perp])  # basis of annihilator
L=7; rows=[]
for i in range(L):
    x=rng.uniform(1,2,8); e=0.1*rng.normal(size=len(E))@E
    rows.append(np.concatenate([x,x+e]))
Wd=np.array(rows)
pts=[np.concatenate([np.array(s),-np.array(s)]) for s in C]
vals=[Wd@p for p in pts]; print('doubled: spread of moment vectors over the 10 points',max(np.abs(v-vals[0]).max() for v in vals))
cols=Wd.T; print('min column separation',min(np.linalg.norm(a-b) for a,b in itertools.combinations(cols,2)),
     ' min |col+col\'|',min(np.linalg.norm(a+b) for a,b in itertools.combinations(cols,2)), ' first row min',Wd[0].min())
a=np.concatenate([np.array(h,float).ravel(),-np.array(h,float).ravel()])
print('doubled lattice vector (1_H,-1_H): |Wd a|=',np.abs(Wd@a).max(),' l1=',np.abs(a).sum(),' vs L+1=',L+1,
      ' a+=',a[a>0].sum(),' vs L/2=',L/2)
# Prop 2.5 'any code' failure: code with unit vector in span(C-C)
C2=[(1,1,1,1),(-1,-1,1,1),(-1,1,-1,1),(1,-1,-1,1)]
d2=sp.Matrix([[ (p[j]-C2[0][j])//2 for j in range(4)] for p in C2[1:]])
print('C2 rank',d2.rank(),' annihilator',[list(v) for v in d2.nullspace()],' -> columns k, D0+k coincide for k=0,1,2')
