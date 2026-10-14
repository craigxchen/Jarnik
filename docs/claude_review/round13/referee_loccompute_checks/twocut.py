import numpy as np, itertools, math, sys
from flint import nmod_mat
from reng import REng, PR
a,b,c,L1,L2=map(int,sys.argv[1:6])
d=a+b+c
S1=np.array([1]*a+[1]*b+[0]*c); S2=np.array([0]*a+[1]*b+[1]*c)
K=np.array(list(itertools.product([-1,0,1],repeat=d)),dtype=np.int64)
A=K[(np.abs(K@S1)<=L1)&(np.abs(K@S2)<=L2)]; print("|A|",len(A))
Dm=int(sys.argv[6])
E=REng(d,Dm)
Zr=(A==0).astype(np.int64); Nr=(A==-1).astype(np.int64)
blocks=[]
for D in range(Dm+1):
    M=np.array(E.ms[D],dtype=np.int64)
    zero=(Zr@(M>0).T.astype(np.int64))>0; par=(Nr@(M%2==1).T.astype(np.int64))%2
    blocks.append(np.where(zero,0,1-2*par))
JA=[];prev=0;cum=0
for D in range(Dm+1):
    S=np.concatenate(blocks[:D+1],axis=1)
    r=nmod_mat(S.shape[0],S.shape[1],[int(x)%PR for x in S.flatten()],PR).rank()
    JA.append(len(E.ms[D])-(r-prev)); prev=r
print("J(A) ub",JA,"sumH",prev)
dirs=[v for v in itertools.product(range(-2,3),repeat=d) if any(v) and math.gcd(*[abs(x) for x in v])==1 and next(x for x in v if x)>0]
W=[]
for v in dirs:
    s=len(np.unique(A@np.array(v)))
    if s<=Dm: W.append((list(v),s))
out=[]
for D in range(Dm+1):
    rows=[]
    for v,s in W:
        if s<=D:
            g=E.pullback([1],[(s,)],[v])
            if g.any(): rows+=E.ideal_rows(g,s,D)
    r=nmod_mat(len(rows),len(E.ms[D]),[int(x) for row in rows for x in row],PR).rank() if rows else 0
    out.append(r)
print("J_w(R=2)",out,"equal" if out==JA else "DIFF")
