import sys
sys.argv=[sys.argv[0],'2','8']
exec(open('m8.py').read())
# width ideal: all primitive dirs in {-2..2}^7
V=np.array(list(itertools.product(range(-2,3),repeat=7)),dtype=np.int64)
V=V[[ (np.gcd.reduce(np.abs(v))==1) and next(x for x in v if x)>0 for v in V if v.any()] if False else None] if False else V
W=[]
import math
for v in V:
    if not v.any(): continue
    if math.gcd(*[int(abs(x)) for x in v])!=1: continue
    if next(x for x in v if x)<0: continue
    s=len(np.unique(A@v))
    if s<=DT: W.append((E.pullback([1],[(s,)],[v]),s))
print("width gens",len(W))
for Dt in [5,6,7,8]:
    acc=Acc(Dt); acc.add(rows_deg(W,Dt)); base=acc.r
    print("D",Dt,"J_w",base,"J(A)",JA[Dt])
lines=[]
for i,j in itertools.combinations(range(7),2):
    a=np.array(pts[i]); b=np.array(pts[j]); c=tuple((a+b)%2); k=pts.index(c)
    l=tuple(sorted((i,j,k)))
    if l not in lines: lines.append(l)
print("lines",lines)
# pairs within line spans: kernel directions n in Z^3 with |n_i|<=3
def prim(v): return math.gcd(*[abs(x) for x in v])==1
pairgens=[]; linegens=[]
for l in lines:
    Y=[COLS[i] for i in l]
    g,sz,grid=fam_gens(Y); linegens+=g
    for n in itertools.product(range(-3,4),repeat=3):
        if not any(n) or not prim(n) or next(x for x in n if x)<0: continue
        import sympy as sp
        ns=sp.Matrix([list(n)]).nullspace(); B=[]
        for v in ns:
            den=1
            for x in v: den=sp.ilcm(den,sp.fraction(x)[1])
            B.append([int(x*den) for x in v])
        U=[sum(int(c)*y for c,y in zip(bb,Y)) for bb in B]
        g,sz,grid=fam_gens(U)
        if not grid: pairgens+=g
for Dt in [6,7,8]:
    acc=Acc(Dt); acc.add(rows_deg(W,Dt)); a0=acc.r
    acc.add(rows_deg(pairgens,Dt)); a1=acc.r
    acc2=Acc(Dt); acc2.add(rows_deg(W,Dt)); acc2.add(rows_deg(linegens,Dt))
    print("D",Dt,"J_w",a0,"J_w+pairs in line spans",a1,"J_w+lines",acc2.r,"J(A)",JA[Dt])
