# adversarial search: planes P (2 integer directions) whose pullback ideal raises dim(J_w + Fano-pencil planes)_6 above 630
from rlib import *
import pickle, random, sys, time, numpy as np
A,H,dimJ,W=pickle.load(open('m8_a.pkl','rb'))
Aa=np.array(A,dtype=np.int64)
p=P; random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 5)
deg=6; N=math.comb(12,6)
pts=[tuple(random.randrange(p) for _ in range(7)) for _ in range(N+8)]
mv={e:[[evmono(m,pt)%p for pt in pts] for m in monos(7,e)] for e in range(0,7)}
def lin_at(v): return [sum(a*b for a,b in zip(v,pt))%p for pt in pts]
rows=[]
for v,t,h in W:
    if h+1<=deg:
        lv=lin_at(v); fv=[pow(x,h+1,p) for x in lv]
        for m in mv[deg-h-1]: rows.append([x*y%p for x,y in zip(fv,m)])
def plane_rows(V, maxe=6):
    Y=sorted(set(map(tuple,(Aa@np.array(V,dtype=np.int64).T).tolist())))
    g=leading_ideal_gens(Y,2,maxe)
    out=[]
    l1=lin_at(V[0]); l2=lin_at(V[1])
    for e in range(maxe+1):
        for f in g[e]:
            fv=[0]*len(pts)
            for ex,c in f.items():
                for i in range(len(pts)):
                    fv[i]=(fv[i]+c*pow(l1[i],ex[0],p)*pow(l2[i],ex[1],p))%p
            for m in mv[deg-e]: out.append([x*y%p for x,y in zip(fv,m)])
    return out,len(Y)
# baseline: Fano pencil planes span(z_l, c) for the three columns of each line
lines=[(a,b,a^b) for a in range(1,8) for b in range(a+1,8) if (a^b)>b]
for (a,b,c) in lines:
    z=tuple((x+y+w)//2 for x,y,w in zip(col(a),col(b),col(c)))
    for cc in (a,b,c):
        r,_=plane_rows([z,col(cc)]); rows+=r
M=nmod_mat(rows,p); base=M.rank()
print("baseline J_w + Fano pencil planes (mod p, lower bound):",base,flush=True)
# reduce baseline to a basis
R=M.rref()[0]; basis=[[int(R[i,j]) for j in range(R.ncols())] for i in range(base)]
cur=basis; r0=base
t0=time.time(); tested=0; hits=0; seen=set()
pal=[(0,0),(1,0),(0,1),(1,1),(1,-1),(2,1),(1,2),(2,-1),(1,-2),(2,0),(0,2),(2,2),(2,-2),(-1,0),(0,-1)]
while time.time()-t0 < float(sys.argv[2]) if len(sys.argv)>2 else 500:
    k=random.choice([2,3,3,4,4,5,6])
    pal_s=random.sample(pal,k)
    rowsP=[random.choice(pal_s) for _ in range(7)]
    V=[tuple(r[0] for r in rowsP), tuple(r[1] for r in rowsP)]
    if np.linalg.matrix_rank(np.array(V))<2: continue
    key=tuple(np.round(np.linalg.qr(np.array(V,dtype=float).T)[0]@np.linalg.qr(np.array(V,dtype=float).T)[0].T,6).flatten())
    if key in seen: continue
    seen.add(key)
    Y=set(map(tuple,(Aa@np.array(V,dtype=np.int64).T).tolist()))
    if len(Y)>36: continue  # admissible needs |Y|<=36 (rho<8)
    tested+=1
    r,ny=plane_rows(V)
    if not r: continue
    trial=nmod_mat(cur+r,p); rk=trial.rank()
    if rk>r0:
        Yl=sorted(Y); rh=rho(Yl,2)
        hits+=1
        print("HIT plane",V,"|Y|",ny,"rho",rh,"rank",r0,"->",rk,flush=True)
        if rh<8:
            Rr=trial.rref()[0]; cur=[[int(Rr[i,j]) for j in range(Rr.ncols())] for i in range(rk)]; r0=rk
print("tested admissible-size planes",tested,"distinct planes",len(seen),"final rank",r0,"hits",hits)
