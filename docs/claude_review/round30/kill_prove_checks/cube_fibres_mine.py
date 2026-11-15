# Independent fibre computation on the cube surface directly from the 8 roots rho in {r_k,-r_k} (normalised Im>0).
import pickle, itertools, mpmath as mp, numpy as np, random
from fractions import Fraction
mp.mp.dps=50
eps=[1,-1,-1,1]
def roots(y):
    R=[]
    for k in range(4):
        r=mp.mpc(eps[k]/y[k],y[k])
        for z in (r,-r):
            R.append(z if z.imag>0 else mp.conj(z))
    return R
def nondeg(R):
    allr=R+[mp.conj(z) for z in R]
    return min(abs(a-b) for a,b in itertools.combinations(allr,2))
def fibres(R):
    mom=[[mp.im(z**m) for m in (1,2,3)] for z in R]
    pts=list(itertools.product([1,-1],repeat=8))
    vecs=[[sum(s[j]*mom[j][m] for j in range(8)) for m in range(3)] for s in pts]
    # cluster by tolerance; rigorous separation check: report min gap between distinct groups vs max spread within
    groups=[]; 
    for i,v in enumerate(vecs):
        for g in groups:
            if max(abs(v[m]-vecs[g[0]][m]) for m in range(3))<mp.mpf(10)**-30: g.append(i); break
        else: groups.append([i])
    within=max([max(abs(vecs[i][m]-vecs[g[0]][m]) for i in g for m in range(3)) for g in groups])
    reps=[vecs[g[0]] for g in groups]
    between=min(max(abs(a[m]-b[m]) for m in range(3)) for a,b in itertools.combinations(reps,2))
    sizes=sorted((len(g) for g in groups),reverse=True)
    big=[g for g in groups if len(g)==sizes[0]]
    ranks=[np.linalg.matrix_rank(np.array([[ (pts[i][j]-pts[g[0]][j])//2 for j in range(8)] for i in g])) for g in big]
    return sizes,within,between,ranks
def onlocus(y0,y1):
    s=-(y0+y1); p=1/(y0*y1); disc=s*s-4*p
    if disc<0: return None
    return [y0,y1,(s+mp.sqrt(disc))/2,(s-mp.sqrt(disc))/2]
# generic members
random.seed(1); out=[]
cnt=0
while cnt<30:
    y=onlocus(mp.mpf(random.uniform(-4,4)),mp.mpf(random.uniform(-4,4)))
    if y is None: continue
    R=roots(y)
    if nondeg(R)<1e-3: continue
    sz,w,b,rk=fibres(R); cnt+=1
    out.append((sz[:4],rk,float(w),float(b)))
prof=set(tuple(o[0]) for o in out); print('generic profiles (top4):',prof,' ranks of max fibres:',set(tuple(o[1]) for o in out),
      ' max within',max(o[2] for o in out),' min between',min(o[3] for o in out))
# special members: orbit under S_4 and y->-y of my elimination solutions
sp_=pickle.load(open('cube_special_mine.pkl','rb'))
orb=set(); pts=[]
for d,yv,nd in sp_:
    if not nd: continue
    yv=[mp.mpf(x) for x in yv]
    for perm in itertools.permutations(range(4)):
        for sg in (1,-1):
            y=[sg*yv[p] for p in perm]
            key=tuple(mp.nstr(t,20) for t in y)
            if key in orb: continue
            orb.add(key); pts.append(y)
print('special orbit size',len(pts))
res=set(); worst=(1e9,0)
for y in pts:
    assert abs(sum(y))<1e-40 and abs(y[0]*y[1]*y[2]*y[3]-1)<1e-40
    R=roots(y); sz,w,b,rk=fibres(R)
    res.add((tuple(sz[:6]),tuple(rk))); worst=(min(worst[0],float(b)),max(worst[1],float(w)))
    assert nondeg(R)>1e-3
print('special profiles:',res,' min between-gap',worst[0],' max within',worst[1])
