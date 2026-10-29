# For one-sided near-tight profiles x_F (0/1 vertices of Pi_nt), volume of Q_F(1) = {m: <m,1_B> >= -1 for B in F}
# Lower bound tau_e(Q) >= (4! vol Q)^{1/4}; sections beat 15/4 at x_F if 4!vol > (15/4)^4 = 197.75.
import itertools, numpy as np
from scipy.spatial import HalfspaceIntersection, ConvexHull
from scipy.optimize import linprog
n=5
Bs=[frozenset(c) for r in range(1,n) for c in itertools.combinations(range(n),r)]
idx={B:i for i,B in enumerate(Bs)}
V=np.load('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks/pi_nt_vertices.npy')
ones=[v for v in V if all(abs(t-round(t))<1e-9 for t in v)]
print("0/1 vertices:", len(ones))
def poly_vol(F, p=None):
    # coordinates m1..m4, m5 = -(m1+..+m4); constraint sum_{a in B} m_a >= -p_B  <=>  -<m,1_B> - p_B <= 0
    hs=[]
    for j,B in enumerate(F):
        coef=np.zeros(4)
        for a in B:
            if a<4: coef[a]+=1
            else: coef-=1
        pb=1.0 if p is None else p[j]
        hs.append(np.concatenate([-coef,[-pb]]))
    hs=np.array(hs)
    # boundedness check: maximize/minimize each coordinate
    for i in range(4):
        for s in (1,-1):
            c=np.zeros(4); c[i]=s
            r=linprog(c,A_ub=hs[:,:4],b_ub=-hs[:,4],bounds=[(None,None)]*4,method='highs')
            if r.status!=0: return None
    hi=HalfspaceIntersection(hs,np.zeros(4))
    return ConvexHull(hi.intersections).volume*24
canon={}
for v in ones:
    F=[B for B in Bs if v[idx[B]]>0.5]
    # canonical form under S_5
    keys=[]
    for perm in itertools.permutations(range(5)):
        keys.append(tuple(sorted(tuple(sorted(perm[a] for a in B)) for B in F)))
    k=min(keys)
    canon.setdefault(k,F)
print("one-sided vertex classes under S_5:", len(canon))
W=[{0,3,4},{0,1,3,4},{3,4},{1,3,4},{0,2,4},{0,1,2,4},{2,4},{1,2,4},{0,2,3},{0,1,2,3},{2,3},{1,2,3},{1},{0},{0,1}]
Wk=min(tuple(sorted(tuple(sorted(perm[a] for a in B)) for B in W)) for perm in itertools.permutations(range(5)))
res=[]
for k,F in canon.items():
    vol=poly_vol(F)
    sizes=sorted(len(B) for B in F)
    res.append((vol,sizes,k==Wk,k))
res.sort(key=lambda r:(r[0] if r[0] is not None else 1e9))
for vol,sizes,isW,k in res:
    print("4!vol=",None if vol is None else round(vol,4), " (vol^(1/4)=", None if vol is None else round(vol**0.25,4),")  sizes",sizes, " W-class" if isW else "")
print("threshold (15/4)^4 =", (15/4)**4)
