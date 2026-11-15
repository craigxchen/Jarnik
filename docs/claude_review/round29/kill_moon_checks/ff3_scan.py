"""Independent F_3(t) circle scan (kill_moon). F_9=F_3(i). Points of squarefree split N: one of lambda, conj(lambda)
for each conjugate pair of monic irreducibles of F_9[t]; ball clusters = equal coefficients in degrees > s,
s=floor(D/2+c). Reports max cluster size per (D,c). Also checks Prop 4.2 and forced collisions on every cluster."""
import itertools, sys
from collections import defaultdict
q=3
E=[(a,b) for a in range(3) for b in range(3)]  # a+bi
idx={e:k for k,e in enumerate(E)}
def m9(x,y): a,b=E[x]; c,d=E[y]; return idx[((a*c-b*d)%3,(a*d+b*c)%3)]
def a9(x,y): a,b=E[x]; c,d=E[y]; return idx[((a+c)%3,(b+d)%3)]
def s9(x,y): a,b=E[x]; c,d=E[y]; return idx[((a-c)%3,(b-d)%3)]
def c9(x): a,b=E[x]; return idx[(a,(-b)%3)]
MUL=[[m9(x,y) for y in range(9)] for x in range(9)]; ADD=[[a9(x,y) for y in range(9)] for x in range(9)]
SUB=[[s9(x,y) for y in range(9)] for x in range(9)]; CJ=[c9(x) for x in range(9)]
ZERO=idx[(0,0)]; ONE=idx[(1,0)]
# polys: tuple of coeffs low->high
def pmul(f,g):
    r=[ZERO]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        if a==ZERO: continue
        for j,b in enumerate(g): r[i+j]=ADD[r[i+j]][MUL[a][b]]
    return tuple(r)
def pconj(f): return tuple(CJ[a] for a in f)
def deg(f):
    for k in range(len(f)-1,-1,-1):
        if f[k]!=ZERO: return k
    return -1
def psub(f,g):
    n=max(len(f),len(g)); f=f+(ZERO,)*(n-len(f)); g=g+(ZERO,)*(n-len(g)); return tuple(SUB[a][b] for a,b in zip(f,g))
def monics(d):
    for c in itertools.product(range(9),repeat=d): yield tuple(c)+(ONE,)
def is_irred(f,lowirr):
    d=deg(f)
    for g in lowirr:
        if 2*deg(g)>d: break
        # trial division
        r=list(f)
        dg=deg(g)
        for k in range(d-dg,-1,-1):
            co=r[k+dg]
            if co!=ZERO:
                for j in range(dg+1): r[k+j]=SUB[r[k+j]][MUL[co][g[j]]]
        if all(x==ZERO for x in r[:dg]): return False
    return True
DMAX=int(sys.argv[1]) if len(sys.argv)>1 else 6
irr=[]; pairs=defaultdict(list)
for d in range(1,DMAX+1):
    if d>3: break
    for f in monics(d):
        if is_irred(f,irr):
            irr.append(f)
            g=pconj(f)
            if g!=f and g>f: pairs[d].append((f,g))
print({d:len(v) for d,v in pairs.items()})
# enumerate sets of pairs with sum of degrees D (multiset-free: distinct pairs)
allp=[(d,p) for d in sorted(pairs) for p in pairs[d]]
best=defaultdict(int); nN=defaultdict(int); viol=0; nclus=0
def rec(start,D,chosen):
    global viol,nclus
    if D>0:
        pass
    for k in range(start,len(allp)):
        d,p=allp[k]
        if D+d>DMAX: continue
        ch=chosen+[(d,p)]; DD=D+d
        evaluate(ch,DD)
        rec(k+1,DD,ch)
def evaluate(ch,D):
    global viol,nclus
    nN[D]+=1
    om=len(ch)
    pts=[]
    for eps in itertools.product((0,1),repeat=om):
        if eps[0]==1: continue  # conjugation symmetry: cluster sizes same for conj set; fix first
        z=(ONE,)
        for (d,(f,g)),e in zip(ch,eps): z=pmul(z,f if e==0 else g)
        pts.append((eps,z))
    # include conj: z and zbar both points; but zbar ball is conj of z ball -> same sizes; need full set though
    full=[]
    for eps in itertools.product((0,1),repeat=om):
        z=(ONE,)
        for (d,(f,g)),e in zip(ch,eps): z=pmul(z,f if e==0 else g)
        full.append((eps,z))
    for c in (-1,0,1):
        s=(D//2) + c if True else None
        s=int((D/2+c)//1)
        groups=defaultdict(list)
        for eps,z in full: groups[z[s+1:]].append((eps,z))
        mx=max(len(v) for v in groups.values())
        if mx>best[(D,c)]: best[(D,c)]=mx
        if c==0:
            for v in groups.values():
                if len(v)<2: continue
                nclus+=1; M=len(v)
                # Prop 4.2 check: sum deg c_xy <= binom(M,2)c + M D/4 - cut;  deg c = deg(z_x-z_y)-deg g
                tot=0
                for (e1,z1),(e2,z2) in itertools.combinations(v,2):
                    dg=sum(d for (d,_),a,b in zip(ch,e1,e2) if a==b)
                    tot+=deg(psub(z1,z2))-dg
                    assert deg(psub(z1,z2))>=dg
                cut=0
                for j,(d,_) in enumerate(ch):
                    S=sum(1 for e,_ in v if e[j]==0); cut+=d*(M-2*S)**2/4
                if tot> 0*(M*(M-1)/2)+M*D/4-cut+1e-9: viol+=1
rec(0,0,[])
print("N counts by D:",dict(nN))
for D in range(1,DMAX+1): print("D=",D," max cluster c=-1,0,1:",[best[(D,c)] for c in (-1,0,1)])
print("clusters checked at c=0:",nclus,"Prop4.2 violations:",viol)
