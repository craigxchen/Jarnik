"""Referee's independent check of lift.md Theorem 4.2(a), by a different route.
V = H^0(O_E(3H)) (x) H^0(O(2)) (36-dim), restricted to forms vanishing on aux_j x P^1.
W = those vanishing at the 25 points (B_A, v_A).  Claim to test: dim W = 2 and W = Pi*Pencil.
Then non-special q: no solution unless Pi(q)=0 or degenerate member. Special q: direct systems.
All exact over Q (Fractions, own Gaussian elimination)."""
from fractions import Fraction as Fr
import sys, json, itertools

def rank(M):
    M=[row[:] for row in M]; r=0; ncol=len(M[0]) if M else 0
    for c in range(ncol):
        piv=None
        for i in range(r,len(M)):
            if M[i][c]!=0: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        pv=M[r][c]
        M[r]=[x/pv for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                f=M[i][c]; M[i]=[x-f*y for x,y in zip(M[i],M[r])]
        r+=1
    return r

def nullspace(M,n):
    M=[row[:] for row in M]; r=0; pivs=[]
    for c in range(n):
        piv=None
        for i in range(r,len(M)):
            if M[i][c]!=0: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; pv=M[r][c]; M[r]=[x/pv for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                f=M[i][c]; M[i]=[x-f*y for x,y in zip(M[i],M[r])]
        pivs.append(c); r+=1
    free=[c for c in range(n) if c not in pivs]; basis=[]
    for fcol in free:
        v=[Fr(0)]*n; v[fcol]=Fr(1)
        for i,pc in enumerate(pivs): v[pc]=-M[i][fcol]
        basis.append(v)
    return basis

def setup(t):
    x=[Fr(z) for z in t]
    e1=sum(x); e2=sum(x[i]*x[j] for i in range(4) for j in range(i+1,4))
    e3=sum(x[i]*x[j]*x[k] for i,j,k in itertools.combinations(range(4),3)); e4=x[0]*x[1]*x[2]*x[3]
    lb=4*e4/(e3**2-4*e2*e4); lc=Fr(-4)/e1**2
    def A(lam,T): # P_lam - T Q_lam = A2 z^2 + A1 z + A0
        return (1+lam*e2+lam*T, -lam*(e3+e1*T), lam*e4-T)
    def Delta_coeffs(lam): # Delta(T)=d2 T^2 + d1 T + d0
        # A1^2 - 4 A2 A0 ; A2 = (1+lam e2) + lam T ; A1 = -lam e3 - lam e1 T ; A0 = lam e4 - T
        p2=(1+lam*e2, lam); p1=(-lam*e3, -lam*e1); p0=(lam*e4, Fr(-1))
        def mul(p,q): return (p[0]*q[0], p[0]*q[1]+p[1]*q[0], p[1]*q[1])
        a=mul(p1,p1); b=mul(p2,p0)
        return tuple(a[i]-4*b[i] for i in range(3))
    db=Delta_coeffs(lb); dc=Delta_coeffs(lc)
    assert db[0]==0 and dc[2]==0
    kb=db[2]; s=-db[1]/kb; kc=dc[1]; tp=-dc[0]/kc
    D=e3**2-e1*e2*e3+e1**2*e4
    return dict(x=x,e=(e1,e2,e3,e4),lb=lb,lc=lc,A=A,kb=kb,s=s,kc=kc,tp=tp,D=D)

def points(S):
    x=S['x']; A=S['A']; lb=S['lb']; lc=S['lc']; kb=S['kb']; s=S['s']; kc=S['kc']; tp=S['tp']
    pts=[]  # (j, Sset, (a,u,w,h))
    for j in range(4):
        y=x[j]**2
        A2b,A1b,A0b=A(lb,y); A2c,A1c,A0c=A(lc,y)
        assert A2b*x[j]**2+A1b*x[j]+A0b==0 and A2c*x[j]**2+A1c*x[j]+A0c==0
        broots=[x[j], -A1b/A2b-x[j]]; croots=[x[j], -A1c/A2c-x[j]]
        for a in (x[j],-x[j]):
            for bb in broots:
                for cc in croots:
                    u=(2*A2b*bb+A1b)/a; w=2*A2c*cc+A1c
                    assert u*u==kb*(a*a-s) and w*w==kc*(a*a-tp)
                    Sset=frozenset(m for m,val in (('a',a),('b',bb),('c',cc)) if val==x[j])
                    pts.append((j,Sset,(a,u,w,Fr(1))))
    return pts

# cubic basis on E (normal monomials)
MONS=[(3,0,0,0),(2,0,0,1),(1,0,0,2),(0,0,0,3),(2,1,0,0),(1,1,0,1),(0,1,0,2),(2,0,1,0),(1,0,1,1),(0,0,1,2),(1,1,1,0),(0,1,1,1)]
def mono(m,p):
    r=Fr(1)
    for e,c in zip(m,p): r*=c**e
    return r
def dmono(m,p,tau):
    tot=Fr(0)
    for k in range(4):
        if m[k]==0: continue
        mm=list(m); mm[k]-=1
        tot+=m[k]*mono(mm,p)*tau[k]
    return tot

def tangent(S,p):
    a,u,w,h=p; kb,s,kc,tp=S['kb'],S['s'],S['kc'],S['tp']
    g1=[-2*kb*a, 2*u, Fr(0), 2*kb*s*h]; g2=[-2*kc*a, Fr(0), 2*w, 2*kc*tp*h]
    ker=nullspace([g1,g2],4); assert len(ker)==2
    for v in ker:
        if rank([list(p),v])==2: return v
    raise Exception

NU=36 # index m*3+r
def row_point(p,v):
    return [mono(MONS[m],p)*v**r for m in range(12) for r in range(3)]
def rows_fiber(p):
    out=[]
    for r0 in range(3):
        out.append([mono(MONS[m],p) if r==r0 else Fr(0) for m in range(12) for r in range(3)])
    return out
def row_dpoint(p,tau,v):
    return [dmono(MONS[m],p,tau)*v**r for m in range(12) for r in range(3)]
def rows_dfiber(p,tau):
    return [[dmono(MONS[m],p,tau) if r==r0 else Fr(0) for m in range(12) for r in range(3)] for r0 in range(3)]

def run(t):
    S=setup(t); x=S['x']; e1,e2,e3,e4=S['e']
    pts=points(S)
    aux=[p for (j,Ss,p) in pts if j<3 and len(Ss)==0]
    contacts=[]
    for (j,Ss,p) in pts:
        if j<3 and len(Ss)>0: contacts.append((('A%d'%(j+1),)+tuple(sorted(Ss)),p,x[j]))
        if j==3 and len(Ss)>=2: contacts.append((tuple(sorted(Ss)),p,x[3]))
    assert len(aux)==3 and len(contacts)==25
    allp=[p for (_,_,p) in pts]; assert len(set(allp))==32
    base=[]
    for p in aux: base+=rows_fiber(p)
    cond=[row_point(p,v) for (_,p,v) in contacts]
    Wb=nullspace(base+cond,NU)
    res={'t':[str(z) for z in t],'dimW':len(Wb),'D':str(S['D'])}
    # plane through aux
    Pi=nullspace([list(p) for p in aux],4); assert len(Pi)==1; Pi=Pi[0]
    def Piv(p): return sum(c*z for c,z in zip(Pi,p))
    on_plane=[name for (name,p,v) in contacts if Piv(p)==0]
    res['contacts_on_Pi']=on_plane
    # Pi*sigma_0 and Pi*sigma_inf as coefficient vectors
    def lin_times(quadpoly):
        # quadpoly: dict over monomials in (a,u,w,h) degree 2 -> dict r->coef
        vec=[Fr(0)]*NU
        lin=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
        for li,c in zip(lin,Pi):
            if c==0: continue
            for qm,rc in quadpoly.items():
                m=tuple(i+j for i,j in zip(li,qm)); k=MONS.index(m)
                for r,cc in rc.items(): vec[k*3+r]+=c*cc
        return vec
    sig0={(0,0,0,2):{2:Fr(1)},(2,0,0,0):{0:Fr(-1)}}
    siginf={(0,0,0,2):{2:e2,1:-e3,0:e4},(2,0,0,0):{1:-e1,2:Fr(1)}}
    v0=lin_times(sig0); vi=lin_times(siginf)
    M=base+cond
    res['Pi_sig0_in_W']=all(sum(a*b for a,b in zip(row,v0))==0 for row in M)
    res['Pi_siginf_in_W']=all(sum(a*b for a,b in zip(row,vi))==0 for row in M)
    res['Pi_pencil_rank']=rank([v0,vi])
    # special q
    spec={}
    for (name,p,vA) in contacts:
        tau=tangent(S,p)
        rows=base+rows_fiber(p)+[row_dpoint(p,tau,vA)]+[row_point(pp,vv) for (nn,pp,vv) in contacts if nn!=name]
        spec[' '.join(name)]=NU-rank(rows)
    for k,p in enumerate(aux):
        tau=tangent(S,p)
        rows=base+rows_dfiber(p,tau)+cond
        spec['aux%d'%(k+1)]=NU-rank(rows)
    res['special_dims']=spec
    return res

if __name__=='__main__':
    T=[(2,3,-5,1),(Fr(1,3),3,2,Fr(5,3)),(-2,Fr(7,5),-3,1),(2,Fr(9,5),Fr(6,5),1),(Fr(5,2),Fr(-9,5),-7,1),
       (Fr(-3,4),Fr(-1,3),Fr(1,2),1),(-3,Fr(1,2),4,1),(Fr(3,5),2,Fr(-3,2),1),(Fr(7,5),Fr(3,2),Fr(7,3),1)]
    if len(sys.argv)>1: T=T[:int(sys.argv[1])]
    out=[]
    for t in T:
        r=run(t); out.append(r)
        nz={k:v for k,v in r['special_dims'].items() if v!=0}
        print(r['t'],'dimW',r['dimW'],'D',r['D'],'onPi',r['contacts_on_Pi'],'Pi*pencil in W',r['Pi_sig0_in_W'],r['Pi_siginf_in_W'],r['Pi_pencil_rank'],'nonzero special',nz,flush=True)
    json.dump(out,open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round18/referee_lift_checks/ref_beta7.json','w'),indent=1)
