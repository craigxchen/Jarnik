"""Referee: independent recomputation of the contact integers d_C(N) of Phi_N(t0) (orbit Q1+N(Q2-Q1) on E'),
their balance (eta_min of Phi_N itself), and (optionally) window-splittability. Exact integers."""
from fractions import Fraction as Fr
from math import gcd, log, isqrt
import itertools, sys, json
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round18/referee_lift_checks')
from ref_beta7 import setup
ALL8=True
SIGMA=[2,3,5,7,11,13,17,19,23,29,37,41,43,61,67,71,73,97,113,311,521,571,599,1321,2281]
S=setup((2,3,-5,1)); x=S['x']; A=S['A']; lb,lc=S['lb'],S['lc']; kb,s,kc,tp=S['kb'],S['s'],S['kc'],S['tp']
def Delta(lam,T):
    a2,a1,a0=A(lam,T); return a1*a1-4*a2*a0
kap=kb*kc
# E': W^2 = Delta_b(T) Delta_c(T) = kap T (T-s)(T-tp) ; X=kap T, Y=kap W : Y^2 = X(X-kap s)(X-kap tp)
c2=-(kap*s+kap*tp); c1=kap*s*kap*tp
def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    (x1,y1),(x2,y2)=P,Q
    if x1==x2 and y1==-y2: return None
    if P==Q: m=(3*x1*x1+2*c2*x1+c1)/(2*y1)
    else: m=(y2-y1)/(x2-x1)
    x3=m*m-c2-x1-x2; y3=-(y1+m*(x3-x1)); return (x3,y3)
def neg(P): return None if P is None else (P[0],-P[1])
def mul(n,P):
    R=None; Q=P
    if n<0: n=-n; Q=neg(P)
    while n:
        if n&1: R=add(R,Q)
        Q=add(Q,Q); n>>=1
    return R
def Qj(j):
    y=x[j]**2
    a2b,a1b,a0b=A(lb,y); a2c,a1c,a0c=A(lc,y)
    bp=-a1b/a2b-x[j]; cp=-a1c/a2c-x[j]
    W=a2b*(x[j]-bp)*a2c*(x[j]-cp)
    assert W*W==Delta(lb,y)*Delta(lc,y)
    return (kap*y,kap*W)
def ratsqrt(q):
    q=Fr(q)
    if q<0: return None
    n,d=q.numerator,q.denominator; rn,rd=isqrt(n),isqrt(d)
    return Fr(rn,rd) if rn*rn==n and rd*rd==d else None
def lifts(P):
    X,Y=P; T=X/kap; W=Y/kap; a=ratsqrt(T)
    if a is None: return []
    out=[]
    a2b,a1b,a0b=A(lb,T); a2c,a1c,a0c=A(lc,T)
    rb=ratsqrt(Delta(lb,T)); rc=ratsqrt(Delta(lc,T))
    if rb is None or rc is None: return []
    for sa in (a,-a):
        for sb in (rb,-rb):
            for sc in (rc,-rc):
                b=(-a1b+sb)/(2*a2b); bb=(-a1b-sb)/(2*a2b); c=(-a1c+sc)/(2*a2c); cb=(-a1c-sc)/(2*a2c)
                if a2b*(b-bb)*a2c*(c-cb)==W or ALL8: out.append((sa,b,c))
    return out
def prim(q): return (q.numerator,q.denominator)
def sigfree(n):
    n=abs(n); sp=1
    for p in SIGMA:
        while n%p==0: n//=p; sp*=p
    return n,sp
LAB=['A1','A2','A3','a','b','c']
sides=[]
for r in (2,3,4):
    for C in itertools.combinations(range(6),r):
        if sum(1 for k in C if k<3)<=1 and sum(1 for k in range(3) if k not in C)>=2: sides.append(C)
assert len(sides)==25
def profile(pts):
    v=[prim(Fr(p)) for p in pts]
    det={}
    for i,j in itertools.combinations(range(6),2):
        det[(i,j)]=abs(v[i][0]*v[j][1]-v[j][0]*v[i][1])
    g={}
    for C in sides:
        gg=0
        for pr in det:
            if pr[0] in C and pr[1] in C: gg=gcd(gg,det[pr])
        g[C]=sigfree(gg)[0]
    d={}
    for C in sorted(sides,key=len,reverse=True):
        q=g[C]
        for C2 in d:
            if set(C)<set(C2):
                assert q%d[C2]==0; q//=d[C2]
        d[C]=q
    # check product formula: det_ij / prod d_C(same side) is a Sigma-unit
    ok=True
    for (i,j),D in det.items():
        pr=1
        for C,val in d.items():
            same=(i in C and j in C)
            if same: pr*=val
        if D%pr!=0 or sigfree(D//pr)[0]!=1: ok=False
    cop=all(gcd(d[C1],d[C2])==1 for C1,C2 in itertools.combinations(sides,2))
    return d,ok,cop
if __name__=='__main__':
    Q1,Q2=Qj(0),Qj(1); D=add(Q2,neg(Q1))
    res=[]
    for N in range(2,int(sys.argv[1]) if len(sys.argv)>1 else 9):
        R=add(Q1,mul(N,D)); L=lifts(R)
        for k,(a,b,c) in enumerate(L):
            d,ok,cop=profile([x[0],x[1],x[2],a,b,c])
            logs=[log(d[C]) if d[C]>1 else 0.0 for C in sides]
            w6=sum(logs)/25 or 1e-9; mn,mx=min(logs),max(logs)
            etamin=(mx-mn)/(mx+mn) if mx+mn>0 else -1
            print(N,k,'w6=%.3f min/w6=%.3f max/w6=%.3f eta_min(Phi_N)=%.3f prodformula=%s coprime=%s'%(w6,mn/w6,mx/w6,etamin,ok,cop),flush=True)
            res.append(dict(N=N,lift=k,w6=w6,etamin=etamin,d={'|'.join(map(str,C)):str(d[C]) for C in sides}))
    json.dump(res,open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round18/referee_lift_checks/ref_phiN.json','w'))
