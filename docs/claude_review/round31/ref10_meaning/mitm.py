# Meet in the middle: min over c in {0,+-1}^s \ 0 of dist(sum c_j theta_j, (period) Z),
# period = pi/4 (all units) or pi/2 (k even). Angles as uint64 fixed point: period <-> 2^64 (exact wraparound).
import sys, math, numpy as np, mpmath
from sympy import primerange
mpmath.mp.dps=80
def gauss(p):
    for b in range(1,math.isqrt(p)+1):
        a=math.isqrt(p-b*b)
        if a*a==p-b*b and a>=b: return (a,b)
def sums(th):  # all 3^k signed sums, uint64, and code arrays implicitly by index (base 3: digit 0,1,2 -> 0,+,-)
    arr=np.zeros(1,dtype=np.uint64)
    for t in th:
        t=np.uint64(t)
        arr=np.concatenate((arr,arr+t,arr-t))
    return arr
def decode(idx,k):
    c=[]
    for _ in range(k):
        d=idx%3; idx//=3; c.append(0 if d==0 else (1 if d==1 else -1))
    return c
def run(P,period_name):
    import os; fam=os.environ.get("FAM","split"); S=[p for p in primerange(P,2*P+1) if (p%4==1 if fam=="split" else p%8==5)]; s=len(S); pis=[gauss(p) for p in S]
    per=mpmath.pi/4 if period_name=='pi/4' else mpmath.pi/2
    exact=[mpmath.atan2(b,a) for a,b in pis]
    th=[int(mpmath.nint((x/per)*2**64)) % 2**64 for x in exact]
    s1=min(s//2,15); s2=s-s1
    A=sums(th[:s1]); order=np.argsort(A,kind="stable"); As=A[order]; del A; nA=len(As)
    best=[]  # (dist, ia, ib_global)
    # c2 = 0 case
    nz=As[As!=0]; dd=np.minimum(nz, np.uint64(0)-nz); 
    j=int(np.argmin(dd)); best.append((int(dd[j]), int(order[np.nonzero(As!=0)[0][j]]), 0))
    # B in chunks: fixed prefix over first (s2-c) coords of second half, full sums over the rest
    c=min(s2,15); pre=s2-c
    Bt=sums(th[s1+pre:])     # inner part, 3^c (coordinates s1+pre .. s-1)
    Pre=sums(th[s1:s1+pre])  # 3^pre prefixes
    for ip,pv in enumerate(Pre):
        B=Bt+pv
        t=np.uint64(0)-B
        i=np.searchsorted(As,t)
        i1=i % nA; i0=(i-1) % nA
        d1=As[i1]-t; d0=t-As[i0]; del t
        d=np.minimum(np.minimum(d1,np.uint64(0)-d1),np.minimum(d0,np.uint64(0)-d0))
        if ip==0: d[0]=np.uint64(2**63)  # exclude c2=0
        # keep a few best
        ks=np.argpartition(d,5)[:5]
        for kk in ks:
            dk=int(d[kk])
            ia = int(order[i1[kk]]) if int(np.minimum(d1[kk],np.uint64(0)-d1[kk]))==dk else int(order[i0[kk]])
            ib = ip + (3**pre)*int(kk)  # global index in second half: prefix digits low, inner digits high
            best.append((dk,ia,ib))
        best=sorted(best)[:10]
    # exact verification of the top candidates
    q=1
    for p in S: q*=p
    q4=mpmath.mpf(q)**0.25
    out=[]
    for dk,ia,ib in best:
        c1=decode(ia,s1)
        # second half index: prefix digits for coords s1..s1+pre-1 = ib % 3^pre, inner = ib // 3^pre
        cp=decode(ib % 3**pre, pre); ci=decode(ib//3**pre, s2-pre)
        cc=c1+cp+ci
        tot=sum(ci_*e for ci_,e in zip(cc,exact))
        r=tot-mpmath.nint(tot/per)*per
        out.append((float(abs(r)), float(dk/2**64*per), cc))
    out.sort(key=lambda x:x[0])
    dmin,_,cc=out[0]
    # arc distance between the two points = R * 2*dmin ; normalised by q^{1/4}: 2*dmin*q^{1/4}
    D=[S[j] for j in range(s) if cc[j]!=0]
    # exact Gaussian integer of the minimiser
    zr,zi=1,0
    for j in range(s):
        if cc[j]==0: continue
        a,b=pis[j]; 
        if cc[j]<0: b=-b
        zr,zi=zr*a-zi*b, zr*b+zi*a
    print({'P':P,'s':s,'period':period_name,'min_dist':dmin,'arc_over_q14':float(2*dmin*q4),
           '|D|':len(D),'t_candidates_exact':[round(o[0]/ (o[1] if o[1] else 1),6) for o in out[:3]],
           'gauss':(zr,zi)},flush=True)
if __name__=='__main__':
    per=sys.argv[1]
    for P in map(int,sys.argv[2:]): run(P,per)
