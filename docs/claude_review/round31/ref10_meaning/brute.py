# Exact-ish brute force: all lattice points of norm q_P = prod of split primes in [P,2P].
# Angles phi(eps)=sum eps_j theta_j (mod pi/2 gives all points up to units).
import sys, math, numpy as np
from sympy import primerange
from sympy.ntheory import factorint
import mpmath
mpmath.mp.dps=60
def gauss(p):
    # p = a^2+b^2, a>b>0
    for b in range(1,int(math.isqrt(p))+1):
        a2=p-b*b; a=math.isqrt(a2)
        if a*a==a2 and a>=b: return (a,b)
def run(P, Cs=(1,4,8), mod=math.pi/2, tag='mod pi/2'):
    S=[p for p in primerange(P,2*P+1) if p%4==1]
    s=len(S); pis=[gauss(p) for p in S]
    th=np.array([math.atan2(b,a) for a,b in pis])
    q=1
    for p in S: q*=p
    arr=np.zeros(1)
    for t in th:
        arr=np.concatenate((arr+t,arr-t))
    arr=np.mod(arr,mod); idx=np.argsort(arr,kind='stable'); a=arr[idx]
    ext=np.concatenate((a,a[:1]+mod))
    gaps=np.diff(ext); i=int(np.argmin(gaps)); g=gaps[i]
    q4=mpmath.mpf(q)**0.25
    # exact check of the closest pair
    e1=idx[i]; e2=idx[(i+1)%len(a)]
    def z(e):
        w=mpmath.mpc(1)
        for j in range(s):
            a_,b_=pis[s-1-j]  # bit j of index corresponds to theta index s-1-j? handle below
        return None
    # recompute exactly the angle via mpmath from index bits
    def ang(e):
        tot=mpmath.mpf(0)
        # index construction: after processing theta_0..theta_{s-1}, index bit (s-1-k) ... derive by simulation
        for k in range(s):
            bit=(e>>k)&1
            sgn=1 if bit==0 else -1
            tot+=sgn*mpmath.atan2(pis[k][1],pis[k][0])
        return tot
    d=ang(e2)-ang(e1)
    d=d-mpmath.floor(d/mod)*mod
    if d>mod/2: d=mod-d
    res={'P':P,'s':s,'log10q':float(mpmath.log10(q)),'mingap_float':float(g),'pair':(int(e1),int(e2)),'mingap_exact':float(d),
         'mingap_arc_over_q14':float(d*mpmath.sqrt(q)/q4)}
    for C in Cs:
        w=float(C/q4)
        e2a=np.concatenate((a,a+mod))
        hi=np.searchsorted(e2a,a+w,side='right')
        cnt=hi-np.arange(len(a))
        res['M_max_C%d'%C]=int(cnt.max())
        res['npairs_C%d'%C]=int((cnt-1).sum())
    res['tag']=tag
    return res
if __name__=='__main__':
    for P in map(int,sys.argv[1:]):
        for mod,tag in ((math.pi/2,'all units (mod pi/2)'),(math.pi,'k even (mod pi)')):
            print(run(P,mod=mod,tag=tag),flush=True)
