# Exact asymptotic template constant: C_inf = max|gamma| / sqrt( prod_l max_{t0 in Q_l} |g(t0)|_l ).
# Local maxima are independent (CRT + archimedean freedom), so the product is attained along a progression of t0.
import sys, itertools, math
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round30/bridge_checks')
from glib import *
import sympy as sp
from fractions import Fraction as Fr

def split_prime(l):
    # Gaussian prime of norm l, l = 1 mod 4
    for a in range(1,int(math.isqrt(l))+1):
        b2=l-a*a; b=math.isqrt(b2)
        if b*b==b2: return (a,b)
def vpi_int(x,pi):
    # x Gaussian integer tuple, pi tuple; valuation
    if x==(0,0): return 10**9
    v=0
    while True:
        q,r=gi_divmod(x,pi)
        if r!=(0,0): return v
        x=q; v+=1
def vpi(z,pi,l):
    # z Gaussian rational G; valuation at pi (pi above l)
    den=1
    for q in (z.a,z.b): den=den*q.denominator//math.gcd(den,q.denominator)
    x=(int(z.a*den),int(z.b*den))
    vd=0; d=den
    while d%l==0: d//=l; vd+=1
    # v_pi(den) : l inert -> vd ; l=2 -> 2vd ; split -> vd
    vden = 2*vd if l==2 else vd
    return vpi_int(x,pi)-vden

def gauss_primes_above(l):
    if l==2: return [((1,1),2)]
    if l%4==3: return [((l,0),l*l)]
    a,b=split_prime(l); return [((a,b),l),((a,-b),l)]

def sqrt_m1(l,K):
    # x with x^2=-1 mod l^K
    x=next(a for a in range(2,l) if (a*a+1)%l==0)
    for _ in range(K.bit_length()+2):
        m=l**K
        x=(x-(x*x+1)*pow(2*x,-1,m))%m
    return x

def fval(t0,pts_roots,l):
    lg=0.0
    for pi,Npi in gauss_primes_above(l):
        mv=min(sum(vpi(t0-r,pi,l) for r in rr) for rr in pts_roots)
        lg+=mv*math.log(Npi)/2
    return lg

def local_logmax(pts_roots, l, J=3, extra=3):
    allr=set()
    for rr in pts_roots: allr.update(rr)
    allr=list(allr)
    vmax=0
    for pi,_ in gauss_primes_above(l):
        for a,b in itertools.combinations(allr,2):
            vmax=max(vmax,vpi(a-b,pi,l))
        for a in allr: vmax=max(vmax,-vpi(a,pi,l))
    k=vmax+extra
    cands=[Fr(0)]
    if l%4==1:
        K=k+6
        for iota in (sqrt_m1(l,K), l**K-sqrt_m1(l,K)):
            for r in allr: cands.append(r.a+r.b*iota)
    elif l%4==3:
        for r in allr: cands.append(r.a)     # nearest-candidate argument: exact for inert l
    else:
        for r in allr: cands+= [r.a, r.a+r.b, r.a-r.b]
        jd=max([0]+[ (q.denominator & -q.denominator).bit_length()-1 for r in allr for q in (r.a,r.b)])
        J2=max(J,jd+1)
        for j in range(0,J2+1):
            mod=2**min(k+j,14)
            for a in range(mod):
                if j>0 and a%2==0: continue
                cands.append(Fr(a,2**j))
    best=None
    for t0 in cands:
        lg=fval(t0,pts_roots,l)
        if best is None or lg>best[0]+1e-12: best=(lg,t0)
    return best

def relevant_primes(pts_roots):
    allr=set()
    for rr in pts_roots: allr.update(rr)
    allr=list(allr); P=set()
    for a in allr:
        for q in (a.a,a.b): P.update(sp.factorint(q.denominator).keys())
    for a,b in itertools.combinations(allr,2):
        n=(a-b).norm(); P.update(sp.factorint(n.numerator).keys()); P.update(sp.factorint(n.denominator).keys())
    return sorted(int(x) for x in P)

def C_inf(pts_roots, D, verbose=False):
    zs=[prod_lin(rr) for rr in pts_roots]
    gm=0
    for a,b in itertools.combinations(zs,2):
        d=psub(a,b)
        if deg(d)>D//2: return None
        if len(d)>D//2: gm=max(gm,float(d[D//2].norm())**0.5)
    tot=0.0; info=[]
    for l in relevant_primes(pts_roots):
        lg,t0=local_logmax(pts_roots,l)
        tot+=lg; info.append((l,round(math.exp(lg),6),t0))
    if verbose: print(" local factors |g|_l max:",info)
    return gm/math.exp(tot/2), gm, info
