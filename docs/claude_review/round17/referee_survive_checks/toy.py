# Independent toy test of Thm 3.1 / (3.1): exact integral c in V(A,m), alpha_i=<c,M_i> exact in Z[x]/(x^q-1),
# conjugate moduli (mpmath 40 digits), LP Chebyshev (HiGHS, evidence) vs bound; direct check of item 2 on a theta grid.
import sys, itertools, math, cmath
import numpy as np, flint
from mpmath import mp, mpf, mpc, exp as mexp, pi as mpi, binomial
from scipy.optimize import linprog
mp.dps=40
q=int(sys.argv[1]); L=int(sys.argv[2]); ms=[int(x) for x in sys.argv[3].split(',')]
def leg(a):
    a%=q; return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
pats=[]
for a in range(q):
    pats.append([a]+[x for x in range(q) if leg(a-x)==1]); pats.append([a]+[x for x in range(q) if leg(a-x)==-1])
inc=[[i for i,P in enumerate(pats) if x in P] for x in range(q)]
rem=[[sum(1 for y in P if y>=x) for P in pats] for x in range(q+1)]
A=[]; k=[0]*q; sums=[0]*len(pats)
def rec(x):
    if x==q: A.append(tuple(k)); return
    for v in (-1,0,1):
        for i in inc[x]: sums[i]+=v
        if all(abs(sums[i])-rem[x+1][i]<=L for i in range(len(pats))):
            k[x]=v; rec(x+1)
        for i in inc[x]: sums[i]-=v
    k[x]=0
rec(0)
U=len(A); print('q',q,'L',L,'|A|',U,flush=True)
Aidx={a:i for i,a in enumerate(A)}
def monos(deg):  # exponent vectors in {0,1,2}^q of total degree < deg
    out=[]
    def r(x,left,cur):
        if x==q: out.append(tuple(cur)); return
        for e in (0,1,2):
            if e<=left: cur.append(e); r(x+1,left-e,cur); cur.pop()
    r(0,deg-1,[]); return out
def mval(k,al):
    v=1
    for kx,e in zip(k,al):
        if e: v*= kx**e
    return v
# M_i(k) as element of Z[x]/(x^q-1): l_1 = sum k_x x^x, l_{-1} = sum k_x x^{-x}
def polmul(a,b):
    c=[0]*q
    for i,ai in enumerate(a):
        if ai:
            for j,bj in enumerate(b):
                if bj: c[(i+j)%q]+=ai*bj
    return c
def polpow(a,e):
    r=[1]+[0]*(q-1)
    for _ in range(e): r=polmul(r,a)
    return r
zeta=[mexp(2j*mpi*t/q) for t in range(q)]
def ev(poly,s):  # sigma_s applied: x -> zeta^s
    return sum(poly[t]*zeta[(s*t)%q] for t in range(q) if poly[t])
import random
random.seed(1)
for m in ms:
    mon=monos(m)
    E=flint.fmpz_mat([[mval(a,al) for al in mon] for a in A])  # U x N
    # left kernel: c^T E = 0  <=> E^T c = 0
    ET=E.transpose()
    X,nul=ET.nullspace()
    if nul==0: print('m',m,'dim V = 0'); continue
    B=flint.fmpz_mat([[X[i,j] for j in range(nul)] for i in range(U)]).transpose()  # nul x U, rows = kernel vectors
    # make primitive & reduce
    B=B.lll()
    rows=[[int(B[i,j]) for j in range(U)] for i in range(nul)]
    rows.sort(key=lambda r: sum(abs(x) for x in r))
    print('m',m,'dim V',nul,flush=True)
    # precompute l1 poly per k
    l1=[[0]*q for _ in range(U)]
    for ii,a in enumerate(A):
        for x,kx in enumerate(a): l1[ii][x]+=kx
    lm1=[[l[(-t)%q] for t in range(q)] for l in l1]
    # LP Cheb for theta list (evidence)
    mon_arr=np.array([[mval(a,al) for al in mon] for a in A],dtype=float)
    thetas=[0.0,0.3,math.pi/q,1.1]
    cheb={}
    for th in thetas:
        u=np.array([math.cos(th+2*math.pi*x/q) for x in range(q)])
        f=np.array(A,dtype=float)@u; f=f**m
        N=len(mon)
        # min t: -t <= f - M p <= t
        c_obj=np.zeros(N+1); c_obj[-1]=1
        Aub=np.vstack([np.hstack([-mon_arr,-np.ones((U,1))]), np.hstack([mon_arr,-np.ones((U,1))])])
        bub=np.concatenate([-f,f])
        res=linprog(c_obj,A_ub=Aub,b_ub=bub,bounds=[(None,None)]*N+[(0,None)],method='highs')
        cheb[th]=res.fun
    found=False
    for c in rows[:40]:
        cn=sum(abs(x) for x in c)
        # alpha_i exact
        best=None
        Mi_cache={}
        for i in range(m+1):
            al=[0]*q
            for ii,ck in enumerate(c):
                if ck:
                    p=polmul(polpow(l1[ii],i),polpow(lm1[ii],m-i))
                    for t in range(q): al[t]+=ck*p[t]
            # reduce mod Phi_q: alpha=0 iff all coefficients equal
            if len(set(al))==1: continue
            mods=[abs(ev(al,s)) for s in range(1,q)]
            normabs=1
            for z in mods: normabs*=z
            b=mpf(2)**(-m)*binomial(m,i)*max(mods)/cn
            if best is None or b>best[0]: best=(b,i,max(mods),normabs,al)
        if best is None: continue
        found=True
        b,i,mx,nrm,al=best
        print(f'  witness |c|_1={cn} |c|_inf={max(abs(x) for x in c)} i={i} max_s|sigma_s alpha|={float(mx):.6g} |N(alpha)|={float(nrm):.6g}')
        print(f'  bound(before max>=1)={float(b):.6g}  bound(3.1, |c|_1^-1 2^-m binom)={float(mpf(2)**(-m)*binomial(m,i)/cn):.6g}')
        for th in thetas: print(f'    theta={th:.4f} LP Cheb={cheb[th]:.6g}  ratio LP/bound={cheb[th]/float(b):.4g}  OK={cheb[th]>=float(b)*(1-1e-9)}')
        # direct check of item 2 on theta grid: max_a |Phi_{tau^a mu_s c}(u_theta)| >= 2^-m binom |sigma_s alpha|
        s_best=max(range(1,q),key=lambda s: abs(ev(al,s)))
        Aarr=np.array(A,dtype=float); carr=np.array(c,dtype=float)
        worst=1e300
        for th in np.linspace(0,2*np.pi,61):
            vals=[]
            for a_ in range(q):
                # Phi_{tau^a mu_s c}(u) = sum_k c_k <g k,u>^m with g = tau^a mu_s : (gk)_x = k_{s^{-1}(x-a)}
                # <gk,u> = sum_y k_y u_{s y + a}
                u=np.array([math.cos(th+2*math.pi*((s_best*y+a_)%q)/q) for y in range(q)])
                vals.append(abs(carr@((Aarr@u)**m)))
            worst=min(worst,max(vals)/float(mpf(2)**(-m)*binomial(m,i)*abs(ev(al,s_best))))
        print(f'    item-2 check: min over 61 thetas of max_a|Phi|/(2^-m binom |sigma_s alpha|) = {worst:.6f} (must be >=1)',flush=True)
        break
    if not found: print('  no witness among 40 shortest LLL vectors (circle may be in TC_V or witness not found)')
