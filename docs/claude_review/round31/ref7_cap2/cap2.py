# Independent construction of the VII.1 cap-two cluster (constant model eta=i, p=3), explicit expansion.
import sys, itertools, time, numpy as np
import gpoly as G, f3
from f27 import ADD, MUL, CHI as CHI27, neg as neg27, pw as pw27
def run(q,k,forms,signs,L,Tkind='s',verbose=True,gcdcheck=False):
    t0=time.time()
    p=3
    # field F_q: q=3 -> ints mod 3 ; q=27 -> f27 tables
    if q==3:
        els=[0,1,2]; add=lambda x,y:(x+y)%3; mulf=lambda x,y:(x*y)%3; negf=lambda x:(-x)%3
        chi=lambda x:[0,1,-1][x%3]; digits=lambda x:[x%3]; frob=lambda x:x; ndig=1
    else:
        els=list(range(27)); add=lambda x,y:ADD[x][y]; mulf=lambda x,y:MUL[x][y]; negf=neg27
        chi=lambda x:CHI27[x]; digits=lambda x:[x%3,(x//3)%3,(x//9)%3]; frob=lambda x:pw27(x,9); ndig=3
    n=q**k; h=n//q
    # C = prod_{j=1}^{k-1} (t^{q^j}-t)
    C=np.ones(1,complex)
    for j in range(1,k):
        C=G.mul(C,G.sub(G.tpow(q**j),G.tpow(1)))
    degC=G.deg(C); C3=G.powr(C,3)
    b=p*(degC+k-1)
    # T
    if Tkind=='s':
        s=L-p*degC; assert s>=1
        T=G.add(np.ones(1,complex), G.mul(C3, G.add(G.tpow(s), G.tpow(1)) if s%3==0 else G.tpow(s)))
    elif Tkind=='astra':   # T = 1 + t*C^(3u), L = 1+3u*degC
        u=(L-1)//(3*degC); assert 1+3*u*degC==L
        T=G.add(np.ones(1,complex), G.mul(G.tpow(1), G.powr(C3,u)))
    assert G.deg(T)==L and np.all(T.imag==0)
    labs=list(itertools.product(els,repeat=k))
    def lf(f,x):
        r=0
        for c,xi in zip(f,x):
            for _ in range(c): r=add(r,xi)
        return r
    E=[[1+chi(lf(f,x) if s_==1 else negf(lf(f,x))) for x in labs] for f,s_ in zip(forms,signs)]
    gv=[min(c) for c in zip(*E)]; g=sum(gv)
    idx={x:i for i,x in enumerate(labs)}
    neglab=lambda x:tuple(negf(c) for c in x)
    # orbits under 9-power on coordinates
    seen=set(); orbits=[]
    for x in labs:
        if x in seen: continue
        orb=[x]; y=tuple(frob(c) for c in x)
        while y!=x: orb.append(y); y=tuple(frob(c) for c in y)
        seen.update(orb); orbits.append(orb)
    # descent check: exponents constant on orbits
    for orb in orbits:
        for row in E: assert len({row[idx[y]] for y in orb})==1
    # label polynomial components: v = -i C^3 sum_j x_j t^{3j}; component r over basis a^r
    def labelcomp(x):
        comps=[]
        for r in range(ndig):
            w=np.zeros(3*(k-1)+1,complex)
            for j,xj in enumerate(x): w[3*j]=digits(xj)[r]
            comps.append(G.mul(-1j*C3, G.norm(w)) if np.any(w) else np.zeros(1,complex))
        return comps
    def amul(A,B):  # multiply in F_9[t][a]/(a^3+2a+1)
        if ndig==1: return [G.mul(A[0],B[0])]
        c=[np.zeros(1,complex)]*5
        for i in range(3):
            for j in range(3): c[i+j]=G.add(c[i+j],G.mul(A[i],B[j]))
        for d in (4,3):   # a^d = a^{d-3}(a+2)
            c[d-2]=G.add(c[d-2],c[d]); c[d-3]=G.add(c[d-3],2*c[d]); c[d]=np.zeros(1,complex)
        return c[:3]
    orbpoly=[]
    for orb in orbits:
        P=[np.ones(1,complex)]+[np.zeros(1,complex)]*(ndig-1)
        for x in orb:
            lc=labelcomp(x)
            fac=[G.sub(T,lc[0])]+[G.norm(-lc[r]) for r in range(1,ndig)]
            P=amul(P,fac)
        for r in range(1,ndig): assert G.deg(P[r])==-1, "descent failed"
        orbpoly.append(G.norm(P[0]))
    rep=[idx[o[0]] for o in orbits]
    m=len(forms)
    Z=[]
    for j in range(m):
        Z.append(G.prod([G.powr(orbpoly[oi], E[j][rep[oi]]-gv[rep[oi]]) for oi in range(len(orbits)) if E[j][rep[oi]]-gv[rep[oi]]>0]))
    D=G.deg(Z[0])
    res={'q':q,'k':k,'m':m,'n':n,'h':h,'g':g,'b':b,'L':L,'degC':degC,'D':D,'D_pred':(n-g)*L}
    assert all(G.deg(z)==D for z in Z) and all(z[-1]==1 for z in Z)
    Ns=[G.mul(z,G.conj(z)) for z in Z]
    assert all(np.all(N.imag==0) for N in Ns)
    res['norms_equal']=all(len(N)==len(Ns[0]) and np.all(N==Ns[0]) for N in Ns)
    N=Ns[0]
    res['distinct']=len({tuple(z) for z in Z})==m
    pd=[G.deg(G.sub(Z[i],Z[j2])) for i in range(m) for j2 in range(i+1,m)]
    res['pair_degs']=pd; res['max_pair']=max(pd); res['D/2']=D/2
    res['bound_pair']=((n+h)//2-g)*L+b*(n-h)//2
    res['critical']=2*max(pd)<=D
    # R and Q with N*Q = R^2
    R=G.prod(orbpoly); assert np.all(R.imag==0); res['degR']=G.deg(R)
    Q=G.prod([G.powr(orbpoly[oi], gv[rep[oi]]+gv[idx[neglab(orbits[oi][0])]]) for oi in range(len(orbits))])
    res['NQ=R^2']=bool(np.array_equal(G.mul(N,Q),G.mul(R,R)))
    # R' = a0 * T'  with a0 = prod_{v != 0} v  (up to sign)
    Rp=G.deriv(R); Tp=G.deriv(T)
    a0=np.ones(1,complex)
    for oi,orb in enumerate(orbits):
        if all(c==0 for c in orb[0]): continue
        P=[np.ones(1,complex)]+[np.zeros(1,complex)]*(ndig-1)
        for x in orb: P=amul(P,labelcomp(x))
        a0=G.mul(a0,G.norm(P[0]))
    a0Tp=G.mul(a0,Tp)
    res["R'=+-a0*T'"]=bool(np.array_equal(Rp,a0Tp) or np.array_equal(Rp,G.norm(-a0Tp)))
    # R mod C (C = prod (t^{q^j}-t)) via exact division by f3
    Rr=np.rint(R.real).astype(np.int64); Cr=np.rint(C.real).astype(np.int64)
    res['R mod C']=f3.pmod(Rr,Cr).tolist() if len(Rr)<200000 else 'skipped(see fold)'
    if k==2:  # C=t^q-t: fold exponents
        e=np.arange(len(Rr)); ee=np.where(e>=q,(e-1)%(q-1)+1,e)
        fold=np.zeros(q,np.int64); np.add.at(fold,ee,Rr); res['R mod (t^q-t) by fold']=tr=(np.mod(fold,3)).tolist()
    # exponent census of N: label v contributes primes of T-v with exponent 2-g(v)-g(-v)
    ex=[2-gv[i]-gv[idx[neglab(x)]] for i,x in enumerate(labs)]
    res['labels_exp2']=ex.count(2); res['labels_exp1']=ex.count(1); res['labels_exp0']=ex.count(0)
    res['frac_degN_exp2']=2*ex.count(2)/(2*ex.count(2)+ex.count(1))
    res['frac_degN_exp1']=ex.count(1)/(2*ex.count(2)+ex.count(1))
    res['degN']=G.deg(N); res['degN_pred']=L*(2*ex.count(2)+ex.count(1))
    if gcdcheck:
        Nr=np.rint(N.real).astype(np.int64)
        res['gcd(R,R\')']=f3.gcd(Rr,f3.deriv(Rr)).tolist()
        S=f3.gcd(Nr,Rr); res['deg rad N']=len(S)-1
        q2,r2=f3.pdivmod(Nr,S); res['deg N/rad N (exp-2 part)']=len(f3.tr(q2))-1
        res['gcd(N/radN, R/radN)... exp<=2 check: gcd(N/rad, N/rad^2?)']=None
    res['fft_maxerr']=G.MAXERR[0]; res['secs']=round(time.time()-t0,1)
    if verbose: print(res, flush=True)
    return res,Z,N
if __name__=='__main__':
    pass
