# general p = 3 mod 4, q = p (r=1), k=2, constant model, explicit expansion over F_{p^2}=Z[i]/(p)
import sys, itertools, numpy as np
import gpolyp as G
def run(p,forms,signs,L):
    G.P=p; q=p; k=2; n=q*q; h=q
    chi=lambda x: 0 if x%p==0 else (1 if pow(x,(p-1)//2,p)==1 else -1)
    C=G.sub(G.tpow(q),G.tpow(1)); degC=q; Cp=G.powr(C,p); b=p*(degC+k-1)
    s=L-p*degC; assert s>=1
    T=G.add(np.ones(1,complex), G.mul(Cp, G.add(G.tpow(s),G.tpow(1)) if s%p==0 else G.tpow(s)))
    labs=list(itertools.product(range(p),repeat=2))
    E=[[1+chi(sg*(f[0]*x[0]+f[1]*x[1])) for x in labs] for f,sg in zip(forms,signs)]
    gv=[min(c) for c in zip(*E)]; g=sum(gv)
    ip=(1j)**(p%4)  # i^p
    fac=[]
    for x in labs:
        w=np.zeros(p+1,complex); w[0]=x[0]; w[p]=x[1]
        v=G.mul(ip*Cp,G.norm(w)) if any(x) else np.zeros(1,complex)
        fac.append(G.sub(T,v))
    Z=[G.prod([G.powr(fac[i],E[j][i]-gv[i]) for i in range(n) if E[j][i]>gv[i]]) for j in range(len(forms))]
    D=G.deg(Z[0]); assert all(G.deg(z)==D for z in Z)
    Ns=[G.mul(z,G.conj(z)) for z in Z]
    R=G.prod(fac)
    neg=lambda i: labs.index(tuple((-c)%p for c in labs[i]))
    Q=G.prod([G.powr(fac[i],gv[i]+gv[neg(i)]) for i in range(n)])
    Rp=G.deriv(R)
    # R mod C by folding t^q = t
    Rr=np.rint(R.real).astype(np.int64); e=np.arange(len(Rr)); ee=np.where(e>=q,(e-1)%(q-1)+1,e)
    fold=np.zeros(q,np.int64); np.add.at(fold,ee,Rr)
    a0=G.prod([G.sub(np.zeros(1,complex),G.sub(f,T)) for i,f in enumerate(fac) if any(labs[i])])
    Tp=G.deriv(T); a0Tp=G.mul(a0,Tp)
    pd=[G.deg(G.sub(Z[i],Z[j])) for i in range(len(Z)) for j in range(i+1,len(Z))]
    print(dict(p=p,m=len(forms),g=g,h=h,b=b,L=L,D=D,Dpred=(n-g)*L,norms_equal=all(np.array_equal(N,Ns[0]) for N in Ns),
      N_real=bool(np.all(Ns[0].imag==0)),distinct=len({tuple(z) for z in Z})==len(Z),pair=pd,half=D/2,crit=2*max(pd)<=D,
      bound=((n+h)//2-g)*L+b*(n-h)//2, NQR2=bool(np.array_equal(G.mul(Ns[0],Q),G.mul(R,R))),
      Rmod=np.mod(fold,p).tolist()[:3], Rp_eq=bool(np.array_equal(Rp,a0Tp) or np.array_equal(Rp,G.norm(-a0Tp))),
      Lmin=-(-b*(n-h)//(g-h)) if g>h else None, maxerr=G.MAXERR[0]))
p=int(sys.argv[1])
for fs,sg in [([(1,0),(0,1)],[1,1]),([(1,0),(0,1),(1,1)],[1,1,1]),([(1,0),(1,1),(2,1)],[1,1,-1])]:
    import itertools as it
    # choose L = Lmin computed inside: first pass to get g
    chi=lambda x: 0 if x%p==0 else (1 if pow(x,(p-1)//2,p)==1 else -1)
    labs=list(it.product(range(p),repeat=2))
    g=sum(min(1+chi(s*(f[0]*x[0]+f[1]*x[1])) for f,s in zip(fs,sg)) for x in labs)
    b=p*(p+1); n=p*p; h=p
    if g<=h: print('g<=h',fs,sg,g); continue
    Lmin=max(-(-b*(n-h)//(g-h)), b+1, p*p+2)
    for L in (Lmin-1,Lmin):
        try: run(p,fs,sg,L)
        except AssertionError as e: print('assert at L',L)
