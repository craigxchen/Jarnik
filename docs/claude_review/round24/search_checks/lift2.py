# lift doubled-k=4 V-null curves to k=5: choose one place per pair (B -> B u {4}), solve for rho_4, then polish full system
import itertools, numpy as np, sys
from scipy.optimize import least_squares
import dbl4
from fastK import Prob, diag, isgood, SUBS, IDX
K=[1,2,3,5,6,7,8]
def defect(S):
    # S: 7 complex points; find rho4 minimising || Im e_k(S u {rho4}) ||, k in K. linear in rho4=(a,b)
    P=np.poly(S)  # length 8, P[k] = (-1)^k e_k(S)
    # coefficient of s^{8-k} in (s - rho4) Pi_S = P[k] - rho4*P[k-1]
    A=[];y=[]
    for k in K:
        pk=P[k] if k<8 else 0; pk1=P[k-1]
        # Im(pk - (a+ib) pk1) = Im pk - a Im pk1 - b Re pk1 = 0
        A.append([pk1.imag, pk1.real]); y.append(pk.imag)
    A=np.array(A); y=np.array(y)
    sol,*_=np.linalg.lstsq(A,y,rcond=None)
    res=np.linalg.norm(A@sol-y)/(1+np.linalg.norm(y))
    return res, sol[0]+1j*sol[1]
if __name__=='__main__':
    N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=open(sys.argv[3],'w')
    pats=[e for e in itertools.product([1,-1],repeat=14) if all(sum(e[m] for m in dbl4.MEM[j])==0 for j in range(3))]
    found=0; tried=0
    for it in range(N):
        eps=pats[rng.integers(len(pats))]; pr=dbl4.P(eps,rng,6)
        x0=np.concatenate([rng.normal(size=14)*rng.choice([0.5,1,2]),rng.normal(size=14)*0.7])
        try: sol=least_squares(pr.F,x0,jac=pr.J,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=1500)
        except Exception: continue
        r=pr.rho(sol.x)
        if np.linalg.norm(pr.F(sol.x)[:21])>1e-10 or np.min(np.abs(r.imag))<1e-3: continue
        # candidate choices
        cands=[]
        for ch in itertools.product([0,1],repeat=7):
            S=np.array([r[2*i+ch[i]] for i in range(7)])
            d,r4=defect(S); cands.append((d,ch,r4))
        cands.sort(key=lambda t:t[0])
        for d,ch,r4 in cands[:6]:
            tried+=1
            # build full k=5 rho_T: T subset [4]
            rho=np.zeros(15,complex)
            for i,b in enumerate(dbl4.B):
                rho[IDX[frozenset(b)]]=r[2*i+1-ch[i]]; rho[IDX[frozenset(b|{4})]]=r[2*i+ch[i]]
            rho[IDX[frozenset({4})]]=r4 if abs(r4.imag)>1e-6 else r4+1e-3j
            eps5=np.sign(rho.imag); eps5[eps5==0]=1
            P5=Prob(K,eps5,1,rng)
            x=np.concatenate([rho.real,np.log(np.abs(rho.imag))])
            def F(y): return P5.F(y)[:29]
            def J(y): return P5.J(y)[:29]
            try: s2=least_squares(F,x,jac=J,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
            except Exception: continue
            rr=P5.rho(s2.x); res=np.linalg.norm(P5.F(s2.x)[:28]); dg=diag(rr)
            ok=isgood(rr,res)
            FULL=frozenset(range(1,6)); PL=list(SUBS)+[FULL-T for T in SUBS]; allv=np.concatenate([rr,np.conj(rr)])
            lab=lambda A:''.join(str(a) for a in sorted(A))
            coinc=[(lab(PL[a]),lab(PL[b])) for a in range(30) for b in range(a+1,30) if abs(allv[a]-allv[b])<1e-6]
            print(it,'coinc',coinc,file=out)
            print(it,f"defect0={d:.2e} final_res={res:.1e} imin={dg['imin']:.1e} sep={dg['sep']:.1e} rmax={dg['rmax']:.1e} rmin={dg['rmin']:.1e}",'FOUND' if ok else '',list(rr) if ok else '',file=out,flush=True)
            found+=ok
    print('tried',tried,'found',found,file=out)
