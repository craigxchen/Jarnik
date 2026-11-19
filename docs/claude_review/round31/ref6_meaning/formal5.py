"""Own formal (X-level) analysis of the VI.4 clusters restricted to subtuples.
F_q = F_3[x]/(f), q=3^r; chi via a log table. V = F_q^k (coefficient vectors of v/i).
e_j(v) = 1 + chi(eps_j * lambda_j(v)), lambda_j with F_3 coefficients.
For a subtuple S: g_S = sum_v min_{j in S} e_j(v); per conjugate pair {v,-v} two layers
A1={j: e_j(v)>=1}, A2={j: e_j(v)=2}; nontrivial cuts get weight 1 (X-degree; x L after X=t^L).
Formal pair difference X-degree (n+h)/2 (top term c^{(q-1)/2} A^{(q+1)/2}, generic), so
chord c_ab has X-degree (n+h)/2 - g_ab; reduced exponent delta_S = (g_S-h)/(2(n-g_S)),
alpha_S = 4/(1+2 delta_S) [definition as in ref5_reconcile]."""
import numpy as np, itertools, sys
from collections import Counter

IRR={3:[1,2,0,1],5:[1,2,0,0,0,1]}   # x^3+2x+1, x^5+2x+1 (coeffs low->high, monic)
def field(r):
    f=IRR[r]; q=3**r
    # element = int base 3, digits low->high; multiply by x
    def mulx(d):
        d=[0]+d; top=d.pop()     # degree r term
        return [ (d[i]-top*f[i])%3 for i in range(r)]
    # find primitive element by brute force over candidates; build powers
    def toint(d): return sum(c*3**i for i,c in enumerate(d))
    def todig(a): return [(a//3**i)%3 for i in range(r)]
    def mul(a,b):
        A=todig(a); B=todig(b); res=[0]*r; cur=A[:]
        for i in range(r):
            if B[i]: res=[(res[t]+B[i]*cur[t])%3 for t in range(r)]
            cur=mulx(cur)
        return toint(res)
    for gcand in range(2,q):
        x=1; seen=set(); logt={}
        for e in range(q-1):
            logt[x]=e; x=mul(x,gcand)
        if len(logt)==q-1: break
    chi=np.zeros(q,np.int64)
    for a,e in logt.items(): chi[a]= 1 if e%2==0 else -1
    return q,chi

def run(r,k,forms,eps,subsize):
    q,chi=field(r)
    n=q**k; h=n//q
    # digits of each coordinate: v is k elements of F_q, each r base-3 digits
    idx=np.arange(n)
    coords=[(idx//q**b)%q for b in range(k)]
    digs=[np.stack([(c//3**i)%3 for i in range(r)]) for c in coords]   # k arrays (r,n)
    def lam(coef):
        acc=np.zeros((r,n),np.int64)
        for b,a in enumerate(coef):
            if a: acc=(acc+a*digs[b])%3
        return (acc*(3**np.arange(r))[:,None]).sum(0)
    # negation map index for v -> -v
    negd=[(3-d)%3 for d in digs]
    negidx=np.zeros(n,np.int64)
    for b in range(k):
        cval=(negd[b]*(3**np.arange(r))[:,None]).sum(0)
        negidx+=cval*q**b
    E=[]
    for coef,s in zip(forms,eps):
        val=lam(coef)
        if s==-1: val=lam([(-a)%3 for a in coef])
        E.append(1+chi[val])
    E=np.array(E)
    m=len(forms)
    assert all((E[j]+E[j][negidx]==2).all() for j in range(m))
    gall=E.min(0).sum()
    print(f"q={q} k={k} n={n} h={h} m={m} formal g(all)={gall}  slack (g-h)/2={(gall-h)/2}")
    rep=idx<negidx   # one representative per conjugate pair (v=0 is self-paired)
    out=[]
    for S in itertools.combinations(range(m),subsize):
        ES=E[list(S)]
        gS=ES.min(0).sum()
        gab={(a,b):np.minimum(E[a],E[b]).sum() for a,b in itertools.combinations(S,2)}
        cuts=Counter()
        R=ES[:,rep]
        for lev in (1,2):
            A=(R>=lev)
            key=(A*(1<<np.arange(subsize))[:,None]).sum(0)
            full=(1<<subsize)-1
            for kk,c in zip(*np.unique(key,return_counts=True)):
                if kk in (0,full): continue
                kk=min(kk,full^kk)
                cuts[kk]+=c
        ncut=2**(subsize-1)-1
        w=[cuts.get(c,0) for c in range(1,1<<subsize) if c==min(c,((1<<subsize)-1)^c)]
        delta=(gS-h)/(2*(n-gS)); alpha=4/(1+2*delta)
        Hs=sum((n+h)/2-g for g in gab.values()); wbar=sum(w)/len(w)
        out.append((S,gS,delta,alpha,max(w)/max(min(w),1e-9),Hs/wbar,len([x for x in w if x>0])))
    lim=1/(2**subsize-2)
    print(f"  {subsize}-subtuples: {len(out)}; Liouville balanced value 1/(2^{subsize}-2) = {lim:.5f}, alpha_lim = {4/(1+2*lim):.4f}")
    ds=[o[2] for o in out]; als=[o[3] for o in out]; br=[o[4] for o in out]; hw=[o[5] for o in out]; nc=[o[6] for o in out]
    print(f"  delta_S in [{min(ds):.5f},{max(ds):.5f}]  alpha in [{min(als):.4f},{max(als):.4f}]  balance max/min in [{min(br):.3f},{max(br):.3f}]  H/wbar in [{min(hw):.3f},{max(hw):.3f}]  nonzero cuts {min(nc)}..{max(nc)} of {2**(subsize-1)-1}")
    return out

if __name__=='__main__':
    x0,x1,x2=[1,0,0],[0,1,0],[0,0,1]
    print("== four-mark (q=27,k=2): x,y,x+y,x-y, eps=+")
    run(3,2,[[1,0],[0,1],[1,1],[1,2]],[1,1,1,1],4)
    print("== Astra five-mark (q=27,k=3): x,y,z,x+y,x+z, eps=+")
    run(3,3,[x0,x1,x2,[1,1,0],[1,0,1]],[1]*5,5)
    print("== referee best five (q=27,k=3): x2,x1,x1+x2,x0,x0+x1-x2 ; eps (+,+,-,+,+)")
    run(3,3,[x2,x1,[0,1,1],x0,[1,1,2]],[1,1,-1,1,1],5)
    print("== referee six (q=27,k=3), 5-subtuples")
    run(3,3,[x2,x1,[0,1,1],x0,[1,1,2],[1,2,1]],[1,1,-1,1,1,1],5)
