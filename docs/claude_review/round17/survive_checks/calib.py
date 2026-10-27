# Toy calibration of Theorem 3.1 (exact witnesses in Z[zeta]; LP Chebyshev values = evidence).
# usage: python3 calib.py q L mmax
import sys, time, math, itertools, numpy as np
from flint import fmpz_mat, nmod_mat
from scipy.optimize import linprog
from slib import box
q,L,mmax=map(int,sys.argv[1:4]); mmin=int(sys.argv[4]) if len(sys.argv)>4 else 1; MAXROWS=60
def monos(q,m):
    out=[]
    def rec(x,rem,cur):
        if x==q: out.append(tuple(cur)); return
        for e in (0,1,2):
            if e<=rem: cur.append(e); rec(x+1,rem-e,cur); cur.pop()
    rec(0,m-1,[]); return out
def cyc_mul(a,b):  # product in Z[x]/(x^q-1)
    c=np.zeros(q,dtype=object)
    for i in range(q):
        if a[i]!=0:
            c=c+a[i]*np.roll(b,i)
    return c
def is_zero_zeta(v):  # element of Z[x]/(x^q-1) is 0 in Z[zeta] iff all coeffs equal
    return all(v[t]==v[0] for t in range(q))
def to_complex(v,s=1):
    return sum(complex(int(v[t]))*np.exp(2j*np.pi*s*t/q) for t in range(q))
A=box(q,L); U=len(A); Aarr=np.array(A,dtype=np.int64)
print(f"q={q} L={L} |A|={U}",flush=True)
# l_1(k), l_{-1}(k) as elements of Z[x]/(x^q-1)
l1=[];lm1=[]
for k in A:
    a=np.zeros(q,dtype=object); b=np.zeros(q,dtype=object)
    for x in range(q):
        a[x%q]+=k[x]; b[(-x)%q]+=k[x]
    l1.append(a); lm1.append(b)
for m in range(mmin,mmax+1):
    if m>=q: break
    t0=time.time()
    al=monos(q,m)
    E=np.array([[int(np.prod([k[x]**e for x,e in enumerate(a)])) for a in al] for k in A],dtype=np.int64)
    p1=2**31-1
    r0=nmod_mat(E.T.tolist(),p1).rank() if E.size else 0
    D=U-r0
    if D==0:
        print(f" m={m}: r0={r0} dimV=0",flush=True); continue
    # exact integer basis of V = {c: E^T c = 0}
    X,nul=fmpz_mat(E.T.tolist()).nullspace()
    B=fmpz_mat([[X[r,cidx] for r in range(U)] for cidx in range(nul)])  # rows = basis vectors
    B=B.lll()
    rows=[[int(B[i,j]) for j in range(U)] for i in range(B.nrows())]
    # powers of l1, lm1 per point
    P1=[[None]*(m+1) for _ in range(U)]; Pm=[[None]*(m+1) for _ in range(U)]
    one=np.zeros(q,dtype=object); one[0]=1
    for idx in range(U):
        P1[idx][0]=one; Pm[idx][0]=one
        for e in range(1,m+1):
            P1[idx][e]=cyc_mul(P1[idx][e-1],l1[idx]); Pm[idx][e]=cyc_mul(Pm[idx][e-1],lm1[idx])
    best=None; inTC=True; nz_i=set()
    rows=sorted(rows,key=lambda c:sum(abs(v) for v in c))[:MAXROWS]
    for c in rows:
        l1n=sum(abs(v) for v in c)
        for i in range(m+1):
            alpha=np.zeros(q,dtype=object)
            for idx,cv in enumerate(c):
                if cv: alpha=alpha+cv*cyc_mul(P1[idx][i],Pm[idx][m-i])
            if not is_zero_zeta(alpha):
                inTC=False; nz_i.add(i)
                conj=[abs(to_complex(alpha,s)) for s in range(1,q)]
                bound=2.0**(-m)*math.comb(m,i)*max(conj)/l1n
                if best is None or bound>best[0]:
                    best=(bound,l1n,i,max(conj),min(conj),max(abs(v) for v in c))
    msg=f" m={m}: r0={r0} dimV={D} circle_in_TC_V={inTC} nonzero_i={sorted(nz_i)}"
    if best:
        bound,l1n,i,cmax,cmin,cinf=best
        # LP Chebyshev at a few theta (evidence)
        chebs=[]
        for th in [0.0,0.3,np.pi/q,1.1]:
            u=np.cos(th+2*np.pi*np.arange(q)/q)
            f=(Aarr@u)**m
            nv=len(al)
            # vars: coeffs (nv), t ; minimize t ; f - E a <= t ; -(f - E a) <= t
            cobj=np.zeros(nv+1); cobj[-1]=1
            Aub=np.vstack([np.hstack([-E,-np.ones((U,1))]),np.hstack([E,-np.ones((U,1))])])
            bub=np.concatenate([-f,f])
            res=linprog(cobj,A_ub=Aub,b_ub=bub,bounds=[(None,None)]*nv+[(0,None)],method='highs')
            chebs.append(res.fun if res.status==0 else float('nan'))
        msg+=f" | best witness: |c|_1={l1n} |c|_inf={cinf} i={i} max|conj|={cmax:.4g} min|conj|={cmin:.4g} Thm3.1 bound (with max|conj|) on Cheb={bound:.4e} LP Cheb(theta)={[f'{x:.4e}' for x in chebs]}"
    print(msg+f" ({time.time()-t0:.1f}s)",flush=True)
