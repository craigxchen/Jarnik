# sdc7, L=2, m=7,8: exact V(A,m), AGL-invariance of the certificate, <c,M_i> in Z[zeta], and LP Chebyshev (evidence).
import numpy as np, math
from flint import fmpz_mat
from scipy.optimize import linprog
from slib import box
q,L=7,2
A=box(q,L); U=len(A); idx={k:i for i,k in enumerate(A)}; Aarr=np.array(A)
def monos(m):
    out=[]
    def rec(x,rem,cur):
        if x==q: out.append(tuple(cur)); return
        for e in (0,1,2):
            if e<=rem: cur.append(e); rec(x+1,rem-e,cur); cur.pop()
    rec(0,m-1,[]); return out
zeta=np.exp(2j*np.pi/q)
for m in (7,8):
    al=monos(m)
    E=[[int(np.prod([k[x]**e for x,e in enumerate(a)])) for a in al] for k in A]
    X,nul=fmpz_mat(np.array(E).T.tolist()).nullspace()
    print(f"m={m}: |A|={U} dim V={nul}")
    if nul==0: continue
    B=fmpz_mat([[X[r,c] for r in range(U)] for c in range(nul)]).lll()
    for row in range(min(nul,3)):
        c=[int(B[row,t]) for t in range(U)]
        # AGL invariance check (up to sign)
        def act(k,s,b): return tuple(k[(pow(s,-1,q)*(y-b))%q] for y in range(q))
        inv={}
        for s in range(1,q):
            for b in range(q):
                cg=[0]*U
                for i,k in enumerate(A): cg[idx[act(k,s,b)]]=c[i]
                inv[(s,b)]= 1 if cg==c else (-1 if cg==[-v for v in c] else 0)
        l1=np.array([sum(k[x]*zeta**x for x in range(q)) for k in A]); lm=np.conj(l1)
        alphas=[complex(np.dot(c,l1**i*lm**(m-i))) for i in range(m+1)]
        nz=[i for i,a in enumerate(alphas) if abs(a)>1e-6]
        cl1=sum(abs(v) for v in c)
        print(f"  basis vec {row}: |c|_1={cl1} |c|_inf={max(abs(v) for v in c)} AGL character values={sorted(set(inv.values()))} "
              f"nonzero <c,M_i> at i={nz} values={[round(alphas[i].real,6)+1j*round(alphas[i].imag,6) for i in nz]}")
        # Phi on circle and Chebyshev LP
        for th in (0.0,0.7):
            u=np.cos(th+2*np.pi*np.arange(q)/q); f=(Aarr@u)**m
            Phi=float(np.dot(c,f))
            Em=np.array(E,dtype=float); nv=len(al)
            cobj=np.zeros(nv+1); cobj[-1]=1
            Aub=np.vstack([np.hstack([-Em,-np.ones((U,1))]),np.hstack([Em,-np.ones((U,1))])])
            res=linprog(cobj,A_ub=Aub,b_ub=np.concatenate([-f,f]),bounds=[(None,None)]*nv+[(0,None)],method='highs')
            print(f"    theta={th}: Phi_c={Phi:.6f}  |Phi_c|/|c|_1={abs(Phi)/cl1:.6e}  LP Cheb_A(u_theta)={res.fun:.6e}")
