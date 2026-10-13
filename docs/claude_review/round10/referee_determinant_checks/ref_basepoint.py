# Referee: Prop 5.1's "65% of points coincide" needs the coinciding set to contain the base point n.
# Full box B_T (no cut), V = V(B_T,m): residual of f_delta outside Poly_{<m}|_B for
#  (a) one point x<n moved (zeta_x != 1, others 1)           -> predicted in base locus if 2T < m
#  (b) point n moved alone (all zeta_x equal != 1, x<n)      -> NOT in base locus unless 2T(n-1) < m
import numpy as np, itertools, sys
out=[]
for (n,T,m) in [(4,1,3),(5,1,3),(5,2,5),(6,1,4)]:
    d=n-1
    B=np.array(list(itertools.product(range(-T,T+1),repeat=d)),dtype=float)
    mons=[a for t in range(m) for a in itertools.product(range(t+1),repeat=d) if sum(a)==t]
    P=np.stack([np.prod(B**np.array(a),axis=1) for a in mons],axis=1)
    Uu,sv,_=np.linalg.svd(P,full_matrices=False); rk=int((sv>sv[0]*1e-10).sum()); Q=Uu[:,:rk]
    def res(dl):
        f=np.exp(1j*(B@dl)); r=f-Q@(Q.T@f); return np.linalg.norm(r)/np.linalg.norm(f)
    a=res(np.array([0.9]+[0.0]*(d-1))); b=res(np.full(d,0.9))
    s=f"n={n} T={T} m={m}: (a) one point x<n moved: residual {a:.1e};  (b) base point n moved alone: residual {b:.1e}"
    print(s); out.append(s)
open(sys.argv[1],"w").write("\n".join(out)+"\n")
