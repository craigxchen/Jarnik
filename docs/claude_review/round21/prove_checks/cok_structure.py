"""For the record clusters: exponent profile a_x, B = rows 2(a_x - a_1), rank, s, and the group Cok = Z^s / B^T Z^(M-1)
(dual of H, Section 4). Reports the order D_j of e_j in Cok (inf = free coordinate) and whether coordinates are pairwise
independent. Exact integer linear algebra (sympy Smith normal form via HNF-based order computation)."""
import math
from sympy import Matrix, ZZ, factorint
from sympy.matrices.normalforms import smith_normal_form, hermite_normal_form
from glib import windows, ggcd, gnorm, gprime, gmul, gconj, gdivmod_round
def val(z, pi):
    k=0
    while True:
        q,r=gdivmod_round(z,pi)
        if r!=(0,0): return k
        z=q; k+=1
def order_in_cok(BT, v):
    # order of v in Z^s / L, L = column lattice of BT (s x (M-1)). inf if v not in Q-span.
    s=BT.rows
    Mx=BT.row_join(v)
    if Mx.rank()>BT.rank(): return math.inf
    for d in range(1,10000):
        # test d*v in L: solve BT x = d v over Z via HNF of [BT | d v] lattice equality
        L1=hermite_normal_form(BT); L2=hermite_normal_form(BT.row_join(d*v))
        if L1.shape==L2.shape and L1==L2: return d
    return None
for n,k,label in [(1176852625,10,'10-pt C=8.24'),(537606725,8,'golden 8-pt'),(235370525,7,'sporadic 7-pt'),(2163838625,5,'5-pt C=2.26')]:
    res,f=windows(n,k,9.0)
    k_,C,W=min([r for r in res if r[0]==k], key=lambda r:r[1])
    g=W[0]
    for z in W[1:]: g=ggcd(g,z)
    us=[gdivmod_round(z,g)[0] for z in W]
    Nn=n//gnorm(g); fn=factorint(Nn); ps=sorted(fn)
    A=[[val(u,gprime(p)) for p in ps] for u in us]
    B=Matrix([[2*(A[x][j]-A[0][j]) for j in range(len(ps))] for x in range(1,k)])
    BT=B.T
    D=[order_in_cok(BT, Matrix([1 if i==j else 0 for i in range(len(ps))])) for j in range(len(ps))]
    snf=smith_normal_form(BT, domain=ZZ)
    inv=[snf[i,i] for i in range(min(snf.shape))]
    print(label,"n=%d k=%d C=%.4f N'=%d s=%d rank B=%d"%(n,k,float(C),Nn,len(ps),B.rank()))
    print("   primes",ps,"exponents",[fn[p] for p in ps])
    print("   Smith invariants of B^T:",inv,"  orders D_j of e_j in Cok:",D)

# Pairwise-independence test of Corollary 4.3 (all D_j = inf case): V cap span(e_i,e_j) != 0 iff rank([B; e_i; e_j]) < rank B + 2
print("\nCorollary 4.3 hypothesis (no row-space vector supported on <= 2 coordinates):")
for n,k,label in [(537606725,8,'golden 8-pt'),(2163838625,5,'5-pt C=2.26')]:
    res,f=windows(n,k,9.0)
    k_,C,W=min([r for r in res if r[0]==k], key=lambda r:r[1])
    g=W[0]
    for z in W[1:]: g=ggcd(g,z)
    us=[gdivmod_round(z,g)[0] for z in W]
    Nn=n//gnorm(g); fn=factorint(Nn); ps=sorted(fn); s=len(ps)
    A=[[val(u,gprime(p)) for p in ps] for u in us]
    B=Matrix([[2*(A[x][j]-A[0][j]) for j in range(s)] for x in range(1,k)])
    r=B.rank(); bad=[]
    for i in range(s):
        for j in range(i+1,s):
            E=Matrix([[1 if t==i else 0 for t in range(s)],[1 if t==j else 0 for t in range(s)]])
            if B.col_join(E).rank()<r+2: bad.append((ps[i],ps[j]))
    print("  ",label,"rank",r,"s",s,"dependent pairs:",len(bad),"of",s*(s-1)//2, bad[:6])
