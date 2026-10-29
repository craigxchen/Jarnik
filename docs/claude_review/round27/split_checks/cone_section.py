# Cost of the quadric-cone section G = (u2-u1)(u3-u1)-(u4-u1)(u5-u1) (ord_e = 2) over Pi_nt; beats 15/4 iff min cost < 8.
import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks')
from reglib import *
from scipy.optimize import linprog
import sympy as sp
u=sp.symbols('u0:5')
G=sp.expand((u[1]-u[0])*(u[2]-u[0])-(u[3]-u[0])*(u[4]-u[0]))
mon=[sp.Poly(t,*u).monoms()[0] for t in G.as_ordered_terms()]
print("monomials of G:",mon)
# Laurent f = G/u^m with any m (deg 2); support S = {mono - m}; cost_x is translation invariant on Pi_nt
S=[tuple(np.array(mm)-np.array([2,0,0,0,0])) for mm in mon]
S4=[s[:4] for s in S]
print("check sums zero:", all(sum(s)==0 for s in S))
r=reg(S4); print("reg(S) =",r, "(ord_e G should be 2)")
pi=pis(S4)
A_eq=[];b_eq=[];seen=set()
for B in Bs:
    C=comp(B)
    if C in seen: continue
    seen.add(B); row=[0]*30; row[idx[B]]=1; row[idx[C]]=1; A_eq.append(row); b_eq.append(1)
for a in range(1,5):
    row=[0]*30
    for B in Bs:
        if 0 in B: row[idx[B]]+=1
        if a in B: row[idx[B]]-=1
    A_eq.append(row); b_eq.append(0)
res=linprog(pi.astype(float),A_eq=np.array(A_eq,float),b_eq=np.array(b_eq,float),bounds=[(0,1)]*30,method='highs')
print("min over Pi_nt of cost =",res.fun," (beat iff < 4*reg =",4*r,")")
print("cost at x=1/2:", pi.sum()/2)
