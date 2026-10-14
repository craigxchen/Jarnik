from rlib import *
import pickle, sympy as sp, random
A,H,dimJ,W=pickle.load(open('m8_a.pkl','rb'))
a,b=1,2; c=a^b
Ta,Tb=set(T[a]),set(T[b])
def ind(S): return tuple(1 if x in S else 0 for x in range(1,8))
p1,p2,p3=ind(Ta&Tb),ind(Ta-Tb),ind(Tb-Ta)
def width(v):
    vals=[sum(x*y for x,y in zip(k,v)) for k in A]; return max(vals)-min(vals), len(set(vals))-1
add=lambda *vs: tuple(sum(t) for t in zip(*vs))
pts={'P1':p1,'P2':p2,'P3':p3,'Q12':add(p1,p2),'Q13':add(p1,p3),'Q23':add(p2,p3),'Z':add(p1,p2,p3)}
print({k:width(v) for k,v in pts.items()})
print("c_a=p1+p2:",add(p1,p2)==col(a)," c_b=p1+p3:",add(p1,p3)==col(b)," c_ab=p2+p3:",add(p2,p3)==col(c))
# ternary octics with multiplicity >= 8-h at the 7 points (coords t1,t2,t3 wrt p1,p2,p3)
t=sp.symbols('t1:4'); mons=list(sp.itermonomials(t,8)); mons=[m for m in mons if sp.Poly(m,*t).total_degree()==8]
cs=sp.symbols('c0:%d'%len(mons)); F=sum(ci*m for ci,m in zip(cs,mons))
coords={'P1':(1,0,0),'P2':(0,1,0),'P3':(0,0,1),'Q12':(1,1,0),'Q13':(1,0,1),'Q23':(0,1,1),'Z':(1,1,1)}
eqs=[]
for name,pt in coords.items():
    mult=8-width(pts[name])[0]
    # all partial derivatives of order mult-1 vanish at pt  <=> multiplicity>=mult (Euler, homogeneous)
    for al in sp.itermonomials(t,mult-1):
        if sp.Poly(al,*t).total_degree()!=mult-1 and al!=1: continue
        if mult-1==0 and al!=1: continue
        d=F
        for var in t:
            e=sp.degree(al,var)
            if e: d=sp.diff(d,var,e)
        eqs.append(d.subs(dict(zip(t,pt))))
M=sp.Matrix([[sp.diff(e,ci) for ci in cs] for e in eqs])
print("octic monomials",len(mons),"conditions",len(eqs),"rank",M.rank(),"=> dim of system",len(mons)-M.rank())
# direct TC check: for random integer u in U_l, l_u^8 in J(A)_8 iff ev(l_u^8) in span ev(R_<=7) on A (exact)
random.seed(7)
ms=[m for m in monos_le(7,7) if max(m)<=2]
base=[[evmono(m,k) for m in ms] for k in A]
r0=fmpz_mat(base).rank()
for trial in range(3):
    co=[random.randint(-5,5) for _ in range(3)]
    u=tuple(co[0]*x+co[1]*y+co[2]*z for x,y,z in zip(p1,p2,p3))
    M2=[row+[sum(x*y for x,y in zip(k,u))**8] for row,k in zip(base,A)]
    print("u",u,"in TC(A,8):",fmpz_mat(M2).rank()==r0)
# control: a random u outside
u=(1,2,0,-1,3,0,1)
M2=[row+[sum(x*y for x,y in zip(k,u))**8] for row,k in zip(base,A)]
print("control u",u,"in TC:",fmpz_mat(M2).rank()==r0)
