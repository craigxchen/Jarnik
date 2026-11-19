# Exact count of V_t ∩ (generic affine 6-plane) using linear eqs + Pfaffians + ONE triangle sextic,
# then numeric check that every slice point satisfies ALL 30 triangle identities (i.e. lies on V_t).
import sympy as sp, random, itertools, sys
from sympy import Rational as R
seed=int(sys.argv[1]); random.seed(seed)
E=[(i,j) for i in range(5) for j in range(i+1,5)]; EI={e:k for k,e in enumerate(E)}
eid=lambda i,j:EI[(min(i,j),max(i,j))]
tlist=[int(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [1]*10
T=lambda i,j: tlist[eid(i,j)]*(1 if i<j else -1)
c=[R(1,sp.prod([T(i,j) for j in range(5) if j!=i])) for i in range(5)]
nn=sp.symbols('n0:10'); u,v=sp.symbols('u v')
lin=[sum(T(i,j)*c[j]*nn[eid(i,j)] for j in range(5) if j!=i) for i in range(5)]
A=sp.Matrix([[sp.Rational(c_) for c_ in [sp.Poly(l,*nn).coeff_monomial(x) for x in nn]] for l in lin])
null=A.nullspace(); assert len(null)==6
rnd=lambda: R(random.randint(-30,30),random.randint(1,7))
base=sum((rnd()*w for w in null),sp.zeros(10,1)); w1=sum((rnd()*w for w in null),sp.zeros(10,1)); w2=sum((rnd()*w for w in null),sp.zeros(10,1))
nsub=[sp.expand(base[k]+u*w1[k]+v*w2[k]) for k in range(10)]
M=sp.Matrix(5,5,lambda i,j: 0 if i==j else T(i,j)*nsub[eid(i,j)])
pfs=[]
for d in range(5):
    idx=[a for a in range(5) if a!=d]; a,b,cc,dd=idx
    pfs.append(sp.expand(M[a,b]*M[cc,dd]-M[a,cc]*M[b,dd]+M[a,dd]*M[b,cc]))
# ratios of pfaffians are constant?
q0=[p for p in pfs if p!=0][0]
props=[sp.simplify(sp.cancel(p/q0)) for p in pfs]
print("pfaffian ratios constant:",all(pr.free_symbols==set() for pr in props))
N=lambda a,b: nsub[eid(a,b)]
def Aijk(i,j,k):
    l,m=[a for a in range(5) if a not in (i,j,k)]
    return T(j,k)**2*N(j,k)*N(i,l)*N(i,m)+T(i,k)**2*N(i,k)*N(j,l)*N(j,m)-T(i,j)**2*N(i,j)*N(k,l)*N(k,m)
def Dij(i,j): return sp.prod([N(i,k)*N(j,k) for k in range(5) if k not in (i,j)])
def sext(i,j,k): return sp.expand(4*T(i,k)**2*T(j,k)**2*Dij(i,j)-Aijk(i,j,k)**2-4*(T(i,j)*T(i,k)*T(j,k))**2)
Psi=sext(0,1,2)
res=sp.Poly(sp.resultant(q0,Psi,v),u)
print("deg resultant in u:",res.degree(), " squarefree:", sp.degree(sp.gcd(res,res.diff(u)),u)==0)
# also check no solutions at infinity: top-degree forms
q2=sp.Poly(q0,u,v).homogenize(sp.Symbol('w')); 
ps=sp.Poly(Psi,u,v)
top_q=sum(cf*u**m[0]*v**m[1] for m,cf in sp.Poly(q0,u,v).terms() if sum(m)==2)
top_p=sum(cf*u**m[0]*v**m[1] for m,cf in ps.terms() if sum(m)==6)
print("common root at infinity:", sp.degree(sp.gcd(sp.Poly(top_q.subs(v,1),u),sp.Poly(top_p.subs(v,1),u)),u)>0)
roots=sp.Poly(res,u).nroots(n=40,maxsteps=200)
pts=[]
for r in roots:
    qv=sp.Poly(q0.subs(u,r),v); cand=qv.nroots(n=40)
    best=min(cand,key=lambda vv: abs(complex(Psi.subs({u:r,v:vv}))))
    pts.append((r,best))
allsext=[(i,j,k) for (i,j) in E for k in range(5) if k not in (i,j)]
worst=0; minn=1e99
for (r,vv) in pts:
    vals=[complex(x.subs({u:r,v:vv})) for x in nsub]; minn=min(minn,min(abs(x) for x in vals))
    scal=max(abs(x) for x in vals)**6
    for tri in allsext:
        worst=max(worst,abs(complex(sext(*tri).subs({u:r,v:vv})))/scal)
print("points:",len(pts)," max rel residual of all 30 triangle identities:",worst," min |n_e|:",minn)
