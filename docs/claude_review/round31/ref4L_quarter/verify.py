import sys, gmpy2
from itertools import combinations
from math import prod, gcd, isqrt
from fractions import Fraction
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from mpmath import mp, mpf, atan2, log, mpc
mp.dps=200
P=int(sys.argv[1]); q=int(sys.argv[2])
x=(-27,-270,135,120,40)
s=[sum(prod(c) for c in combinations(x,k)) for k in range(6)]
s1,s2,s3,s4,s5=s[1:]
assert s4==0
l=isqrt(-s1*s5); assert l*l==-s1*s5 and l==97200
assert (3*s5)%2==0; b0=3*s5//2
E=list(combinations(range(5),2))
b={e:x[e[0]]*x[e[1]]*(s1-x[e[0]]-x[e[1]]) for e in E}
J={e:P+q*b[e] for e in E}
for e in E: assert gmpy2.is_prime(J[e],50) and J[e]%4==1
assert len(set(J.values()))==10
# q divisible by primes of l*b0*prod x
for p in (2,3,5): assert q%p==0
h=[xi**3-s1*xi**2+s2*xi for xi in x]
U=[P+q*hi for hi in h]; B=b0*q*q; L=l*q
V=[(P*U[i]+B*x[i], L*U[i]) for i in range(5)]
for v in V: assert gcd(*v)==1
n=[prod(J[e] for e in E if i in e) for i in range(5)]
N=prod(J.values())
det=lambda u,v:u[0]*v[1]-u[1]*v[0]
cof={}
for (i,j) in E:
    g=gcd(n[i],n[j]); assert g==J[(i,j)]
    d=det(V[i],V[j]); assert d==L*B*(x[i]-x[j])*J[(i,j)]
    assert d%g==0; cof[(i,j)]=d//g
print("ten primes J_e = P+q b_e, all = 1 mod 4, distinct; columns primitive; all ten minors = LB(x_i-x_j)J_ij exactly")
print("cofactors t_ij = l*b0*q^3*(x_i-x_j); l*b0*q^3 =",L*B, " max|t| =",max(abs(c) for c in cof.values()))
# index one: valuation matrix rows = marks, columns = primes (endpoint marking), anchored (subtract row 5)
M=Matrix([[1 if i in e else 0 for e in E] for i in range(5)])
A=Matrix([[M[i,k]-M[4,k] for k in range(10)] for i in range(4)])
print("Smith form of anchored valuation matrix:",[smith_normal_form(A)[i,i] for i in range(4)])
# defects
for i in range(5):
    nv=V[i][0]**2+V[i][1]**2; dl=nv-n[i]
    lead=P*q**3*s5*x[i]*(s2*x[i]-s1*s2-s3)
    # exact q^4 coefficient
    others=[b[e] for e in E if i in e]
    e4=prod(others)
    c4=(b0**2*x[i]**2+l*l*h[i]**2-e4)
    assert dl==lead+c4*q**4, i
    print(f"mark {i}: delta={dl:.3e}  log|delta|/log n_i = {float(log(abs(mpf(dl)))/log(mpf(n[i]))):.4f}  sign {'+' if dl>0 else '-'}  log n_i/log N={float(log(mpf(n[i]))/log(mpf(N))):.4f}")
print("exact residual = leading P q^3 term + (const) q^4 verified for all five marks")
ang=[atan2(mpf(v[1]),mpf(v[0])) for v in V]
diam=max(ang)-min(ang)
print("diam of star angles =",mp.nstr(diam,6)," log diam/log N =",mp.nstr(log(diam)/log(mpf(N)),6))
print("diam*N^(3/10) =",mp.nstr(diam*mpf(N)**(mpf(3)/10),8)," predicted l*b0*q^3*max|x_i-x_j| =",L*B*405)
print("log(diam*N^(3/10))/log N =",mp.nstr(log(diam*mpf(N)**(mpf(3)/10))/log(mpf(N)),6))
# Liouville check: sin(angle_ij) vs |t_ij|/sqrt(n_i n_j/g^2)
worst=max(abs(abs(mp.sin(ang[i]-ang[j]))*mp.sqrt(mpf(n[i])*n[j])/J[(i,j)]/abs(cof[(i,j)])-1) for i,j in E)
print("max | sin(angle)/( |t|/sqrt(n_i n_j/g^2) ) - 1 | =",mp.nstr(worst,5))
# doubled ratios leading direction
xi=[mpc(v[0],v[1])/mpc(v[0],-v[1]) for v in V]
u=[xi[i]/xi[4]-1 for i in range(4)]
w=[ui/u[0] for ui in u]
print("normalized (xi_i/xi_5 -1) :",[mp.nstr(wi,8) for wi in w], " target",[Fraction(x[i]-x[4],x[0]-x[4]) for i in range(4)])
# actual height of the torus point xi_i/xi_5 = beta/conj(beta), beta = conj(V_5) V_i
print("actual heights of the model's torus points vs profile height n_i n_5/g^2:")
for i in range(4):
    re=V[i][0]*V[4][0]+V[i][1]*V[4][1]; im=det(V[4],V[i])
    c=gcd(re,im); a,bb=re//c,im//c
    H=a*a+bb*bb
    if (a-bb)%2==0: H//=2
    prof=n[i]*n[4]//J[(min(i,4),max(i,4))]**2
    print(f"  pair ({i},4): log H(xi)/log N = {float(log(mpf(H))/log(mpf(N))):.4f}; profile log(n_i n_5/g^2)/log N = {float(log(mpf(prof))/log(mpf(N))):.4f}; content c ~ N^{float(log(mpf(c))/log(mpf(N))):.4f}")
# isotropy: does J_e divide |V_i|^2 (equivalently delta_i) for marked e at i ? and dot-product divisibility
inc=sum(1 for i in range(5) for e in E if i in e and (V[i][0]**2+V[i][1]**2)%J[e]==0)
dot=sum(1 for (i,j) in E if (V[i][0]*V[j][0]+V[i][1]*V[j][1])%J[(i,j)]==0)
print("isotropic marked incidences:",inc,"of 20; pairs with J_ij | V_i.V_j:",dot,"of 10")
print("min block prime =",min(J.values())," max|delta| / min J =",max(abs(V[i][0]**2+V[i][1]**2-n[i]) for i in range(5))/min(J.values()))
