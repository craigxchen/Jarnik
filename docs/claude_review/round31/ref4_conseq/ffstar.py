# Independent check: genuine C(t)/Q(i)(t) clusters map to the Part IV C[t] data
# (n_i | N, gcd(n_i,n_j) | D_ij, coprime real column entries, deg D_ij <= (a_i+a_j)/2 - L/4 + delta
#  with delta = max deg(z_i-z_j) - D/2), and (IV.2) is then evaluated on them.
import itertools, random, sympy as sp
from fractions import Fraction as Fr
t=sp.symbols('t'); I=sp.I
def conjp(p): return sp.Poly([sp.conjugate(c) for c in p.all_coeffs()],t)
def re_im(p):
    c=p.all_coeffs(); X=sp.Poly([sp.re(a) for a in c],t); Y=sp.Poly([sp.im(a) for a in c],t); return X,Y
def deg(p): return -10**9 if p.is_zero else p.degree()
def gcdp(a,b): return sp.gcd(a,b)   # over Q(i) since coefficients gaussian rationals
def run(roots, eps_list, label):
    # roots rho_k = x_k + i y_k (y_k != 0), z_eps = prod (t - x_k - i eps_k y_k)
    D=len(roots); L=2*D
    f=[sp.Poly(t-(x+I*y),t) for x,y in roots]          # orientation: f_k
    N=sp.Poly(1,t)
    for fk in f: N=N*fk*conjp(fk)
    zs=[]; Ps=[]; alpha=[]; marks=[]
    for eps in eps_list:
        z=sp.Poly(1,t); P=sp.Poly(1,t); mk=[]
        for k,(x,y) in enumerate(roots):
            z=z*sp.Poly(t-(x+I*eps[k]*y),t)
            if eps[k]==1: P=P*f[k]; mk.append(k)
        zs.append(z); Ps.append(P); alpha.append(2*deg(P)); marks.append(set(mk))
    m=len(zs)
    # star identity check: z_i = Fbar * P_i / Pbar_i, Fbar = prod conj f_k  (all monic, so phase c_i = 1)
    Fb=sp.Poly(1,t)
    for fk in f: Fb=Fb*conjp(fk)
    for z,P in zip(zs,Ps):
        assert (z*conjp(P)-Fb*P).is_zero
    cols=[re_im(P) for P in Ps]
    for X,Y in cols: assert deg(gcdp(X,Y))==0, "column not primitive"
    ns=[P*conjp(P) for P in Ps]
    for n in ns: assert sp.rem(N,n).is_zero
    delta=max(deg(zs[i]-zs[j]) for i in range(m) for j in range(i+1,m))-Fr(D,2)
    E=0
    for i,j in itertools.combinations(range(m),2):
        X1,Y1=cols[i]; X2,Y2=cols[j]; Dij=X1*Y2-X2*Y1
        assert not Dij.is_zero
        g=gcdp(ns[i],ns[j])
        Eij=deg(g)-deg(gcdp(g,Dij)); E+=Eij
        # degree identity deg D_ij = deg(z_i-z_j) - D + (a_i+a_j)/2
        assert deg(Dij)==deg(zs[i]-zs[j])-D+Fr(alpha[i]+alpha[j],2)
        assert deg(Dij)<=Fr(alpha[i]+alpha[j],2)-Fr(L,4)+delta
    # layers: roots of N, i.e. f_k and conj f_k; s for f_k = #marks containing k, for conj = m - that
    S=0
    for k in range(D):
        s=sum(1 for mk in marks if k in mk)
        S+=(Fr(s)-Fr(m,2))**2+(Fr(m-s)-Fr(m,2))**2
    lhs=Fr(3,4)*(m-24)*L+S; rhs=m*(m+3)*delta+4*E
    print(f"{label}: D={D} m={m} delta={delta} sumE={E} IV.2: lhs={lhs} rhs={rhs} ok={lhs<=rhs}")
    return E,delta
random.seed(1)
# D=4, y_k=1: choices with sum eps = 0 give Im p_1 equal -> deg(z_i-z_j)<=2=D/2 (c=0)
roots=[(sp.Rational(x,3),1) for x in random.sample(range(-20,20),4)]
eps=[e for e in itertools.product([1,-1],repeat=4) if sum(e)==0]
run(roots,eps,"D=4 balanced (c=0)")
# D=4 all 16 choices: c>0
run(roots,list(itertools.product([1,-1],repeat=4)),"D=4 all (c>0)")
# D=6, y=(1,1,1,1,2,...) random gaussian roots, take a balanced subfamily
roots6=[(sp.Rational(x,2),1) for x in random.sample(range(-20,20),6)]
eps6=[e for e in itertools.product([1,-1],repeat=6) if sum(e)==0]
run(roots6,eps6,"D=6 Im p1 equal (deg diff<=4, delta=1)")
# non-squarefree parent: repeated root pair (prime powers); dedupe equal z
rootsR=[(sp.Rational(1,3),1),(sp.Rational(1,3),1),(sp.Rational(-5,2),1),(sp.Rational(7,3),1),(sp.Rational(-4,3),1),(sp.Rational(-4,3),1)]
seen=set(); epsR=[]
for e in itertools.product([1,-1],repeat=6):
    if sum(e)!=0: continue
    key=(e[0]+e[1],e[2],e[3],e[4]+e[5])
    if key in seen: continue
    seen.add(key); epsR.append(e)
run(rootsR,epsR,"D=6 non-squarefree, Im p1 equal")
