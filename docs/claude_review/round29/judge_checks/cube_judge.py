# Judge's independent check of ff Theorem C (exact, sympy over Q(i)).
import itertools, sympy as sp
t=sp.symbols('t'); I=sp.I
y=[sp.Rational(-9,4),sp.Rational(-2,3),sp.Rational(1,4),sp.Rational(8,3)]
assert sum(y)==0 and sp.prod(y)==1
eps=[1,-1,-1,1]
r=[sp.nsimplify(eps[k]/y[k])+I*y[k] for k in range(4)]
rho=[r[0],r[1],r[2],r[3],-r[0],-r[1],-r[2],-r[3]]   # roots of Q: +-r_k (8 of them); conj roots are bar(rho)
# make all roots in the upper half plane: choose for each the representative with Im>0, track if conj taken
up=[]; 
for z in rho:
    up.append(z if sp.im(z)>0 else sp.conjugate(z))
assert len(set(up))==8
D=8
def zpoly(a):  # a in {0,1}^8: take root up[k] if a_k=1 else conj
    return sp.expand(sp.prod([(t-(up[k] if a[k] else sp.conjugate(up[k]))) for k in range(D)]))
Npoly=sp.expand(zpoly([1]*8)*sp.expand(sp.prod([t-sp.conjugate(u) for u in up])))
polys={a:sp.Poly(zpoly(a),t) for a in itertools.product([0,1],repeat=D)}
base=(1,)*8
# Q itself corresponds to taking rho (not up); find a with zpoly==Q
Q=sp.Poly(sp.expand(sp.prod([t-z for z in rho])),t)
a0=[a for a,p in polys.items() if p==Q]; assert len(a0)==1; a0=a0[0]
clus=[a for a,p in polys.items() if (p-Q).is_zero or (p-Q).degree()<=4]
print("cluster through Q:",len(clus))
mx=max(((polys[a]-polys[b]).degree() for a in clus for b in clus if a!=b))
print("max pairwise diff degree:",mx)
print("N squarefree:", sp.gcd(sp.Poly(Npoly,t),sp.Poly(Npoly,t).diff(t)).degree()==0, "N real:", all(sp.im(c)==0 for c in sp.Poly(Npoly,t).all_coeffs()))
