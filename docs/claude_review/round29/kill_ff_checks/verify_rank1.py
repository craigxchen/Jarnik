import sympy as sp, itertools, sys
t=sp.symbols('t')
def check(y,K,x0=sp.Rational(1,3)):
    D=len(y); roots=[x0+sp.I*v for v in y]
    Q=sp.expand(sp.prod([t-r for r in roots]))
    zs=[]
    for F in itertools.product([0,1],repeat=D):
        z=sp.expand(sp.prod([(t-sp.conjugate(roots[k])) if F[k] else (t-roots[k]) for k in range(D)]))
        dz=sp.expand(z-Q)
        if dz==0 or sp.degree(dz,t)<=K: zs.append(z)
    pw=max(sp.degree(sp.expand(a-b),t) for a,b in itertools.combinations(zs,2))
    N=sp.Poly(sp.expand(Q*sp.prod([t-sp.conjugate(r) for r in roots])),t)
    sqf=sp.degree(sp.gcd(N,N.diff(t)))==0
    return len(zs),pw,sqf
print("D=7,c=1 (K=4):",check([1,2,3,4,5,-6,-7],4))
print("D=8,c=1 (K=5):",check([1,2,3,4,-5,-6,-7,8],5))
