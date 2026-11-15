# independent exact check of Theorem C
import sympy as sp, itertools
a,b,c=sp.symbols('a b c')
d=-(a+b+c)
ys=[a,b,c,d]
expr=sp.together(sum(3/y-y**3 for y in ys))
num=sp.factor(sp.numer(expr)); den=sp.factor(sp.denom(expr))
print("num",num); print("den",den)
# claimed: abcd*sum = -3(a+b)(a+c)(b+c)(1-abcd)
print("identity:", sp.simplify(a*b*c*d*sum(3/y-y**3 for y in ys) + 3*(a+b)*(a+c)*(b+c)*(1-a*b*c*d)))
I=sp.I; t=sp.symbols('t')
def cluster(yv, eps=(1,-1,-1,1)):
    r=[sp.Rational(eps[k])/yv[k]+I*yv[k] for k in range(4)]
    roots=r+[-x for x in r]  # 8 roots of Q
    Q=sp.expand(sp.prod([t-x for x in roots]))
    cnt=0; sets=[]
    maxdeg=-1
    for S in itertools.product([0,1],repeat=8):
        z=sp.expand(sp.prod([(t-sp.conjugate(roots[k])) if S[k] else (t-roots[k]) for k in range(8)]))
        dd=sp.Poly(sp.expand(z-Q),t).degree() if sp.expand(z-Q)!=0 else -1
        if dd<=4:
            cnt+=1; sets.append(S); maxdeg=max(maxdeg,dd)
    # pairwise check
    zs=[sp.expand(sp.prod([(t-sp.conjugate(roots[k])) if S[k] else (t-roots[k]) for k in range(8)])) for S in sets]
    pw=max((sp.Poly(sp.expand(x-y),t).degree() for x,y in itertools.combinations(zs,2)),default=-1)
    N=sp.expand(Q*sp.expand(sp.prod([t-sp.conjugate(x) for x in roots])))
    sqf=sp.degree(sp.gcd(sp.Poly(N,t,domain='QQ<I>'),sp.Poly(sp.diff(N,t),t,domain='QQ<I>')))==0
    realN=all(sp.im(co)==0 for co in sp.Poly(N,t).all_coeffs())
    return cnt,sets,pw,sqf,realN
for yv in [(sp.Rational(-9,4),sp.Rational(-2,3),sp.Rational(1,4),sp.Rational(8,3))]:
    print(sum(yv),sp.prod(yv))
    cnt,sets,pw,sqf,realN=cluster(yv)
    print("cluster size",cnt,"max pairwise deg",pw,"N squarefree",sqf,"N real",realN)
    for S in sets: print(S)
