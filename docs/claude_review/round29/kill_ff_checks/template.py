import sympy as sp, itertools
I=sp.I; u=sp.symbols('u')
y=[sp.Rational(-9,4),sp.Rational(-2,3),sp.Rational(1,4),sp.Rational(8,3)]; eps=[1,-1,-1,1]
lam=72
r=[lam*(sp.Rational(eps[k])/y[k]+I*y[k]) for k in range(4)]
roots=r+[-x for x in r]
print("roots",roots)
sets=[(0,)*8,(0,0,0,0,1,1,1,1),(0,0,1,1,0,0,1,1),(0,1,0,1,0,1,0,1),(1,0,1,0,1,0,1,0),(1,1,0,0,1,1,0,0),(1,1,1,1,0,0,0,0),(1,)*8]
zs=[sp.Poly(sp.expand(sp.prod([(u-sp.conjugate(roots[k])) if S[k] else (u-roots[k]) for k in range(8)])),u) for S in sets]
print("coeffs Gaussian integers:",all(all(sp.re(c).is_integer and sp.im(c).is_integer for c in z.all_coeffs()) for z in zs))
N=sp.expand(zs[0].as_expr()*sp.conjugate(zs[0].as_expr()))
print("same norm:",all(sp.expand(z.as_expr()*sp.conjugate(z.as_expr())-N)==0 for z in zs))
c4=[z.coeff_monomial(u**4) for z in zs]
print("deg diffs:",max(sp.degree(sp.expand(a.as_expr()-b.as_expr()),u) for a,b in itertools.combinations(zs,2)))
spread=max(abs(complex(a-b)) for a,b in itertools.combinations(c4,2))
print("max |c4 diff| (= limit chord/u^4 = chord/sqrt R):",float(spread))
# common gcd of the 8 polynomials as polys? and numeric gcd at a sample u
g=zs[0]
for z in zs[1:]: g=sp.gcd(g,z)
print("poly gcd",g)
# Gaussian integer gcd of values at u=10**6
import math
def ggcd(a,b):
    while b!=0:
        q=a/b; q=complex(round(q.real),round(q.imag)) if False else None
        break
uval=10**6+1
vals=[complex(0)]
def conjp(z): return sp.Poly([sp.conjugate(c) for c in z.all_coeffs()],u)
Ns=[sp.expand((z*conjp(z)).as_expr()) for z in zs]
print("same norm (coeff conj):",all(sp.expand(n-Ns[0])==0 for n in Ns))
# check at integer u: 8 distinct Gaussian integers with same norm; ratio chord/sqrt R
for uv in [10**6, 10**9+7]:
    vals=[complex(z.eval(uv)) for z in zs]
    Zs=[(int(sp.re(z.eval(uv))),int(sp.im(z.eval(uv)))) for z in zs]
    norms={a*a+b*b for a,b in Zs}
    R=sp.sqrt(list(norms)[0])
    import math
    ch=max(math.hypot(a[0]-b[0],a[1]-b[1]) for a,b in itertools.combinations(Zs,2))
    print(uv,"distinct",len(set(Zs)),"one norm",len(norms)==1,"chord/sqrtR",ch/math.sqrt(float(R)))
