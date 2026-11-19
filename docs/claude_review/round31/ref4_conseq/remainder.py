# One-parameter remainder lemma sanity check, and the two-parameter failure of the same argument.
import sympy as sp, math
t=sp.symbols('t')
G=sp.Poly(t**2+1,t); F=sp.Poly(3*t**5+t+7,t)
q,r=sp.div(F,G); den=sp.ilcm(*[sp.fraction(c)[1] for c in r.all_coeffs()+q.all_coeffs()])
print("remainder",r.as_expr(),"denominator",den)
hits=[T for T in range(-10**5,10**5) if F.eval(T)%G.eval(T)==0]
print("integers T with G(T)|F(T):",hits, " bound |R'(T)|<|G(T)| for |T|>",
      max(abs(T) for T in range(1,10**4) if abs((den*r).eval(T))>=G.eval(T)) if any(abs((den*r).eval(T))>=G.eval(T) for T in range(1,10**4)) else 0)
# Bezout constant for a coprime column (X,Y): gcd(X(T),Y(T)) divides Res(X,Y)
X=sp.Poly(t**3+2*t+5,t); Y=sp.Poly(2*t**2-t+3,t); R=sp.resultant(X.as_expr(),Y.as_expr(),t)
gs=set(math.gcd(int(X.eval(T)),int(Y.eval(T))) for T in range(-3000,3000))
print("Res=",R," observed gcds:",sorted(gs)," all divide Res:",all(R%g==0 for g in gs))
# two parameters: G=T1, F=T2^2+1. G | F on {(a, b0+ka): b0^2=-1 mod a}, Zariski dense, but T1 does not divide T2^2+1.
pts=[(a,b) for a in range(2,60) for b in range(0,200) if (b*b+1)%a==0]
print("two-parameter points with T1 | T2^2+1 (sample):",pts[:8],"... count",len(pts))
