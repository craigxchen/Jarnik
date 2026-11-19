# Re-run Astra's accounting for an inequality  kappa*e_orb >= (K+B)^2  with general kappa, section weight a.
import sympy as sp
k,a,b0,b,K,m,s,e=sp.symbols('k a b0 b kappa m s e')
F=K*(1-2*b+b*b0)+K*a*k*(b-b0)+(1+b-b0-a*k)**2-2*a**2*sp.binomial(k,2)
unit=(K-2)*a**2                      # contribution of one transverse unit of intersection
# baseline remainder (should be (2a-a^2) S): kappa*(4-2am) - (8-4am) plus 4 half-fibres
rem=sp.expand(K*(4-2*a*m)-(8-4*a*m)+4*(K*(-1+a*m/2)-(-2+a*m)))
print("baseline remainder besides S-term:",rem)
Fk=sp.expand(sp.expand_func(F.subs({b0:0,b:0})))
phi=sp.expand(unit*k*(k-1)/2-Fk)
print("large centre F(k,0,0) =",Fk,";  deduction phi(k) =",sp.factor(phi))
# per layer, delta=E=0, using I=(m-1)S/2, S<=sum(s_h-m/4):
per=sp.expand(unit*((m-1)*(s-m/4)/2-s*(s-1)/2)+(2*a-a**2)*(s-m/4)+K+1-2*a*s+a**2*s)
target=K+1+a*(K*a-4)*m/8-unit*(s-m/2)**2/2
print("per-layer bound identity holds:",sp.simplify(per-target)==0)
F1=sp.expand(F.subs({k:2,b0:0,b:a})); F2=sp.expand(F.subs({k:2,b0:a,b:a}))
def check(Kv,av):
    f1=F1.subs({K:Kv,a:av}); f2=F2.subs({K:Kv,a:av}); u=unit.subs({K:Kv,a:av})
    chain_ok=all(f1+(ee-1)*f2<=ee*u for ee in range(2,200))
    coef=sp.Rational(av)*(Kv*av-4)/8
    mmax=sp.floor((Kv+1)/(-coef)) if coef<0 else None
    print(f"kappa={Kv} a={av}: unit={u}, F1={f1}, F2={f2}, chains ok={chain_ok}, m-coefficient={coef}, bound m<={mmax}")
for Kv,av in [(3,sp.Rational(2,3)),(3,1),(4,sp.Rational(1,2)),(4,1),(5,sp.Rational(2,5)),(6,sp.Rational(1,3)),(10,sp.Rational(1,5)),(sp.Rational(7,2),sp.Rational(4,7))]:
    check(Kv,av)
