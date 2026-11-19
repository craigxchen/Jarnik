# Optimise the section weight a in Astra's ledger (IV.2 proof), k>=3 centres weight 0.
import sympy as sp
a,s,m,b0,b,k=sp.symbols('a s m b0 b k',real=True)
# per-blowup change of 3e_orb-(K+B)^2 (before re-charging a^2 per pair), centre on k sections
G=3*(1-2*b+b*b0)+3*a*k*(b-b0)+(1+b-b0-a*k)**2-3*a**2*k*(k-1)/2
Gk=sp.expand(G.subs({b0:0,b:0}))
cS=2*a-a**2+a**2*(m-1)/2      # coefficient of S after using I=(m-1)S/2
H=sp.expand(cS*(s-m/4)+Gk.subs(k,s))   # per-layer contribution after (IV.3)
print('G(k,0,0)=',sp.factor(Gk))
print('H coeffs in s:',sp.Poly(H,s).all_coeffs())
smax=sp.solve(sp.diff(H,s),s)[0]; Hmax=sp.simplify(H.subs(s,smax))
print('argmax s=',sp.simplify(smax),' max H=',sp.expand(Hmax))
amax=sp.solve(sp.diff(Hmax,a),a); print('optimal a',amax,'value',[sp.expand(Hmax.subs(a,x)) for x in amax])
for av in [sp.Rational(1,2),sp.Rational(2,3),sp.Rational(3,4),1]:
    print('a=',av,' max H =',sp.expand(Hmax.subs(a,av)),' => m <=',sp.floor(sp.solve(sp.Eq(Hmax.subs(a,av),0),m)[0]))
# check a=2/3 identity
print(sp.simplify(H.subs(a,sp.Rational(2,3))-(4-m/6-sp.Rational(2,9)*(s-m/2)**2)))
# chain values at a=2/3 and a=1 (b=a on chain)
for av in [sp.Rational(2,3),1]:
    g1=G.subs({k:2,b0:0,b:av,a:av}); g2=G.subs({k:2,b0:av,b:av,a:av})
    print('a=',av,'chain first',g1,' later',g2)
