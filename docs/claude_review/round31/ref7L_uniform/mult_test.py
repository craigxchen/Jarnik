"""Covering inequality with multiplicities and a common place: q_i=((t-a_i)/(t-c))^{m_i}, m=(1,2,3,1).
Params a_2..a_5,c (dominant). h=1, W=10, charged places=5. Log-det of (d/dt, d/da_2, d/da_3, d/dc)."""
import sympy as sp, itertools
t=sp.Symbol('t'); a=sp.symbols('a2:6'); c=sp.Symbol('c'); m=[1,2,3,1]
q=[((t-a[i])/(t-c))**m[i] for i in range(4)]
params=list(a)+[c]
vals={a[0]:1,a[1]:3,a[2]:7,a[3]:-2,c:5}
for ch in itertools.combinations(params,3):
    cols=[[sp.diff(sp.log(x),v) for x in q] for v in [t]+list(ch)]
    D=sp.cancel(sp.Matrix(4,4,lambda i,k:cols[k][i]).det().subs(vals))
    if D==0: print(ch,'det 0'); continue
    n,d=sp.fraction(D); n=sp.Poly(n,t); d=sp.Poly(d,t)
    print(ch,'ord_inf',d.degree()-n.degree(),'den',sp.factor(d.as_expr()))
