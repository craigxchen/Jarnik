# Positive control: the algebra of nonlin_ref.py on an H built from actual rational roots.
import sympy as sp, itertools, random
random.seed(1)
a,b=sp.symbols('a b')
for trial in range(3):
    xs=[sp.Rational(random.randint(-9,9),random.randint(1,5)) for _ in range(8)]
    ys=[sp.Rational(random.randint(1,9),random.randint(1,5)) for _ in range(8)]
    rho=[x+sp.I*y for x,y in zip(xs,ys)]
    H=sp.Matrix([[sp.im(sp.expand(r**k)) for k in (1,2,3)] for r in rho])  # u=e1, 2xy column -> v=e2/2, w=e3
    g=7; hg=H.row(g)
    Gi=sp.Matrix.vstack(sp.Matrix([[0,1,0]]),sp.Matrix([[0,0,1]]),hg); G=Gi.inv(); H2=H*G
    unew=Gi*sp.Matrix([1,0,0]); unew=unew/unew[2]
    sub={a:unew[0],b:unew[1]}
    u=sp.Matrix([a,b,1]); y=[(H2.row(j)*u)[0] for j in range(8)]
    mon=lambda h:[h[0]**2,h[1]**2,h[2]**2,2*h[0]*h[1],2*h[0]*h[2],2*h[1]*h[2]]
    A=sp.Matrix([mon(list(H2.row(j))) for j in range(8)])
    lk=A.T.nullspace(); comp=[sum(z[j]*y[j]**4 for j in range(8)) for z in lk]
    sel=next(T for T in itertools.combinations(range(8),6) if A.extract(list(T),list(range(6))).rank()==6)
    q=A.extract(list(sel),list(range(6))).LUsolve(sp.Matrix([y[j]**4 for j in sel]))
    Q=sp.Matrix([[q[0],q[3],q[4]],[q[3],q[1],q[5]],[q[4],q[5],q[2]]])
    F=sp.Matrix([[1,0,-a],[0,1,-b]]); det=(F*Q*F.T).det()
    j,k=0,1; z=(y[k]*H2.row(j)-y[j]*H2.row(k)).T
    Bjk=(y[j]-y[k])**2*y[j]**2*y[k]**2+sp.Rational(1,3)*(z.T*Q*z)[0]
    lam=1/hg.dot(sp.Matrix([1,0,0]))  # scaling: y_g normalised to 1
    true=lam**6*ys[0]**2*ys[1]**2*abs(rho[0]-rho[1])**2
    print([sp.simplify(c.subs(sub)) for c in comp], sp.simplify(det.subs(sub)), sp.simplify(Bjk.subs(sub)-true), [sp.simplify(y[i].subs(sub)-lam*ys[i]) for i in range(3)])
