import sympy as sp
t,s,u=sp.symbols('t s u',real=True)
A=sp.Matrix([[t**2+1,t],[t,1]])
for n in [1,2]:
    An=A**n
    def Gb(seed,sign):
        v=An*sp.Matrix([1,seed]); return sp.expand(v[0]+sign*sp.I*v[1])
    R=sp.resultant(Gb(s,1),Gb(u,-1),t)
    R=sp.expand(R); re,im=sp.expand(sp.re(R)),sp.expand(sp.im(R))
    print("n",n,"Re:",sp.factor(re)); print("   Im:",sp.factor(im))
    sol=sp.solve([re,im],[s,u],dict=True)
    print("   real solutions:",[x for x in sol if all(v.is_real for v in x.values())][:10])
    print("   all:",sol[:10])
