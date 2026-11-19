import sympy as sp
t,s,u=sp.symbols('t s u',real=True)
A=sp.Matrix([[t**2+1,t],[t,1]])
for n in [1,2,3]:
    An=A**n
    def Gb(seed,sign):
        v=An*sp.Matrix([1,seed]); return sp.expand(v[0]+sign*sp.I*v[1])
    R=sp.expand(sp.resultant(Gb(s,1),Gb(u,-1),t)); re,im=sp.expand(sp.re(R)),sp.expand(sp.im(R))
    print("n",n,"Im factors:",sp.factor_list(im))
    print("  Re at s=u:",sp.factor(re.subs(u,s)))
    print("  Re at u=1/s:",sp.factor(sp.together(re.subs(u,1/s))))
