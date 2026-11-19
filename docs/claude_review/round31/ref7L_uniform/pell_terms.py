"""Term budget of Astra's five-mark Pell ratios at fixed n (in t): b(n). Composition with t^k keeps b, multiplies W,h by k."""
import sympy as sp
t=sp.Symbol('t'); I=sp.I
for n in [1,2,3,4]:
    A=sp.Matrix([[t**2+1,t],[t,1]])**n
    G=[];Gb=[]
    for s in [1,3,5,8]:
        v=A*sp.Matrix([1,s]); G.append(sp.expand(v[0]+I*v[1])); Gb.append(sp.expand(v[0]-I*v[1]))
    beta=t+2*I; betab=t-2*I
    b=0
    for j,k in [(0,1),(0,2),(1,2),(0,3)]:
        num=sp.expand(betab*G[j]*G[k]); den=sp.expand(beta*Gb[j]*Gb[k])
        b=max(b,len(sp.Add.make_args(num)),len(sp.Add.make_args(den)))
    print('n',n,'b(n)',b,'h/W',sp.Rational(4*n,16*n+2))
