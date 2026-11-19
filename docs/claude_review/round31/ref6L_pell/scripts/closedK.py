import sympy as sp
c,s=sp.cos(sp.pi/8),sp.sin(sp.pi/8)
seeds=[(1,0),(0,1),(2,-1)]
al=[x*c+y*s for x,y in seeds]; be=[-x*s+y*c for x,y in seeds]
u=[sp.nsimplify(b/a) for a,b in zip(al,be)]
u=[sp.simplify(b/a) for a,b in zip(al,be)]
S=[0,u[0]+u[1],u[0]+u[2],u[1]+u[2]]
vals=[sp.N(t,30) for t in S]; print(vals)
sp_=max(vals)-min(vals)
K3=8*sp.N((al[0]*al[1]*al[2])**2,30)*sp_**3
print("K^3=",K3,"K=",K3**(sp.Rational(1,3)))
print("nsimplify K^3:",sp.nsimplify(K3,[sp.sqrt(2)],tolerance=1e-20))
print("CG bounds: (40+40/3 sqrt10)^(1/3)=",sp.N((40+sp.Rational(40,3)*sp.sqrt(10))**sp.Rational(1,3)),
 " (40+20 sqrt5)^(1/3)=",sp.N((40+20*sp.sqrt(5))**sp.Rational(1,3)), " sqrt2=",sp.N(sp.sqrt(2)))
