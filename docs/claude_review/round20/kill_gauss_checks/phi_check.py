# Independent check of the k=4 fourth-root formula and Phi (gauss.md 4.1)
import sympy as sp
u1,v1,u2,v2,u3,v3=sp.symbols('u1 v1 u2 v2 u3 v3',real=True)
r=[u1+sp.I*v1,u2+sp.I*v2,u3+sp.I*v3]
s1=sum(r); s2=r[0]*r[1]+r[0]*r[2]+r[1]*r[2]; s3=r[0]*r[1]*r[2]
re1,im1=sp.re(sp.expand(s1)),sp.im(sp.expand(s1))
re2,im2=sp.re(sp.expand(s2)),sp.im(sp.expand(s2))
r4=(re1*im1-im2)/im1 - sp.I*im1
E1=sp.expand(s1+r4); E2=sp.expand(s2+s1*r4); E3=sp.expand(s3+s2*r4); E4=sp.expand(s3*r4)
print('Im e1 =',sp.simplify(sp.im(E1)))
print('Im e2 =',sp.simplify(sp.im(E2)))
im3=sp.together(sp.im(E3))
num,den=sp.fraction(sp.factor(im3))
Phi=sum(vv[0]*vv[1]*(uu[0]-uu[1])**2 for vv,uu in [((v1,v2),(u1,u2)),((v1,v3),(u1,u3)),((v2,v3),(u2,u3))]) + (v1+v2)*(v1+v3)*(v2+v3)*(v1+v2+v3)
print('Im e3 numerator / Phi =',sp.factor(sp.cancel(num/Phi)),' denominator=',den)
