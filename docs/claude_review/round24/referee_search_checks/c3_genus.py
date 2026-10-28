import sympy as sp
m,Y,W=sp.symbols('m Y W')
c=m**3-11*m**2+39*m-13
u=-(m+3)*c/(16*(m-3)**2*(m-1)**2); v=-(m+1)*c/(4*(m-3)**3*(m-1))
# Y-relation check: resultant of W-quartic and Y + W^2/v
Wq=W**4+v*W**2-v**2*W+u*v**2
G=(Y**2-Y+u)**2+v*Y
R=sp.resultant(sp.numer(sp.together(Wq)),sp.numer(sp.together(Y+W**2/v)),W)
print('Res_W(Wquartic, Y+W^2/v) proportional to G^2? ', sp.simplify(sp.factor(R)/sp.factor(sp.numer(sp.together(G)))**2))
Gn=sp.Poly(sp.numer(sp.together(G)),Y)
print('leading coeff', sp.factor(Gn.LC()))
D=sp.factor(Gn.discriminant()); print('disc_Y =',D)
# orders at roots
fl=sp.factor_list(Gn.discriminant())
lc=Gn.LC()
for f,e in fl[1]:
    print(' factor',f,' mult',e,' lc vanishes on it:', sp.rem(sp.Poly(lc,m),sp.Poly(f,m)).is_zero)
quart=m**4+70*m**2-48*m+9
print('squarefree & coprime check: disc(c)=',sp.discriminant(c,m),' disc(quart)=',sp.discriminant(quart,m),' res(c,quart)=',sp.resultant(c,quart,m))
print('real roots of quartic:',sp.Poly(quart,m).count_roots(), ' of c:',sp.Poly(c,m).count_roots())
# behaviour at m=infinity: substitute m=1/w
w=sp.symbols('w')
Gw=sp.Poly(sp.numer(sp.together(G.subs(m,1/w))),Y)
Dw=sp.factor(Gw.discriminant()); print('disc at infinity chart:',Dw)
