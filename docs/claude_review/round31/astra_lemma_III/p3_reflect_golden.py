from p3_reflection import *
import sympy as sp
ctx=fmpz_mpoly_ctx.get(['k','v']);kf,vf=ctx.gens();H,U,W,raw,Im=data(ctx)
R=[canonical(raw[j],ctx)[0] for j in (0,1,4,8)]
cr,g=cramer(R,(0,1,2),ctx)
K=sp.QQ.algebraic_field(sp.sqrt(5)); phi=(1+sp.sqrt(5))/2
x=sp.Symbol('x')
f=sp.Poly(x**2-x-1,x)
def at(poly):
    a=sp.Poly(0,x,domain=K)
    for (i,j),c in (poly % (kf**2-kf-1)).to_dict().items():
        e=sp.rem(sp.Poly(x**i,x),f).as_expr().subs(x,phi)
        a+=sp.Poly.from_dict({(j,):int(c)*e},(x,),domain=K)
    return a
ex=extra(ctx,U,W,cr)
J=[at(det(R)),at(cr[1]*cr[3]-cr[2]**2)]
for p in ex:J.append(at(p))
g=J[0]
for h in J[1:]:g=g.gcd(h)
print("good v",sp.factor(g.as_expr(),extension=sp.sqrt(5)),flush=True)
for h in cr:print('cr',sp.rem(at(h),g).as_expr(),flush=True)
