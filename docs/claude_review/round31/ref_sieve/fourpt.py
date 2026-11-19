# Independent check of section 6 (four-point family): identities and exact Jacobian over Q(i)
import sympy as sp
x,m,s=sp.symbols('x m s')
I=sp.I
y=x*(1-x)/(m**2+x); z=1-x-y; u=m*y
# abc=u^2 symmetric
print('xyz-u^2 =', sp.simplify(x*y*z-u**2))
A,B,C=sp.symbols('A B C')
# cubic identity with c=1-a-b
a,b=sp.symbols('a b'); c=1-a-b
print('cubic identity:', sp.expand(1+(a-1)**3+(b-1)**3+c**3-3*(1+c)*a*b))
def roots(a,b,c):
    return [I, u+I*(a-1), u+I*(b-1), u*(1+c)/c+I*c]
perms=[(x,y,z),(x,z,y),(y,z,x)]
t=sp.symbols('t',real=True)
for P in perms:
    rs=roots(*P)
    for k in (1,2,3):
        pk=sum(r**k for r in rs)
        print('perm',P[0],'k',k,'Im p_k =',sp.simplify(sp.im(sp.expand(pk.subs({x:sp.Symbol('X',real=True),m:sp.Symbol('M',real=True)})))))
# Jacobian of log F_j, F_j=P_j(s)/conj P_j(s)
def conjr(r):  # params real
    return sp.re(r)-I*sp.im(r)
X_,M_=sp.symbols('X M',real=True)
vars_=(X_,M_,s)
rows=[]
for P in perms:
    rs=[sp.expand_complex(r.subs({x:X_,m:M_})) for r in roots(*P)]
    logF=sum(sp.log(s-r)-sp.log(s-sp.conjugate(r)) for r in rs)
    rows.append([sp.diff(logF,v) for v in vars_])
J=sp.Matrix(rows)
pt={X_:-3,M_:3,s:7}
Jv=J.subs(pt).applyfunc(lambda e: sp.nsimplify(sp.simplify(e)))
print('exact Jacobian at (-3,3,7):'); sp.pprint(Jv)
d=sp.simplify(Jv.det()); print('det =',d, ' numeric', sp.N(d))
# reduce mod 101 with i->10
mod=101
def red(e):
    e=sp.nsimplify(sp.expand(e)); re_,im_=sp.re(e),sp.im(e)
    def rr(q):
        q=sp.Rational(q); return int(q.p)*pow(int(q.q),-1,mod)%mod
    return (rr(re_)+10*rr(im_))%mod
print('mod 101 matrix:',[[red(Jv[i,j]) for j in range(3)] for i in range(3)], 'det mod 101:', red(d))
# seven roots at the point
allr=set()
for P in perms:
    for r in roots(*P): allr.add(sp.nsimplify(sp.expand(r.subs({x:-3,m:3}))))
print('distinct roots:',sorted(allr,key=str), len(allr))
