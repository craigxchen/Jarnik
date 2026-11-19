# Independent test of the log-determinant mechanism of III.6 on the record's k=4 family (gauss.md 4.1).
import sympy as sp
s,X,m=sp.symbols('s X m')
I=sp.I
def conjf(e): return e.subs(I,-I)  # parameters real
xa=X; xb=X*(1-X)/(m**2+X); xc=1-xa-xb; u=m*xb
r={ (1,2,3):I, (1,2):u+I*(xa-1),(1,3):u+I*(xb-1),(2,3):u+I*(xc-1)}
def fourth(roots):
    s1=sum(roots); s2=roots[0]*roots[1]+roots[0]*roots[2]+roots[1]*roots[2]
    re1,im1=sp.re(sp.expand(s1)),sp.im(sp.expand(s1))
    return None
# use explicit Re/Im via real symbols
Xr,mr=sp.symbols('X m',real=True)
sub={X:Xr,m:mr}
for k in list(r): r[k]=sp.simplify(r[k].subs(sub))
def fourth(roots):
    s1=sp.expand(sum(roots)); s2=sp.expand(roots[0]*roots[1]+roots[0]*roots[2]+roots[1]*roots[2])
    re1,im1=sp.re(s1),sp.im(s1); im2=sp.im(s2)
    return sp.simplify((re1*im1-im2)/im1 - I*im1)
for x in (1,2,3):
    others=[r[T] for T in r if x in T and len(T)>=2]
    r[(x,)]=fourth(others)
blocks=list(r.keys())
# u_x = P_x/Pbar_x, P_x = prod_{T∋x}(s-r_T); check Im P_x constant
for x in (1,2,3):
    P=sp.expand(sp.prod([s-r[T] for T in blocks if x in T]))
    Pb=sp.expand(sp.prod([s-sp.conjugate(r[T]) for T in blocks if x in T]))
    D=sp.simplify(sp.expand(P-Pb))
    print('P-Pbar for x=',x,':',sp.factor(D))
def logder(var):
    v=[]
    for x in (1,2,3):
        e=0
        for T in blocks:
            if x in T:
                e+= sp.diff(sp.log(s-r[T]),var)-sp.diff(sp.log(s-sp.conjugate(r[T])),var)
        v.append(e)
    return v
# numeric rational point to keep it light
pt={Xr:sp.Rational(-3),mr:sp.Rational(3)}
cols=[ [sp.together(e) for e in logder(var)] for var in (s,Xr,mr)]
cols=[[sp.simplify(e.subs(pt)) for e in c] for c in cols]
M=sp.Matrix(3,3,lambda i,j: cols[j][i])
d=sp.factor(sp.together(M.det()))
num,den=sp.fraction(d)
print('3x3 det: deg num',sp.degree(num,s),'deg den',sp.degree(den,s),' order at inf',sp.degree(den,s)-sp.degree(num,s))
# projected modulo b_T for each T: 2x2 minors of [b_T | L_s | L_p] = 3x3 det with b_T column
for T in blocks:
    b=sp.Matrix([1 if x in T else 0 for x in (1,2,3)])
    for p in (1,2):
        Mp=sp.Matrix.hstack(b,M[:,0],M[:,p])
        dd=sp.factor(sp.together(Mp.det())); n2,d2=sp.fraction(dd)
        print('proj',T,'param',p,': deg num',sp.degree(n2,s) if n2!=0 else None,'deg den',sp.degree(d2,s),'order',(sp.degree(d2,s)-sp.degree(n2,s)) if n2!=0 else 'inf')
