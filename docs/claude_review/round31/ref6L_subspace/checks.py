# Independent checks for referee report ref6L_subspace (VI.1)
import itertools, random
from fractions import Fraction
import sympy as sp
from mpmath import mp, mpf, sqrt, log, fabs
mp.dps=60

# (1) quotient rule: deg_w D(f/g) <= deg_w f - deg_w g for monomial valuations, random Laurent polys
X,Y,t=sp.symbols('X Y t')
def degw(expr,w):
    expr=sp.expand(expr)
    if expr==0: return -sp.oo
    P=sp.Poly(sp.expand(expr*X**20*Y**20),X,Y)
    return max(w[0]*(a-20)+w[1]*(b-20) for (a,b) in P.monoms())
def rnd():
    return sum(random.choice([-3,-1,1,2,t,t+1])*X**random.randint(-3,3)*Y**random.randint(-3,3) for _ in range(4))
random.seed(1); bad=0
for trial in range(60):
    w=(sp.Rational(random.randint(-5,5),3),sp.Rational(random.randint(-5,5),2))
    f,g=rnd(),rnd()
    if sp.expand(f)==0 or sp.expand(g)==0: continue
    for D in [lambda h:X*sp.diff(h,X), lambda h:Y*sp.diff(h,Y), lambda h:sp.diff(h,t)]:
        num=sp.expand(g*D(f)-f*D(g))   # D(f/g)=num/g^2
        if num==0: continue
        lhs=degw(num,w)-2*degw(g,w); rhs=degw(f,w)-degw(g,w)
        if lhs>rhs: bad+=1
print("quotient-rule violations:",bad)

# (2) varying-orders example v_n
for n in [10,100,1000]:
    vn=(n+sqrt(n*n-4))/2; vm=(n+1+sqrt((n+1)**2-4))/2
    terms=[vm,1/vm,-vn,-1/vn,mpf(-1)]
    print(n,"relation residue",sp.N(sum(terms),5) if False else fabs(sum(terms)))
    # proper zero subsums (numerically, all conjugate sign choices of the two sqrt's)
    zs=[]
    for s1 in [1,-1]:
      for s2 in [1,-1]:
        a=(n+1+s1*sqrt((n+1)**2-4))/2; b=(n+s2*sqrt(n*n-4))/2
        T=[a,1/a,-b,-1/b,mpf(-1)]
        for r in range(1,5):
          for c in itertools.combinations(range(5),r):
            if fabs(sum(T[i] for i in c))<mpf(10)**-40: zs.append((s1,s2,c))
    print("  proper zero subsums:",zs)
    # height of ratio v_{n+1}/v_n over the degree-4 compositum: average over 4 embeddings of log max(1,|.|)
    hs=[]
    for s1 in [1,-1]:
      for s2 in [1,-1]:
        a=(n+1+s1*sqrt((n+1)**2-4))/2; b=(n+s2*sqrt(n*n-4))/2
        hs.append(max(0,log(fabs(a/b))))
    print("  h(v_{n+1}/v_n)=",sum(hs)/4," log n=",log(n), " cancellation |v_{n+1}-v_n|/v_n=",fabs(vm-vn)/vn)

# (3) one-chord parity: in Im(A*prod(x_j+i y_j)) the all-x and all-y terms (6 factors) carry the same A-component
x=sp.symbols('x1:7'); y=sp.symbols('y1:7'); a,b=sp.symbols('a b',real=True)
for signs in itertools.product([1,-1],repeat=3):
    sg=list(signs)+[1,1,-1]
    prod=sp.expand((a+sp.I*b)*sp.prod([x[j]+sg[j]*sp.I*y[j] for j in range(6)]))
    im=sp.expand(sp.im(prod).subs({**{v:sp.Symbol(str(v),real=True) for v in x+y}}))
    # pick big set = some coordinates; check for each choice of 'large' coordinate per block
    for pick in [(0,)*6,(1,0,0,0,0,0),(1,1,0,1,0,0)]:
        big=[ (y[j] if pick[j] else x[j]) for j in range(6)]
        small=[ (x[j] if pick[j] else y[j]) for j in range(6)]
        P=sp.Poly(prod.as_real_imag()[1] if False else sp.expand(sp.im(sp.expand((a+sp.I*b)*sp.prod([sp.Symbol('X%d'%j,real=True)+sg[j]*sp.I*sp.Symbol('Y%d'%j,real=True) for j in range(6)])))))
        break
    break
# simpler direct check
Xs=[sp.Symbol('X%d'%j,real=True) for j in range(6)]; Ys=[sp.Symbol('Y%d'%j,real=True) for j in range(6)]
ok=True
for sg in itertools.product([1,-1],repeat=6):
  e=sp.expand(sp.im(sp.expand((a+sp.I*b)*sp.prod([Xs[j]+sg[j]*sp.I*Ys[j] for j in range(6)]))))
  for pick in itertools.product([0,1],repeat=6):
    big=sp.prod([Ys[j] if pick[j] else Xs[j] for j in range(6)])
    small=sp.prod([Xs[j] if pick[j] else Ys[j] for j in range(6)])
    cb=e.coeff(big) if True else 0
    cb=sp.Poly(e,*Xs,*Ys).coeff_monomial(big); cs=sp.Poly(e,*Xs,*Ys).coeff_monomial(small)
    # same A-component: both multiples of a, or both of b
    if not ((cb.has(a) and cs.has(a)) or (cb.has(b) and cs.has(b))) or (cb.has(a) and cb.has(b)): ok=False
print("one-chord parity claim holds for all sign patterns/picks:",ok)
