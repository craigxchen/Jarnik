# Checks for the non-existence theorem at x_[5] (exact rational arithmetic on random instances).
import sympy as sp, random, itertools
random.seed(3)
y,w,s=sp.symbols('y w s')
def R(): return sp.Rational(random.randint(-30,30),random.randint(1,9))
for trial in range(5):
    rho=[R() for _ in range(5)]; kap=[R() or 1 for _ in range(5)]
    L=[w-rho[c]*y-kap[c] for c in range(5)]
    # nodes
    nodes=[]
    for a,b in itertools.combinations(range(5),2):
        sol=sp.solve([L[a],L[b]],[y,w],dict=True)[0]; nodes.append((sol[y],sol[w]))
    # general quartic: dim of quartics through 10 nodes
    mons=[y**i*w**j for i in range(5) for j in range(5-i)]
    A=sp.Matrix([[m.subs({y:p[0],w:p[1]}) for m in mons] for p in nodes])
    dimV=len(mons)-A.rank()
    # quartics through nodes with order>=3 at origin: kill monomials of degree<=2
    keep=[m for m in mons if sp.Poly(m,y,w).total_degree()>=3]
    A3=sp.Matrix([[m.subs({y:p[0],w:p[1]}) for m in keep] for p in nodes])
    dim3=len(keep)-A3.rank()
    # Veronese rank of v_c=(rho_c,-1)/kappa_c
    V=sp.Matrix([[1,rho[c]/kap[c],-1/kap[c],(rho[c]/kap[c])**2,-rho[c]/kap[c]**2,1/kap[c]**2] for c in range(5)])
    print(f"trial {trial}: dim quartics thru 10 nodes = {dimV} (expect 5); dim with triple pt at O = {dim3} (expect 0); Veronese rank = {V.rank()}")
# converse sanity: 4 concurrent lines -> a triple-point quartic exists (degenerate)
rho=[R() for _ in range(5)]; P=(R(),R())
kap=[P[1]-rho[c]*P[0] for c in range(4)]+[R()]
nodes=[]
L=[w-rho[c]*y-kap[c] for c in range(5)]
for a,b in itertools.combinations(range(5),2):
    sol=sp.solve([L[a],L[b]],[y,w],dict=True)[0]; nodes.append((sol[y],sol[w]))
mons=[y**i*w**j for i in range(5) for j in range(5-i)]
keep=[m for m in mons if sp.Poly(m,y,w).total_degree()>=3]
A3=sp.Matrix([[m.subs({y:p[0],w:p[1]}) for m in keep] for p in nodes])
print("4 concurrent lines: dim triple-point quartics through the (now 5 distinct) nodes =",len(keep)-A3.rank())
