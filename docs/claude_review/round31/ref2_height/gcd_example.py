# Bounded-gcd squareclass example u_j = t + c_j + i, c=(0,1,3,7); five of six 2-subsets A.
# Checks: equal norms, distinct points, chord*N^{1/4}/R bounded (critical arc), exact relation,
# Gaussian gcd bounded over t, and the resulting primitive tuple has H_S = N_prim.
from itertools import combinations
from math import gcd, isqrt
import sympy as sp
c=(0,1,3,7)
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def ggcd(x,y):
    # Euclid in Z[i]
    while y!=(0,0):
        n=y[0]**2+y[1]**2
        a=x[0]*y[0]+x[1]*y[1]; b=x[1]*y[0]-x[0]*y[1]   # x*conj(y)
        q=((2*a+n)//(2*n),(2*b+n)//(2*n))
        r=(x[0]-(q[0]*y[0]-q[1]*y[1]), x[1]-(q[0]*y[1]+q[1]*y[0]))
        x,y=y,r
    return x
sets=list(combinations(range(4),2))
for drop in range(6):
  use=[A for k,A in enumerate(sets) if k!=drop]
  maxg=0; maxratio=0
  for t in list(range(1,3000))+[10**6+17,10**9+7,10**15+3]:
    u=[(t+cj,1) for cj in c]; ub=[(t+cj,-1) for cj in c]
    z=[]
    for A in use:
        p=(1,0)
        for j in range(4): p=gmul(p,u[j] if j in A else ub[j])
        z.append(p)
    N=[a*a+b*b for a,b in z]; assert len(set(N))==1
    assert len(set(z))==5
    G=z[0]
    for w in z[1:]: G=ggcd(G,w)
    g2=G[0]**2+G[1]**2; maxg=max(maxg,g2)
    # chord/R * N^{1/4}: max |z_a-z_b|^2 / N * sqrt(N) ... use exact floats via big ints
    Nn=N[0]
    dmax=max((a[0]-b[0])**2+(a[1]-b[1])**2 for a in z for b in z)
    ratio=(dmax/ Nn) * Nn**0.5   # (diam/R)^2 * N^{1/2} = (Delta N^{1/4})^2
    if t>50: maxratio=max(maxratio,ratio)
  # relation: complementary pairs present
  comp=[(A,tuple(sorted(set(range(4))-set(A)))) for A in use if tuple(sorted(set(range(4))-set(A))) in use and A<tuple(sorted(set(range(4))-set(A)))]
  print("drop",sets[drop],"max |gcd|^2 =",maxg," max (Delta N^{1/4})^2 for t>50:",round(maxratio,3)," complementary pairs:",comp)
# symbolic gcd and resultant bound
t=sp.symbols('t')
u=[t+cj+sp.I for cj in c]; ub=[t+cj-sp.I for cj in c]
polys=[sp.expand(sp.prod([u[j] if j in A else ub[j] for j in range(4)])) for A in sets[:5]]
gg=polys[0]
for p in polys[1:]: gg=sp.gcd(gg,p,extension=sp.I)
print("polynomial gcd over Q(i):",gg)
# after division: odd primitive norm, gcd 1
use=sets[:5]; bad=0; gvals=set()
for tt in range(1,3000):
    u=[(tt+cj,1) for cj in c]; ub=[(tt+cj,-1) for cj in c]
    z=[]
    for A in use:
        p=(1,0)
        for j in range(4): p=gmul(p,u[j] if j in A else ub[j])
        z.append(p)
    G=z[0]
    for w in z[1:]: G=ggcd(G,w)
    g2=G[0]**2+G[1]**2; gvals.add(g2)
    Np=(z[0][0]**2+z[0][1]**2)//g2
    if Np%2==0: bad+=1
print("gcd norms seen:",sorted(gvals),"; even primitive norms:",bad)
