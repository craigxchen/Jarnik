# Exact indefinite rank-one models: genuine Z[sqrt2] stars, Q* = X^2-2Y^2 = I - 3*Y^2.
import itertools, math
from sympy import primerange
from fractions import Fraction as Fr
def elt(p):  # a+b*sqrt2 of norm p (p = +-1 mod 8), search
    for b in range(1,10**5):
        a2=p+2*b*b; a=math.isqrt(a2)
        if a*a==a2: return (a,b)
ps=[p for p in primerange(7,400) if p%8 in (1,7)][:10]
E=list(itertools.combinations(range(5),2)); pi={e:elt(p) for e,p in zip(E,ps)}
def mul(u,v): return (u[0]*v[0]+2*u[1]*v[1], u[0]*v[1]+u[1]*v[0])
V=[]
for i in range(5):
    z=(1,0)
    for e in E:
        if i in e: z=mul(z,pi[e])
    V.append(z)
n=[math.prod(p for e,p in zip(E,ps) if i in e) for i in range(5)]
print('primes',ps)
for i in range(5): print('V',V[i],' X^2-2Y^2==n_i:',V[i][0]**2-2*V[i][1]**2==n[i],' defect |V|^2-n =',V[i][0]**2+V[i][1]**2-n[i],' =3Y^2:',V[i][0]**2+V[i][1]**2-n[i]==3*V[i][1]**2)
ok=all((V[i][0]*V[j][1]-V[i][1]*V[j][0])%p==0 for (i,j),p in zip(E,ps))
print('all ten divisibilities p_ij | det(V_i,V_j):',ok, ' primitive:',all(math.gcd(*v)==1 for v in V))
