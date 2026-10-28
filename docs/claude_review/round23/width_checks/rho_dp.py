# Recursive peeling bound rho_n(a,b) >= reg(T_n(a,b)), T_n(a,b) = lattice points of a*Delta - b*Delta in {sum y = a-b} in Z^n.
# Slices y_1 = c: c>=0 -> T_{n-1}(a-c,b); c<0 -> T_{n-1}(a,b-|c|).  rho_2(a,b) = a+b (exact).
from functools import lru_cache
import sys
def peel(vals):
    v=sorted(vals,reverse=True); return max(x+i for i,x in enumerate(v))
@lru_cache(None)
def rho(n,a,b):
    if n==1: return 0
    if n==2: return a+b
    vals=[rho(n-1,a-c,b) for c in range(a+1)]+[rho(n-1,a,b-c) for c in range(1,b+1)]
    return peel(vals)
N=int(sys.argv[1]) if len(sys.argv)>1 else 30
for n in (3,4,5):
    print('n',n,[ (m,rho(n,m,m)) for m in range(1,N+1)])
print('rho4 table a,b<=6:')
for a in range(7): print([rho(4,a,b) for b in range(7)])
print('rho5(m,m) - floor(10m/3):',[rho(5,m,m)-(10*m)//3 for m in range(1,N+1)])
