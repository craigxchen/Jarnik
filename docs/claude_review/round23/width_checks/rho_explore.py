import sys; sys.setrecursionlimit(1000000)
from functools import lru_cache
def peel(vals):
    v=sorted(vals,reverse=True); return max(x+i for i,x in enumerate(v))
@lru_cache(None)
def rho(n,a,b):
    if n==1: return 0
    if n==2: return a+b
    vals=[rho(n-1,a-c,b) for c in range(a+1)]+[rho(n-1,a,b-c) for c in range(1,b+1)]
    return peel(vals)
# symmetry check
print('sym4', all(rho(4,a,b)==rho(4,b,a) for a in range(15) for b in range(15)))
print('sym5', all(rho(5,a,b)==rho(5,b,a) for a in range(12) for b in range(12)))
for a in range(0,13):
    print(a,[rho(4,a,b) for b in range(a,a+25)])
print('rho5 minus (a+b), rows a, cols b=a..a+30')
for a in range(0,13):
    print(a,[rho(5,a,b)-a-b for b in range(a,a+31)])
def g4(a,b): return a+b+min((a+b)//2, 2*min(a,b))
print('rho4 formula ok', all(rho(4,a,b)==g4(a,b) for a in range(40) for b in range(40)))
def g3(a,b): return a+b+min(a,b)
print('rho3 formula ok', all(rho(3,a,b)==g3(a,b) for a in range(60) for b in range(60)))
