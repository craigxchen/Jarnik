# Peeling bound for S_m = mR_5 ∩ M sliced by a functional u (exact reg of each slice via rank mod p, an upper bound
# certified over Q). Peeling lemma: reg(S) <= max_i (reg(S_(i)) + i - 1), slices ordered by decreasing reg.
import itertools, sys
from reglib import reg
def Sm(m):
    pts=[]
    for x in itertools.product(range(-m,m+1),repeat=4):
        y=list(x)+[-sum(x)]
        if sum(v for v in y if v>0)<=m: pts.append(tuple(y))
    return pts
def coords(y):  # lattice coords on M: drop last coordinate
    return y[:4]
u=tuple(int(t) for t in sys.argv[2].split(','))
for m in range(1,int(sys.argv[1])+1):
    S=Sm(m); sl={}
    for y in S:
        c=sum(a*b for a,b in zip(u,y)); sl.setdefault(c,[]).append(coords(y))
    regs=sorted(((reg(v),c) for c,v in sl.items()),reverse=True)
    pb=max(r+i for i,(r,c) in enumerate(regs))
    print(m, 'slices', [(c,r) for r,c in sorted(regs,key=lambda t:t[1])], 'peel', pb, '10m/3=',10*m/3, flush=True)
