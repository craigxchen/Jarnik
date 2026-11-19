import itertools
from f27 import *
labels=list(itertools.product(range(27),repeat=2))
forms=[(1,0),(0,1),(1,1),(1,2)]
def lin(f,x): return fadd(fmul(from_f3(f[0]),x[0]),fmul(from_f3(f[1]),x[1]))
CH={ (fi,x): chi(lin(f,x)) for fi,f in enumerate(forms) for x in labels}
for m in (3,4):
  for fs in itertools.combinations(range(4),m):
    res=[]
    for signs in itertools.product((1,-1),repeat=m):
        g=0
        for x in labels:
            g+=min(1+s*CH[(fi,x)] for fi,s in zip(fs,signs))
        res.append((g,signs))
    print(m,fs,max(res))
