import json, sys, time, itertools
from math import log
import flint
R=json.load(open('ref_phiN.json'))
for N in map(int,sys.argv[1:]):
    rec=[r for r in R if r['N']==N][0]   # first point = author's choice (w6 matches)
    t0=time.time(); uns=0
    for C,ds in rec['d'].items():
        d=int(ds)
        if d<=1: uns+=1; continue
        fac=flint.fmpz(d).factor()
        ps=[(int(p),e) for p,e in fac]
        L=log(d); divs=[0.0]
        for p,e in ps: divs=[x+k*log(p) for x in divs for k in range(e+1)]
        ok=any(0.25*L<=x<=0.75*L for x in divs)
        if not ok: uns+=1
    print(N,'w6=%.3f'%rec['w6'],'unsplittable',uns,'time %.1f'%(time.time()-t0),flush=True)
