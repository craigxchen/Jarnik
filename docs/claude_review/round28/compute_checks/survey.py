import sys; sys.path.insert(0,'.')
import slib, models
from run_box import parse
specs=[]
for L in (1,2,3):
    specs += ['rsyl:4:0,1,2,3,4,5,6,7,8:%d'%L, 'rsyl:4:0,1,2,3,4,5,6,7,8,9:%d'%L, 'rsyl:4:0,1,2,4,8,15,3,5,6:%d'%L,
              'rsdc:11:0,1,2,3,4,5,6,7,8:%d'%L, 'rsdc:11:0,1,2,3,4,5,6,7,8,9:%d'%L,
              'rhad:11:0,1,2,3,4,5,6,7,8:%d'%L, 'rhad:11:0,1,2,3,4,5,6,7,8,9:%d'%L, 'sdc:11:%d'%L, 'had:11:%d'%L]
for s in specs:
    n,pats,L=parse(s)
    if n-1>9: 
        print(s,'n',n,'r',len(pats),'(d=%d, box not built)'%(n-1)); continue
    A=slib.build_box(n,[(S,L) for S in pats])
    print(s,'n',n,'r',len(pats),'|A|',len(A),'fullrank',models.fullrank(n,pats))
for n in (9,10):
    sz=[]
    for seed in range(1,41):
        r,L,pats=models.random_profile(n,seed)
        A=slib.build_box(n,[(S,L) for S in pats]); sz.append((seed,r,L,len(A)))
    print(n, sz)
