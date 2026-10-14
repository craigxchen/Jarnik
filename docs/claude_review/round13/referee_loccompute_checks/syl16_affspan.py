# Each width-7 affine-2-flat indicator of Syl(16,2) lies in the span of some plane (7-column) family.
import numpy as np
from syl16_base import COLS, pts
A=np.load("syl16_A.npy").astype(np.int64)
W7=np.load("syl16_w7dirs.npy"); aff=[v for v in W7 if np.abs(v).sum()==4]
planes=[[i for i,p in enumerate(pts) if sum(x*y for x,y in zip(a,p))%2==0] for a in pts]
cnt=0
for u in aff:
    hits=[k for k,V in enumerate(planes) if np.linalg.matrix_rank(np.vstack([COLS[V],u]).astype(float))==7]
    cnt+= len(hits)>0
print("affine 4-sums in the span of some plane family:",cnt,"of",len(aff))
u=aff[0]; hits=[k for k,V in enumerate(planes) if np.linalg.matrix_rank(np.vstack([COLS[V],u]).astype(float))==7]
print("example",u.tolist(),"in spans of planes",hits,"width on A",len(np.unique(A@u)))
# is it in the span of a line (3 columns)? or of any <=6 columns?
import itertools
lines=[]
for i,j in itertools.combinations(range(15),2):
    c=tuple((a+b)%2 for a,b in zip(pts[i],pts[j])); k=pts.index(c); l=tuple(sorted((i,j,k)))
    if l not in lines: lines.append(l)
print("in some line span:",any(np.linalg.matrix_rank(np.vstack([COLS[list(l)],u]).astype(float))==3 for l in lines))
mins=min(len(S) for b in range(1,8) for S in itertools.combinations(range(15),b) if b<=7 and np.linalg.matrix_rank(np.vstack([COLS[list(S)],u]).astype(float))==np.linalg.matrix_rank(COLS[list(S)].astype(float))) if False else None
for b in range(1,8):
    found=any(np.linalg.matrix_rank(np.vstack([COLS[list(S)],u]).astype(float))==np.linalg.matrix_rank(COLS[list(S)].astype(float)) for S in itertools.combinations(range(15),b))
    if found: print("smallest column family whose span contains it: b =",b); break
