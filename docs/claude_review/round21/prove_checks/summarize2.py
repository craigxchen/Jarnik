import sys, math
from collections import defaultdict
fn=sys.argv[1]
B=defaultdict(dict)  # (bin,k) -> {p: (C,n)}
for line in open(fn):
    t=line.split()
    if t[0]!='best': continue
    b,k,j,p,C,n=int(t[1]),int(t[2]),int(t[3]),int(t[4]),float(t[5]),int(t[6])
    B[(b,k)][p]=(C,n)
bins=sorted(set(b for b,k in B))
print("best C_k per decade of n: all (B-smooth) | P_+<=k^2 | P_+<=k^1.5   [n, P_+ of the overall best]")
for k in range(3,13):
    print("k=%d"%k)
    for b in bins:
        d=B.get((b,k))
        if not d: continue
        allb=min((v[0],v[1],p) for p,v in d.items())
        w2=[v[0] for p,v in d.items() if p<=k*k]
        w15=[v[0] for p,v in d.items() if p<=k**1.5]
        print("  1e%-2d %9.4f | %9s | %9s   n=%d P+=%d"%(b,allb[0], "%.4f"%min(w2) if w2 else '-', "%.4f"%min(w15) if w15 else '-', allb[1], allb[2]))
