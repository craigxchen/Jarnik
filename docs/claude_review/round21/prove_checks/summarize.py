import sys, math
from sympy import factorint
fn=sys.argv[1]
best={}; byj={}
first={}
for line in open(fn):
    t=line.split()
    if t[0]=='best':
        k,j,p,C,n=int(t[1]),int(t[2]),int(t[3]),float(t[4]),int(t[5])
        byj[(k,p)]=(C,n)
        if k not in best or C<best[k][0]: best[k]=(C,n,p)
    elif t[0]=='first':
        C0,k,j,p,n=float(t[1]),int(t[2]),int(t[3]),int(t[4]),int(t[5])
        key=(C0,k)
        first.setdefault(key,[]).append((n,p))
print("k  bestC  n  P+  P+/k^2  factor")
for k in sorted(best):
    C,n,p=best[k]
    print(k,round(C,6),n,p,round(p/k**2,3),factorint(n))
print("\nbest C_k restricted to P_+ <= k^a:")
ps=sorted(set(p for (k,p) in byj))
for a in [1.0,1.5,2.0]:
    row=[]
    for k in sorted(best):
        cands=[byj[(k,p)] for p in ps if (k,p) in byj and p<=k**a]
        row.append((k, round(min(cands)[0],4) if cands else None))
    print(a,row)
print("\nleast n with C_k<=C0 (all P_+ <= B), and least n with P_+<=k^2:")
for (C0,k),L in sorted(first.items()):
    n,p=min(L); Lw=[x for x in L if x[1]<=k*k]
    nw=min(Lw) if Lw else None
    print(C0,k,n,p, 'W/(k log k)=%.3f'%(math.log(n)/(k*math.log(k))), 'window:',nw)
