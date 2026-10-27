# forced check with first coordinate reduced modulo unit scaling (collision-freeness is invariant under c -> u c, u in (Z/m)^*).
import sys, json, rlib, time, subprocess, os
from math import gcd
P=json.loads(sys.argv[1]); K=int(sys.argv[2]); QMIN=int(sys.argv[3]); QMAX=int(sys.argv[4])
pis=[rlib.gpi(p) for p in P]; r=len(P)
V=[v for v in rlib.VK(r,K) if next(t for t in v if t)>0]
F=[]; t0=time.time()
for q in range(QMIN,QMAX+1):
    if not rlib.isprime(q) or q in P: continue
    m,L=rlib.labels(pis,q); par=[x%2 for x in L]
    reps=sorted(set(gcd(x,m)%m for x in range(m) if x%2==par[0]))
    lines=["%d %d 0 0"%(r,m)," ".join(map(str,par)),"%d "%len(reps)+" ".join(map(str,reps)),str(len(V))]
    for v in V: lines.append(" ".join(map(str,v))+" 15")
    lines.append("0"); lines.append("")
    out=subprocess.run([os.path.join(os.path.dirname(os.path.abspath(__file__)),'enum2')],input="\n".join(lines)+"\n",capture_output=True,text=True).stdout
    if 'NONE' in out: F.append(q)
    elif 'FOUND' not in out: raise RuntimeError(out)
    print("q",q,"forced" if 'NONE' in out else "free",round(time.time()-t0,1),file=sys.stderr,flush=True)
print(json.dumps({'P':P,'K':K,'qmin':QMIN,'qmax':QMAX,'F':F}))
