import sys, json, rlib, time
P=json.loads(sys.argv[1]); K=int(sys.argv[2]); QMAX=int(sys.argv[3])
pis=[rlib.gpi(p) for p in P]; r=len(P)
V=[v for v in rlib.VK(r,K) if next(t for t in v if t)>0]
F=[]; t0=time.time()
for q in range(3,QMAX+1):
    if not rlib.isprime(q) or q in P: continue
    m,L=rlib.labels(pis,q)
    par=[L[j]%2 for j in range(r)]
    out=rlib.run_enum(r,m,0,0,par,[(v,15) for v in V],[])
    if "NONE" in out: F.append(q)
    elif 'FOUND' not in out: raise RuntimeError(out)
    print("q",q,"forced" if "NONE" in out else "free",round(time.time()-t0,1),file=sys.stderr,flush=True)
print(json.dumps({'P':P,'K':K,'qmax':QMAX,'F':F,'n':len(F),'max':max(F)}))
print('time',time.time()-t0,file=sys.stderr)
