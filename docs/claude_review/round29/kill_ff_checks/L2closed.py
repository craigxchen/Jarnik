# exhaustive closed-set enumeration: M_L(D) for L=1 (cond: e_k notin span) and L=2 (also e_k +- e_l notin span)
import itertools, numpy as np, sys, time
D=int(sys.argv[1]); L=int(sys.argv[2])
allv=np.array(list(itertools.product([0,1],repeat=D)),dtype=float)  # index = binary number
bad=[np.eye(D)[k] for k in range(D)]
if L==2:
    for k,l in itertools.combinations(range(D),2):
        bad.append(np.eye(D)[k]-np.eye(D)[l]); bad.append(np.eye(D)[k]+np.eye(D)[l])
bad=np.array(bad)
def closure(idx):
    B=allv[list(idx)]
    u,s,vt=np.linalg.svd(B)
    rk=int((s>1e-9).sum())
    Nn=vt[rk:]          # complement basis rows
    if Nn.shape[0]==0: return None
    inside=np.abs(allv@Nn.T).max(axis=1)<1e-9
    if (np.abs(bad@Nn.T).max(axis=1)<1e-9).any(): return None
    return frozenset(np.nonzero(inside)[0].tolist())
start=frozenset([0])
seen={start}; frontier=[start]; best=1; bestset=start
t0=time.time()
while frontier:
    nf=[]
    for S in frontier:
        for v in range(1,2**D):
            if v in S: continue
            T=closure(S|{v})
            if T is None or T in seen: continue
            seen.add(T); nf.append(T)
            if len(T)>best: best=len(T); bestset=T
    frontier=nf
    print("level done, closed sets so far",len(seen),"best",best,round(time.time()-t0,1),flush=True)
print("D",D,"L",L,"M =",best, [format(i,'0%db'%D) for i in sorted(bestset)])
