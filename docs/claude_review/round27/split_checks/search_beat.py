# Search: finite S in M with 4*reg(S) > min_{x in Pi_nt} sum_B x_B pi_B(S)  (=> sections beat 15/4 at that profile x)
import sys, time; sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks')
from reglib import *
from scipy.optimize import linprog
A_eq=[];b_eq=[];seen=set()
for B in Bs:
    C=comp(B)
    if C in seen: continue
    seen.add(B); row=[0]*30; row[idx[B]]=1; row[idx[C]]=1; A_eq.append(row); b_eq.append(1)
for a in range(1,5):
    row=[0]*30
    for B in Bs:
        if 0 in B: row[idx[B]]+=1
        if a in B: row[idx[B]]-=1
    A_eq.append(row); b_eq.append(0)
A_eq=np.array(A_eq,float); b_eq=np.array(b_eq,float)
def mincost(pi):
    r=linprog(pi.astype(float),A_eq=A_eq,b_eq=b_eq,bounds=[(0,1)]*30,method='highs')
    return r.fun, r.x
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
best=0; t0=time.time(); cnt=0
mode=sys.argv[2] if len(sys.argv)>2 else 'rand'
while time.time()-t0<float(sys.argv[3]) if len(sys.argv)>3 else 300:
    if mode=='rand':
        p=list(rng.integers(0,4,size=30))
    else:  # lopsided: small on a random one-sided family, larger on complements
        p=[None]*30
        for B in Bs:
            if p[idx[B]] is not None: continue
            C=comp(B); s=rng.integers(0,2)
            small=int(rng.integers(0,3)); big=int(rng.integers(2,7))
            p[idx[B]],p[idx[C]]=(small,big) if s==0 else (big,small)
    S=lattice_points(p)
    if S is None or len(S)<2 or len(S)>1500: continue
    r=reg(S)
    pi=pis(S)
    mc,x=mincost(pi)
    cnt+=1
    ratio=4*r/mc if mc>1e-9 else float('inf')
    if ratio>best-1e-12:
        best=max(best,ratio)
        print(f"#{cnt} |S|={len(S)} reg={r} mincost={mc:.4f} ratio={ratio:.4f} p={p}",flush=True)
print("done",cnt,"best",best)
