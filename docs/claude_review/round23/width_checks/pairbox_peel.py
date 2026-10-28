# For random pair-box polytopes P(l,u) = {x in M_5 : l_A <= x_a+x_b <= u_A}, compare
#   (i) reg(S) (rank mod p; upper bound certified), (ii) the W bound (1/6) sum_{|A|=2} width_A(S),
#   (iii) the iterated coordinate-peeling bound (best over orders of the first two slicing coordinates; bottom level =
#         exact-mod-p reg of 2-dim slices).
import itertools, random, sys
from reglib import reg, rank_at
pairs=list(itertools.combinations(range(5),2))
def pts_box(l,u,R):
    out=[]
    for x in itertools.product(range(-R,R+1),repeat=4):
        y=x+(-sum(x),)
        if all(l[i]<=y[a]+y[b]<=u[i] for i,(a,b) in enumerate(pairs)): out.append(y)
    return out
def peelvals(v):
    v=sorted(v,reverse=True); return max(x+i for i,x in enumerate(v))
cache={}
def reg2(pts,drop):  # pts in Z^5 lying in a 2-dim affine plane; use the coords not in drop (keep 3, then one dependent)
    keep=[i for i in range(5) if i not in drop][:2]
    key=tuple(sorted((p[keep[0]],p[keep[1]]) for p in pts))
    if key not in cache: cache[key]=reg(list(key))
    return cache[key]
def iter_peel(S,i,j):
    sl={}
    for y in S: sl.setdefault(y[i],[]).append(y)
    vals=[]
    for c,T in sl.items():
        sl2={}
        for y in T: sl2.setdefault(y[j],[]).append(y)
        vals.append(peelvals([reg2(U,(i,j)) for U in sl2.values()]))
    return peelvals(vals)
random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 1)
N=int(sys.argv[1]); worst=0
for t in range(N):
    R=random.randint(2,4)
    # random set of points -> its pair-box hull
    k=random.randint(2,6)
    gen=[]
    for _ in range(k):
        x=[random.randint(-R,R) for _ in range(4)]; gen.append(tuple(x)+(-sum(x),))
    l=[min(g[a]+g[b] for g in gen) for a,b in pairs]; u=[max(g[a]+g[b] for g in gen) for a,b in pairs]
    S=pts_box(l,u,3*R+2)
    W=sum(u[i]-l[i] for i in range(10))/6
    ip=min(iter_peel(S,i,j) for i in range(5) for j in range(5) if i!=j)
    r=reg([y[:4] for y in S]) if len(S)<=700 else None
    worst=max(worst,ip/W)
    print(t,'|S|',len(S),'reg',r,'iterpeel',ip,'W',round(W,3),'ratio_ip/W',round(ip/W,3),flush=True)
print('max iterpeel/W',worst)
