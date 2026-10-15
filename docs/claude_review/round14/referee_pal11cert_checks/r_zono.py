# Referee: exact zonotope test (some translate k0 + {subset sums} inside A) for each leaf.
import numpy as np, itertools, sys
pre=sys.argv[1]
A=np.load(pre+'_A.npy'); Aset=set(map(tuple,A))
pool=[list(map(int,l.split())) for l in open(pre+'_pool.txt').read().split('\n')[1:] if l.strip()]
nz=0; tot=0
for ln in open(pre+'_leaves.txt'):
    zs=[np.array(pool[i][1:]) for i in map(int,ln.split())]; tot+=1
    sums=set()
    for S in itertools.product((0,1),repeat=len(zs)):
        sums.add(tuple(sum(s*z for s,z in zip(S,zs))))
    sums=[np.array(s) for s in sums]
    if any(all(tuple(k0+s) in Aset for s in sums) for k0 in A): nz+=1
print(pre,'leaves',tot,'zonotopes',nz)
