import itertools
Ts=[T for r in range(1,5) for T in itertools.combinations(range(1,5),r)]
pats=[]
for bits in range(1<<15):
    eps=[1 if bits>>i&1 else -1 for i in range(15)]
    if all(sum(eps[i] for i,T in enumerate(Ts) if j in T)==0 for j in range(1,5)): pats.append(tuple(eps))
print('patterns',len(pats))
idx={T:i for i,T in enumerate(Ts)}
seen=set(); classes=0
for p in pats:
    if p in seen: continue
    classes+=1
    for perm in itertools.permutations(range(1,5)):
        for s in (1,-1):
            q=[0]*15
            for i,T in enumerate(Ts):
                q[idx[tuple(sorted(perm[t-1] for t in T))]]=s*p[i]
            seen.add(tuple(q))
print('classes under S4 x flip',classes)
