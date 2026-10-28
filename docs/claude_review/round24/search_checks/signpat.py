import itertools
SUBS=[frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,5),k)]
pats=[]
for bits in itertools.product([1,-1],repeat=15):
    if all(sum(b for b,T in zip(bits,SUBS) if j in T)==0 for j in range(1,5)): pats.append(bits)
print(len(pats))
# orbits under S4 x global flip
def act(p,perm,fl):
    d={frozenset(perm[x-1] for x in T):b*fl for b,T in zip(p,SUBS)}
    return tuple(d[T] for T in SUBS)
seen=set(); reps=[]
for p in pats:
    if p in seen: continue
    orb={act(p,perm,fl) for perm in itertools.permutations(range(1,5)) for fl in (1,-1)}
    seen|=orb; reps.append((p,len(orb)))
print(len(reps))
for p,n in reps: print(n, ''.join('+' if b>0 else '-' for b in p))
