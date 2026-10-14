# Thm 3.4(iii) check: J(P_C) = grid + sum over circuits of pullbacks, and whether each circuit is needed.
import itertools, sys
from rlib import *
from echelon import Ech
gens=eval(sys.argv[1]); j=int(sys.argv[2])
b=len(gens[0]); cap=[2*j+1]*b
def span(gens,b):
    C={tuple([0]*b)}
    for g in gens: C|={tuple((a+x)%2 for a,x in zip(c,g)) for c in C}
    return sorted(C)
C=span(gens,b)
circ=[s for s in C if any(s) and not any(t!=s and any(t) and all(ti<=si for ti,si in zip(t,s)) for t in C)]
print("circuits",circ)
G=list(itertools.product(range(-j,j+1),repeat=b))
P=[y for y in G if all(sum(yy for yy,ss in zip(y,s) if ss)%2==0 for s in C)]
Dmax=2*j*b
hs=[h_func_exact(P,e,cap) for e in range(Dmax+1)]
H=[hs[0]]+[hs[i]-hs[i-1] for i in range(1,len(hs))]
target=[len(monos_eq(b,D,cap))-H[D] for D in range(Dmax+1)]
# local generators of each circuit: J(P_s)_e for e up to Dmax, in the circuit's variables
loc={}
for s in circ:
    idx=[q for q in range(b) if s[q]]; c=len(idx)
    Gs=list(itertools.product(range(-j,j+1),repeat=c)); Ps=[y for y in Gs if sum(y)%2==0]
    gl={}
    for e in range(1,min(Dmax,2*j*c)+1):
        sp=J_space(Ps,e,[2*j+1]*c)
        if sp:
            gl[e]=[{tuple(ex[idx.index(q)] if q in idx else 0 for q in range(b)):v for ex,v in g.items()} for g in sp]
    loc[s]=gl
def dims(use):
    out=[]
    for D in range(Dmax+1):
        E=Ech()
        for s in use:
            for r in ideal_degree_rows(loc[s],b,D,cap): E.add(r)
        out.append(E.rank())
    return out
full=dims(circ)
print("target ",target); print("all circ",full,"equal:",full==target)
for s in circ:
    d=dims([t for t in circ if t!=s])
    print("without",s,d,"deficit",[a-b_ for a,b_ in zip(target,d)])
