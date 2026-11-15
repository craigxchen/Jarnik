# combinatorial fibres when only coincidences are alpha-alpha' in Z(1,1,1,1) (plus exact m=2 condition)
import itertools, numpy as np
eps=[1,-1,-1,1]
states=[(1,0),(-1,0),(0,1),(0,-1)]
pts=list(itertools.product(range(4),repeat=4))
def key(p):
    al=[states[s][0] for s in p]; be=[states[s][1] for s in p]
    # alpha class modulo Z*1: normalise by subtracting min? alpha in {-1,0,1}^4; alpha~alpha' iff alpha-alpha' const vector
    a0=al[0]; rel=tuple(x-a0 for x in al)
    return (rel, sum(b*e for b,e in zip(be,eps)))
from collections import defaultdict
F=defaultdict(list)
for p in pts: F[key(p)].append(p)
sizes=sorted((len(v) for v in F.values()),reverse=True); print(sizes[:6])
def svec(p):  # back to s in {+-1}^8: s+=al+be? s_k+ = alpha+beta, s_k- = beta-alpha
    out=[]
    for s in p:
        a,b=states[s]; out+= [a+b, b-a]
    return out
for k,v in F.items():
    if len(v)==8:
        S=np.array([svec(p) for p in v]); print('8-fibre key',k,'rank',np.linalg.matrix_rank(S[1:]-S[0]))
