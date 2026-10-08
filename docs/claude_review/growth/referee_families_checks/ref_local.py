# Independent local solubility test for symmetric clusters.
# AXIS type: exist integers n, a with n-(a-d)^2 a square for all d in D (0 in D).
# DIAG type (reduced cluster on odd n, reflection z -> i conj z): with u=x+y, v=x-y, u^2+v^2=2n,
#   u,v odd; u-values a-d with d in D (D subset of 2Z, 0 in D); need 2n-(a-d)^2 = odd square.
import itertools
MODS=[16,32,64,128,9,27,81,5,25,125,7,49,11,13,17,19,23,29,31,37,41]
SQ={M:{x*x%M for x in range(M)} for M in MODS}
def ok_mod(D,M,diag):
    sq=SQ[M]
    for a in range(M):
        if diag and M%2==0 and a%2==0: continue
        for N in range(M):
            if diag and M%2==0 and N%4!=2 % M: 
                if M>=4 and N%4!=2: continue
            if all((N-(a-d)**2)%M in sq for d in D): return True
    return False
def obstruction(D,diag=False):
    for M in MODS:
        if not ok_mod(D,M,diag): return M
    return 0
def table(m,diag,maxspan):
    step=2 if diag else 1
    res={}
    for span in range((m-1)*step, maxspan+1, step):
        sols=[]
        for mid in itertools.combinations(range(step,span,step), m-2):
            D=(0,)+mid+(span,)
            if obstruction(D,diag)==0: sols.append(D)
        res[span]=sols
    return res
import math
print("AXIS type (C >= sqrt(8*dmax) - o(1)):")
for m,maxspan in ((3,4),(4,6),(5,10),(6,11)):
    t=table(m,False,maxspan)
    first=min([s for s in t if t[s]],default=None)
    print("  m=%d (2m=%d pts): smallest locally soluble dmax=%s sets=%s -> C>=%.4f"%(m,2*m,first,t.get(first),math.sqrt(8*first) if first else float('nan')))
print("DIAG type on odd n (C >= sqrt(8*dmax_uv)/2^(1/4) - o(1)):")
for m,maxspan in ((3,8),(4,10),(5,16),(6,18)):
    t=table(m,True,maxspan)
    first=min([s for s in t if t[s]],default=None)
    print("  m=%d (2m=%d pts): smallest locally soluble dmax_uv=%s sets=%s -> C>=%.4f"%(m,2*m,first,(t.get(first) or [])[:6],math.sqrt(8*first)/2**0.25 if first else float('nan')))
