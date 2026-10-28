# for each homotopy endpoint at tau=1: cluster coincident places, sum valuation vectors, test 'all multiples of +-1_a'
import sys, ast, numpy as np
def vec(lab):
    v=np.zeros(5,int)
    for ch in lab: v[int(ch)-1]=1
    return v
def reduce(v):  # mod the all-ones vector: normalise min to 0
    return tuple(v-v.min())
for fn in sys.argv[1:]:
    res={}
    for line in open(fn):
        if ' coinc ' not in line or 'reached' not in line: continue
        it=line.split()[0]; i=line.index('coinc ')+6; j=line.index(' real ')
        pairs=ast.literal_eval(line[i:j]); par={}
        def f(x):
            while par.get(x,x)!=x: x=par[x]
            return x
        for a,b in pairs:
            ra,rb=f(a),f(b)
            if ra!=rb: par[ra]=rb
        cl={}
        for a,b in pairs:
            for x in (a,b): cl.setdefault(f(x),set()).add(x)
        vs=[reduce(sum(vec(x) for x in c)) for c in cl.values()]
        # all clusters' vectors supported on a single coordinate after reduction (or its complement)?
        def single(v):
            v=np.array(v); nz=np.nonzero(v)[0]
            if len(nz)==1: return ('+',nz[0]+1)
            z=np.nonzero(v==0)[0]
            if len(z)==1 and len(set(v[v>0]))==1: return ('-',z[0]+1)
            return None
        vs=[v for v in vs if any(v)]  # zero-vector clusters are not places
        tags=[single(v) for v in vs]
        idx={t[1] for t in tags if t}
        res[it]=('collapse to point %s'%idx if all(tags) and len(idx)==1 else 'other', len(cl))
    for k,v in res.items(): print(fn,k,v)
