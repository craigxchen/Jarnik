# non-squarefree types: group cut columns of a squarefree type into nested chains (root pairs with multiplicity)
import itertools, pickle, sys
from cuttypes import CUTS
ALL=frozenset(range(5))
perms=list(itertools.permutations(range(5)))
def chain_vector(cuts):
    # cuts: list of frozensets; find orientation making them nested; return exponent vector (a, e) canonical mod flip
    for orient in itertools.product([0,1],repeat=len(cuts)):
        S=[c if o==0 else ALL-c for c,o in zip(cuts,orient)]
        S.sort(key=len,reverse=True)
        if all(S[i+1]<=S[i] for i in range(len(S)-1)):
            a=tuple(sum(1 for s in S if i in s) for i in range(5))
            e=len(cuts)
            b=tuple(e-x for x in a)
            return (min(a,b),e)
    return None
def canon(cols):
    best=None
    for p in perms:
        cc=[]
        for a,e in cols:
            pa=[0]*5
            for i in range(5): pa[p[i]]=a[i]
            pa=tuple(pa); pb=tuple(e-x for x in pa)
            cc.append((min(pa,pb),e))
        t=tuple(sorted(cc))
        if best is None or t<best: best=t
    return best
def partitions_into_chains(cuts):
    n=len(cuts); res=set()
    def rec(i,chains):
        if i==n:
            if any(len(c)>1 for c in chains):
                cols=[chain_vector([cuts[j] for j in c]) for c in chains]
                res.add(canon(cols))
            return
        for c in chains:
            c.append(i)
            if chain_vector([cuts[j] for j in c]) is not None:
                rec(i+1,chains)
            c.pop()
        chains.append([i]); rec(i+1,chains); chains.pop()
    rec(0,[])
    return res
if __name__=="__main__":
    D,h=int(sys.argv[1]),int(sys.argv[2])
    T=pickle.load(open(f'types_{D}_{h}.pkl','rb'))
    allc=set()
    for ty in T:
        cuts=[]
        for k,c in enumerate(ty): cuts+= [CUTS[k]]*c
        allc|=partitions_into_chains(cuts)
    allc=sorted(allc)
    print(D,h,len(allc))
    from collections import Counter
    print(Counter(len(c) for c in allc))
    pickle.dump(allc,open(f'nsf_{D}_{h}.pkl','wb'))
