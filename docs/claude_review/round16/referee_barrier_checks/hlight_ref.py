# Exhaustive min of sigma_lambda over integral sum-zero lambda of l1-norm 4, 6 (and 8 for small q), sdc_q.
# sigma = sum_S |lambda(S)| - q over the 2q relative patterns. Translation symmetry: support contains 0 (or only base+..., impossible: sum zero needs >=2 coords; if support avoids F_q entirely impossible).
import sys, itertools, numpy as np, time
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/referee_barrier_checks')
from sdc_sym import check
def types(norm):
    # multisets of nonzero ints, sum 0, abs sum = norm, as sorted tuples
    out=set()
    def parts(k,maxp):
        if k==0: yield (); return
        for p in range(min(k,maxp),0,-1):
            for rest in parts(k-p,p): yield (p,)+rest
    h=norm//2
    for pos in parts(h,h):
        for neg in parts(h,h):
            out.add(tuple(sorted(pos+tuple(-x for x in neg))))
    return sorted(out)
def run(q,norms,tlimit=240):
    _,_,_,M=check(q); n=q+1; r=2*q
    Mc=M.astype(np.int64)  # 2q x n ; column q = base (zeros)
    for norm in norms:
        for vals in types(norm):
            s=len(vals)
            if s==2: continue  # pair multiples
            perms=sorted(set(itertools.permutations(vals)))
            best=None; t0=time.time(); cnt=0
            for rest in itertools.combinations(range(1,n),s-1):
                pass
                break
            combos=np.array([(0,)+c for c in itertools.combinations(range(1,n),s-1)],dtype=np.int64)
            for p in perms:
                pv=np.array(p,dtype=np.int64)
                for st in range(0,len(combos),200000):
                    C=combos[st:st+200000]
                    sums=np.zeros((len(C),2*q),dtype=np.int64)
                    for i in range(s):
                        sums+=pv[i]*Mc[:,C[:,i]].T
                    sig=np.abs(sums).sum(1)-q
                    mi=int(sig.min()); cnt+=len(C)
                    if best is None or mi<best[0]:
                        j=int(sig.argmin()); best=(mi,tuple(C[j]),p)
                if time.time()-t0>tlimit: print("  timeout",flush=True); break
            print(f"q={q} norm={norm} type={vals} count={cnt} min sigma={best[0]} min sigma/r={best[0]/r:.4f} witness coords={best[1]} vals={best[2]}",flush=True)
if __name__=="__main__":
    for q in [7,11,19,23]:
        run(q,[4,6])
    for q in [7,11]:
        run(q,[8])
