from itertools import permutations, combinations
from collections import Counter
def run(k):
    full=(1<<k)-1
    Sp=[m for m in range(1<<k) if not m&1 and 2<=bin(m).count('1')<=k-2]
    def canon(m): return full^m if m&1 else m
    def compat(a,b):
        return (a&b)==0 or (a&~b&full)==0 or (b&~a&full)==0 or ((full^a)&(full^b))==0
    best=None; bytype=Counter(); mins={}
    for s in permutations(range(k)):
        if s==tuple(range(k)): continue
        def img(m): 
            r=0
            for i in range(k):
                if m>>i&1: r|=1<<s[i]
            return canon(r)
        c=sum(1 for T in Sp if not compat(T,img(T)))
        # cycle type
        seen=set(); ct=[]
        for i in range(k):
            if i in seen: continue
            j=i; L=0
            while j not in seen: seen.add(j); j=s[j]; L+=1
            ct.append(L)
        ct=tuple(sorted(ct,reverse=True))
        mins[ct]=min(mins.get(ct,10**9),c)
        best=c if best is None else min(best,c)
    return best, sorted(mins.items(), key=lambda t:t[1])[:4]
for k in [5,6,7,8]:
    print(k, run(k), flush=True)
