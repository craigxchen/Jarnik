"""(1) uniqueness: pattern with agreement product h > C^2 has at most one admissible magnitude.
(2) random sign-free tables on maximum cliques: largest line-embeddable subset; p-sides embed.
"""
import random, sys
from itertools import product, combinations
sys.path.insert(0,'.')
from stmt1 import run

def embed(S, d):
    """S list of vertices, d(i,j) magnitude; return coords or None (try all signs relative to S[0])."""
    if len(S)<=1: return [0]*len(S)
    r=S[0]
    # backtracking
    y={r:0}
    order=S[1:]
    def rec(t):
        if t==len(order): return True
        v=order[t]
        for sgn in (1,-1):
            yv=sgn*d(r,v)
            if all(abs(yv-y[w])==d(v,w) for w in y if w!=r):
                y[v]=yv
                if rec(t+1): return True
                del y[v]
        return False
    return [y[v] for v in S] if rec(0) else None

def max_embeddable(S,d):
    n=len(S)
    for r in range(n,0,-1):
        for T in combinations(S,r):
            if embed(list(T),d) is not None: return r,T
    return 0,()

if __name__=='__main__':
  random.seed(7)
  sets=[[5,7,11],[7,11,13],[11,13,17],[13,17,19],[5,7,11,13],[7,11,13,17],[17,19,23,29],[11,13,17,19,23]]
  tot_unique_viol=0
  for primes in sets:
      for C2 in (1,2,4):
          if min(primes)<=C2: continue
          q=1
          for p in primes:q*=p
          k=len(primes)
          alld=list(product(*[range(1,(p-1)//2+1) for p in primes]))
          alld=random.sample(alld,min(len(alld),60))
          hist={}
          for dl in alld:
              bf,bi,v,cl,al=run(primes,C2,dl,None)
              for pat,lst in al.items():
                  h=1
                  for t,p in enumerate(primes):
                      if not (pat>>t)&1: h*=p
                  if h>C2 and len(lst)>1: tot_unique_viol+=1; print("non-unique",primes,dl,pat,lst)
              if bf<5: continue
              for trial in range(5):
                  ch={}
                  def d(i,j):
                      key=(min(i,j),max(i,j))
                      if key not in ch: ch[key]=random.choice(al[i^j])
                      return ch[key]
                  r,T=max_embeddable(cl,d)
                  ncomp=sum(1 for i,j in combinations(cl,2) if (i^j)==(1<<k)-1)
                  hist[(bf,bf-r,ncomp)]=hist.get((bf,bf-r,ncomp),0)+1
                  # p-sides
                  for t,p in enumerate(primes):
                      if p<=C2: continue
                      for b in (0,1):
                          side=[x for x in cl if ((x>>t)&1)==b]
                          if len(side)>=5 and embed(side,d) is None:
                              print("SIDE FAILS",primes,C2,dl,side)
          print(primes,"C2",C2,"(B_free, deletions needed, #complementary pairs):count",sorted(hist.items()),flush=True)
  print("nonunique violations",tot_unique_viol)
