# Sections at the even-set vertex x_[5] = 1[|B| even]: reg of P(alpha,beta) = {m in M: m_a+m_b >= -alpha, m_c <= beta}
# vs cost_{x5}(S) = sum_{B even} pi_B(S). reg computed as max r with rank[m^alpha]_{|alpha|<r} < |S| modulo two primes.
import itertools, sys, flint, time
from math import comb
Bs=[frozenset(s) for r in range(1,5) for s in itertools.combinations(range(5),r)]
def lattice(alpha,beta):
    pts=[]
    R=range(-alpha-beta-1, beta+1)
    for m in itertools.product(R,repeat=4):
        m5=-sum(m); mm=list(m)+[m5]
        if max(mm)>beta: continue
        if any(mm[a]+mm[b]< -alpha for a,b in itertools.combinations(range(5),2)): continue
        pts.append(tuple(mm))
    return pts
def pi(S,B): return max(-sum(m[a] for a in B) for m in S)
def monos(r):
    out=[]
    for d in range(r):
        for c in itertools.combinations_with_replacement(range(4),d): out.append(c)
    return out
def rank_mod(S,r,p):
    ms=monos(r); M=flint.nmod_mat(len(ms),len(S),p)
    for j,m in enumerate(S):
        for i,c in enumerate(ms):
            v=1
            for k in c: v=v*m[k]%p
            M[i,j]=v%p
    return M.rank()
def reg(S,p):
    r=1
    while True:
        if comb(r+4,4)<len(S): r+=1; continue   # #conditions for ord>=r+1 is comb(r+4,4); if < |S| a section with ord>=r+1 exists
        if rank_mod(S,r+1,p)<len(S): r+=1
        else: return r
P1,P2=2**31-1,1000000007
for alpha in range(0,7):
    for beta in range(0,7):
        if alpha+beta==0: continue
        S=lattice(alpha,beta)
        if len(S)>2500 or len(S)<2: continue
        cost=sum(pi(S,B) for B in Bs if len(B)%2==0)
        costodd=sum(pi(S,B) for B in Bs if len(B)%2==1)
        r1=reg(S,P1); r2=reg(S,P2)
        flag="  <-- BEATS (4reg>cost)" if 4*min(r1,r2)>cost else ""
        print(f"alpha={alpha} beta={beta} |S|={len(S)} reg={r1},{r2} cost_even={cost} 4reg/cost={4*min(r1,r2)/cost:.4f} cost_odd={costodd}{flag}",flush=True)
