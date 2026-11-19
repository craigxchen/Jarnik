# For V (F_2-space of F_2[t] containing 1, given by basis) and evaluation Y=T+theta, enumerate ALL sub-transversal
# products z = prod_{w in S}(Y+w), S a transversal of a union W of pairs {w,w+1}; group by (W, top coefficients above
# floor(D_W/2)); same W gives same norm prod_{w in W}(Y+w)(Y+w+1).  Report, for each cluster size m, the least D found.
import sys, itertools, random
from f2lib import *
def span(basis):
    S={0}
    for b in basis: S|={x^b for x in S}
    return sorted(S)
def pairs_of(V):
    P=[];seen=set()
    for v in V:
        if v in seen: continue
        seen|={v,v^1}; P.append(v)  # representative v, partner v^1
    return P
def search(V,T,best,tag):
    P=pairs_of(V); n=len(P)
    groups={}
    def dfs(k,mask,prod):
        if k==n:
            if mask:
                D=deg4(prod); K=D//2
                key=(mask,prod[0]>>(K+1),prod[1]>>(K+1))
                groups.setdefault(key,[]).append(prod)
            return
        dfs(k+1,mask,prod)
        w=P[k]
        dfs(k+1,mask|1<<k,m4(prod,(T^w,1)))
        dfs(k+1,mask|1<<k,m4(prod,(T^w^1,1)))
    dfs(0,0,(1,0))
    for key,pts in groups.items():
        m=len(pts)
        if m<3: continue
        D=deg4(pts[0])
        if m not in best or D<best[m][0]: best[m]=(D,tag,bin(key[0]))
    return best
if __name__=="__main__":
    best={}
    mode=sys.argv[1]
    if mode=="full4":
        V=span([1,2,4,8])
        for T in list(range(0,1<<9,16))+[1<<L for L in range(4,21)]:
            search(V,T,best,("V=F2[t]<4","T=%s"%bin(T)))
    elif mode=="rand4":
        random.seed(int(sys.argv[2]))
        for it in range(int(sys.argv[3])):
            e=random.randint(4,8)
            basis=[1]+[random.randrange(2,1<<e) for _ in range(3)]
            V=span(basis)
            if len(V)<16: continue
            T=random.choice([random.randrange(0,1<<random.randint(1,10)), 1<<random.randint(e,16)])
            search(V,T,best,("basis="+str([bin(b) for b in basis]),"T=%s"%bin(T)))
    for m in sorted(best): print(m,best[m])
