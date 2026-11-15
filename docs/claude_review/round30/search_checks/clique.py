# max size of a set of sign vectors in {+-1}^D (normalised: +1 = upper root) with pairwise
# pos>=h, neg>=h, pos+neg>=L+1  (balance (S2, proved in search.md section 1) + Lemma 1.2)
import igraph as ig, itertools, sys, time
D=int(sys.argv[1]); c=int(sys.argv[2]); bal=sys.argv[3]=='1'
K=(D+2*c)//2; L=D-1-K; h=(L+1)//2
V=list(range(2**D)); pc=[bin(v).count('1') for v in range(2**D)]
E=[]
for x in V:
    for y in range(x+1,2**D):
        p=pc[x&~y]; n=pc[y&~x]
        if p+n>=L+1 and (not bal or (p>=h and n>=h)): E.append((x,y))
g=ig.Graph(n=2**D,edges=E)
t=time.time(); w=g.clique_number()
cl=g.largest_cliques()[0] if D<=9 else None
print("D",D,"c",c,"L",L,"h",h if bal else 0,"max clique",w,"time %.1f"%(time.time()-t), [format(v,'0%db'%D) for v in cl] if cl else "")
