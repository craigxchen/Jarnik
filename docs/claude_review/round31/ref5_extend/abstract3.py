# Abstract char-3 (and char-2 control) combinatorics of the large-L product design.
# p=3: V=F_3^n, pairs {v,-v} (v!=0). A transversal S picks one element per pair. Two transversals S,S' are "coset-adjacent" if
# A=S\S' is an affine coset J+u (u notin J) of a subspace J (then g_A(X)-g_{-A}(X)=2P_J(u) is constant in X).
# Large-L cluster condition (necessary and, for L large, sufficient within coset-type pairs): |A_ij| > |S\C|/2, C = common part.
# We search max cliques.
import itertools, sys, networkx as nx
p=int(sys.argv[1]); n=int(sys.argv[2])
V=list(itertools.product(range(p),repeat=n))
def add(a,b): return tuple((x+y)%p for x,y in zip(a,b))
def neg(a): return tuple((-x)%p for x in a)
def subspaces():
    seen=set()
    for k in range(0,n+1):
        for gens in itertools.combinations(V,k):
            S={tuple([0]*n)}
            for g in gens:
                S={add(s,tuple(c*x%p for x in g)) for s in S for c in range(p)}
            fs=frozenset(S)
            if fs not in seen: seen.add(fs); yield fs
SUB=list(subspaces())
COS=set()
for J in SUB:
    for u in V:
        if u in J: continue
        COS.add(frozenset(add(u,j) for j in J))
if p==3:
    pairs=[];seen=set()
    for v in V:
        if any(v) and v not in seen: seen|={v,neg(v)}; pairs.append((v,neg(v)))
else:  # p=2 control: pairs {v, v+e0}
    e0=tuple([1]+[0]*(n-1)); pairs=[];seen=set()
    for v in V:
        if v not in seen: seen|={v,add(v,e0)}; pairs.append((v,add(v,e0)))
np_=len(pairs)
trans=[frozenset(pairs[k][(x>>k)&1] for k in range(np_)) for x in range(1<<np_)]
G=nx.Graph(); G.add_nodes_from(range(len(trans)))
for i in range(len(trans)):
    for j in range(i+1,len(trans)):
        A=trans[i]-trans[j]
        if A in COS: G.add_edge(i,j)
print("p",p,"n",n,"pairs",np_,"transversals",len(trans),"coset edges",G.number_of_edges())
best=0;bestc=None
for cl in nx.find_cliques(G):
    # check size condition on all subsets? condition depends on C; check maximal clique and greedily its subsets of size k
    for k in range(len(cl),1,-1):
        ok=None
        for sub in itertools.combinations(cl,k):
            C=frozenset.intersection(*[trans[i] for i in sub]); nprime=len(trans[sub[0]]-C)
            if all(2*len(trans[a]-trans[b])>nprime for a,b in itertools.combinations(sub,2)): ok=sub;break
        if ok:
            if k>best: best=k;bestc=ok
            break
print("max large-L coset cluster size",best)
if bestc:
    C=frozenset.intersection(*[trans[i] for i in bestc]); print("n'=",len(trans[bestc[0]]-C),"pair sizes",sorted(set(len(trans[a]-trans[b]) for a,b in itertools.combinations(bestc,2))))
