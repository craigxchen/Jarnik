# kill_search: circle-locus bound: delta needs >= L sign changes in angular order, + balance + distance.
import sys, itertools, numpy as np, networkx as nx, time
D, c = int(sys.argv[1]), int(sys.argv[2])
K = (D + 2*c)//2; L = D-1-K; h = -(-L//2)
V = np.array(list(itertools.product([1,-1], repeat=D)))
n = len(V)
def sc(d):
    s = [x for x in d if x != 0]
    return sum(1 for i in range(len(s)-1) if s[i] != s[i+1])
G = nx.Graph(); G.add_nodes_from(range(n))
for i in range(n):
    for j in range(i+1, n):
        d = (V[i]-V[j])//2
        p = int((d == 1).sum()); m = int((d == -1).sum())
        if min(p, m) >= h and p+m >= L+1 and sc(d) >= L: G.add_edge(i, j)
t0 = time.time()
cl, _ = nx.max_weight_clique(G, weight=None)
print("D", D, "circle-locus clique number", len(cl), "time %.1f" % (time.time()-t0))
