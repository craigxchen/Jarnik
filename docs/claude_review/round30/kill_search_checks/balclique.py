# kill_search: exact max size of balanced codes (pairwise p,n >= h, distance >= L+1); networkx, neighbourhood split by S_D.
import sys, itertools, numpy as np, networkx as nx, time
D, c = int(sys.argv[1]), int(sys.argv[2])
K = (D + 2*c)//2; L = D-1-K; h = -(-L//2)
V = np.array(list(itertools.product([1,-1], repeat=D)))
n = len(V)
P = ((V[:,None,:] == 1) & (V[None,:,:] == -1)).sum(2)
Nn = P.T
adj = (np.minimum(P, Nn) >= h) & (P + Nn >= L+1)
best = 0
t0 = time.time()
for w in range(D+1):
    v = np.array([1]*w + [-1]*(D-w)); vi = int(np.where((V == v).all(1))[0][0])
    nb = np.where(adj[vi])[0]
    G = nx.Graph(); G.add_nodes_from(nb)
    sub = adj[np.ix_(nb, nb)]
    ii, jj = np.where(np.triu(sub, 1))
    G.add_edges_from(zip(nb[ii], nb[jj]))
    cl, _ = nx.max_weight_clique(G, weight=None)
    best = max(best, 1 + len(cl))
print("D", D, "c", c, "L", L, "h", h, "balanced clique number", best, "time %.1f" % (time.time()-t0))
