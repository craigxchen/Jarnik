# kill_search: independent enumeration of M-point balanced types, up to S_D x S_M x global flip; Gordan; saturated balance.
import sys, itertools, numpy as np, networkx as nx
from scipy.optimize import linprog
D, M, c = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
K = (D + 2*c)//2; L = D-1-K; h = -(-L//2)
V = np.array(list(itertools.product([1,-1], repeat=D)))
P = ((V[:,None,:] == 1) & (V[None,:,:] == -1)).sum(2)
adj = (np.minimum(P, P.T) >= h) & (P + P.T >= L+1)
G = nx.Graph(); ii, jj = np.where(np.triu(adj, 1)); G.add_edges_from(zip(ii, jj))
def canon(rows):
    best = None
    for perm in itertools.permutations(range(M)):
        for f in (1, -1):
            Mx = f*rows[list(perm)]
            cols = tuple(sorted(tuple(col) for col in Mx.T))
            if best is None or cols < best: best = cols
    return best
types = {}
for cl in nx.enumerate_all_cliques(G):
    if len(cl) < M: continue
    if len(cl) > M: break
    rows = V[cl]
    types.setdefault(canon(rows), rows)
print("D", D, "M", M, "c", c, "types:", len(types))
def gordan(rows):
    Dl = np.array([(rows[i]-rows[0])//2 for i in range(1, M)])
    # need y>0 (y>=1 wlog) with Dl y = 0
    r = linprog(np.zeros(D), A_eq=Dl, b_eq=np.zeros(M-1), bounds=[(1, None)]*D)
    return r.status == 0
alive = []
for key, rows in types.items():
    Dl = np.array([(rows[i]-rows[0])//2 for i in range(1, M)], float)
    rk = np.linalg.matrix_rank(Dl)
    g = gordan(rows)
    # saturated balance: all a in {-2..2}^D in the real span of Dl
    viol = None
    if g:
        Q, _ = np.linalg.qr(Dl.T); Q = Q[:, :rk]
        Aall = np.array(list(itertools.product(range(-2, 3), repeat=D)), float)
        resid = np.linalg.norm(Aall - (Aall @ Q) @ Q.T, axis=1)
        In = Aall[(resid < 1e-9) & (np.abs(Aall).sum(1) > 0)]
        ap = np.where(In > 0, In, 0).sum(1); am = np.where(In < 0, -In, 0).sum(1)
        bad = (ap + am < L+1) | (np.minimum(ap, am) < L/2)
        if bad.any(): viol = In[np.where(bad)[0][0]].astype(int)
    if g and viol is None: alive.append(rows)
    print(" rank", rk, "gordan", g, "sat-balance violation", None if viol is None else viol.tolist())
print("alive after Gordan + saturated balance (coeffs in [-2,2]):", len(alive))
for rows in alive:
    if len(alive) <= 6: print(rows)
