# kill_search: random nondegenerate points on the cube locus; full cluster structure of N (all 256 points).
import numpy as np, itertools, networkx as nx
rng = np.random.default_rng(1)
eps = np.array([1,-1,-1,1]); D=8; L=3
W = np.array(list(itertools.product([-1,0,1], repeat=D)), float)
subs = list(itertools.product([0,1], repeat=D))
stats = {}
done = 0
while done < 60:
    a, b = rng.normal(size=2)*1.5
    disc = (a*b*(a+b))**2 - 4*a*b
    if disc < 0: continue
    c = (-a*b*(a+b) + rng.choice([-1,1])*np.sqrt(disc))/(2*a*b)
    y = np.array([a, b, c, -(a+b+c)])
    assert abs(np.prod(y)-1) < 1e-9 and abs(y.sum()) < 1e-12
    r = eps/y + 1j*y
    q = np.concatenate([r, -r])
    allr = np.concatenate([q, q.conj()])
    sep = min(abs(allr[i]-allr[j]) for i in range(16) for j in range(i))
    if sep < 1e-2: continue
    A = np.array([np.imag(q**m) for m in range(1, L+1)])
    rr = np.linalg.norm(W @ A.T, axis=1)
    ker = W[rr < 1e-9].astype(int)
    zsQ = sum(1 for s in subs if np.linalg.norm(np.array(s) @ A.T) < 1e-9)
    G = nx.Graph(); G.add_nodes_from(range(256))
    for i, s in enumerate(subs):
        for k in ker:
            g = np.array(s) - k
            if np.all((g == 0) | (g == 1)):
                j = int(''.join(map(str, g)), 2)
                if j > i: G.add_edge(i, j)
    w, _ = nx.max_weight_clique(G, weight=None)
    key = (zsQ, len(w), len(ker), round(np.min(rr[rr >= 1e-9]), 3) > 0)
    stats[key] = stats.get(key, 0) + 1
    done += 1
print("(cluster through Q, max cluster on N, #kernel vectors, gap>0): count", stats)
