import time, pickle, itertools
from m8lib import *
A, cols, H, JA, W = pickle.load(open('m8_step1.pkl','rb'))
cache = {}
gw = {}
for v, tau, h in W:
    gw.setdefault(tau+1, []).append(ppow(lin(v), tau+1, d))
print({e: len(g) for e, g in gw.items()})
res = []
for i in range(9):
    t = time.time()
    rows = ideal_degree_rows(gw, i)
    de = dim_exact(rows, i) if rows else 0
    res.append(de)
    print(f"deg {i}: rows {len(rows)} dim_Q (J_w)_i = {de}  dim J(A)_i = {JA[i]}  [{time.time()-t:.1f}s]", flush=True)
print("J_w dims exact:", res)
pickle.dump(gw, open('m8_gw.pkl','wb'))
