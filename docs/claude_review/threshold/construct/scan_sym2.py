import sys, itertools
from fractions import Fraction as Fr
from clib import reg2, proj, widths, sym_adversary, popc
from gens import SUsym
M = int(sys.argv[1]); Umax = int(sys.argv[2]); maxpts = int(sys.argv[3])
res = []; seen = set()
for s in range(0, M):
    for U in itertools.product(range(-Umax, Umax+1), repeat=M-1):
        S = SUsym(M, s, U)
        if len(S) < 2 or len(S) > maxpts: continue
        key = tuple(sorted(S))
        if key in seen: continue
        seen.add(key)
        w = widths(S, M)
        ws = {}
        for P, v in w.items():
            ws.setdefault(min(popc(P), M-popc(P)), set()).add(v)
        assert all(len(v)==1 for v in ws.values())
        ws = {j: v.pop() for j, v in ws.items()}
        adv, mu = sym_adversary(ws, M)
        r = reg2(proj(S))
        res.append((Fr(r)/(2*adv), r, adv, s, U, len(S), ws, mu))
res.sort(key=lambda x: -x[0])
for x in res[:12]:
    print(f"ratio={float(x[0]):.4f} reg={x[1]} adv={x[2]} s={x[3]} U={x[4]} |S|={x[5]} widths={x[6]} mu={ {j:str(v) for j,v in x[7].items()} }")
print("distinct sets:", len(res), " max ratio:", res[0][0])
