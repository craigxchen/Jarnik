"""Scan S_M-symmetric cut-polytopes S(U) = {sum k = s, k(P) <= U_|P|} and compare reg(S) with
2 * (best half-separating adversary payoff).  Prints the best ratios found."""
import sys, itertools
from fractions import Fraction as Fr
from clib import reg2, proj, widths, sym_adversary, popc, adversary
from gens import SU
M = int(sys.argv[1]); Umax = int(sys.argv[2]); box = int(sys.argv[3])
res = []
seen = set()
for s in range(0, M):
    for U in itertools.product(range(-Umax, Umax+1), repeat=M-1):
        S = SU(M, s, lambda P: U[len(P)-1], box)
        if len(S) < 2: continue
        key = (s, tuple(sorted(S)))
        if key in seen: continue
        seen.add(key)
        if max(max(abs(x) for x in k) for k in S) >= box: continue   # avoid truncation
        w = widths(S, M)
        ws = {}
        for P, v in w.items():
            j = min(popc(P), M - popc(P))
            ws.setdefault(j, set()).add(v)
        assert all(len(v)==1 for v in ws.values())
        ws = {j: v.pop() for j, v in ws.items()}
        adv, mu = sym_adversary(ws, M)
        r = reg2(proj(S))
        res.append((Fr(r) / (2*adv), r, adv, s, U, len(S), ws))
res.sort(key=lambda x: -x[0])
for x in res[:15]:
    print(f"ratio={float(x[0]):.4f} reg={x[1]} adv={x[2]} s={x[3]} U={x[4]} |S|={x[5]} widths={x[6]}")
print("number of distinct sets:", len(res))
