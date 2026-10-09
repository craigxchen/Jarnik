"""Centrally symmetric S_M-symmetric cut polytopes: s = 0, U_j = U_{M-j} = u_j (j <= M/2).
For each shape u and dilation lam, compute reg exactly (2 primes) and the adversary value."""
import sys, itertools, time
from fractions import Fraction as Fr
from clib import reg2, proj, widths, sym_adversary, popc
from gens import SUsym
M = int(sys.argv[1]); umax = int(sys.argv[2]); maxpts = int(sys.argv[3])
h = M // 2
rows = []
for u in itertools.product(range(1, umax+1), repeat=h):
    U = [0]*(M-1)
    for j in range(1, M):
        U[j-1] = u[min(j, M-j)-1]
    S = SUsym(M, 0, U)
    if len(S) < 2 or len(S) > maxpts: continue
    w = widths(S, M)
    ws = {}
    for P, v in w.items():
        ws.setdefault(min(popc(P), M-popc(P)), set()).add(v)
    ws = {j: v.pop() for j, v in ws.items()}
    # skip non-tight profiles (some constraint redundant): require width_j == 2 u_j
    if any(ws[j] != 2*u[j-1] for j in ws): continue
    adv, mu = sym_adversary(ws, M)
    t0 = time.time()
    r = reg2(proj(S))
    rows.append((Fr(r)/(2*adv), r, adv, u, len(S), ws))
    print(f"u={u} |S|={len(S)} reg={r} widths={ws} adv={adv} ratio={float(Fr(r)/(2*adv)):.4f} ({time.time()-t0:.1f}s)", flush=True)
rows.sort(key=lambda x: -x[0])
print("TOP:")
for x in rows[:8]:
    print(f"  ratio={x[0]} = {float(x[0]):.4f} reg={x[1]} u={x[3]} |S|={x[4]} widths={x[5]} adv={x[2]}")
