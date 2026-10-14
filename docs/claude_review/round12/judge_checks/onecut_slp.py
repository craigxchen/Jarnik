"""J6d: one cut. Check J(A_S(L)) = (x_i^3 (all i), sigma_S^{2L+1}) in ALL degrees (cumulative sandwich mod p), i.e.
the one-cut leading-form ideal is the monomial complete intersection plus one power of sigma_S, and its Hilbert function
is the Lefschetz truncation h_i - h_{i-2L-1} of (1+t+t^2)^s times (1+t+t^2)^{d-s}."""
import itertools, math, sys
import numpy as np
sys.path.insert(0, '.')
from jw_test import PR, eval_rank_profile, onecut
from jproj_test import powl, jdims

def poly_series(s, d, L):
    h = [1]
    for _ in range(s):
        h = [sum(h[i-t] for t in range(3) if 0 <= i-t < len(h)) for i in range(len(h)+2)]
    hs = [max(h[i] - (h[i-2*L-1] if i-2*L-1 >= 0 else 0), 0) if i <= s+L else 0 for i in range(len(h))]
    hs = [x for x in hs]
    g = [1]
    for _ in range(d-s):
        g = [sum(g[i-t] for t in range(3) if 0 <= i-t < len(g)) for i in range(len(g)+2)]
    out = [0]*(len(hs)+len(g))
    for a, x in enumerate(hs):
        for b, y in enumerate(g): out[a+b] += x*y
    while out and out[-1] == 0: out.pop()
    return out

out = []
for d, S, L in [(4,[0,1,2,3],1),(5,[0,1,2,3,4],1),(5,[0,1,2,3,4],2),(6,[0,1],1),(6,[0,1,2,3],1),(6,[0,1,2,3,4,5],2)]:
    A = onecut(d, S, L)
    rho = L + len(S) + 2*(d-len(S))
    ranks = eval_rank_profile(A, d, rho)
    H = [ranks[0]] + [ranks[i]-ranks[i-1] for i in range(1, rho+1)]
    cache = {}
    g = {3: [powl(tuple(1 if y == x else 0 for y in range(d)), 3, d, cache) for x in range(d)]}
    sig = tuple(1 if x in S else 0 for x in range(d))
    g.setdefault(2*L+1, []).append(powl(sig, 2*L+1, d, cache))
    J0 = jdims(g, d, rho)
    N = [math.comb(i+d-1, d-1) for i in range(rho+1)]
    HJ0 = [N[i]-J0[i] for i in range(rho+1)]
    pred = poly_series(len(S), d, L)
    s = (f"d={d} S={S} L={L} |A|={len(A)} sum H_A={sum(H)} rho<= {rho}\n  H_A      ={H}\n  H(S/J0)  ={HJ0}\n  Lefschetz={pred[:rho+1]}\n"
         f"  J(A)=J0 in all degrees (cum sandwich): {all(a==b for a,b in zip(H,HJ0)) and sum(H)==len(A)}")
    print(s, flush=True); out.append(s)
open('onecut_slp_out.txt','w').write("\n".join(out)+"\n")
