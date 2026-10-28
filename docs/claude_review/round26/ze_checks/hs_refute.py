"""Prop 4.5 (Hilbert-sharp variant HS is false on T=1 regime boxes): exact big-integer verification that
NO connected subgroup H != G' = G_a^d x G_m satisfies  card(phi(A) mod H) * h_H(D) <= h_G'(D),  D = (m-1, 1),
at the cheap-subtorus point Y (one column direction a_p, L = h_A(a_p) in [3, m-1]).
Lower bounds used (proved in ze.md Prop 4.5):
  card for V x G_m (dim V = e <= d-1):  >= LB(e);   card for V x {1} (dim V = e <= d):  >= max(L+1, LB(e));
  LB(e) = max(3^{min(d-e, floor(w/2))}, ceil(ceil(3^d/2)/3^e));
  h_{V x G_m} = 2 binom(m-1+e, e), h_{V x 1} = binom(m-1+e, e), h_{G'} = 2 binom(m-1+d, d).
L is taken = 3 (worst case allowed)."""
import math
from mpmath import mp, sqrt, log, ceil, floor
mp.dps=30
def m_max_exact(d):
    lo,hi=0,2*d; t=3**d
    while lo<hi:
        mid=(lo+hi+1)//2
        if 8*math.comb(mid+d-1,d)<=t: lo=mid
        else: hi=mid-1
    return lo
for n in [2966,4096]:
    d=n-1; w=sqrt(8*d*log(4*n)); hw=int(floor(w/2)); L=3
    A_lb=(3**d+1)//2
    for m in (int(ceil(4*w)), m_max_exact(d)):
        hG=2*math.comb(m-1+d,d)
        worst_ratio=None; ok=True
        pw3=[1]*(d+1)
        for i in range(1,d+1): pw3[i]=pw3[i-1]*3
        b=1  # binom(m-1+e,e) iteratively
        for e in range(0,d+1):
            if e>0: b=b*(m-1+e)//e
            LB=max(pw3[min(d-e,hw)], -(-A_lb//pw3[e]))
            # V x {1}
            c1=max(L+1,LB)*b
            if c1<=hG: ok=False; print("  V x 1 obstruction at e",e)
            if e<=d-1:
                c2=LB*2*b
                if c2<=hG: ok=False; print("  V x Gm obstruction at e",e)
            r=min(c1, c2 if e<=d-1 else c1)
            # track min of log(lhs/hG)
            lr=math.log(r)-math.log(hG) if r>0 else -1e9
            if worst_ratio is None or lr<worst_ratio[0]: worst_ratio=(lr,e)
        print(f"n={n} m={m} floor(w/2)={hw}: HS refuted (no obstruction subgroup) = {ok}; min ln(LHS/RHS) = {worst_ratio[0]:.3f} at dim V = {worst_ratio[1]}")
