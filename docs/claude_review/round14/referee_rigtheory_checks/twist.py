"""Associate twist pi_1 -> i*pi_1: labels change by m/2 on column 1 at every modulus; check strictness (author's notion:
l' = u (eps o L) mod m, u unit, u = 1 mod 4) of the twisted genuine system, modulus by modulus."""
import itertools, math
from glib import *
from forced import forced_set
ps = [5, 13, 17, 29]; pis = [two_sq(p) for p in ps]
a, b = pis[0]; pis2 = [(-b, a)] + pis[1:]          # i*pi_1
def strict_forms(m, L):
    S = set()
    for eps in itertools.product([1, -1], repeat=len(L)):
        for u in range(1, m, 4):
            if math.gcd(u, m) == 1: S.add(tuple((u*e*x) % m for e, x in zip(eps, L)))
    return S
for K in (3, 4, 5):
    F, _ = forced_set(ps, K, 400)
    bad = []; good = []
    for q in F:
        T = Tor(q); m = T.m; n = m//4
        L = [T.ratio_log(z) for z in pis]; L2 = [T.ratio_log(z) for z in pis2]
        assert L2[0] == (L[0] + m//2) % m and L2[1:] == L[1:]
        (good if tuple(L2) in strict_forms(m, L) else bad).append(q)
    pred = [q for q in F if ((q-1 if q % 4 == 1 else q+1)//4) % 2 == 1 and nu(ps, q)[0] == 0]
    print(f"K={K}: twisted genuine system NON-strict (no eps,u=1 mod 4) at {len(bad)}/{len(F)} forced moduli: {bad}")
    print(f"   predicted by lemma (n odd and (5/q)=1): {pred}; subset: {set(pred) <= set(bad)}")
