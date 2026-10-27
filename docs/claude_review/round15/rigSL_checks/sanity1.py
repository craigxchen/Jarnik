# sanity1.py NAME K Qsub : (i) cpc with genuine budgets and NO coherence requirement must be feasible;
# (ii) the genuine labels' level-1 signature at every q is an option of the PU table and fits the genuine budgets.
import sys
from sl import *
name, K, Qsub = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
S = setup(name, K); qs = [q for q in S.F if q <= Qsub]; tabs = tables(S, qs)
Bud = budgets_genuine(S)
st, nodes, sols = run_cpc(S, tabs, set(), 0, Bud, 10 ** 6)
print("no-coherence-requirement run:", st, nodes)
load = [0.0] * len(Bud)
for q in qs:
    m, n, L, iota = modulus_data(S, q)
    sig = []
    for i, v in enumerate(S.V):
        d = sum(a * b for a, b in zip(v, L)) % m
        if d % n == 0:
            s = [u for u in range(4) if (u * iota) % m == d][0]; sig.append(4 * i + s)
            assert S.D[i][s] % q == 0
    sig = tuple(sorted(sig)); T = tabs[q]
    present = any(tuple(el) == sig for el, mu, mks in T["sigs"])
    for e in sig: load[e] += math.log(q)
    assert present, q
print("genuine signatures present in all tables; max load-budget =", max(l - b for l, b in zip(load, Bud)))
