import sys, time
from sl import *
name, K, Qsub = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
S = setup(name, K); qs = [q for q in S.F if q <= Qsub]
t0=time.time(); tabs = tables(S, qs); print("tables", time.time()-t0, {q:(len(tabs[q]['sigs']),len(tabs[q]['min'])) for q in qs}, flush=True)
full = full_mask(S)
for lab, Bud in (("PU", budgets_PU(S)), ("GEN", budgets_genuine(S))):
    for Q in [q for q in qs if q <= 31]:
        C = set(q for q in qs if q <= Q)
        t0=time.time(); st, nodes, sols = run_cpc(S, tabs, C, full, Bud, nodecap=2*10**6)
        print(lab, Q, st, nodes, f"{time.time()-t0:.1f}s", flush=True)
        if st != "INFEASIBLE": break
