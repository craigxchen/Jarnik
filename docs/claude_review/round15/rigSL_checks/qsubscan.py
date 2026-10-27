import sys
from sl import *
name, K, Q = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
S = setup(name, K); Bud = budgets_genuine(S); full = full_mask(S)
for Qsub in [q for q in S.F if q >= Q]:
    qs = [q for q in S.F if q <= Qsub]; tabs = tables(S, qs); C = set(q for q in qs if q <= Q)
    st, nodes, sols = run_cpc(S, tabs, C, full, Bud, 10 ** 7)
    print(Qsub, st, nodes, flush=True)
    if st == "INFEASIBLE": break
