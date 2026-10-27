# drv_repair.py NAME K Qsub [seeds] : GB witness at Q0 = 11 by cpc (level 1 + deep) + extension, then exact repair of
# the remaining violations by local search at full depth (ls.py), then independent check (slverify).
import sys, time, pickle
from sl import *
import ls, slverify
name, K, Qsub = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); Q = 11
nseeds = int(sys.argv[4]) if len(sys.argv) > 4 else 5
S = setup(name, K); t0 = time.time(); lg = lambda s: print(s, f"[{time.time()-t0:.0f}s]", flush=True)
capm = caps_genuine(S); set_caps(S, capm)
qs = [q for q in S.F if q <= Qsub]; tabs = tables_tag(S, qs, "gen")
Bud = budgets_from_caps(S, capm); full = full_mask(S); C = set(q for q in qs if q <= Q)
for ds in (((3, 2), (5, 2), (7, 2)), ((3, 3), (5, 3), (7, 2))):
    td = dict(tabs)
    for qd, a in ds:
        if qd in td: td[qd] = deep_tag(S, qd, a, "gen")
    st, nodes, sols = run_cpc(S, td, C, full, Bud, 10 ** 6, maxsol=nseeds)
    lg(f"deep {ds}: {st} nodes={nodes} sols={len(sols)}")
    for k, sol in enumerate(sols):
        ch, vt = fj.to_verify(S, td, sol)
        ext = jointcore.extend(S, dict(tabs=vt), ch)
        if ext is None: lg("  extend failed"); continue
        chf, tabsf = ext
        init = {q: (chf[q]["rep"], tabsf[q].get("M", tabsf[q]["m"])) for q in chf}
        S2 = setup(name, K)
        W = ls.search(S2, Q, seed=k + 1, maxit=20000, init=init, log=lg, C=1500)
        if W is None: lg("  repair failed"); continue
        fn = f"W_rep_{name}_K{K}_Q{Q}.pkl"; pickle.dump(W, open(fn, "wb")); lg("WITNESS saved " + fn)
        slverify.check(W); sys.exit(0)
lg("RESULT: no witness")
