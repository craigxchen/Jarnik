# drv_ladder.py NAME K Qsub Qstart Qend [nodecap] : genuine-angle-budget ladder (= Rig^SL when phi is pinned, pin.py)
# For each forced Q in [Qstart, Qend]: level-1 cpc on F ∩ [3,Qsub] with budgets log|D_s(Z_v)| (PU tables: valid
# relaxation); INFEASIBLE = certificate; otherwise build_witness (deep certificates or explicit verified witness).
import sys, time, pickle
from sl import *
from conv import convert
import slverify
name, K, Qsub, Qa, Qb = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
nodecap = int(sys.argv[6]) if len(sys.argv) > 6 else 10 ** 7
S = setup(name, K); t0 = time.time()
lg = lambda s: print(s, f"[{time.time()-t0:.0f}s]", flush=True)
qs = [q for q in S.F if q <= Qsub]; set_caps(S, caps_genuine(S)); tabs = tables_tag(S, qs, "gen"); lg("GB tables: " + " ".join(f"{q}:{len(tabs[q]['sigs'])}" for q in qs))
Bud = budgets_genuine(S); full = full_mask(S)
for Q in [q for q in qs if Qa <= q <= Qb]:
    C = set(q for q in qs if q <= Q)
    st, nodes, sols = run_cpc(S, tabs, C, full, Bud, nodecap)
    lg(f"Q={Q}: level-1 genuine-budget cpc on F<={Qsub}: {st} nodes={nodes}")
    if st == "INFEASIBLE": lg(f"CERT Q={Q} level1"); continue
    if st == "ABORTED": lg(f"UNDECIDED Q={Q} (aborted)"); break
    res = build_witness(S, caps_genuine(S), "gen", Q, Qsub, log=lg, nodecap=nodecap)
    if res["status"] == "INFEASIBLE_DEEP": lg(f"CERT Q={Q} deep {res['dset']}"); continue
    if res["status"] == "WITNESS":
        Wc = convert(S, res["labels"], res["deep"], "genuine"); Wc["Q"] = Q
        fn = f"W_{name}_K{K}_gen_Q{Q}.pkl"; pickle.dump(Wc, open(fn, "wb")); lg(f"WITNESS Q={Q} saved {fn}")
        slverify.check(Wc)
    else: lg(f"UNDECIDED Q={Q} ({res['status']})")
    break
