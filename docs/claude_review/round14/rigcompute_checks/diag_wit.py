# diag_wit.py -- why do level-1/deep solutions fail to become witnesses?  (diagnostic only)
import sys, json, pickle
from collections import Counter
from fj import *
name, K, Qsub, Q = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); maxsol = int(sys.argv[5])
dset = eval(sys.argv[6])
thr = json.load(open(f"thr_{name}_K{K}.json")); P, F = thr["P"], thr["forced"]
S = Setup(P, K, F); full = (1 << (2 ** (len(P) - 1))) - 1
Fs = [q for q in F if q <= Qsub]; tabs = {}
for q in Fs:
    cf = f"l1cache_{name}_K{K}_{q}_None.pkl"
    if os.path.exists(cf): tabs[q] = pickle.load(open(cf, "rb"))
    else: tabs[q] = tab_q(S, q); tabs[q]["min"] = minimal_ids(tabs[q]["sigs"]); pickle.dump(tabs[q], open(cf, "wb"))
for qd, a in dset:
    cf = f"deepcache_{name}_K{K}_{qd}_{a}.pkl"
    if os.path.exists(cf): tabs[qd] = pickle.load(open(cf, "rb"))
    else:
        tabs[qd] = tab_q(S, qd, None, deep=a); tabs[qd]["min"] = list(range(len(tabs[qd]["sigs"]))); pickle.dump(tabs[qd], open(cf, "wb"))
C = set(q for q in Fs if q <= Q)
margin = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
st, nodes, sols = run_cpc(S, tabs, C, full, 10 ** 6, maxsol=maxsol, margin=margin)
print(st, nodes, len(sols), flush=True)
reasons = Counter(); qbad = Counter()
for sol in sols:
    ch, vt = to_verify(S, tabs, sol)
    ext = extend(S, dict(tabs=vt), ch)
    if ext is None: reasons["extend failed"] += 1; continue
    chf, tabsf = ext
    ver = verify(P, S.V, tabsf, chf, None, "mix")
    if ver["ok"]:
        reasons["OK"] += 1
        pickle.dump(dict(q0=Q, Q=Q, target="mix", P=P, K=K, choice=chf, labels=ver["labels"],
                         deep={q: tabsf[q]["deep"] for q in tabsf if tabsf[q].get("deep")}), open(f"jwmix_{name}_K{K}_{Q}.pkl", "wb"))
        break
    reasons[ver["why"].split()[-1] if "persist" in ver["why"] else "budget"] += 1
    for b in ver.get("bad", []) or []:
        for (qq, lev) in (b[2] or []): qbad[(qq, lev)] += 1
print(reasons); print(qbad.most_common(8))
