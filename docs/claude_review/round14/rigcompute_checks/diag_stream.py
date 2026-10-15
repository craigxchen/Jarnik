# diag_stream.py -- reasons why relaxation solutions do not become witnesses (streaming tables; diagnostic only)
import sys, json, pickle
from collections import Counter
from fj import *
name, K, Qsub, Q, maxsol = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
deepmap = dict(eval(sys.argv[6])) if len(sys.argv) > 6 else {}
thr = json.load(open(f"thr_{name}_K{K}.json")); P, F = thr["P"], thr["forced"]
S = Setup(P, K, F); full = (1 << (2 ** (len(P) - 1))) - 1
qs = [q for q in F if q <= Qsub]; C = set(q for q in qs if q <= Q)
st, n, sols = run_cpc_stream(S, name, K, qs, C, full, 10 ** 6, maxsol=maxsol, deepmap=deepmap, tmp=f"cpcin_diag_{name}_{K}.txt")
print(st, n, len(sols), flush=True)
R = Counter(); B = Counter()
for sol in sols:
    comp, vt = resolve(S, name, K, sol, C, deepmap); ch, vt2 = to_verify(S, vt, comp)
    ext = extend(S, dict(tabs=vt2), ch)
    if ext is None: R["extend failed"] += 1; continue
    chf, tabsf = ext; ver = verify(P, S.V, tabsf, chf, None, "mix")
    R["ok" if ver["ok"] else ver["why"].split(" at ")[0]] += 1
    if ver["ok"]:
        pickle.dump(dict(q0=Q, Q=Q, target="mix", P=P, K=K, choice=chf, labels=ver["labels"],
                         deep={q: tabsf[q]["deep"] for q in tabsf if tabsf[q].get("deep")}), open(f"jwmix_{name}_K{K}_{Q}.pkl", "wb"))
        print("WITNESS SAVED", flush=True); break
    for b in ver.get("bad") or []:
        for x in (b[2] or []): B[x] += 1
print(R); print(B.most_common(6))
