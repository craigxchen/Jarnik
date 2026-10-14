# witness_table.py -- re-check every joint witness jwmix_*.pkl with check_system.check and tabulate it.
import glob, pickle, json
from check_system import check
rows = []
for fn in sorted(glob.glob("jwmix_*.pkl")):
    W = pickle.load(open(fn, "rb"))
    name = fn[len("jwmix_"):].rsplit("_K", 1)[0]
    thr = f"thr_{name}_K{W['K']}.json"
    st, epsq = check(fn, thr)
    F = sorted(st); common = None; seg = None
    for q in F:
        c2 = epsq[q] if common is None else common & epsq[q]
        if not c2: break
        common = c2; seg = q
    none = [q for q in F if st[q] == "NONE"]
    rows.append((name, W["K"], W["Q"], len(F), seg, len(none), min(none) if none else None,
                 sum(1 for q in F if q > W["Q"] and st[q] == "NONE"), sum(1 for q in F if q > W["Q"])))
print("| witness file | set | K | non-coherent on F ∩ [3,Q], Q = | |F| | coherent (one eps) on forced q <= | moduli with label outside every eps-class | of the forced moduli above Q: outside every eps-class |")
print("|---|---|---|---|---|---|---|---|")
for name, K, Q, nF, seg, nn, mn, na, tot in rows:
    print(f"| `jwmix_{name}_K{K}_{Q}.pkl` | {name} | {K} | {Q} | {nF} | {seg} | {nn} (smallest {mn}) | {na}/{tot} |")
