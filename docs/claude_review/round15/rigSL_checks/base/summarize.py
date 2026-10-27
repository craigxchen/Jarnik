# summarize.py -- table of one-modulus results per (set, K) from thr_*.json and def_wit_*.json
import json, glob, os, sys, math
from rl import mq

def summary(name, K):
    thr = json.load(open(f"thr_{name}_K{K}.json"))
    fn = f"def_wit_{name}_K{K}.json"
    if not os.path.exists(fn): return None
    d = json.load(open(fn)); F = thr["forced"]
    ver = {o["q0"]: o["verdict"] for o in d}
    done = [q for q in F if q in ver]
    complete = len(done) == len(F)
    nonvac = [q for q in done if ver[q] != "VACUOUS"]
    flex = [q for q in done if ver[q].startswith("FLEX")]
    und = [q for q in done if ver[q] == "UNDECIDED"]
    rigid = [q for q in done if ver[q] == "RIGID"]
    # one-modulus rigid initial segment: all forced q < qf are RIGID/VACUOUS, qf = first FLEX/UNDECIDED
    first_bad = next((q for q in done if ver[q] in ("FLEX", "FLEX-L1-ONLY", "UNDECIDED")), None)
    seg = [q for q in done if first_bad is None or q < first_bad]
    q_seg = max(seg) if seg else None
    rig_above = [q for q in rigid if first_bad and q > first_bad]
    return dict(name=name, K=K, P=thr["P"], nF=len(F), qstar=max(F), complete=complete, ndone=len(done),
                first_flex=first_bad, q_seg=q_seg, n_rigid=len(rigid), rig_above=rig_above, n_flex=len(flex),
                undecided=und, nstar=mq(max(F)) // 4,
                nseg=(mq(q_seg) // 4 if q_seg else None))

if __name__ == '__main__':
    rows = []
    for fn in sorted(glob.glob("thr_*_K*.json")):
        b = os.path.basename(fn)[4:-5]; name, k = b.rsplit("_K", 1)
        if name.startswith(("v4", "v5", "t4")): continue
        s = summary(name, int(k))
        if s: rows.append(s)
    rows.sort(key=lambda s: (len(s["P"]), s["name"], s["K"]))
    print("| set | r | P | K | |F| | q* (n*) | one-mod rigid segment: last q (n) | first flexible q | rigid above it | #rigid/#flex | complete |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for s in rows:
        print(f"| {s['name']} | {len(s['P'])} | {s['P']} | {s['K']} | {s['nF']} | {s['qstar']} ({s['nstar']}) | {s['q_seg']} ({s['nseg']}) | "
              f"{s['first_flex']} | {s['rig_above']} | {s['n_rigid']}/{s['n_flex']} | {'yes' if s['complete'] else s['ndone']} |")
