# mech_joint.py -- what makes the JOINT core rigid (round 14).  For (NAME, K, Q = Q_core, Qsub) re-run the exact
# coherence test on C = F ∩ [3,Q] under relaxations, each of which can only make it MORE feasible:
#   (a) drop one outer modulus q' in (Q, Qsub]           -> which outer forced moduli are needed;
#   (b) lift the budgets of all monomials of norm k      -> which norms are needed;
#   (c) lift the budgets of all (v,s) except those whose norm is k (keep only one norm).
# 'needed' = the relaxation becomes feasible (no certificate any more).
# Usage: python3 mech_joint.py NAME K Q Qsub [depthlevel]
import sys, json, os, pickle, math, subprocess
from fj import *

def run(S, tabs, C, full, drop_norms=(), keep_norms=None, nodecap=10 ** 6):
    opts, inp = cpc_input(S, tabs, C, full)
    lines = inp.split("\n")
    pairs = sorted(set(e for q in sorted(tabs) for o in opts[q] for e in o[0]))
    Bs = lines[1].split()
    for i, e in enumerate(pairs):
        nv = sum(abs(t) for t in S.V[e // 4])
        if nv in drop_norms or (keep_norms is not None and nv not in keep_norms): Bs[i] = "1e9"
    lines[1] = " ".join(Bs)
    out = subprocess.run([os.path.join(HERE, "cpc"), str(nodecap), "1"], input="\n".join(lines), capture_output=True, text=True).stdout.splitlines()
    return out[-1].split()[0]

if __name__ == '__main__':
    name, K, Q, Qsub = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    dl = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    thr = json.load(open(f"thr_{name}_K{K}.json")); P, F = thr["P"], thr["forced"]
    S = Setup(P, K, F); full = (1 << (2 ** (len(P) - 1))) - 1
    tabs = {q: pickle.load(open(f"l1cache_{name}_K{K}_{q}_None.pkl", "rb")) for q in F if q <= Qsub}
    for qd, a in (((3, 3), (5, 3), (7, 2)), ((3, 3), (5, 3), (7, 2), (11, 2), (13, 2)))[dl - 1] if dl else ():
        cf = f"deepcache_{name}_K{K}_{qd}_{a}.pkl"
        if qd in tabs: tabs[qd] = pickle.load(open(cf, "rb"))
    C = set(q for q in tabs if q <= Q)
    print(f"{name} K={K} C=F∩[3,{Q}] moduli<= {Qsub} depth {dl}: baseline {run(S, tabs, C, full)}")
    for qo in sorted(q for q in tabs if q > Q):
        t2 = {q: t for q, t in tabs.items() if q != qo}
        print(f"  drop outer modulus {qo}: {run(S, t2, C, full)}")
    for qc in sorted(q for q in tabs if q <= Q):
        t2 = {q: t for q, t in tabs.items() if q != qc}
        print(f"  drop core modulus {qc}: {run(S, t2, C - {qc}, full)}")
    for k in range(1, K + 1):
        print(f"  lift budgets of norm {k}: {run(S, tabs, C, full, drop_norms=(k,))}")
    for k in range(1, K + 1):
        print(f"  keep only norm-{k} budgets: {run(S, tabs, C, full, keep_norms=(k,))}")
