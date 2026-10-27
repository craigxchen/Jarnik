# ladder.py -- drive fj.py with growing Q_sub (dropping moduli above Q_sub is a relaxation; tables are cached).
# Usage: python3 ladder.py NAME K [frac0=0.6]
# Final record: ladder_<NAME>_K<K>.json  (core_lower = largest Q certified coherent, witness_Q = first Q with a
# verified non-coherent good system; Q_core = core_lower exactly when witness_Q is the next forced modulus).
import sys, json, subprocess, math
from summarize import summary
name, K = sys.argv[1], int(sys.argv[2]); frac0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.6
F = json.load(open(f"thr_{name}_K{K}.json"))["forced"]; s = summary(name, K)
qseg = s["q_seg"]; qflex = s["first_flex"]
Qsub = next(q for q in F if q >= frac0 * qseg); Qstart = 0; lower = None; hist = []
while True:
    log = f"fj_{name}_K{K}_{Qsub}.log"
    subprocess.run(f"python3 -u fj.py {name} {K} {Qsub} 1000000 0 {Qstart} > {log} 2>&1", shell=True)
    r = json.load(open(f"jointcore_{name}_K{K}.json")); hist.append(dict(Qsub=Qsub, Qstart=Qstart, **{k: r[k] for k in ("core_lower", "first_noncoherent", "log")}))
    if r["core_lower"] is not None: lower = r["core_lower"]
    print(f"{name} K={K} Qsub={Qsub} Qstart={Qstart}: lower={lower} witness={r['first_noncoherent']}", flush=True)
    if r["first_noncoherent"] is not None: break
    last = r["log"][-1][0] if r["log"] else None
    if last is None or (lower is not None and lower >= max(q for q in F if q <= Qsub)):
        # every Q <= Qsub certified: continue with a larger window
        Qstart = (lower or 0) + 1
    else:
        Qstart = last
    if Qsub >= max(F): break
    nxt = [q for q in F if q >= math.ceil(1.3 * Qsub)]; Qsub = nxt[0] if nxt else max(F)
    if qflex and lower is not None and lower >= qseg: break
json.dump(dict(name=name, K=K, core_lower=lower, witness_Q=hist[-1]["first_noncoherent"], qseg=qseg, qflex1=qflex, history=hist),
          open(f"ladder_{name}_K{K}.json", "w"))
print(f"LADDER {name} K={K}: joint core >= {lower}; witness at {hist[-1]['first_noncoherent']}; one-mod segment {qseg}, q_flex1 {qflex}")
