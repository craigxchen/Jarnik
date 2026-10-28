"""Re-verify sequentially built vectors (seq_out.jsonl) at deeper truncations, and report the
Poisson-expected number of failures in the remaining shells (evidence)."""
import sys, json, math
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from mc_gapset import enum_vectors, survivors
cfg = json.loads(sys.argv[1])
rows = [json.loads(l) for l in open(cfg["file"])]
row = [r for r in rows if r["primes"] == cfg["primes"]][0]
ps = row["primes"]; logs = [math.log(p) for p in ps]
pts = [r["angles"] for r in row["results"] if r.get("verified_T2")]
V = enum_vectors(logs, cfg["Tmax"]); w = np.abs(V) @ np.array(logs)
order = np.argsort(w); V = V[order]; w = w[order]
out = []
for T in cfg["Ts"]:
    sel = w <= T
    ok = survivors(np.array(pts), V[sel], w[sel], tol=1e-11)
    tail = 2.2 * np.exp(-w[w > T] / 2).sum()
    out.append(dict(T=T, nV=int(sel.sum()), surviving=int(ok.sum()), of=len(pts),
                    expected_failures_T_to_Tmax=round(float(tail), 4)))
print(json.dumps(dict(primes=ps, Tmax=cfg["Tmax"], checks=out)))
