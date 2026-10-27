# batch: joint cores.  Qsub ladder: start at the one-modulus segment end q_seg (dropping moduli is a relaxation),
# enlarge (x1.3 steps, up to q*) while the result is neither 'core >= q_seg' nor a verified witness.
import sys, json, subprocess, math, os
from summarize import summary
jobs = [x.split(":") for x in sys.argv[1].split(",")]   # name:K
cexp = sys.argv[2] if len(sys.argv) > 2 else "0"
for jb in jobs:
    name, K = jb[0], int(jb[1])
    s = summary(name, K); qseg = s["q_seg"]; F = json.load(open(f"thr_{name}_K{K}.json"))["forced"]
    Qsub = qseg
    while True:
        args = ["python3", "fj.py", name, str(K), str(Qsub), "1000000"] + ([cexp] if cexp != "0" else [])
        out = subprocess.run(args, capture_output=True, text=True)
        open(f"log_joint_{name}_K{K}_{Qsub}.txt", "w").write(out.stdout + out.stderr)
        res = json.load(open(f"jointcore_{name}_K{K}.json"))
        print(name, K, f"Qsub={Qsub}", [l for l in out.stdout.splitlines() if l.startswith("RESULT")], flush=True)
        if res["first_noncoherent"] is not None or (res["core_lower"] or 0) >= qseg or Qsub >= max(F): break
        nxt = [q for q in F if q >= math.ceil(1.3 * Qsub)]
        Qsub = nxt[0] if nxt else max(F)
