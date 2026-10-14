# probe single moduli with freescan, streaming (for timing)
import sys, subprocess, time
from rl import *
P = [int(x) for x in sys.argv[1].split(",")]; K = int(sys.argv[2]); qs = [int(x) for x in sys.argv[3].split(",")]
maxnodes = sys.argv[4] if len(sys.argv) > 4 else "0"
r = len(P); V = monos_norm(r, K)
for q in qs:
    m = mq(q); n = m // 4; need = 1 if n % 2 == 0 else 0; par = nu(P, q)
    inp = [f"{r} {len(V)}"] + [" ".join(map(str, v)) for v in V] + ["1", f"{q} {n} {need} " + " ".join(map(str, par))]
    t0 = time.time()
    out = subprocess.run(["./freescan", maxnodes], input="\n".join(inp) + "\n", capture_output=True, text=True).stdout
    print(out.strip(), f"[{time.time()-t0:.1f}s]", flush=True)
