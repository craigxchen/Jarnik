# batch driver: exact forced sets for all (set, r, K); Qscan = 4 A_K + 5 (Theorem 1(b) with greedy B_h set),
# so every prime q > Qscan is free by Theorem 1(b) and every q <= Qscan is decided exactly here.
import sys, subprocess, os, time
from sets import sets
from bh import greedy_Bh, AK
Ks = [int(x) for x in sys.argv[1].split(",")]; rs = [int(x) for x in sys.argv[2].split(",")]
tags = sys.argv[3].split(",")
MAXN = sys.argv[4] if len(sys.argv) > 4 else "0"
for K in Ks:
    for r in rs:
        h = K // 2; Bs = greedy_Bh(r, h) if h >= 2 else list(range(1, r + 1)); Q = 4 * AK(K, Bs[-1]) + 5
        for tag in tags:
            P = sets(r)[tag]; name = f"{tag}{r}"
            fn = f"thr_{name}_K{K}.json"
            if os.path.exists(fn): continue
            out = subprocess.run(["python3", "thresholds.py", name, str(K), ",".join(map(str, P)), str(Q), MAXN],
                                 capture_output=True, text=True).stdout
            print(out.splitlines()[0], flush=True)
