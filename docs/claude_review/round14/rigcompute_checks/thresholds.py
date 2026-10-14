# thresholds.py -- exact forced/free status of every prime q <= Qscan (q not dividing 2N) for V_K, via ./freescan.
# Usage: python3 thresholds.py TAG K p1,p2,... Qscan [maxnodes]
# Writes thr_<TAG>_K<K>.txt : one line per q (FREE vector / FORCED), and a summary line.
# Every FREE vector is re-verified here in numpy (exact integer arithmetic).
import sys, os, time, subprocess, json
import numpy as np
from rl import *

def scan(P, K, qs, maxnodes=0):
    r = len(P); V = monos_norm(r, K)
    inp = [f"{r} {len(V)}"] + [" ".join(map(str, v)) for v in V] + [str(len(qs))]
    meta = {}
    for q in qs:
        m = mq(q); n = m // 4; need = 1 if n % 2 == 0 else 0
        par = nu(P, q); meta[q] = (n, need, par)
        inp.append(f"{q} {n} {need} " + " ".join(map(str, par)))
    out = subprocess.run([os.path.join(HERE, "freescan"), str(maxnodes)], input="\n".join(inp) + "\n",
                         capture_output=True, text=True).stdout
    Va = np.array(V, dtype=np.int64)
    res = {}
    for ln in out.splitlines():
        t = ln.split(); q = int(t[0]); n = int(t[1]); st = t[2]
        if st == "FREE":
            c = np.array([int(x) for x in t[3:3 + r]], dtype=np.int64)
            nn, need, par = meta[q]
            assert nn == n
            if need: assert all(int(c[j]) % 2 == par[j] for j in range(r))
            assert ((Va @ c) % n != 0).all(), ("bad free vector", q)
            res[q] = ("FREE", [int(x) for x in c])
        else:
            res[q] = (st, None)
    assert set(res) == set(qs)
    return res

if __name__ == '__main__':
    tag, K = sys.argv[1], int(sys.argv[2]); P = [int(x) for x in sys.argv[3].split(',')]
    Qscan = int(sys.argv[4]); maxnodes = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    qs = [q for q in primerange(3, Qscan + 1) if q not in P]
    t0 = time.time(); res = scan(P, K, qs, maxnodes)
    forced = [q for q in qs if res[q][0] == "FORCED"]; ab = [q for q in qs if res[q][0] == "ABORT"]
    fn = os.path.join(HERE, f"thr_{tag}_K{K}.json")
    json.dump(dict(P=P, K=K, Qscan=Qscan, forced=forced, aborted=ab,
                   free={q: res[q][1] for q in qs if res[q][0] == "FREE"}), open(fn, "w"))
    print(f"{tag} K={K} r={len(P)} P={P} Qscan={Qscan}: #forced={len(forced)} q*={max(forced) if forced else None} "
          f"aborted={ab} [{time.time()-t0:.1f}s]")
    print("  forced:", forced)
