# drv_gen.py NAME K Qsub Q [dsetlevel] : SL witness at phi = arg pi_j (genuine angles), non-coherent on F ∩ [3,Q]
import sys, time, pickle
from sl import *
from conv import convert
import slverify
name, K, Qsub, Q = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
S = setup(name, K); t0 = time.time()
lg = lambda s: print(s, f"[{time.time()-t0:.0f}s]", flush=True)
res = build_witness(S, caps_genuine(S), "gen", Q, Qsub, log=lg)
lg(f"RESULT {name} K={K} Q={Q}: {res['status']}")
if res["status"] == "WITNESS":
    Wc = convert(S, res["labels"], res["deep"], "genuine"); Wc["Q"] = Q
    fn = f"W_{name}_K{K}_gen_Q{Q}.pkl"; pickle.dump(Wc, open(fn, "wb")); lg("saved " + fn)
    slverify.check(Wc)
