import sys, time, pickle
from sl import *
name, K, Qsub, Q = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
S = setup(name, K)
t0 = time.time()
res = build_witness(S, caps_genuine(S), "gen", Q, Qsub)
print(res["status"], f"{time.time()-t0:.0f}s")
if res["status"] == "WITNESS":
    pickle.dump(dict(P=S.P, K=K, Q=Q, phi="genuine", labels=res["labels"], deep=res["deep"]), open(f"wit_gen_{name}_K{K}_{Q}.pkl", "wb"))
