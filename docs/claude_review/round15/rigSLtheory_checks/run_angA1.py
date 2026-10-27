"""For every ANG-good survivor outside every eps-class that is not SL at phi = arg pi_j (it has A1 extras),
run the rigorous interval B&B of slab_iv.py.  Usage: python3 run_angA1.py FFILE q0 K dumpfile p1..pr"""
import sys, json
from slab_iv import constraints, bnb
Fpath, q0, K, dump = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]; ps = [int(x) for x in sys.argv[5:]]
labs = json.load(open(dump)); res = {}
for ell in labs:
    r, nodes, wit = bnb(constraints(Fpath, q0, K, ell, ps), len(ps))
    res[r] = res.get(r, 0) + 1
    print(f"  ell={ell}: {r} ({nodes} nodes) witness={wit}", flush=True)
print(f"ps={ps} K={K} q0={q0}: {len(labs)} labels; outcomes {res}")
