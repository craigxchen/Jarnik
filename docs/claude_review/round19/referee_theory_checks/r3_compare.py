"""Compare referee histograms of gain+n (r3_enum_D*.jsonl) with the author's c3_D*.jsonl (read as data only)."""
import json, sys, os
here = os.path.dirname(os.path.abspath(__file__))
auth = os.path.join(here, "..", "theory_checks")
for D in (1, 2, 3):
    mine = {}; tot_m = 0; maxg = -10**9; maxex = -10**9; orbs = 0
    for line in open(os.path.join(here, f"r3_enum_D{D}.jsonl")):
        r = json.loads(line); orbs += 1
        tot_m += r["triples"]; maxg = max(maxg, r["max_gain"]); maxex = max(maxex, r["max_gain_plus_n_minus_Lambda_k"])
        for v, c in r["hist_gain_plus_n"]:
            mine[v] = mine.get(v, 0) + c
    theirs = {}; tot_a = 0; maxga = -10**9
    for line in open(os.path.join(auth, f"c3_D{D}.jsonl")):
        r = json.loads(line)
        if "hist_gain_plus_n_top" not in r:
            continue
        tot_a += r.get("pencils", 0); maxga = max(maxga, r.get("max_gain", -10**9))
        for v, c in r["hist_gain_plus_n_top"]:
            theirs[v] = theirs.get(v, 0) + c
    print(f"D={D}: referee orbits {orbs}, triples {tot_m}, max gain {maxg}, max (gain+n-Lambda_k(f)) {maxex}")
    print(f"      author triples {tot_a}, max gain {maxga}")
    print(f"      hist(gain+n) referee {sorted(mine.items(), reverse=True)}")
    print(f"      hist(gain+n) author  {sorted(theirs.items(), reverse=True)}")
    print(f"      identical: {mine == theirs}")
