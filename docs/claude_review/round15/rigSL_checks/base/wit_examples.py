# first one-modulus FLEX witness per (set, K)
import json, glob, os
rows = []
for fn in sorted(glob.glob("def_wit_*_K*.json")):
    b = os.path.basename(fn)[8:-5]; name, k = b.rsplit("_K", 1)
    if name.startswith(("v4", "v5", "t4")): continue
    d = json.load(open(fn))
    o = next((o for o in d if o["verdict"] == "FLEX"), None)
    if o is None: continue
    w = next(x for x in o["witnesses"] if x["ok"])
    coll = "; ".join(f"{tuple(v)}:{s}" for v, s in w["coll"][:4]) + (" ..." if len(w["coll"]) > 4 else "")
    rows.append(f"| {name} | {k} | {o['q0']} | {o['m']} | {w['ell']} | {w['ncoll']} | {coll} | {w['lift']} |")
print("| set | K | q0 | m | ell (level 1) | #collisions | collisions (v : unit s) | level-2 lift |")
print("|---|---|---|---|---|---|---|---|")
print("\n".join(rows))
