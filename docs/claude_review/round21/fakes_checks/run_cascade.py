import itertools, json, subprocess, sys
base = [(5,(2,1))]
rest = [(13,(3,2)),(17,(4,1)),(29,(5,2)),(37,(6,1)),(41,(5,4)),(53,(7,2)),(61,(6,5)),(73,(8,3))]
target = int(sys.argv[1]); T = sys.argv[2]
idx = [p for p,_ in rest].index(target)
bgs = rest[:idx]
allok = True; n = 0
for signs in itertools.product([1,-1], repeat=len(bgs)):
    bb = [[2,1]] + [[b[0], s*b[1]] for (p,b), s in zip(bgs, signs)]
    cfg = json.dumps(dict(bg=[5]+[p for p,_ in bgs], bg_beta=bb, p=target, T=float(T)))
    out = subprocess.run(["python3", "cascade_cert.py", cfg], capture_output=True, text=True).stdout
    d = json.loads(out); n += 1
    with open("cascade_cert_out.jsonl", "a") as f: f.write(out)
    allok &= d["certified"]
    print(bb, target, "windows", d["windows"], "uncov", d["uncovered_intervals"], "bad", d["bad_gaps"][:2], "cert", d["certified"], flush=True)
print("ALL", target, n, "branches certified:", allok)
