# table of exact forced sets
import json, glob, os, math
from rl import mq
from bh import greedy_Bh, AK
def S_h(r, h): return sum(2 ** i * math.comb(r, i) * math.comb(h, i) for i in range(0, h + 1))
rows = []
for fn in glob.glob("thr_*_K*.json"):
    b = os.path.basename(fn)[4:-5]; name, k = b.rsplit("_K", 1)
    if name.startswith(("v4", "v5", "t4")): continue
    d = json.load(open(fn)); K = int(k); r = len(d["P"]); h = K // 2
    Bs = greedy_Bh(r, h) if h >= 2 else list(range(1, r + 1))
    rows.append((r, name, K, d["P"], len(d["forced"]), max(d["forced"]), mq(max(d["forced"])) // 4, S_h(r, h), AK(K, Bs[-1]), d["Qscan"], d.get("aborted", [])))
rows.sort()
print("| set | r | K | P | |F_K^ex| | q* | n* | pigeonhole |S_h| | Thm 1(b) A_K | scanned to | status-unknown (node cap) |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for r, name, K, P, nF, qs, ns, sh, ak, Q, ab in rows:
    print(f"| {name} | {r} | {K} | {P} | {nF} | {qs} | {ns} | {sh} | {ak} | {Q} | {ab if ab else '-'} |")
