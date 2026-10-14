# joint-core summary table (from joint_results.json, thresholds, one-modulus summaries, local ranges)
import json, math
from rl import mq
from summarize import summary
LOCE = {"first4": 31, "first5": 11, "first6": 11, "first7": 11, "rand4": 7, "bal4": 7, "last4": 7}
def S_h(r, h): return sum(2 ** i * math.comb(r, i) * math.comb(h, i) for i in range(h + 1))
d = json.load(open("joint_results.json"))
print("| set | r | K | Q_core | n(Q_core) | first non-coherent Q (witness) | LOCe | one-mod seg (n) | q_flex1 | n* | pigeonhole |S_h| | n(Q_core)/|S_h| |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
for key in sorted(d, key=lambda k: (len(k), k)):
    name, K = key.split(":"); K = int(K); e = d[key]; s = summary(name, K); r = len(s["P"])
    core = e.get("core", e.get("core_lower")); exact = "core" in e
    n = mq(core) // 4
    print(f"| {name} | {r} | {K} | {'' if exact else '>= '}{core} | {'' if exact else '>= '}{n} | {e['witness'] or '-'} | {LOCE.get(name, '-')} | "
          f"{s['q_seg']} ({s['nseg']}) | {s['first_flex']} | {s['nstar']} | {S_h(r, K // 2)} | {'' if exact else '>= '}{n / S_h(r, K // 2):.2f} |")
