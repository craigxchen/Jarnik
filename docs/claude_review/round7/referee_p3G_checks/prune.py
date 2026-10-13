"""Adaptive pruning against weight-3 relations (evidence).
Start from a greedy (weight<=2-clean) set; repeatedly find ALL level-0 certificates with <= 3 terms
among tight relations of weight <= 3 (brute_k3.search) and delete a hitting set of columns
(greedy: most frequent column first) until none is left. Then run the general closure (gclosure,
k<=2 and k=3, weight<=3) on the survivors and re-check every certificate.
usage: python3 prune.py primes.json M
"""
import sys, json
from collections import Counter
from rcore import gauss, Lprime, tight_list
from eval_set import triple_list
from brute_k3 import search
import gclosure

ps = json.load(open(sys.argv[1])); M = int(sys.argv[2]); L = Lprime(2 * M, set())
it = 0
while True:
    gp = [gauss(p) for p in ps]
    rels = tight_list(gp, L) + triple_list(gp, L)
    certs = search(rels)
    print(json.dumps(dict(iter=it, r=len(ps), tight=len(rels), level0_certs=len(certs))), flush=True)
    if not certs: break
    while certs:
        cnt = Counter(c for cert in certs for (idx, n) in cert for c, _ in rels[idx][0])
        col = cnt.most_common(1)[0][0]
        certs = [cert for cert in certs if all(col not in [c for c, _ in rels[idx][0]] for idx, n in cert)]
        ps = ps[:col] + ps[col + 1:]
        # re-index: rebuild everything next round (columns shift)
        break
    it += 1
lat, log = gclosure.closure(len(gp), rels)
ok, ex, hf = gclosure.check(gp, L, rels, log)
out = dict(file=sys.argv[1], final_r=len(ps), tight_w3=len(rels), certs=len(log), certs_checked=ok, free=len(ps) - hf)
print(json.dumps(out), flush=True)
json.dump(ps, open(sys.argv[1].replace('.json', '_pruned.json'), 'w'))
