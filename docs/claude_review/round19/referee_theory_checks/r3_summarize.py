import json, collections
recs = []; allok = True; lams = set()
for line in open("r3_sample.jsonl"):
    r = json.loads(line); allok &= r["all_ok"]; recs += r["recs"]; lams |= set(r["special_lams"])
print("UNIT SAMPLE: chunks 200, (A,C) pairs", len(recs), " all checks ok:", allok)
byD = collections.defaultdict(list)
for D, k, N, nc, g, ex, lam in recs: byD[D].append((g, ex, nc, N, k, lam))
for D in sorted(byD):
    L = byD[D]
    print(f"  D={D}: pairs {len(L)}, candidate lambdas evaluated {sum(x[2] for x in L)}, max gain {max(x[0] for x in L)}, "
          f"max(gain+n-Lambda_k) {max(x[1] for x in L)}, min N {min(x[3] for x in L)}, k-values {sorted(set(x[4] for x in L))}")
sp = [x for x in recs if x[6] not in ('1', '-1')]
print("  pairs whose max gain is attained at a special lambda != +-1:", len(sp), " gains there:", sorted(collections.Counter(x[4] for x in sp).items(), reverse=True)[:6])
print("  distinct special lambdas attaining a max:", len(lams))
hist = collections.Counter(x[4] for x in recs)
print("  hist of per-pair max gain (top):", sorted(hist.items(), reverse=True)[:8])
nrecs = []; tested = 0
for line in open("r4_nonunit.jsonl"):
    r = json.loads(line); nrecs += r["recs"]; tested += r["pencils"]
print("NON-UNIT SAMPLE: (A,C) pairs", len(nrecs), " pencils (A,C,lambda) evaluated", tested)
print("  max (gain + n - Lambda_{kbound}) =", max(x[3] for x in nrecs), "; max gain =", max(x[4] for x in nrecs))
print("  hist of per-pair max excess (top):", sorted(collections.Counter(x[3] for x in nrecs).items(), reverse=True)[:8])
print("  hist of per-pair max gain (top):", sorted(collections.Counter(max(0, 0) or x[4] for x in nrecs).items(), reverse=True)[:8])
