# Compare the round-19 D=4, lam=+-1 histogram with the referee's round-18 output, entry by entry.
import json, collections, sys
ref = collections.Counter()
for l in open('/home/user/Jarnik/docs/claude_review/round18/referee_theory_checks/ref_gain_D4_chunks.jsonl'):
    r = json.loads(l); ref.update({int(k): v for k, v in r['hist'].items()})
ours = collections.Counter()
for l in open(sys.argv[1]):
    r = json.loads(l); ours.update({int(k): v for k, v in r['hist_pm'].items()})
print('referee total', sum(ref.values()), 'ours total', sum(ours.values()))
print('histograms identical:', ref == ours, ' distinct gain values:', len(ref))
print('full histogram:', sorted(ours.items(), reverse=True))
