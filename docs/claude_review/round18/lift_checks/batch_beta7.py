import random, json, sys
from fractions import Fraction as Fr
import beta7lift2 as B
out = []
ts = [[Fr(2), Fr(3), Fr(-5), Fr(1)], [Fr(1, 3), Fr(3), Fr(2), Fr(5, 3)]]
rnd = random.Random(2026)
while len(ts) < 12:
    t = [Fr(rnd.randint(-9, 9), rnd.randint(1, 5)) for _ in range(3)] + [Fr(1)]
    ts.append(t)
for t in ts:
    try:
        r = B.run(t, verbose=False)
        ok = (r['gcd_deg'] == 0 and not r['common_zero_at_infinity'] and r['qstar_nullity'] == 2
              and r['special_nullities'] == [0] * 27 + [2] and r['qstar_special'])
        print([str(x) for x in t], 'gcd_deg', r['gcd_deg'], 'inf', r['common_zero_at_infinity'],
              'q*null', r['qstar_nullity'], 'specials', r['special_nullities'][-3:], 'PASS' if ok else 'FAIL', flush=True)
        out.append(dict(t=[str(x) for x in t], ok=ok, res={k: v for k, v in r.items() if k != 'minor_info'}))
    except AssertionError as e:
        print([str(x) for x in t], 'skipped (not in U / degenerate):', str(e)[:80], flush=True)
json.dump(out, open('batch_beta7.json', 'w'), indent=1)
