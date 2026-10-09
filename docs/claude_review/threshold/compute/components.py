"""Statistics on the components attaining reg(A)."""
import json, glob, sys
from collections import Counter, defaultdict
from modla import P1
kind = sys.argv[1] if len(sys.argv) > 1 else 'shell'
def is_hook(l): return all(x == 1 for x in l[1:])
stats = defaultdict(Counter)
for fn in sorted(glob.glob(f'data/{kind}_M*.jsonl')):
    for line in open(fn):
        r = json.loads(line)
        if r['prime'] != P1: continue
        M = r['M']; am = [tuple(l) for l in r['argmax']]
        st = stats[M]
        st['shells'] += 1
        hasign = any(tuple(c['lam']) == tuple([1]*M) for c in r['comps'])
        st['sign present'] += hasign
        st['sign attains'] += (tuple([1]*M) in am)
        st['some hook attains'] += any(is_hook(l) for l in am)
        st['unique argmax'] += (len(am) == 1)
        st['trivial attains'] += ((M,) in am)
        # if the shell max is attained by a non-hook only
        if not any(is_hook(l) for l in am): st['only non-hooks attain'] += 1
for M in sorted(stats):
    print(M, dict(stats[M]))
