"""excess(M,t) = max_s reg(S_{s,t}) - t, and the s attaining it; also per fixed s."""
import json, glob, sys
from collections import defaultdict
from modla import P1
kind = sys.argv[1] if len(sys.argv) > 1 else 'shell'
D = defaultdict(dict)
for fn in glob.glob(f'data/{kind}_M*.jsonl'):
    for line in open(fn):
        r = json.loads(line)
        if r['prime'] == P1:
            D[r['M']][(r['s'], r['t'])] = r['reg']
for M in sorted(D):
    ts = sorted(set(t for (s, t) in D[M]))
    row = []
    for t in ts:
        ex = max(D[M][(s, tt)] - t for (s, tt) in D[M] if tt == t)
        row.append(f'{t}:{ex}')
    print(f'M={M} excess max_s(reg)-t by t: ' + ' '.join(row))
    for s in range(0, 5):
        row = [f'{t}:{D[M][(s,t)]-t}' for t in ts if (s, t) in D[M]]
        print(f'   s={s}: ' + ' '.join(row))
