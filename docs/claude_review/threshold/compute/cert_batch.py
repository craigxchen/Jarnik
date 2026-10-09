"""Exact lower-bound certificates for selected shells (see cert.py).
usage: python3 cert_batch.py selection out.jsonl   selection in {top, trend, small}"""
import sys, json, glob, time
from fractions import Fraction
from collections import defaultdict
from modla import P1
from cert import certify2
sel, out = sys.argv[1], sys.argv[2]
kind = sys.argv[3] if len(sys.argv) > 3 else 'shell'
Mmax = int(sys.argv[4]) if len(sys.argv) > 4 else 9
recs = {}
for fn in glob.glob(f'data/{kind}_M*.jsonl'):
    for line in open(fn):
        r = json.loads(line)
        if r['prime'] == P1:
            recs[(r['M'], r['s'], r['t'])] = r
def vand(M):
    k = (M + 1) // 2
    return Fraction(2 * k - 1, k)
todo = []
byM = defaultdict(list)
for k, r in recs.items():
    byM[k[0]].append(r)
for M, L in byM.items():
    if M > Mmax:
        continue
    if sel == 'top':
        best = max(Fraction(r['reg'], r['t']) for r in L)
        todo += [r for r in L if Fraction(r['reg'], r['t']) == best]
    elif sel == 'trend':
        for t in sorted(set(r['t'] for r in L)):
            Lt = [r for r in L if r['t'] == t]
            b = max(r['reg'] for r in Lt)
            todo += [r for r in Lt if r['reg'] == b][:1]
    elif sel == 'small':
        todo += [r for r in L if r['t'] <= 12]
done = set()
try:
    for line in open(out):
        c = json.loads(line); done.add((c['M'], c['s'], c['t'], tuple(c['lam'])))
except FileNotFoundError:
    pass
for r in sorted(todo, key=lambda r: (r['M'], r['t'], r['s'])):
    lam = tuple(r['argmax'][-1])        # the argmax component closest to the sign rep
    key = (r['M'], r['s'], r['t'], lam)
    if key in done:
        continue
    T = time.time()
    c = certify2(r['M'], r['p'], r['n'], lam, r['reg'], kind)
    c.update(M=r['M'], s=r['s'], t=r['t'], p=r['p'], n=r['n'], lam=list(lam), reg=r['reg'], kind=kind,
             time=round(time.time() - T, 2))
    with open(out, 'a') as f:
        f.write(json.dumps(c) + '\n')
    print(c, flush=True)
