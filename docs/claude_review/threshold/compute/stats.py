"""Small statistics for compute.md.  usage: python3 stats.py key"""
import json, glob, sys, re
from collections import Counter, defaultdict
from modla import P1, P2
key = sys.argv[1]


def load(kind):
    recs = defaultdict(dict)
    for fn in glob.glob(f'data/{kind}_M*.jsonl'):
        for line in open(fn):
            r = json.loads(line)
            recs[(r['M'], r['s'], r['t'])][r['prime']] = r
    return recs


if key in ('nshell', 'nslice'):
    recs = load('shell' if key == 'nshell' else 'slice')
    both = sum(1 for d in recs.values() if P1 in d and P2 in d)
    print(f'{len(recs)} ({both} with both primes)')

if key == 'crosscheck':
    n = ok = 0
    for fn in ('test_iso_log4.txt', 'data/log_crosscheck_big.txt'):
        try:
            for line in open(fn):
                if 'profiles equal' in line:
                    n += 1
                    ok += line.strip().endswith('True')
        except FileNotFoundError:
            pass
    print(f'{ok} of {n} shells agree')

if key == 'certs':
    seen = {}
    for fn in glob.glob('data/certs_*.jsonl'):
        for line in open(fn):
            c = json.loads(line)
            k = (c['M'], c['s'], c['t'], tuple(c['lam']))
            seen[k] = c
    byM = defaultdict(Counter)
    for k, c in seen.items():
        lvl = 'elementary' if c.get('level', '').startswith('elem') else ('ATY' if c.get('level', '').startswith('ATY') else 'other')
        byM[c['M']][(lvl, c['ok'])] += 1
    print('| M | certified (elementary) | certified (ATY level) | failed |')
    print('|---|---|---|---|')
    for M in sorted(byM):
        C = byM[M]
        fail = sum(v for (l, ok), v in C.items() if ok is not True)
        print(f"| {M} | {C[('elementary', True)]} | {C[('ATY', True)]} | {fail} |")

if key == 'ranges':
    recs = load(sys.argv[2])
    byM = defaultdict(list)
    for (M, s, t), d in recs.items():
        byM[M].append((t, s))
    for M in sorted(byM):
        ts = sorted(set(t for t, s in byM[M]))
        full = [t for t in ts if all((t, s) in byM[M] for s in range(t % 2, t + 1, 2))]
        part = [t for t in ts if t not in full]
        desc = f'all s for t <= {max(full)}' if full else ''
        if part:
            desc += '; partial s for t = ' + ', '.join(
                f"{t} (s = {','.join(str(s) for (tt, s) in sorted(byM[M]) if tt == t)})" for t in part)
        print(f'* M = {M}: {len(byM[M])} shells/slices, {desc}')
