"""Merge the two-prime sweep logs and produce the tables of compute.md.
usage: python3 summarize.py [kind=shell] > summary.txt"""
import json, glob, sys
from fractions import Fraction
from math import comb, factorial
from collections import defaultdict
from modla import P1, P2

kind = sys.argv[1] if len(sys.argv) > 1 else 'shell'


def shell_size(M, p, n):
    tot = 0
    for a in range(0, M + 1):
        for b in range(0, M + 1 - a):
            ca = (1 if p == 0 else 0) if a == 0 else (comb(p - 1, a - 1) if p >= a else 0)
            cb = (1 if n == 0 else 0) if b == 0 else (comb(n - 1, b - 1) if n >= b else 0)
            tot += factorial(M) // (factorial(a) * factorial(b) * factorial(M - a - b)) * ca * cb
    return tot


def vand(M):
    k = (M + 1) // 2
    return Fraction(2 * k - 1, k)      # 2 - 1/ceil(M/2)


def lamstr(l):
    # compact notation, e.g. (3,1,1) -> 31^2
    out = []
    i = 0
    while i < len(l):
        j = i
        while j < len(l) and l[j] == l[i]:
            j += 1
        out.append(str(l[i]) + ('^%d' % (j - i) if j - i > 1 else ''))
        i = j
    return ''.join(out) if all(x < 10 for x in l) else ','.join(map(str, l))


recs = defaultdict(dict)
for fn in sorted(glob.glob(f'data/{kind}_M*.jsonl')):
    for line in open(fn):
        r = json.loads(line)
        recs[(r['M'], r['s'], r['t'])][r['prime']] = r

problems = []
rows = []
for key in sorted(recs):
    M, s, t = key
    d = recs[key]
    both = P1 in d and P2 in d
    r = d.get(P1) or d.get(P2)
    if both:
        a, b = d[P1], d[P2]
        ca = {tuple(c['lam']): (c['N'], c['reg'], c['profile']) for c in a['comps']}
        cb = {tuple(c['lam']): (c['N'], c['reg'], c['profile']) for c in b['comps']}
        if ca != cb:
            problems.append(('prime disagreement', key))
    if kind == 'shell' and r['size'] != shell_size(M, r['p'], r['n']):
        problems.append(('size mismatch', key, r['size']))
    if r['reg'] > r['phi_floor']:
        problems.append(('exceeds floor(phi)', key, r['reg'], r['phi_floor']))
    if any(c.get('beyond_phi') for c in r['comps']):
        problems.append(('component beyond phi', key))
    rows.append((M, s, t, r, both))

print('# kind =', kind)
print('# problems:', problems if problems else 'none')
print()
# per-M overview
byM = defaultdict(list)
for row in rows:
    byM[row[0]].append(row)
print('## per-M maxima')
print('| M | t range | shells | both primes | max reg/t | attained at (s,t) | components attaining | V(M) | any > V(M)? | any > 2? |')
print('|---|---|---|---|---|---|---|---|---|---|')
for M in sorted(byM):
    L = byM[M]
    best = max(Fraction(r['reg'], t) for (_, s, t, r, _) in L)
    at = [(s, t) for (_, s, t, r, _) in L if Fraction(r['reg'], t) == best]
    comps = set()
    for (_, s, t, r, _) in L:
        if Fraction(r['reg'], t) == best:
            for l in r['argmax']:
                comps.add(lamstr(l))
    nb = sum(1 for x in L if x[4])
    tr = (min(x[2] for x in L), max(x[2] for x in L))
    over = any(Fraction(r['reg'], t) > vand(M) for (_, s, t, r, _) in L)
    over2 = any(Fraction(r['reg'], t) > 2 for (_, s, t, r, _) in L)
    print(f'| {M} | {tr[0]}..{tr[1]} | {len(L)} | {nb} | {best} = {float(best):.4f} | {at} | {sorted(comps)} | '
          f'{vand(M)} = {float(vand(M)):.4f} | {"YES" if over else "no"} | {"YES" if over2 else "no"} |')
print()

# trend: for each M and t, max over s of reg/t and the argmax s
print('## trend in t (max over s of reg/t)')
for M in sorted(byM):
    L = byM[M]
    ts = sorted(set(x[2] for x in L))
    line = []
    for t in ts:
        cand = [(Fraction(r['reg'], t), s, r['reg']) for (_, s, tt, r, _) in L if tt == t]
        bst = max(cand)
        line.append(f't={t}: {bst[2]}/{t}={float(bst[0]):.3f} (s={",".join(str(c[1]) for c in cand if c[0]==bst[0])})')
    print(f'M={M}: ' + '; '.join(line))
print()

# full table per M: s, t, |S|, reg, ratio, argmax, phi, per-component regs
print('## full table')
for M in sorted(byM):
    print(f'### M = {M}')
    print('| s | t | (p,n) | |S| | reg | reg/t | attained by | floor(phi) | both primes | component regs (lam:reg) |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    for (_, s, t, r, both) in sorted(byM[M], key=lambda x: (x[2], x[1])):
        comps = ' '.join(f"{lamstr(c['lam'])}:{c['reg']}" for c in r['comps'])
        print(f"| {s} | {t} | ({r['p']},{r['n']}) | {r['size']} | {r['reg']} | {r['reg']/t:.4f} | "
              f"{' '.join(lamstr(l) for l in r['argmax'])} | {r['phi_floor']} | {'yes' if both else 'NO'} | {comps} |")
    print()
