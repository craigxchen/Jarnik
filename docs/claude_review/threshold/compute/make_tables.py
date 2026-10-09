"""Markdown tables for compute.md from the two-prime sweep logs.
usage: python3 make_tables.py kind(shell|slice) section(maxima|trend|vand|full|phi)"""
import json, glob, sys
from fractions import Fraction
from collections import defaultdict
from math import comb, factorial
from modla import P1, P2


def vand(M):
    k = (M + 1) // 2
    return Fraction(2 * k - 1, k)


def lamstr(l):
    out = []
    i = 0
    while i < len(l):
        j = i
        while j < len(l) and l[j] == l[i]:
            j += 1
        out.append(str(l[i]) + ('^%d' % (j - i) if j - i > 1 else ''))
        i = j
    return ''.join(out)


kind = sys.argv[1]
section = sys.argv[2]
recs = defaultdict(dict)
for fn in sorted(glob.glob(f'data/{kind}_M*.jsonl')):
    for line in open(fn):
        r = json.loads(line)
        recs[(r['M'], r['s'], r['t'])][r['prime']] = r


def agree(d):
    if P1 not in d or P2 not in d:
        return False
    a = {tuple(c['lam']): (c['N'], c['reg'], c['profile']) for c in d[P1]['comps']}
    b = {tuple(c['lam']): (c['N'], c['reg'], c['profile']) for c in d[P2]['comps']}
    return a == b


rows = defaultdict(list)
for (M, s, t), d in recs.items():
    r = d.get(P1) or d.get(P2)
    rows[M].append((s, t, r, agree(d)))

if section == 'maxima':
    print('| M | V(M) | shells computed (t range) | two primes agree | max reg/t | attained at (s,t): components | any > V(M) | any > 2 |')
    print('|---|---|---|---|---|---|---|---|')
    for M in sorted(rows):
        L = rows[M]
        best = max(Fraction(r['reg'], t) for (s, t, r, _) in L)
        att = []
        for (s, t, r, _) in sorted(L, key=lambda x: (x[1], x[0])):
            if Fraction(r['reg'], t) == best:
                att.append(f"({s},{t}): {' '.join(lamstr(l) for l in r['argmax'])}")
        over = any(Fraction(r['reg'], t) > vand(M) for (s, t, r, _) in L)
        over2 = any(Fraction(r['reg'], t) > 2 for (s, t, r, _) in L)
        ts = sorted(set(t for (s, t, r, _) in L))
        print(f"| {M} | {vand(M)} | {len(L)} (t = {ts[0]}..{ts[-1]}) | {sum(1 for x in L if x[3])}/{len(L)} | "
              f"{best} | {'; '.join(att)} | {'YES' if over else 'no'} | {'YES' if over2 else 'no'} |")

elif section == 'trend':
    for M in sorted(rows):
        L = rows[M]
        print(f'**M = {M}** (V = {vand(M)} = {float(vand(M)):.4f})\n')
        print('| t | max_s reg | ratio | excess reg-t | attained at s | components (at the first such s) |')
        print('|---|---|---|---|---|---|')
        for t in sorted(set(x[1] for x in L)):
            Lt = [x for x in L if x[1] == t]
            b = max(x[2]['reg'] for x in Lt)
            ss = [x[0] for x in sorted(Lt) if x[2]['reg'] == b]
            r0 = [x for x in Lt if x[0] == ss[0]][0][2]
            mark = ' **' if Fraction(b, t) == vand(M) else ''
            print(f"| {t} | {b} | {b/t:.4f}{mark} | {b - t} | {','.join(map(str, ss))} | "
                  f"{' '.join(lamstr(l) for l in r0['argmax'])} |")
        print()

elif section == 'vand':
    # Vandermonde shells: M odd (k(k+1)/2, k(k+1)/2), M even (k(k+1)/2, k(k-1)/2)
    for M in sorted(rows):
        k = (M + 1) // 2
        if M % 2:
            kk = (M - 1) // 2
            p, n = kk * (kk + 1) // 2, kk * (kk + 1) // 2
        else:
            p, n = k * (k + 1) // 2, k * (k - 1) // 2
        s, t = p - n, p + n
        for (ss, tt, r, ag) in rows[M]:
            if ss == s and tt == t:
                comps = ', '.join(f"{lamstr(c['lam'])}: {c['reg']}" for c in r['comps'])
                print(f"* M = {M}, (p,n) = ({p},{n}), t = {t}, |S| = {r['size']}, reg = {r["reg"]}, C({M},2) = {comb(M, 2)}"
                      f" {'(two primes agree)' if ag else '(ONE PRIME)'}; components: {comps}")

elif section == 'phi':
    tot = eq = 0
    for M in sorted(rows):
        L = rows[M]
        e = sum(1 for (s, t, r, _) in L if r['reg'] == r['phi_floor'])
        mx = max(r['reg'] - r['phi_floor'] for (s, t, r, _) in L)
        print(f'M={M}: {len(L)} {kind}s, reg == floor(phi) in {e}, max(reg - floor(phi)) = {mx}')

elif section == 'full':
    for M in sorted(rows):
        print(f'### M = {M}\n')
        print('| s | t | (p,n) | size | reg | reg/t | attained by | floor(phi) | 2 primes | all components lam:reg_lam |')
        print('|---|---|---|---|---|---|---|---|---|---|')
        for (s, t, r, ag) in sorted(rows[M], key=lambda x: (x[1], x[0])):
            comps = ' '.join(f"{lamstr(c['lam'])}:{c['reg']}" for c in r['comps'])
            print(f"| {s} | {t} | ({r['p']},{r['n']}) | {r['size']} | {r['reg']} | {r['reg']/t:.4f} | "
                  f"{' '.join(lamstr(l) for l in r['argmax'])} | {r['phi_floor']} | {'yes' if ag else 'NO'} | {comps} |")
        print()

if section == 'trendline':
    print('| M | t : max_s reg (ratio; `*` = equals V(M)) |')
    print('|---|---|')
    for M in sorted(rows):
        L = rows[M]
        parts = []
        for t in sorted(set(x[1] for x in L)):
            b = max(x[2]['reg'] for x in L if x[1] == t)
            star = '*' if Fraction(b, t) == vand(M) else ''
            parts.append(f'{t}:{b} ({b/t:.3f}{star})')
        print(f'| {M} | ' + ', '.join(parts) + ' |')

if section == 'excess':
    print('| M | max_s reg(S_{s,t}) - t for t = 1, 2, ... |')
    print('|---|---|')
    for M in sorted(rows):
        L = rows[M]
        parts = []
        for t in sorted(set(x[1] for x in L)):
            b = max(x[2]['reg'] for x in L if x[1] == t)
            parts.append(str(b - t))
        print(f'| {M} | ' + ' '.join(parts) + ' |')
