"""Compare with earlier independent (direct, non-equivariant) computations in scratchpad/tau:
verify.md's e1_M*_shell.txt / e1_M*_ball*.txt and upper.md's upper_checks/log_M*.txt."""
import re, json, glob
from collections import defaultdict
from modla import P1
mine = {}
for kind in ('shell', 'slice'):
    for fn in glob.glob(f'data/{kind}_M*.jsonl'):
        for line in open(fn):
            r = json.loads(line)
            if r['prime'] == P1:
                mine[(kind, r['M'], r['p'], r['n'])] = r['reg']
cmp = defaultdict(lambda: [0, 0, []])
for fn in glob.glob('../e1_M*_shell.txt'):
    for line in open(fn):
        m = re.match(r'M=(\d+) shell D=(\d+) s=(-?\d+) \|A\|=(\d+) reg=(\d+)', line)
        if m:
            M, D, s, A, reg = map(int, m.groups())
            p, n = (D + s) // 2, (D - s) // 2
            key = ('shell', M, max(p, n), min(p, n)) if ('shell', M, max(p, n), min(p, n)) in mine else ('shell', M, p, n)
            if key in mine:
                c = cmp[('verify.md e1 shell', M)]
                c[0] += 1
                if mine[key] != reg: c[2].append((key, reg, mine[key]))
                else: c[1] += 1
for fn in glob.glob('../e1_M*_ball*.txt'):
    for line in open(fn):
        m = re.match(r'M=(\d+) ball D=(\d+) s=(-?\d+) \|A\|=(\d+) reg=(\d+)', line)
        if m:
            M, D, s, A, reg = map(int, m.groups())
            p, n = (D + s) // 2, (D - s) // 2
            for key in (('slice', M, p, n), ('slice', M, n, p)):
                if key in mine:
                    c = cmp[('verify.md e1 ball', M)]
                    c[0] += 1
                    if mine[key] != reg: c[2].append((key, reg, mine[key]))
                    else: c[1] += 1
                    break
for fn in glob.glob('../upper_checks/log_M*.txt'):
    for line in open(fn):
        m = re.match(r'M=(\d+) (shell|ball)\s+\(p,n\)=\((\d+),(\d+)\)\s+\|A\|=\s*(\d+) reg_q1=(\d+) reg_q2=(\d+)', line)
        if m:
            M = int(m.group(1)); kind = 'shell' if m.group(2) == 'shell' else 'slice'
            p, n, A, r1, r2 = map(int, m.groups()[2:])
            for key in ((kind, M, p, n), (kind, M, n, p)):
                if key in mine:
                    c = cmp[('upper.md log ' + kind, M)]
                    c[0] += 1
                    if mine[key] != r1: c[2].append((key, r1, mine[key]))
                    else: c[1] += 1
                    break
for k in sorted(cmp):
    tot, agree, bad = cmp[k]
    print(k, 'compared', tot, 'agree', agree, 'disagree', bad)
