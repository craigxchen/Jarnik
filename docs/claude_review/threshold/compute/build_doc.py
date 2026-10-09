"""Assemble ../compute.md from template.md, substituting generated tables:
   {{TABLE:kind:section}}  -> output of make_tables.py kind section
   {{CMD:...}}             -> output of a python3 command run in this directory
Also writes the appendix tables ../compute_tables_shell.md and ../compute_tables_slice.md."""
import subprocess, re, os, sys
here = os.path.dirname(os.path.abspath(__file__))


def run(args):
    return subprocess.run(['python3'] + args, cwd=here, capture_output=True, text=True).stdout.rstrip()


tpl = open(os.path.join(here, 'template.md')).read()


def sub(m):
    kind, section = m.group(1), m.group(2)
    return run(['make_tables.py', kind, section])


import json as _json
secs = {'METHOD_SECTION': 'draft_method.md', 'TREND_SECTION': 'draft_trend.md'}
for name in ('MAXIMA_NOTES', 'VAND_NOTES', 'TREND_NOTES', 'COMPONENT_NOTES', 'CERT_SECTION',
             'SLICE_SECTION', 'CONSEQ_SECTION', 'FILES_SECTION'):
    secs[name] = 'sec_' + name.lower() + '.md'
for name, fn in secs.items():
    path = os.path.join(here, fn)
    txt = open(path).read().rstrip() if os.path.exists(path) else '(' + name + ' pending)'
    tpl = tpl.replace(name, txt)
vars_path = os.path.join(here, 'doc_vars.json')
if os.path.exists(vars_path):
    for k, v in _json.load(open(vars_path)).items():
        tpl = tpl.replace(k, v)
out = re.sub(r'\{\{TABLE:(\w+):(\w+)\}\}', sub, tpl)
out = re.sub(r'\{\{CMD:([^}]*)\}\}', lambda m: run(m.group(1).split()), out)
open(os.path.join(here, '..', 'compute.md'), 'w').write(out)
for kind in ('shell', 'slice'):
    with open(os.path.join(here, '..', f'compute_tables_{kind}.md'), 'w') as f:
        f.write(f'# Full {kind} tables (appendix to compute.md)\n\n')
        f.write('Columns: s, t, (p,n), number of points, reg (= max ord) mod P1, reg/t, the '
                'components attaining it, floor(phi_M(p,n)) (upper.md bound), whether the two '
                'primes give identical per-component profiles, and reg_lam for every component '
                'with N_lam > 0 (lam in exponent notation, e.g. 21^3 = (2,1,1,1)).\n\n')
        f.write('## Trend tables (max over s for each t)\n\n')
        f.write(run(['make_tables.py', kind, 'trend']) + '\n\n')
        f.write('## Every (M, s, t)\n\n')
        f.write(run(['make_tables.py', kind, 'full']) + '\n')
print('written')
