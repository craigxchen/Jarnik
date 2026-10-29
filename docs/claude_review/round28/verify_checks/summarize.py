"""Collect my reach (Q1, Q2) and zeta results, compare with the primary's boxes.csv; print a markdown table."""
import json, csv
base = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/'
P = {}
for p in (2147483647, 1073741827):
    for l in open(base + 'verify_checks/reach_%d.jsonl' % p):
        r = json.loads(l)
        if r['shifts'] != 'pair' or r['emin'] != 0: continue
        P.setdefault(r['spec'], {})[p] = r
Z = {}
for l in open(base + 'verify_checks/zeta.jsonl'):
    r = json.loads(l); Z[r['spec']] = r
C = {r['spec']: r for r in csv.DictReader(open(base + 'boxes.csv'))}
print('| box | d | L | `|A|` | `rho` | reach `Q1` / `Q2` (fail dims) | `zeta_pair` / `zeta_all` | reach − zeta | primary (reach; zeta_all/zeta_pair; fail) | agree |')
print('| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |')
for spec, d in P.items():
    if 1 not in [1] or len(d) < 2: continue
    a, b = d[2147483647], d[1073741827]
    def fd(r):
        e = str(r['reach'])
        return '%s/%s' % tuple(r['dims'][e]) if e in r['dims'] else 'none below rho'
    z = Z.get(spec, {})
    c = C.get(spec)
    agree = ''
    prim = 'not in primary set'
    if c:
        prim = '%s; %s/%s; %s' % (c['reach_pair_P1'], c['zeta_all'], c['zeta_pair'], c['fail_pair'])
        ok = (str(a['reach']) == c['reach_pair_P1'] == str(b['reach']) and str(a['size']) == c['size']
              and ' '.join(map(str, a['H'])) == c['H'] and a['H'] == b['H']
              and (not z or (str(z.get('zeta_all')) == c['zeta_all'] and str(z.get('zeta_pair')) == c['zeta_pair'])))
        if c['fail_pair'] and '/' in c['fail_pair'] and str(a['reach']) in a['dims']:
            ok = ok and fd(a) == c['fail_pair'] == fd(b)
        agree = 'yes' if ok else '**NO**'
    za = z.get('zeta_all', '?')
    diff = (a['reach'] - za) if isinstance(za, int) else '?'
    print('| %s | %d | %d | %d | %d | %s / %s (%s) | %s / %s | %s | %s | %s |' % (
        spec, a['n'] - 1, a['L'], a['size'], a['rho'], a['reach'], b['reach'], fd(a), z.get('zeta_pair', '?'), za,
        diff, prim, agree))
