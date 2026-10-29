"""Merge results_<P1>.jsonl, results_<P2>.jsonl, zeta.jsonl into ../boxes.csv and print a summary."""
import json, csv, sys
from fractions import Fraction
from collections import OrderedDict
B = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/'
C = B + 'compute_checks/'
P1, P2 = 2147483629, 2147483587
def load(fn):
    out = OrderedDict()
    try:
        for l in open(fn):
            r = json.loads(l); out[r['spec']] = r
    except FileNotFoundError:
        pass
    return out
r1, r2, z = load(C + 'results_%d.jsonl' % P1), load(C + 'results_%d.jsonl' % P2), load(C + 'zeta.jsonl')
def family(spec):
    f = spec.split(':')[0]
    return {'syl': 'Sylvester', 'sdc': 'sdc', 'had': 'Paley-Hadamard', 'rsyl': 'Sylvester16|J', 'rsdc': 'sdc11|J',
            'rhad': 'PaleyHad12|J', 'rand': 'random', 'randL': 'random(L fixed)'}[f]
rows = []
specs = list(OrderedDict.fromkeys(list(r1) + list(r2)))
for s in specs:
    a, b, zz = r1.get(s), r2.get(s), z.get(s)
    base = a or b
    row = OrderedDict(spec=s, family=family(s), d=base['d'], n=base['n'], cuts=base['r'], L=base['L'],
                      size=base['size'], rho=base['rho'], rho_over_L=round(base['rho'] / base['L'], 3),
                      H=' '.join(map(str, base['H'])))
    row['H_agree'] = (a is not None and b is not None and a['H'] == b['H'])
    row['reach_pair_P1'] = a['reach_pair'] if a else ''
    row['reach_pair_P2'] = b['reach_pair'] if b else ''
    row['fail_pair'] = '%d/%d' % tuple(a['fail_pair']) if a and a['fail_pair'] else ''
    row['reach_cube_P1'] = a.get('reach_cube', '') if a else ''
    row['reach_cube_P2'] = b.get('reach_cube', '') if b else ''
    row['cube_complete'] = (a.get('cube_complete', '') if a else '')
    rp = a['reach_pair'] if a else b['reach_pair']
    rc = a.get('reach_cube') if a else None
    if zz:
        za, zp = zz.get('all', {}), zz.get('pair', {})
        row['zeta_all'] = za.get('witness_order', ''); row['zeta_all_LP'] = za.get('lp', '')
        row['zeta_all_exact'] = za.get('exact', '')
        row['zeta_all_milp_bound'] = za.get('milp_dual_bound', '')
        row['zeta_pair'] = zp.get('witness_order', ''); row['zeta_pair_LP'] = zp.get('lp', '')
        row['zeta_pair_exact'] = zp.get('exact', '')
        zeta = za.get('witness_order')
    else:
        zeta = None
    row['reach_minus_zeta_pair'] = (rp - zeta) if zeta is not None else ''
    row['reach_cube_minus_zeta'] = (rc - zeta) if (zeta is not None and rc is not None) else ''
    row['reach_over_rho'] = round(rp / base['rho'], 3)
    lo = max([v for v in (rp, rc, zeta) if v is not None])
    hi = base['rho']
    if zeta is not None and zeta < base['rho'] and base['H'][base['rho']] == 1:
        hi = base['rho'] - 1          # single top form, not a binomial product (lemma in compute.md §1)
    row['truth_lo'] = lo; row['truth_hi'] = hi
    row['H_top'] = base['H'][base['rho']]
    rows.append(row)
keys = list(rows[0].keys()) if rows else []
allk = []
for r in rows:
    for k in r:
        if k not in allk: allk.append(k)
with open(B + 'boxes.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=allk); w.writeheader()
    for r in rows: w.writerow(r)
for r in rows:
    print("%-34s d=%d |A|=%5d L=%d rho=%2d rho/L=%5.2f reach_pair=%s/%s cube=%s/%s zeta=%s/%s(%s,%s) r-z=%s" % (
        r['spec'][:34], r['d'], r['size'], r['L'], r['rho'], r['rho_over_L'], r['reach_pair_P1'], r['reach_pair_P2'],
        r['reach_cube_P1'], r['reach_cube_P2'], r.get('zeta_all', ''), r.get('zeta_pair', ''),
        r.get('zeta_all_LP', ''), 'ex' if r.get('zeta_all_exact') else 'milp', r['reach_minus_zeta_pair']))
