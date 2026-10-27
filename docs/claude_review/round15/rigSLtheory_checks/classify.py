"""Round 15 (rigSLtheory): one-modulus test at q0, with every PU survivor's extra collisions classified.

Setting (as round 14 J2): genuine labels (all levels) on F \\ {q0}; free moduli collision-free; at q0 an arbitrary
parity-respecting level-1 label ell (level-1 relaxation).  For v in V_K and unit s, g_s(v) = F\\{q0}-part of D_s(Z_v)
(all levels).  Under ell, v collides at q0 with unit s' iff <v,ell> = s' n mod m; then m_{s'} = q0 g_{s'}, else m = g.

Per-v classes of a colliding (v, s') (exact integer / rational arithmetic):
  GEN  : the genuine label also collides v at q0 with the same unit s'  (no change of demand)
  C    : q0 g_{s'} <= |D_{s'}(Z_v)|                     (cofactor type: fits at phi = arg pi_j)
  else : PU (per-unit budget) is assumed (survivors are PU-good);
         ANGok  : J0 box meets circle                    (type A1: archimedean, opposite unit not saturated)
         OPP    : J0 fails via opposite pair  a_s + a_{s+2} > 1
         ADJ    : J0 fails via adjacent condition (min |z|^2 > 1), opposite pairs fine
         MAX    : J0 fails via max |z|^2 < 1 only
Here a_s = m_s^2 / (nu_s^-2 N), nu_s^-2 = 1 (s even), 2 (s odd).
Output per survivor outside every eps-class: counts of extras per class; whether ANG-good, SL at genuine phi.
Usage: python3 classify.py FFILE_or_json q0 K p1..pr
"""
import sys, os, math, itertools, json
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sllib import two_sq, genuine, nu, Zv, Ds, VK, vq


def load_F(path):
    if path.endswith('.json'):
        return json.load(open(path))['forced']
    for line in open(path):
        if line.startswith('F '):
            return eval(line[2:])


def j0_reason(ms, N):
    """exact J0 test; returns 'ok' or the failing condition"""
    a = [Fr(ms[s]**2, N * (2 if s % 2 else 1)) for s in range(4)]
    if a[0] + a[2] > 1 or a[1] + a[3] > 1:
        return 'OPP'
    dC = max(Fr(0), 2*a[0]-1, 2*a[2]-1); dS = max(Fr(0), 2*a[1]-1, 2*a[3]-1)
    if dC*dC + dS*dS > 1:
        return 'ADJ'
    MC = max(abs(2*a[2]-1), abs(1-2*a[0])); MS = max(abs(2*a[3]-1), abs(1-2*a[1]))
    if MC*MC + MS*MS < 1:
        return 'MAX'
    return 'ok'


def main():
    Fpath, q0, K = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    ps = [int(x) for x in sys.argv[4:]]
    F = load_F(Fpath); assert q0 in F
    r = len(ps); pis = [two_sq(p) for p in ps]; V = VK(r, K)
    m = q0-1 if q0 % 4 == 1 else q0+1; n = m//4
    Mgen, L = genuine(ps, q0); assert Mgen == m
    data = []  # per v: N, D[4], g[4], genuine unit at q0 (or None), PU-bad flags per s
    for v in V:
        z = Zv(pis, v); N = z[0]**2 + z[1]**2
        D = [Ds(z, s) for s in range(4)]
        g = []
        for s in range(4):
            gg = 1
            for q in F:
                if q != q0: gg *= q**vq(D[s], q)
            g.append(gg)
        dot = sum(a*b for a, b in zip(v, L)) % m
        gu = dot // n if dot % n == 0 else None
        pubad = []
        for s in range(4):
            ms = list(g); ms[s] *= q0
            pubad.append(not all(ms[t]**2 <= N*(2 if t % 2 else 1) for t in range(4)))
        data.append((v, N, D, g, gu, pubad))
    # DFS enumeration of PU survivors: columns in order, check v supported on first k columns
    par = nu(ps, q0)
    order = sorted(range(len(V)), key=lambda i: max([j for j in range(r) if V[i][j] != 0]))
    lastcol = [max(j for j in range(r) if V[i][j] != 0) for i in range(len(V))]
    cur = np.zeros((1, 0), dtype=np.int64)
    for k in range(r):
        vals = np.arange(par[k], m, 2, dtype=np.int64)
        cur = np.concatenate([np.repeat(cur, len(vals), axis=0), np.tile(vals, len(cur))[:, None]], axis=1)
        keep = np.ones(len(cur), dtype=bool)
        for i in range(len(V)):
            if lastcol[i] != k: continue
            v = np.array(V[i][:k+1], dtype=np.int64)
            d = (cur @ v) % m; hit = (d % n == 0); s = (d // n) % 4
            bad = np.array(data[i][5]); keep &= ~(hit & bad[s])
        cur = cur[keep]
        print(f"  col {k}: {len(cur)} partial survivors", file=sys.stderr, flush=True)
    units = [u for u in range(1, n+1) if math.gcd(u, n) == 1] if n > 1 else [1]
    gen = set(tuple((u*x) % n for x in L) for u in units)
    epsc = set(tuple((u*e*x) % n for e, x in zip(eps, L)) for u in units
               for eps in itertools.product([1, -1], repeat=r))
    tot = {'all': len(cur), 'out_gen': 0, 'out_eps': 0, 'ANG': 0, 'ANG_out_eps': 0, 'SLg_out_eps': 0}
    reasons = {}
    profile = {}
    examples = []
    angA1 = []
    slgw = []
    slall = []
    pu_all = [] if os.environ.get('DUMPPU') else None
    for row in cur:
        c = tuple(int(x) % n for x in row)
        oe = c not in epsc
        tot['out_gen'] += c not in gen; tot['out_eps'] += oe
        cnt = {'GEN': 0, 'C': 0, 'A1': 0, 'OPP': 0, 'ADJ': 0, 'MAX': 0}
        oppU = []
        for (v, N, D, g, gu, _) in data:
            d = sum(int(a)*b for a, b in zip(row, v)) % m
            if d % n: continue
            s = d // n
            if gu == s: cnt['GEN'] += 1; continue
            ms = list(g); ms[s] *= q0
            if ms[s] <= abs(D[s]): cnt['C'] += 1; continue
            rs = j0_reason(ms, N)
            if rs == 'ok': cnt['A1'] += 1
            else:
                cnt[rs] += 1
                if rs == 'OPP':
                    so = (s+2) % 4
                    oppU.append(abs(D[so])//g[so])
        angok = cnt['OPP'] + cnt['ADJ'] + cnt['MAX'] == 0
        slg = angok and cnt['A1'] == 0
        tot['ANG'] += angok
        if slg: slall.append(tuple(int(x) for x in row))
        if pu_all is not None: pu_all.append((tuple(int(x) for x in row), slg))
        if oe:
            tot['ANG_out_eps'] += angok; tot['SLg_out_eps'] += slg
            key = tuple(sorted((k2, v2) for k2, v2 in cnt.items() if v2 and k2 != 'GEN'))
            profile[key] = profile.get(key, 0) + 1
            for u in oppU: reasons[u] = reasons.get(u, 0) + 1
            if len(examples) < 3: examples.append((tuple(int(x) for x in row), cnt))
            if angok and not slg: angA1.append(tuple(int(x) for x in row))
            if slg and len(slgw) < 20: slgw.append(tuple(int(x) for x in row))
    print(f"ps={ps} K={K} q0={q0} n={n} |F|={len(F)} genuine L={L}")
    print(f"  totals: {tot}")
    print(f"  extra-type profiles of survivors outside every eps-class (count):")
    for k2, v2 in sorted(profile.items(), key=lambda t: -t[1]):
        print(f"    {dict(k2)}: {v2}")
    print(f"  opposite-unit cofactor U_(s'+2) on OPP-failing extras: {dict(sorted(reasons.items()))}")
    print(f"  examples: {examples}")
    print(f"  ANG-good, not SL at genuine phi, outside every eps-class: {len(angA1)}")
    if os.environ.get('DUMP'): json.dump(angA1, open(os.environ['DUMP'], 'w'))
    if os.environ.get('DUMPPU'): json.dump(pu_all, open(os.environ['DUMPPU'], 'w'))
    if os.environ.get('DUMPSL'): json.dump(slgw, open(os.environ['DUMPSL'], 'w'))


if __name__ == '__main__':
    main()
