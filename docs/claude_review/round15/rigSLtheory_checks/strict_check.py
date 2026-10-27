"""Strictness of PU / SL(phi_0) survivors at q0 (one-modulus): a level-1 label ell is (eps,t)-strict if
ell = u (eps*L + (m/2) t) mod m for some u = 1 mod 4 unit, eps in {+-1}^r, t in {0,1}^r.  Counts, among survivors in
the eps-class (mod n), how many are (eps,t)-strict, for PU survivors and for SL-at-phi_0 survivors.
Usage: python3 strict_check.py dumpPU q0 p1..pr"""
import sys, json, itertools, math
from sllib import genuine
dump, q0 = sys.argv[1], int(sys.argv[2]); ps = [int(x) for x in sys.argv[3:]]; r = len(ps)
m, L = genuine(ps, q0); n = m//4
U1 = [u for u in range(1, m) if math.gcd(u, m) == 1 and u % 4 == 1]
strict = set(); strict0 = set(); strict00 = set()
for eps in itertools.product([1, -1], repeat=r):
    for t in itertools.product([0, 1], repeat=r):
        base = [(e*x + (m//2)*tt) % m for e, x, tt in zip(eps, L, t)]
        for u in U1:
            lab = tuple((u*b) % m for b in base); strict.add(lab)
            if not any(t): strict0.add(lab)
            if not any(t) and all(e == 1 for e in eps): strict00.add(lab)
units = [u for u in range(1, n+1) if math.gcd(u, n) == 1] if n > 1 else [1]
epsc = set(tuple((u*e*x) % n for e, x in zip(eps, L)) for u in units for eps in itertools.product([1, -1], repeat=r))
cnt = {'PU_in': 0, 'PU_in_strict': 0, 'PU_in_strict_t0': 0, 'PU_in_strict_gen': 0, 'SL_in': 0, 'SL_in_strict': 0, 'SL_in_strict_t0': 0, 'SL_in_strict_gen': 0}
for ell, slg in json.load(open(dump)):
    ell = tuple(ell)
    if tuple(x % n for x in ell) not in epsc: continue
    st = ell in strict
    s0 = ell in strict0; s00 = ell in strict00
    cnt['PU_in'] += 1; cnt['PU_in_strict'] += st; cnt['PU_in_strict_t0'] += s0; cnt['PU_in_strict_gen'] += s00
    if slg: cnt['SL_in'] += 1; cnt['SL_in_strict'] += st; cnt['SL_in_strict_t0'] += s0; cnt['SL_in_strict_gen'] += s00
print(f"ps={ps} q0={q0} m={m}: {cnt}")
