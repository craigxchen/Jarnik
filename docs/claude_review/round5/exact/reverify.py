"""Independent exact re-verification of every reported optimal solution (results/*.json) with
verify.py, using the reported windings kappa: checks sigma_lambda >= r(lambda) (> if strict) for all
relevant lambda (integral and rational classes with lambda.k in Z), exact arithmetic.
Also checks the Plotkin lower bound and the rank bound against the value (exactness flag)."""
import json, glob, sys
from fractions import Fraction
from core import pattern_of
from verify import verify

rows = []
for f in sorted(glob.glob('results/c_*.json')) + sorted(glob.glob('results/b_*.json')):
    d = json.load(open(f))
    if d.get('status') not in ('optimal', 'feasible_time_limit'):
        continue
    M = d['M']; rho0 = Fraction(d['rho0']); rho1 = Fraction(d['rho1']); strict = d['strict']
    T = [pattern_of(M, r) for r in d['support']]
    w = [Fraction(x) for x in d['weights']]
    r = verify(M, T, w, d['kappa'], rho0, rho1, strict=strict)
    Om = sum(w)
    pl = (2 * (M - 1) if M % 2 == 0 else 2 * M) * (rho0 + rho1)
    lb = max(Fraction(M - 1), pl)
    exact = (Om == lb)
    print(f'{f:32s} M={M} rho0={rho0} rho1={rho1} strict={strict}: Omega={Om} verify_ok={r["ok"]} '
          f'min_slack={r["min_slack"]} Plotkin/rank bound={lb} -> value {"EXACT (meets proven bound)" if exact else "solver-optimal"}',
          flush=True)
