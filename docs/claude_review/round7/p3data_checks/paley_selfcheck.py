"""Self-consistency of the round-5 Paley profiles with their OWN actual primes (part 1 supplement).

For each trial: consecutive primes = 1 mod 4 from a random start in [1e7, 1e9] (as genuine_stats.py),
random conjugates, genuine unit-1 moduli m_xy at X = 2M and 4M (exact, demands7), and the profile's
exact pair margins mu_xy = sum_j log p_j ([a_xj != a_yj] - 1/2).  With C = sqrt 2 the pair inequality is
prod_sep p >= m^4 prod_same p (checked in floating point with a 1e-9 relative guard; the margin of every
reported decision exceeds 1e-6).  Reports: fraction of trials in which all pair inequalities hold,
and tau_req = max_xy D_xy / nu_xy where nu_xy = mu_xy / tau_bar is the normalised margin
(tau_bar = mean log p_j): the least common weight scale making this profile pair-feasible for these
residues; W_req = r * tau_req.  Evidence (exact integers for D, floating logs for margins).
"""
import json, math, random, sys
import numpy as np
from gauss7 import *
from demands7 import *

def run(M, variant, ntr):
    q0 = M - 1
    b = 8 if variant == 'construct' else 33
    S = profile_construct(q0, b) if variant == 'construct' else profile_lower(q0, b)
    A = ((1 + S) // 2).astype(np.int64)
    r = S.shape[1]
    out = []
    for t in range(ntr):
        rng = random.Random(555 + 1000 * M + t + (7 if variant == 'lower' else 0))
        start = rng.randrange(10 ** 7, 10 ** 9)
        ps = primes_1mod4_window(start, r)
        gA, gB = gaussian_primes(ps, rng)
        w = np.log(np.array(ps, dtype=float))
        sep = (A[:, None, :] != A[None, :, :])
        mu = (sep * w).sum(2) - w.sum() / 2
        tau = w.mean()
        nu = mu / tau
        classes = point_classes(A, gA, gB, ps, 4 * M)
        rec = dict(M=M, variant=variant, t=t, r=r, tau=tau)
        for Xf in (2, 4):
            m1, _ = collision_moduli(classes, Xf * M, M)
            D = demand_matrix(m1)
            iu = np.triu_indices(M, 1)
            slack = (mu - D)[iu]
            rec['ok_X%d' % Xf] = bool(slack.min() > 0)
            rec['minslack_X%d' % Xf] = float(slack.min())
            rec['tau_req_X%d' % Xf] = float((D[iu] / nu[iu]).max())
            rec['nviol_X%d' % Xf] = int((slack < 0).sum())
        out.append(rec)
    return out

if __name__ == '__main__':
    res = []
    for M in (12, 20, 24, 32, 44, 60, 72, 84):
        for variant in ('construct', 'lower'):
            recs = run(M, variant, 60 if variant == 'construct' else 20)
            res += recs
            r = recs[0]['r']
            line = 'M=%2d %-9s r=%4d tau~%.2f |' % (M, variant, r, np.mean([x['tau'] for x in recs]))
            for Xf in (2, 4):
                ok = np.mean([x['ok_X%d' % Xf] for x in recs])
                tr = np.median([x['tau_req_X%d' % Xf] for x in recs])
                nv = np.mean([x['nviol_X%d' % Xf] for x in recs])
                line += ' X=%dM: P(feasible)=%.2f tau_req med %.2f (W_req med %.0f) mean#viol %.2f |' % (Xf, ok, tr, r * tr, nv)
            print(line); sys.stdout.flush()
    json.dump(res, open('paley_selfcheck.json', 'w'))
