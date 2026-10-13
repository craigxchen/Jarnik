"""Lemma 6.3 of p3data.md (Perron attainment) on genuine demands: how close W_cut is to 2 lambda_max(D),
and, at the LP optimum, how u-balanced the support cuts are and how tight the margins are.

For each instance: ratio W_cut/(2 lambda_max), the certified dual bound, sum_s lam_s <s,u>^2 (which
equals u^T G u = W - 2 sum_{x != y} u_x u_y mu_xy, Lemma 6.3) and the u-weighted slack
2 sum_{x<y} 2 u_x u_y (mu_xy - D_xy); by Lemma 6.3 their sum is W - 2 lambda_max exactly.
Usage: python3 perron_check.py
"""
import math, json
import numpy as np
from qsolve import lp_full, lp_colgen, surplus
from qcore import pair_index, mask_to_s

out = []
for M, trials in ((8, range(6)), (12, range(6)), (20, range(4)), (24, range(2))):
    iu = pair_index(M)
    for t in trials:
        for Xf in (2, 4):
            if M == 24 and Xf == 4:
                continue
            D = np.load('data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
            ev, V = np.linalg.eigh(D)
            lam = float(ev[-1]); u = np.abs(V[:, -1])
            if M <= 12:
                W, x, masks, nu = lp_full(D)
                lb = float((nu[iu] * D[iu]).sum() / surplus(nu))
            else:
                W, x, masks, nu, lb = lp_colgen(D)
            bal = 0.0
            mu = np.zeros((M, M))
            for m, l in zip(masks, x):
                if l > 1e-12:
                    s = mask_to_s(m, M).astype(float)
                    bal += l * float(s @ u) ** 2
                    mu += l * ((s[:, None] != s[None, :]) - 0.5)
            slack = float(sum(4 * u[a] * u[b] * (mu[a, b] - D[a, b]) for a, b in zip(*iu)))
            rec = dict(M=M, t=t, X=Xf * M, W=W, lb=lb, lam2=2 * lam, ratio=W / (2 * lam), bal=bal, slack=slack,
                       identity_err=W - 2 * lam - bal - slack, min_slack=float((mu - D)[iu].min()))
            out.append(rec)
            print('M=%2d t=%d X=%3d  W_cut/2lam = %.5f  (lb/2lam = %.5f)  sum lam<s,u>^2 = %.4f  u-slack = %.4f  identity err %.1e' % (
                M, t, Xf * M, rec['ratio'], lb / (2 * lam), bal, slack, rec['identity_err']), flush=True)
json.dump(out, open('perron_check.json', 'w'), indent=1)
