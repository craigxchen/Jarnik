"""Theorem 3.3 (support-uncertainty budget) evaluated on exact models and on regime Sylvester sizes.
For a profile containing a Hadamard core (n-1 sign columns which with the all-ones vector form an n x n Hadamard matrix):
a product prod_{i<=m} <z_i,u> in L_m (width system) forces  m^2 <= (2d + tau_0) * T_core / n,
T_core = sum over core columns of tau_A(1_P), tau_0 = tau_A(base direction) (only min(tau_0, m) matters).
"""
import numpy as np, math
from fplib import *

def hadamard_core(M, pats):
    # find n-1 patterns whose +-1 vectors with all-ones form a Hadamard matrix (greedy maximal orthogonal set)
    vecs = []
    for P in pats:
        v = np.ones(M, dtype=int); v[list(P)] = -1
        vecs.append((P, v))
    core = []
    for P, v in vecs:
        if v.sum() != 0: continue
        if all(int(v @ w) == 0 for _, w in core):
            core.append((P, v))
    return core

with open('uncert_budget_out.txt', 'w') as out:
    for name, L in [('sdc7', 2), ('pal7', 2), ('syl8', 2), ('syl8', 3), ('sdc11', 2), ('pal11', 2), ('paley12', 2)]:
        M, pats = MODELS[name](); A = np.array(box(M, pats, L)); d = M - 1
        core = hadamard_core(M, pats)
        Tcore = 0
        for P, v in core:
            vals = A[:, [x - 1 for x in P]].sum(axis=1)
            Tcore += len(np.unique(vals)) - 1
        tau0 = len(np.unique(-A.sum(axis=1))) - 1
        bound = lambda m: m * m <= (2 * d + min(tau0, m)) * Tcore / M
        mmax = max(m for m in range(1, 60) if bound(m))
        line = f'{name} L={L}: n={M} core columns={len(core)} (Hadamard core needs {M-1}) T_core={Tcore} tau_0={tau0} -> certificates in L_m impossible for m > {mmax}'
        print(line); out.write(line + '\n')
    # regime Sylvester: m'^2 <= (2d + m') * 2j(n-1)/n
    for k in range(12, 21):
        n = 2 ** k; d = n - 1; r = d
        w = math.sqrt(8 * d * math.log(4 * r)); j = math.floor(w / 2); mp = math.ceil(4 * w)
        lhs = mp * mp; rhs = (2 * d + mp) * 2 * j * (n - 1) / n
        line = f'Sylvester n=2^{k}: m\'={mp} j={j}: m\'^2={lhs}  (2d+m\')*2j(n-1)/n={rhs:.0f}  -> width certificates {"IMPOSSIBLE" if lhs > rhs else "not excluded"}'
        print(line); out.write(line + '\n')
