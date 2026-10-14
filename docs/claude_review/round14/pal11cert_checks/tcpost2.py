"""EVIDENCE: classify the 200 numerical minimisers of tcnum.py (seed 1).  For each: F, support, whether the support lies in a
certified coordinate subspace of TC (tcpost.py), the number of distinct positions among the 12 points (0, u_1..u_11)
at tolerance 1e-3 (collision stratum), and the minimal gap sep(u) between the 12 positions (normalised |u|_inf = 1)."""
import pickle, numpy as np
from collections import Counter
zeros = pickle.load(open('tcnum_zeros_1.pkl', 'rb'))
maximal = pickle.load(open('tc_coord_maximal.pkl', 'rb'))
rows = []
for F, u in zeros:
    u = u / np.abs(u).max()
    S = tuple(np.nonzero(np.abs(u) > 1e-7)[0].tolist())
    inside = any(set(S) <= set(T) for T in maximal)
    pos = np.sort(np.concatenate([[0.0], u]))
    gaps = np.diff(pos)
    ncl = 1 + int((gaps > 1e-3).sum())
    rows.append((inside, len(S), ncl, gaps.min(), F))
print('total', len(rows))
print('certified coordinate-subspace zeros:', sum(1 for r in rows if r[0]), ' F max among them %.1e' % max(r[4] for r in rows if r[0]))
rest = [r for r in rows if not r[0]]
print('others:', len(rest))
print('  number of position clusters (tol 1e-3) -> count:', sorted(Counter(r[2] for r in rest).items()))
print('  max over others of min gap sep(u): %.2e' % max(r[3] for r in rest))
print('  F range among others: %.1e .. %.1e' % (min(r[4] for r in rest), max(r[4] for r in rest)))
print('  largest cluster count found:', max(r[2] for r in rest), ' (12 = all positions distinct at 1e-3)')
