import random, itertools
from lamA_direct import *
random.seed(3)
pR, secs = lam_data(good_mob(()))
print('degrees: R', max(c.degree() for c in pR), ' sections', [max(c.degree() for c in s) for s in secs])
M = [[meet(secs[a], secs[b]) if a != b else None for b in range(6)] for a in range(6)]
for row in M: print(row)
print('Sigma.s_i:', [meet(pR, s) for s in secs])
# additivity test: M_ab + M_cd == M_ac + M_bd
bad = 0
for a, b, c, d in itertools.permutations(range(6), 4):
    if M[a][b] + M[c][d] != M[a][c] + M[b][d]: bad += 1
print('additivity violations:', bad)
