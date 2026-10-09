"""Test the summed good-pair inequality  E[W(eps) ((sum eps)^2 - M)] <= 0  for widths W of finite sets,
and find the worst ratio E[W (sum eps)^2] / (M E[W])."""
import random, itertools, sys
import numpy as np
from fractions import Fraction as Fr
random.seed(int(sys.argv[3]))
M = int(sys.argv[1]); trials = int(sys.argv[2])
eps = np.array(list(itertools.product([1, -1], repeat=M)), dtype=np.int64)
s2 = eps.sum(axis=1)**2
worst = None; viol = 0
for t in range(trials):
    n = random.choice([2, 2, 3, 4, 6, 10]); R = random.choice([1, 2, 3, 5, 10, 30])
    S = []
    for _ in range(n):
        k = [random.randint(-R, R) for _ in range(M-1)]
        S.append(k + [-sum(k)])
    A = np.array(S, dtype=np.int64); V = eps @ A.T
    W = V.max(axis=1) - V.min(axis=1)
    if W.sum() == 0: continue
    val = Fr(int((W*s2).sum()), M*int(W.sum()))
    if worst is None or val > worst[0]: worst = (val, S)
    if val > 1: viol += 1
print(f"M={M}: violations of sum-GPL: {viol}/{trials}; max E[W s^2]/(M E W) = {worst[0]} = {float(worst[0]):.4f} at S={worst[1]}")
