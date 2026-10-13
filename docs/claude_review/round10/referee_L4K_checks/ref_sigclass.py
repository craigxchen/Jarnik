"""Of the parity-respecting vectors with the genuine (v,s)-signature (K=2), how many are in the genuine class
c = u L mod n (u unit mod n)?"""
import itertools, math
import numpy as np
from rlib import *
ps = [5, 13, 17, 29]; K = 2; r = 4
V = VK_mod_sign(r, K); Varr = np.array(V)
pis = gaussian_primes(ps)
for q in [3, 7, 11, 19, 23, 31, 41]:
    m = q - chi4(q); n = m // 4; nv = nu(ps, q)
    T = Torus(q); io = T.iota()
    L = [T.ratio_log(pi) for pi in pis]
    if io != m // 4: L = [(-x) % m for x in L]
    Lg = np.array(L)
    choices = [[x for x in range(m) if x % 2 == nv[j]] for j in range(r)]
    allv = np.array(list(itertools.product(*choices)), dtype=np.int64)
    t = (allv @ Varr.T) % m
    tg = (Varr @ Lg) % m
    same = np.all(t == tg[None, :], axis=1) | False
    # signature equality: collision pattern and unit at colliding v equal
    cg = (tg % n == 0); sg = (tg // n) % 4
    c = (t % n == 0); s = (t // n) % 4
    sig_eq = np.all(c == cg[None, :], axis=1) & np.all(~cg[None, :] | (s == sg[None, :]), axis=1)
    units = [u for u in range(1, n + 1) if math.gcd(u, n) == 1] if n > 1 else [0]
    cls = set(tuple((u * x) % n for x in L) for u in units) if n > 1 else {tuple([0] * r)}
    incls = np.array([tuple(row % n) in cls for row in allv]) if n > 1 else np.ones(len(allv), bool)
    print(f"q={q} n={n} genuine-signature vectors={int(sig_eq.sum())} of which in genuine class={int((sig_eq & incls).sum())}; class size={int(incls.sum())}")
