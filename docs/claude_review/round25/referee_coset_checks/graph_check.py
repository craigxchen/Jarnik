"""Referee: independent check of Cor 4.2 (step sums in graph coordinates) and of Thm 5.1 (total step mass lower bound).
Coset g + H_0, H_0 = {(theta, C theta)} : phi_J = theta in T^m, phi_{J^c} = C theta (normalised coords x = phi/(pi/2)).
mu_v = v_J + C^T v_{Jc}. last(mu) = largest index k with mu_k != 0."""
import itertools, math, json
import numpy as np
from fractions import Fraction
pJ  = [17, 29, 37]          # free primes j_1, j_2, j_3 (order = coordinate order)
pD  = [13, 41]              # dependent primes
Cm  = np.array([[1, -2, 1], [0, 3, -1]])   # rows c_i (i in J^c), entries over k = 1..3 ; row 2 last entry at k=3, row1 at k=3
m = 3
def y(p, n): return p ** (-abs(n) / 2)
Pi = lambda p: 1 + 2 / (math.sqrt(p) - 1)
# brute force over a box
Bx = 26
res_bf = [0.0] * m
rng = range(-Bx, Bx + 1)
for x in itertools.product(range(-12, 13), repeat=2):
    yx = y(pD[0], x[0]) * y(pD[1], x[1])
    if yx < 1e-14: continue
    f = Cm.T @ np.array(x)
    for vJ in itertools.product(rng, repeat=3):
        yv = yx * y(pJ[0], vJ[0]) * y(pJ[1], vJ[1]) * y(pJ[2], vJ[2])
        if yv < 1e-16: continue
        mu = np.array(vJ) + f
        nz = np.nonzero(mu)[0]
        if len(nz) == 0: continue
        res_bf[nz[-1]] += yv
# formula
res_f = [0.0] * m
for x in itertools.product(range(-40, 41), repeat=2):
    yx = y(pD[0], x[0]) * y(pD[1], x[1])
    f = Cm.T @ np.array(x)
    for k in range(m):
        E = 1.0
        for kp in range(k + 1, m): E *= pJ[kp] ** (-abs(f[kp]) / 2)
        Pl = 1.0
        for kp in range(k): Pl *= Pi(pJ[kp])
        res_f[k] += Pl * yx * E * (Pi(pJ[k]) - pJ[k] ** (-abs(f[k]) / 2))
rel = [abs(a - b) / b for a, b in zip(res_bf, res_f)]
print('brute', res_bf); print('formula', res_f); print('rel err', rel)
# Thm 5.1 check: true U(mu) (measure of union of windows over the fibre) for a random g, truncated fibres (lower bound)
rngg = np.random.default_rng(5)
g = rngg.uniform(0, 1, 5)    # (g_J, g_Jc)
def halfw(yv): return (2 / math.pi) * math.asin(yv), (2 / math.pi) * math.asin(yv / math.sqrt(2))
fib = {}
for x in itertools.product(range(-3, 4), repeat=2):
    f = Cm.T @ np.array(x)
    for vJ in itertools.product(range(-6, 7), repeat=3):
        v = list(vJ) + list(x)
        if not any(v): continue
        w = sum(abs(vi) * math.log(p) for vi, p in zip(v, pJ + pD))
        if w > 16: continue
        mu = tuple(int(t) for t in (np.array(vJ) + f))
        if not any(mu): continue
        nz = [i for i in range(m) if mu[i]][-1]
        if mu[nz] < 0: continue        # keep F_k representative
        c = sum(vi * gi for vi, gi in zip(v, list(g)))   # <v,g>
        a0, a1 = halfw(math.exp(-w / 2))
        fib.setdefault(mu, []).append(((-c) % 1, a0)); fib[mu].append(((0.5 - c) % 1, a1))
def union_len(iv):
    segs = []
    for cen, a in iv:
        lo, hi = cen - a, cen + a
        for s in (-1, 0, 1):
            segs.append((lo + s, hi + s))
    segs.sort(); tot = 0; cur_lo, cur_hi = None, None
    for lo, hi in segs:
        lo, hi = max(lo, 0), min(hi, 1)
        if lo >= hi: continue
        if cur_hi is None or lo > cur_hi:
            if cur_hi is not None: tot += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
        else: cur_hi = max(cur_hi, hi)
    if cur_hi is not None: tot += cur_hi - cur_lo
    return tot
U = [0.0] * m
for mu, iv in fib.items():
    k = [i for i in range(m) if mu[i]][-1]
    U[k] += union_len(iv)
kappa = (4 / math.pi) * (1 + 2 ** -0.5)
PiJ = Pi(17) * Pi(29) * Pi(37)
out = {'brute': res_bf, 'formula': res_f, 'rel_err': rel, 'U_k_truncated_lower': U, 'sum_U': sum(U),
       'barrier_bound_(kappa/2)(PiJ-1)': kappa / 2 * (PiJ - 1)}
print(out)
json.dump(out, open('graph_check_out.json', 'w'), indent=1)
