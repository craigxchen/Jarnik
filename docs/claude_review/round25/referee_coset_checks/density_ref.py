"""Referee: replicate §8.1 (survival of uniform points in cubes around the transplant vector), own code, T = 24."""
import math, json, sys
import numpy as np
sys.path.insert(0, '.')
from cover_check import enum, PS
split = [5, 13, 17, 29, 37, 41, 53, 61, 73]
def gp(p):
    for bb in range(1, int(math.isqrt(p)) + 1):
        a2 = p - bb * bb; a = math.isqrt(a2)
        if a * a == a2 and a > bb: return a, bb
qs = [split[split.index(p) - 1] for p in PS]
xb = np.array([math.atan2(gp(q)[1], gp(q)[0]) / (math.pi / 2) for q in qs])
V, W = enum(24.0)
y = np.exp(-W / 2); a0 = (2 / math.pi) * np.arcsin(y); a1 = (2 / math.pi) * np.arcsin(y / math.sqrt(2))
def ok(X):
    C = X @ V.T
    d0 = np.abs(C - np.round(C)); d1 = np.abs(C - 0.5 - np.round(C - 0.5))
    return np.all((d0 >= a0) & (d1 >= a1), axis=1)
rng = np.random.default_rng(2024)
res = {'n_monomials': int(len(V)), 'transplant_ok': bool(ok(xb[None, :])[0])}
for eps in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
    eps_n = eps / (math.pi / 2)   # coset.md eps is in angle units? report both conventions
    for name, e in (('eps_normalised', eps),):
        cnt = 0; N = 2000
        for _ in range(4):
            X = xb + rng.uniform(-e, e, size=(N // 4, 8)); cnt += int(ok(X).sum())
        res[f'{name}={eps}'] = cnt / N
X = rng.uniform(0, 1, size=(2000, 8)); res['uniform'] = float(ok(X).mean())
print(json.dumps(res, indent=1)); json.dump(res, open('density_ref_out.json', 'w'), indent=1)
