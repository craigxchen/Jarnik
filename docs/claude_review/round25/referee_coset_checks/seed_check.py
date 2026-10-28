"""Referee: Lemma 7.2 (Gaussian seeds on base cosets). psi_0 = g_0 + B arg(beta), g_0 in (pi/4 Z)^r.
Check: every v with mu_v = B^T v != 0 and w(v) <= T satisfies both gap conditions for weights p' = p/tau^2, exactly
(Gaussian integer arithmetic for gamma_mu, and the rotation by g_0 done exactly in Z[i] with (1+i)/sqrt2)."""
import math, itertools, json
from fractions import Fraction
B = [[1, 0], [0, 1], [1, 1], [1, -1]]
q = [5, 13]; beta = [(2, 1), (3, 2)]
tau2 = 1.5
p = [17, 41, 197, 229]
for i, row in enumerate(B):
    assert math.prod(q[k] ** abs(row[k]) for k in range(2)) <= p[i] / (2 * tau2), i
def gmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def gpow(a, n):
    r = (1, 0)
    if n < 0: a = (a[0], -a[1]); n = -n
    for _ in range(n): r = gmul(r, a)
    return r
T = 22.0
logs = [math.log(x) for x in p]
worst = {}
for g0 in itertools.product(range(4), repeat=4):    # g_0 = (pi/4) * g0 (mod pi/2 only parity matters for (pi/4)Z mod pi/2... keep 0..3)
    mins = math.inf; cnt = 0
    rngs = [range(-int(T / l), int(T / l) + 1) for l in logs]
    for v in itertools.product(*rngs):
        w = sum(abs(a) * l for a, l in zip(v, logs))
        if w > T or not any(v): continue
        mu = [sum(B[i][k] * v[i] for i in range(4)) for k in range(2)]
        if not any(mu): continue
        cnt += 1
        gam = gmul(gpow(beta[0], mu[0]), gpow(beta[1], mu[1]))
        s = sum(a * b for a, b in zip(v, g0))          # <v, g0> in units of pi/4
        X, Y = gam
        # rotation by s*pi/4: odd s swaps the two families
        N = X * X + Y * Y
        ax = min(abs(X), abs(Y)); dg = min(abs(X - Y), abs(X + Y))
        if s % 2 == 0:
            d_axis = math.asin(ax / math.sqrt(N)); d_diag = math.asin(dg / math.sqrt(2 * N))
        else:
            d_axis = math.asin(dg / math.sqrt(2 * N)); d_diag = math.asin(ax / math.sqrt(N))
        yp = math.exp(-(w - sum(abs(a) for a in v) * math.log(tau2)) / 2)   # y'_v for weights p/tau^2
        m = min(d_axis - math.asin(min(1, yp)), d_diag - math.asin(min(1, yp / math.sqrt(2))))
        # relative margin
        mins = min(mins, m / yp)
    worst[str(g0)] = (cnt, mins)
vals = [x[1] for x in worst.values()]
out = {'n_g0': len(worst), 'n_noncharacter_monomials': list(worst.values())[0][0], 'min_relative_margin': min(vals)}
print(out); json.dump(out, open('seed_check_out.json', 'w'), indent=1)
