# Exact equal-norm pair counts E_off(theta) = sum_n r_I(n)(r_I(n)-1) over circles n in [N-W,N+W],
# sector I=[theta,theta+Theta], Theta = C sqrt(R)/R. Checks: (i) mean over theta equals the exact
# average A = sum_n sum_{z!=w} (Theta-angle)_+/(2pi); (ii) compares A with kappa*(W/R)*C^2*log R;
# (iii) sup_theta E_off vs sup_theta M_N(M_N-1) on the central circle.
import sys, math, numpy as np
from collections import defaultdict
N = int(sys.argv[1]); C = float(sys.argv[2]); Wfac = float(sys.argv[3])
R = math.sqrt(N); W = int(Wfac * R); Theta = C * math.sqrt(R) / R
pts = defaultdict(list)
r = int(math.isqrt(N + W)) + 1
for x in range(-r, r + 1):
    lo = N - W - x * x; hi = N + W - x * x
    if hi < 0: continue
    ylo = 0 if lo <= 0 else math.isqrt(lo - 1) + 1
    yhi = math.isqrt(hi)
    for y in range(ylo, yhi + 1):
        for yy in ([y] if y == 0 else [y, -y]):
            pts[x * x + yy * yy].append(math.atan2(yy, x) % (2 * math.pi))
# pairs with angular separation < Theta (ordered pairs counted twice later)
events = []  # (start, end) arcs of theta for which both points in sector
A = 0.0; npairs = 0
for n, angs in pts.items():
    a = np.sort(np.array(angs)); m = len(a)
    if m < 2: continue
    ext = np.concatenate([a, a + 2 * math.pi])
    for i in range(m):
        j = i + 1
        while j < i + m and ext[j] - a[i] < Theta:
            psi = ext[j] - a[i]
            A += 2 * (Theta - psi) / (2 * math.pi)   # ordered pairs
            events.append((ext[j] - Theta, a[i], n))  # theta in [ext[j]-Theta, a[i]]
            npairs += 1; j += 1
# sweep max over theta of E_off (ordered pairs = 2 * unordered)
ev = []
for s, e, n in events:
    ev.append((s % (2*math.pi), +2)); ev.append((e % (2*math.pi), -2))
    if (s % (2*math.pi)) > (e % (2*math.pi)): ev.append((0.0, +2))  # wrap
ev.sort(key=lambda t: (t[0], -t[1]))
cur = 0; best = 0; best_off = 0; hist = []
lo_ok, hi_ok = 2*Theta, math.pi/4 - 2*Theta
for t, d in ev:
    cur += d; best = max(best, cur)
    if lo_ok <= t <= hi_ok: best_off = max(best_off, cur)
# central circle cluster max
a = np.sort(np.array(pts.get(N, [])))
Mmax = 0
if len(a):
    ext = np.concatenate([a, a + 2 * math.pi])
    for i in range(len(a)):
        if not (lo_ok <= a[i] <= hi_ok): continue
        Mmax = max(Mmax, int(np.searchsorted(ext, a[i] + Theta, side='right') - i))
logR = math.log(R)
print(f"N={N} R={R:.1f} C={C} W={W} Theta={Theta:.4g} circles={len(pts)} r(N)={len(pts.get(N,[]))}")
print(f"mean_theta E_off (exact A) = {A:.4f}   A/((W/R) C^2 log R) = {A/((W/R)*C*C*logR):.4f}")
print(f"max_theta E_off = {best}  max over theta in [2Theta, pi/4-2Theta] = {best_off}   central-circle max M = {Mmax}  M(M-1) = {Mmax*(Mmax-1)}")
