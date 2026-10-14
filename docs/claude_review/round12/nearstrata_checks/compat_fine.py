"""Step-1 rescan of condition A of compat.py around its threshold (the ratio is not monotone because of floor/ceil).
Scans n in [0.9 a, 1.6 a] with step 1 for each r-choice and reports the last failing n and the max ratio beyond it;
then the exact rational value of the A-ratio at last+1 and the float max over [last+1, 1.6a]. Also checks that beyond 1.6 a
the coarse scan of compat.py found no failure (up to 4*10^6)."""
import math
from compat import conds_float, condA_exact
for name, rf, a in [('r = n', lambda n: n, 473466), ('r = n ln n', lambda n: n * math.log(n), 563900), ('r = n^2', lambda n: n * n, 947150)]:
    last = None; mx = 0.0
    lo, hi = int(0.9 * a), int(1.6 * a)
    vals = []
    for n in range(lo, hi + 1):
        A = conds_float(n, rf(n))[0]
        if A > 1: last = n
        vals.append((n, A))
    beyond = max(A for n, A in vals if n > last)
    print(f"{name}: last failing n in [{lo},{hi}] = {last}; max A-ratio on ({last},{hi}] = {beyond:.6f}; exact A-ratio at {last+1} = {float(condA_exact(last+1, rf(last+1))):.6f}")
