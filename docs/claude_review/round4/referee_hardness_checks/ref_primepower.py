#!/usr/bin/env python3
"""Referee check R12: the elementary pair bound is NOT 'attained up to a constant' on every form.
For u = (2+i)^a (pair norm n = 5^a, gcd(u, ubar) = 1) the best pair (u, eta*ubar) has
|u - eta ubar| = min(2|B|, 2|A|, sqrt2|A-B|, sqrt2|A+B|) where u = A + Bi (exact integers), and the
Liouville ratio |lambda|/(sqrt2 n^(-1/2)) is >= |u - eta ubar|/sqrt2 (chord <= R*angle).  Baker's
theorem forces this ratio to tend to infinity; we print the exact minimum chord for a <= 400."""
import math
u = (1, 0)
rows = []
for a in range(1, 401):
    u = (u[0]*2 - u[1]*1, u[0]*1 + u[1]*2)
    A, B = u
    ch2 = min(4*B*B, 4*A*A, 2*(A-B)**2, 2*(A+B)**2)   # exact squared chord |u - eta ubar|^2
    rows.append((a, ch2))
mins = []
for lo, hi in [(1, 50), (51, 100), (101, 200), (201, 400)]:
    m = min(rows[lo-1:hi], key=lambda t: t[1])
    mins.append((lo, hi, m[0], math.log10(m[1])/2 - math.log10(math.sqrt(2))))
for lo, hi, a, l in mins:
    print(f"a in [{lo},{hi}]: min over a of log10(chord/sqrt2) = {l:.2f} (at a = {a}); "
          f"i.e. ratio >= 10^{l:.1f}")
