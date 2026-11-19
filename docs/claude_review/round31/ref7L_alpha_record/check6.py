from fractions import Fraction as Fr
from check3 import marks, stats
# seed-dependence of the constant Delta*H^(1/4) along fixed-T curves (T=2), n = 40
T = 2
print('denominator growth: seeds (q+1)/q, 2, 3, 5')
for q in (2, 5, 17, 101, 1009, 10007):
    sd = [(q+1, q), (2, 1), (3, 1), (5, 1)]
    vals = [stats(marks(T, n, sd))[2] for n in (39, 40)]
    print(f'  q={q:6d}: Delta*H^(1/4) at n=39,40 = {vals[0]:.4g}, {vals[1]:.4g}')
print('integer seeds growing: (1,2,3,m)')
for m in (5, 50, 500, 5000, 50000):
    sd = [(1, 1), (2, 1), (3, 1), (m, 1)]
    vals = [stats(marks(T, n, sd))[2] for n in (39, 40)]
    print(f'  m={m:6d}: {vals[0]:.4g}, {vals[1]:.4g}')
print('integer seeds spread: (k,2k,3k,5k)')
for k in (1, 10, 100, 1000):
    sd = [(k, 1), (2*k, 1), (3*k, 1), (5*k, 1)]
    vals = [stats(marks(T, n, sd))[2] for n in (39, 40)]
    print(f'  k={k:6d}: {vals[0]:.4g}, {vals[1]:.4g}')
