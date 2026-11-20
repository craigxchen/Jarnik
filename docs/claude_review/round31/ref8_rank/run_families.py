from families import full_detail
from math import isqrt
from sympy import fibonacci

print("== record single clusters")
for N, C, lab, m in [(1176852625, 9, '10-pt record', 10), (1176852625, 8, '9-pt record', 9),
                     (235370525, 6, '7-pt record', 7), (32988780417125, 6, '8-pt record', 8),
                     (537606725 * 4, 8.5, 'golden n=4 (unnormalised, x4)', 8),
                     (9520005625, 1, '3-pt C=.37', 3), (412333805, 1.01, '4-pt C<=1', 4),
                     (1762645, 1.42, '4-pt C<=sqrt2', 4), (683617125425, 2.01, '5-pt C<=2', 5),
                     (10414625, 3, '5-pt C<=3', 5)]:
    full_detail(N, C, lab, m)

print("== golden family (Z_w(n), n = 1 mod 3)")
for n in [4, 7, 10, 13]:
    N = 5 * fibonacci(2*n+1) * fibonacci(2*n+3) * fibonacci(2*n+5) * fibonacci(2*n+7)
    full_detail(int(N), 8.5, f'golden n={n}', 8)

print("== axis Pell family (6 points)")
y2, y1 = 1656, 1171
for k in range(5):
    a = (y1 * y1 - 143) // 2
    N = a * a + 144
    full_detail(N, 4.3, f'Pell a={a}', 6)
    y2, y1 = 3 * y2 + 4 * y1, 2 * y2 + 3 * y1
