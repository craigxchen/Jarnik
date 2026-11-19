import math, sys
sys.path.insert(0, '.')
from gauss import *
c = [(2, -3), (0, -3), (-1, -1), (-3, -1), (-2, -1), (-2, -2), (1, -2), (2, -2), (0, -2), (1, -3)]
for T in [10**8, 10**16]:
    for lab, B in [('T+i*isqrt(T)', (T, math.isqrt(T))), ('T*(7+i)', (7*T, T)), ('T+i (rational line Im B=1)', (T, 1))]:
        G = [(B[0]+a, B[1]+b) for a, b in c]
        z = tuple_from_blocks(G); H, logD = width_and_H(z)
        print(f'T=1e{int(math.log10(T))} B={lab}: nu = {6+20*logD/math.log(H):.3f}')
