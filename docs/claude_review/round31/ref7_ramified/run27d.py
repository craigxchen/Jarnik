import sys, time
from construct_q27 import *
# argv: m L Delta(coeffs low->high, digits)
m = int(sys.argv[1]); L = int(sys.argv[2]); Delta = trim(np.array([int(c) for c in sys.argv[3]]))
C = sub(mono(27), mono(1))
_, r = divmod_(C, Delta) if deg(Delta) > 0 else (None, Z0)
print('Delta', sys.argv[3], 'Delta | C:', len(r) == 0)
forms=[(1,0),(0,1),(1,1),(1,2)][:m]
s = L - 81
T = add(const(1), mul(powp(C,3), mono(s) if s%3 else add(mono(s), mono(1))))
out, zs, R, N0 = build(Delta, C, T, forms, [1]*m, verbose=False)
print(m, L, out, flush=True)
