import sys, time
from construct_q27 import *
Delta = trim([1,0,2]); C = sub(mono(27), mono(1))
m = int(sys.argv[1]); mode = sys.argv[2]; L = int(sys.argv[3])
forms=[(1,0),(0,1),(1,1),(1,2)][:m]
if mode == 'her':   # T = 1 + t C^((L-1)/27)
    e = (L-1)//27; assert 27*e+1 == L
    T = add(const(1), mul(mono(1), powp(C, e)))
else:
    s = L - 81
    T = add(const(1), mul(powp(C,3), mono(s) if s%3 else add(mono(s), mono(1))))
out, zs, R, N0 = build(Delta, C, T, forms, [1]*m, verbose=False)
print(m, mode, L, out, flush=True)
