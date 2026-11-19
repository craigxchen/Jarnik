import time
from construct_q27 import *
Delta = trim([1,0,2]); C = trim([2]+[0]*25+[0]); C = sub(mono(27), mono(1))
forms=[(1,0),(0,1),(1,1)]
s=2
T = add(const(1), mul(powp(C,3), mono(s)))
out, zs, R, N0 = build(Delta, C, T, forms, [1,1,1])
print(out)
t0=time.time()
print('direct gcd(R,R\') degree:', deg(gcd(R, deriv(R))), time.time()-t0)
# check R' = a0 * T' structure: R'/T' has no common factor with R is the gcd above; also check T' | R'
q_, r_ = divmod_(deriv(R), deriv(T)); print("T' | R':", len(r_)==0)
