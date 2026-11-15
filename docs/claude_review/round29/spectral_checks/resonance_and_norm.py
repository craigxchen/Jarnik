# (1) Resonance (Prop 2.4): g=1+2i, A=floor(sqrt(N/5)); z_s=g(s+iA), z'_s=g(-s+iA) equal norm 5(A^2+s^2),
#     in window [N-W,N+W] with W=5R for 1<=s<=S. Verify exactly and count ordered pairs within angle Theta.
# (2) Li-Liu functional norm (7.15) for a bump of relative width delta: Q ~ delta^{-1} (e^{2}/2pi) ||phi'||^2.
import math, sys, numpy as np
from math import isqrt
for N in [10**8+1, 2576450045, 10**12+39, 10**16+61]:
    R = math.isqrt(N); Rf = math.sqrt(N); W = 5*math.ceil(Rf); C = 4.0
    Theta = C/math.sqrt(Rf)
    A = isqrt(N//5)
    assert N - W <= 5*A*A <= N
    cnt = 0; s = 1
    while 5*(A*A+s*s) <= N+W:
        zx, zy = s-2*A, 2*s+A; px, py = -s-2*A, -2*s+A
        assert zx*zx+zy*zy == px*px+py*py == 5*(A*A+s*s)
        z = complex(zx, zy); zp = complex(px, py)
        ang = abs(math.atan2(z.imag, z.real) - math.atan2(zp.imag, zp.real))
        if ang <= Theta: cnt += 2
        s += 1
    th = math.degrees(math.atan2(1, -2))
    print(f"N={N} sqrtR={Rf**0.5:.1f} ordered resonant pairs within Theta (C=4): {cnt}  ratio to sqrtR: {cnt/Rf**0.5:.3f}  direction {th:.2f} deg")
# (2)
from scipy.integrate import quad
def bump(x): return math.exp(-1/(1-x*x)) if abs(x)<1 else 0.0
def dbump(x): return bump(x)*(-2*x/(1-x*x)**2) if abs(x)<1 else 0.0
I1 = quad(lambda x: dbump(x)**2, -1, 1)[0]
for delta in [0.1, 0.03, 0.01, 0.003]:
    # Phi(u)=sqrt(u) phi((u-1)/delta); V = Phi - u Phi'
    def V(u):
        x=(u-1)/delta
        return math.sqrt(u)*bump(x) - u*(0.5/math.sqrt(u)*bump(x) + math.sqrt(u)*dbump(x)/delta)
    Q = quad(lambda u: V(u)**2/u, 1-delta, 1+delta, limit=200)[0]
    print(f"delta={delta}: int |V|^2 du/u = {Q:.5g}; delta*that = {delta*Q:.5g} (predicted ||phi'||^2 = {I1:.5g})")
