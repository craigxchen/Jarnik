# Judge check of IX.7 mechanism (exact integers; cosines in mpmath).
# (A) fake lift r = pi*gamma/conj(gamma) with genuine angle phi = arg pi: find the least n and cutoff at which
#     the unit -1 ray test |cos(n phi)| >= m(n e_j,-1) p^{-n/2} fails (m = exterior odd part of Re of
#     (r/conj r)^n-collisions = exterior part of Re(alpha^n) when the unit is 1).
# (B) genuine lift r = pi with a wrong angle phi = -arg pi + pi/3: least failing n.
# (C) genuine lift and genuine angle: no failure for n <= 200.
from mpmath import mp, mpf, cos, atan2, pi as PI
mp.dps = 60
from sympy import factorint
def gm(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gp(a,n):
    r=(1,0)
    for _ in range(n): r=gm(r,a)
    return r
def test(alpha, unitsign, p, phi, E, nmax=60):
    # r/conj r = unitsign * alpha/conj alpha ; collisions with -1 at q^k: (r/conj r)^n = -1 mod q^k
    # <=> unitsign^n alpha^n + conj(alpha)^n = 0 mod q^k ; for unitsign^n = 1 this is q^k | Re alpha^n.
    for n in range(1, nmax+1):
        an = gp(alpha, n)
        s = unitsign**n
        val = an[0] if s == 1 else an[1]   # s=-1: alpha^n - conj alpha^n = 2i Im
        if val == 0: continue
        f = factorint(abs(val))
        ext = 1; cut = 1
        for qq, e in f.items():
            if qq % 2 == 1 and qq not in E:
                ext *= qq**e; cut = max(cut, qq**e)
        lhs = abs(cos(n*phi)); rhs = mpf(str(ext))/mpf(p)**(mpf(n)/2)
        if lhs < rhs*(1-mpf(10)**-40): return n, ext, cut, float(lhs), float(rhs)
    return None
p, pi_ = 13, (3, 2)
phi_gen = atan2(2, 3)
# (A) gamma = 2+i : alpha = pi*gamma^2
alpha = gm(pi_, gp((2,1),2))
print('A fake lift, genuine angle:', test(alpha, 1, p, phi_gen, {2,5,13}))
# (A') another gamma: 4+i (norm 17)
alpha = gm(pi_, gp((4,1),2))
print("A' fake lift (gamma=4+i):", test(alpha, 1, p, phi_gen, {2,13,17}))
# (B) genuine lift, angle -arg pi + pi/3
print('B genuine lift, wrong angle:', test(pi_, 1, p, -phi_gen + PI/3, {2,13}))
print('B2 genuine lift, angle +arg pi + pi/5:', test(pi_, 1, p, phi_gen + PI/5, {2,13}))
# (C) genuine
print('C genuine:', test(pi_, 1, p, phi_gen, {2,13}, nmax=80))
print('C2 genuine conj angle:', test(pi_, 1, p, -phi_gen, {2,13}, nmax=80))
