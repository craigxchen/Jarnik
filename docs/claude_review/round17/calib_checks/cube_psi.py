"""Psi_cube on the circle C_1 (uncut ternary cube {-1,0,1}^q, probability normalisation), 50-digit arithmetic, from the
exact formula of reduction Thm 5.1(a):  Psi_cube(u) = (m!)^2 [t^m] prod_x (1 + (2/3) u_x^2 t + (1/18) u_x^4 t^2).
Also E_cube l_u^{2m} (exact moment recursion over independent coordinates) for the relative residual.
Prints both at theta = 0 and theta = 0.37 (constancy check, Prop 2.1 applies to the cube).
usage: python3 cube_psi.py q mmax
"""
import sys
import mpmath as mp
mp.mp.dps = 50
q, mmax = int(sys.argv[1]), int(sys.argv[2])


def psi_cube(u, m):
    poly = [mp.mpf(1)]
    for ux in u:
        a, b = mp.mpf(2) / 3 * ux ** 2, ux ** 4 / 18
        new = [mp.mpf(0)] * (len(poly) + 2)
        for i, c in enumerate(poly):
            new[i] += c; new[i + 1] += a * c; new[i + 2] += b * c
        poly = new[:m + 1] + [mp.mpf(0)] * 0
    return mp.factorial(m) ** 2 * (poly[m] if m < len(poly) else 0)


def moment(u, p):
    # E (sum_x u_x k_x)^p, k_x uniform on {-1,0,1}: distribution convolution via moment generating polynomial
    # use cumulant-free direct: E X^p by dynamic programming over coordinates on the moment sequence
    mom = [mp.mpf(1)] + [mp.mpf(0)] * p                 # moments of partial sum
    for ux in u:
        # moments of k ux: E (k ux)^r = (2/3) ux^r for r even > 0, 0 for odd, 1 for r = 0
        mk = [mp.mpf(1)] + [(mp.mpf(2) / 3 * ux ** r if r % 2 == 0 else 0) for r in range(1, p + 1)]
        new = []
        for r in range(p + 1):
            new.append(mp.fsum(mp.binomial(r, s) * mom[s] * mk[r - s] for s in range(r + 1)))
        mom = new
    return mom[p]


for th in (mp.mpf(0), mp.mpf('0.37')):
    u = [mp.cos(th + 2 * mp.pi * x / q) for x in range(q)]
    for m in range(1, mmax + 1):
        P = psi_cube(u, m); E = moment(u, 2 * m)
        print(f"q={q} theta={float(th):.2f} m={m}: Psi_cube = {mp.nstr(P, 20)}  E_cube l^2m = {mp.nstr(E, 20)}  rel = {mp.nstr(P/E, 12)}")
