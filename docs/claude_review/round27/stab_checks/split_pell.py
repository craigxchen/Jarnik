"""Split twin of the Cox form on the round-23 Pell curve (a,b,c)=(1,2,4): z^2 = 56 x^2 + 65 (T3/T5 calibration).
For Pell points (x,z): u_j = F_j(x,z) (j=0..4), made primitive positive; split Cox data u_A, v_A per cut
A subset of the first four indices (side not containing the 5th point); check u_a = prod_{A ni a} u_A prod_{A not ni a} v_A,
u_5 = prod v_A, residues r_a = U_a - V_a, balance of the 15 cut norms n_A = u_A v_A, and alpha-type constant.
"""
from math import gcd, log
from itertools import combinations
from sympy import factorint

def F(x, z):
    return [(x-1)*(x-2)*(x-4)*(x+7), x**4 - 35*x**2 + 56 + 2*x*z, (x*x-4)*(x*x-16), (x*x-1)*(x*x-16), (x*x-1)*(x*x-4)]

CUTS = [frozenset(s) for k in range(1, 5) for s in combinations([0, 1, 2, 3], k)]

def split_cox(u):
    """u: 5 positive coprime ints, index 4 is the base point."""
    uA = {A: 1 for A in CUTS}; vA = {A: 1 for A in CUTS}
    primes = set()
    for x in u:
        primes |= set(factorint(x).keys())
    for p in primes:
        v = [0]*5
        for j in range(5):
            x = u[j]
            while x % p == 0:
                x //= p; v[j] += 1
        for t in range(1, max(v)+1):
            S = frozenset(j for j in range(5) if v[j] >= t)
            if 4 not in S:
                uA[S] *= p
            else:
                vA[frozenset(range(5)) - S] *= p
    return uA, vA

if __name__ == '__main__':
    pts = [(1, 11)]
    for _ in range(7):
        x, z = pts[-1]
        pts.append((15*x + 2*z, 15*z + 112*x))
    for (x, z) in pts:
        assert z*z == 56*x*x + 65
        u = F(x, z)
        if min(u) <= 0:
            u = [-y for y in u] if max(u) < 0 else u
        if min(u) <= 0:
            print('x=%d: nonpositive entries, skip' % x); continue
        g = 0
        for y in u: g = gcd(g, y)
        u = [y // g for y in u]
        # reorder so the base point (index 4) is the smallest
        order = sorted(range(5), key=lambda j: u[j])
        u = [u[j] for j in order[1:]] + [u[order[0]]]
        uA, vA = split_cox(u)
        ok = True
        for a in range(4):
            pa = 1
            for A in CUTS:
                pa *= uA[A] if a in A else vA[A]
            ok &= (pa == u[a])
        p5 = 1
        for A in CUTS: p5 *= vA[A]
        ok &= (p5 == u[4])
        U = [1]*4; V = [1]*4
        for a in range(4):
            for A in CUTS:
                if a in A:
                    U[a] *= uA[A]; V[a] *= vA[A]
        r = [U[a] - V[a] for a in range(4)]
        nA = [uA[A]*vA[A] for A in CUTS]
        lcm = 1
        for A in CUTS: lcm *= uA[A]*vA[A]
        H = max(u) / min(u) * lcm
        dist = (max(u) - min(u)) / min(u)
        alpha = -log(H) / log(dist) if dist > 0 else float('nan')
        print('x=%-12d cox_ok=%s residues U_a-V_a=%s  log n_A range=[%.1f,%.1f] (log x=%.1f)  alpha=%.4f' % (
            x, ok, r, min(log(n) for n in nA), max(log(n) for n in nA), log(x), alpha))
