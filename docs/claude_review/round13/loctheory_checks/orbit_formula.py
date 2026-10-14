"""Check of the orbit formula for binary-code sections (Theorem 3.5 of loctheory.md).

1D: on I = {-j..j}, eps(t) = (-1)^t, F_e = polys of degree <= e.  G_e = F_e cap eps F_{2j-e} is 1-dimensional
(computed exactly here) and eps acts on G_j by lambda_j = +-1.
b-dim, binary code C (list of supports): H_{P_C}(D) = #{C-orbits O on [0,2j]^b, not killed, min degree D},
where s in C acts by e_q -> 2j - e_q (q in s) and an orbit is killed iff some s in Stab(e) has
prod_{q in s} lambda_{e_q} = -1 (only e_q = j coordinates are fixed, so this is lambda_j^{|s|}).
Compared with the Hilbert function computed directly from the point set (exact ranks for small cases).
"""
import itertools, sys
from fractions import Fraction
import sympy as sp
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from lfi import hilbert

def one_d(j):
    t = sp.symbols('t')
    pts = list(range(-j, j + 1))
    N = 2 * j + 1
    # functions as value vectors; F_e = span of t^0..t^e
    def vec(poly):
        return sp.Matrix([poly.subs(t, x) for x in pts])
    eps = sp.diag(*[(-1) ** abs(x) for x in pts])
    basisF = [vec(t ** k) for k in range(N)]
    lam = {}
    dims = {}
    for e in range(N):
        A = sp.Matrix.hstack(*basisF[:e + 1])
        B = eps * sp.Matrix.hstack(*basisF[:2 * j - e + 1])
        M = sp.Matrix.hstack(A, -B)
        ns = M.nullspace()
        dims[e] = len(ns)
        if e == j:
            v = A * ns[0][:e + 1, 0]
            w = eps * v
            ratio = [w[i] / v[i] for i in range(N) if v[i] != 0]
            lam[e] = ratio[0]
            assert all(r == ratio[0] for r in ratio)
    return dims, lam[j]

def orbit_H(b, j, C, lamj):
    """C: list of codeword supports (as tuples), closed under symmetric difference, including ()."""
    seen = set()
    H = {}
    for e in itertools.product(range(2 * j + 1), repeat=b):
        if e in seen:
            continue
        orb = set()
        killed = False
        for s in C:
            f = list(e)
            for q in s:
                f[q] = 2 * j - f[q]
            f = tuple(f)
            orb.add(f)
            if f == e and len(s) > 0 and lamj ** len(s) == -1:
                killed = True
        seen |= orb
        if not killed:
            D = min(sum(x) for x in orb)
            H[D] = H.get(D, 0) + 1
    return [H.get(D, 0) for D in range(max(H) + 1)]

def span_code(gens, b):
    C = {()}
    for g in gens:
        new = set()
        for c in C:
            new.add(tuple(sorted(set(c) ^ set(g))))
        C |= new
    return sorted(C)

def section(b, j, C):
    return [y for y in itertools.product(range(-j, j + 1), repeat=b)
            if all(sum(y[q] for q in s) % 2 == 0 for s in C)]

if __name__ == '__main__':
    for j in range(1, 6):
        dims, lj = one_d(j)
        print(f'1D j={j}: dim G_e = {[dims[e] for e in range(2*j+1)]}, lambda_j = {lj}')
    tests = [
        (2, [(0, 1)]), (3, [(0, 1, 2)]), (4, [(0, 1, 2, 3)]), (4, [(0, 1, 2)]),
        (5, [(0, 1, 2), (0, 3, 4)]),             # Sylvester F=[1..5] in F_2^3 type
        (4, [(0, 1), (2, 3)]), (4, [(0, 1, 2), (1, 2, 3)]),
        (5, [(0, 1, 2, 3, 4)]), (6, [(0, 1, 2), (0, 3, 4), (1, 3, 5)]),
    ]
    for b, gens in tests:
        C = span_code(gens, b)
        for j in (1, 2, 3):
            if (2 * j + 1) ** b > 4000:
                continue
            _, lj = one_d(j)
            Horb = orbit_H(b, j, C, lj)
            P = section(b, j, C)
            Hdir = hilbert(P, 10 ** 6, cap=2 * j, modp=False)
            print(f'b={b} code gens={gens} j={j}: |P|={len(P)} orbit formula {Horb} direct {Hdir} ',
                  'EQUAL' if Horb == Hdir else 'DIFF')
            sys.stdout.flush()
