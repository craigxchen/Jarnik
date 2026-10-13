# Exact (over Z) version of lamA_frame.py: lambda_A(C) for integer Mobius test curves, by local frames
# of T_X in the 35 cross-ratio factors. Same formulas; all polynomial arithmetic in Z[t] (python-flint
# fmpz_poly). Usage: python3 lamA_exact.py [seed] [forgetful]   ('forgetful' = validate on F_i instead).
import sys, random, itertools
from fractions import Fraction as Fr
import flint
from star import point_of, quadric_through, LAB
from arr import PLANES, CHOSEN

X = flint.fmpz_poly([0, 1])
def Z(c): return flint.fmpz_poly([int(c)])
ZERO = Z(0)
def prim(p):
    if p.is_zero(): return p
    c = p.content()
    q = flint.fmpz_poly([x // c for x in p.coeffs()])
    return -q if q.leading_coefficient() < 0 else q
def gcd_list(L):
    g = None
    for k in L:
        if k.is_zero(): continue
        g = prim(k) if g is None else prim(g.gcd(k))
    return g
def content_div(v):
    g = gcd_list(v)
    return [prim(c) * 0 + (c // g) if not c.is_zero() else c for c in v]
def ints(row):
    den = 1
    for c in row: den = den * Fr(c).denominator // __import__('math').gcd(den, Fr(c).denominator)
    return [int(Fr(c) * den) for c in row]
pl = [ints(list(H)) for H in PLANES]
qb = []
for Q in quadric_through([point_of([PLANES[i], PLANES[j], PLANES[k]]) for (i, j, k) in CHOSEN]):
    flat = ints([Q[i, j] for i in range(4) for j in range(4)])
    qb.append([[Z(flat[4 * i + j]) for j in range(4)] for i in range(4)])
PAIRS = [(i, j) for i in range(4) for j in range(i + 1, 4)]
PI = {(a, b): [pl[a][i] * pl[b][j] - pl[a][j] * pl[b][i] for (i, j) in PAIRS]
      for a in range(7) for b in range(7) if a != b}

def bil(Q, u, v):
    s = ZERO
    for i in range(4):
        for j in range(4):
            s = s + u[i] * Q[i][j] * v[j]
    return s
def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
def wedge(u, v):
    return [u[i] * v[j] - u[j] * v[i] for (i, j) in PAIRS]

def interp_int(xs, ys):
    # Lagrange interpolation over Q; result must be integral
    n = len(xs); coef = [Fr(0)] * n
    for i in range(n):
        num = [Fr(1)]; den = Fr(1)
        for j in range(n):
            if j == i: continue
            num = [Fr(0)] + num
            for k in range(len(num) - 1):
                num[k] -= xs[j] * num[k + 1]
            den *= (xs[i] - xs[j])
        for k in range(n):
            coef[k] += ys[i] * num[k] / den
    assert all(c.denominator == 1 for c in coef)
    return flint.fmpz_poly([int(c) for c in coef])
def kernel_poly(mob):
    xs = list(range(-5, 6))
    vals = []
    for t0 in xs:
        rows = []
        for a in range(7):
            al, be, ga, de = mob[a]
            Nv = al * t0 + be; Dv = ga * t0 + de
            rows.append([Nv * c for c in pl[a]] + [Dv * c for c in pl[a]])
        vals.append([(-1) ** j * int(flint.fmpz_mat([[r[k] for k in range(8) if k != j] for r in rows]).det())
                     for j in range(8)])
    K = [interp_int(xs, [v[j] for v in vals]) for j in range(8)]
    return content_div(K)

def leaf_brackets(mob):
    K = kernel_poly(mob)
    p, pp = K[:4], K[4:]
    rows = [[bil(Qk, u, v) for Qk in qb] for (u, v) in [(p, p), (p, pp), (pp, pp)]]
    cs = []
    for j in range(4):
        d = det3([[r[k] for k in range(4) if k != j] for r in rows])
        cs.append(d if j % 2 == 0 else -d)
    cs = content_div(cs)
    Q = [[sum((cs[k] * qb[k][i][j] for k in range(4)), ZERO) for j in range(4)] for i in range(4)]
    def other_pt(a, b):
        tau = [sum((a[i] * Q[i][j] for i in range(4)), ZERO) for j in range(4)]
        z = [tau[1], -tau[0], ZERO, ZERO]
        Qz = bil(Q, z, z); Bbz = bil(Q, b, z)
        return content_div([Qz * b[i] - 2 * Bbz * z[i] for i in range(4)])
    w = other_pt(p, pp); w2 = other_pt(pp, p)
    pQw2 = bil(Q, p, w2); wQw2 = bil(Q, w, w2); wQpp = bil(Q, w, pp)
    q20 = [pQw2 * pp[i] for i in range(4)]
    q21 = [wQw2 * pp[i] - wQpp * w2[i] for i in range(4)]
    P0 = content_div(wedge(p, q20))
    P1 = content_div([x + y for x, y in zip(wedge(p, q21), wedge(w, q20))])
    lin = lambda coef, Pv: sum((x * c for c, x in zip(coef, Pv) if c), ZERO)
    return {k: lin(PI[k], P0) for k in PI}, {k: lin(PI[k], P1) for k in PI}

def forget_brackets(mob, i):
    N = [Z(m[1]) + Z(m[0]) * X for m in mob]; D = [Z(m[3]) + Z(m[2]) * X for m in mob]
    dN = [ZERO] * 7; dD = [ZERO] * 7; dN[i] = Z(1)
    br0 = {}; br1 = {}
    for a in range(7):
        for b in range(7):
            if a != b:
                br0[(a, b)] = N[a] * D[b] - N[b] * D[a]
                br1[(a, b)] = dN[a] * D[b] + N[a] * dD[b] - dN[b] * D[a] - N[b] * dD[a]
    return br0, br1

def lam_from_br(br0, br1):
    # brackets independent of the 4-subset (G(2,4) model, forgetful)
    K = list(itertools.combinations(range(7), 4))
    b0 = {(k, x, y): br0[(x, y)] for k in K for x in k for y in k if x != y}
    b1 = {(k, x, y): br1[(x, y)] for k in K for x in k for y in k if x != y}
    return lam_from_br_k(b0, b1)

def lam_from_br_k(br0, br1):
    # br0[(k, x, y)]: bracket used for the cross-ratio of the 4-subset k (k-dependent normalisations allowed)
    comps = []
    for k in itertools.combinations(range(7), 4):
        a, b, c, d = k
        N0 = br0[(k, a, c)] * br0[(k, b, d)]; D0 = br0[(k, a, d)] * br0[(k, b, c)]
        Nu = br1[(k, a, c)] * br0[(k, b, d)] + br0[(k, a, c)] * br1[(k, b, d)]
        Du = br1[(k, a, d)] * br0[(k, b, c)] + br0[(k, a, d)] * br1[(k, b, c)]
        W = Nu * D0 - N0 * Du
        assert not N0.is_zero() and not D0.is_zero()
        if W.is_zero(): continue
        g = prim(N0.gcd(D0))
        degh = max((N0 // g).degree(), (D0 // g).degree())
        num, den = prim(W), prim(g * g)
        h = prim(num.gcd(den)); num = num // h; den = den // h
        comps.append((num, den, degh))
    L = Z(1)
    for (num, den, _) in comps:
        L = prim(L * den // prim(L.gcd(den)))
    G = gcd_list([num * (L // den) for (num, den, _) in comps])
    fin = G.degree() - L.degree()
    inf = min(2 * dh - (num.degree() - den.degree()) for (num, den, dh) in comps)
    return fin + inf

def int_mob(T):
    while True:
        r = lambda: random.choice([k for k in range(-9, 10) if k])
        c = [[r() for _ in range(4)] for a in range(7)]
        if T:
            num, den = r(), r()
            for b in T:
                k = r(); c[b][1] = num * k; c[b][3] = den * k
        if any(ci[0] * ci[3] - ci[1] * ci[2] == 0 for ci in c): continue
        # distinct points at t = 0 except inside T, and distinct at t = infinity
        v0 = [Fr(ci[1], ci[3]) if ci[3] else None for ci in c]
        vinf = [Fr(ci[0], ci[2]) if ci[2] else None for ci in c]
        ok = all(v0[a] != v0[b] for a, b in itertools.combinations(range(7), 2) if not (T and a in T and b in T))
        ok = ok and len(set(vinf)) == 7
        if not ok: continue
        # transversality / no codim-2 strata: the collision polynomials of all pairs, with the event
        # factor t removed once per pair of T, must be coprime and squarefree (and nonzero at 0)
        F = Z(1)
        for a, b in itertools.combinations(range(7), 2):
            f = (Z(c[a][1]) + Z(c[a][0]) * X) * (Z(c[b][3]) + Z(c[b][2]) * X) - \
                (Z(c[b][1]) + Z(c[b][0]) * X) * (Z(c[a][3]) + Z(c[a][2]) * X)
            if T and a in T and b in T:
                f = f // X
            F = F * f
        if F.degree() < 1 or prim(F.gcd(F.derivative())).degree() > 0: continue
        if T and int(F.coeffs()[0]) == 0: continue
        return c

if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
    mode = sys.argv[2] if len(sys.argv) > 2 else 'A'
    labels = range(7) if mode == 'forgetful' else [None]
    for i in labels:
        br = (lambda m: forget_brackets(m, i)) if i is not None else leaf_brackets
        lm = lam_from_br(*br(int_mob(())))
        tot = 0; diffs = {}
        for T in itertools.combinations(range(7), 3):
            v = lam_from_br(*br(int_mob(T)))
            diffs[T] = v - lm; tot += v
        tag = 'F_' + LAB[i] if i is not None else 'V_A'
        print(tag, 'lambda(mu) =', lm, ' lambda(beta) =', tot - 32 * lm, flush=True)
        print('   lambda(mu_T) - lambda(mu):', {''.join(LAB[a] for a in T): d for T, d in diffs.items() if d}, flush=True)
