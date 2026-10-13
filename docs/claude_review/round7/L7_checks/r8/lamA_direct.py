# Intrinsic lambda_A(C) = deg F_A|_C for Mobius test curves C, mod p = 2^61-1.
# U = universal conic over B; along C the members form a P^1-bundle over C with section Sigma (R_t)
# and event sections s_i (the rho_t-line through P_i). Sigma_U^2 = 2 Sigma.s_a - s_a.s_b - s_a.s_c + s_b.s_c.
# lambda_X(C) = Sigma_U^2 - k(C), k = number of contacts with the 42 divisors not met by gamma_A.
import sys, random, itertools
from fractions import Fraction as Fr
import flint
from star import point_of, quadric_through, LAB
from arr import PLANES, CHOSEN
from elamp import P, md, det_mod, newton_interp, pl, rnd

X = flint.nmod_poly([0, 1], P)
def C(c): return flint.nmod_poly([c % P], P)
pts = [[md(c) for c in point_of([PLANES[i], PLANES[j], PLANES[k]])] for (i, j, k) in CHOSEN]
qb = [[[md(Q[i, j]) for j in range(4)] for i in range(4)] for Q in quadric_through(
      [point_of([PLANES[i], PLANES[j], PLANES[k]]) for (i, j, k) in CHOSEN])]

def gcd_list(L):
    g = L[0]
    for k in L[1:]:
        g = g.gcd(k)
    return g
def prim(v):
    g = gcd_list([c for c in v if not c.is_zero()])
    return [c // g for c in v]
def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
def bil(Q, u, v):
    s = C(0)
    for i in range(4):
        for j in range(4):
            if Q[i][j] is not None and not Q[i][j].is_zero():
                s = s + u[i] * Q[i][j] * v[j]
    return s
def kernel_poly(mob):
    xs = list(range(3, 3 + 12))
    from elamp import kernel_at
    kers = [kernel_at(mob, x) for x in xs]
    return prim([newton_interp(xs, [k[j] for k in kers]) for j in range(8)])
def plucker(u, v):
    return prim([u[i] * v[j] - u[j] * v[i] for i in range(4) for j in range(i + 1, 4)])
def meet(p, q):
    # number of t (with multiplicity, finite t) where Plucker vectors p, q are proportional
    mins = [p[i] * q[j] - p[j] * q[i] for i in range(6) for j in range(i + 1, 6)]
    mins = [m for m in mins if not m.is_zero()]
    return gcd_list(mins).degree()
def lam_data(mob):
    K = kernel_poly(mob)
    L1, L2 = K[:4], K[4:]
    pR = plucker(L1, L2)
    # quadric Q_t = sum c_k qb_k with c = 3x3 minors of the condition matrix
    qbp = [[[C(Qk[i][j]) for j in range(4)] for i in range(4)] for Qk in qb]
    rows = [[bil(Qk, u, v) for Qk in qbp] for (u, v) in [(L1, L1), (L1, L2), (L2, L2)]]
    cs = []
    for j in range(4):
        sub = [[r[k] for k in range(4) if k != j] for r in rows]
        d = det3(sub)
        cs.append(d if j % 2 == 0 else -d)
    cs = prim(cs)
    Qt = [[sum((cs[k] * qbp[k][i][j] for k in range(4)), C(0)) for j in range(4)] for i in range(4)]
    secs = []
    for e in range(6):
        Pe = [C(c) for c in pts[e]]
        tau = [sum((Pe[i] * Qt[i][j] for i in range(4)), C(0)) for j in range(4)]   # tangent plane form
        tl1 = sum((tau[j] * L1[j] for j in range(4)), C(0)); tl2 = sum((tau[j] * L2[j] for j in range(4)), C(0))
        r = [L1[j] * tl2 - L2[j] * tl1 for j in range(4)]                             # R_t ∩ T
        # basis of T / <P_e>: u1 = tau1 E0 - tau0 E1, u2 = tau2 E0 - tau0 E2
        u1 = [tau[1], -tau[0], C(0), C(0)]; u2 = [tau[2], C(0), -tau[0], C(0)]
        # coordinates of r mod P_e in (u1,u2): (x1:y1) = (det[P,r,u2] : det[P,u1,r]) on rows 0..2 (generic)
        rowsel = [0, 1, 3]
        def d3(a, b, c):
            return det3([[a[i], b[i], c[i]] for i in rowsel])
        x1 = d3(Pe, r, u2); y1 = d3(Pe, u1, r)
        al = bil(Qt, u1, u1); ga = bil(Qt, u2, u2)
        w = [ga * y1 * u1[i] + al * x1 * u2[i] for i in range(4)]
        secs.append(plucker(Pe, w))
    return pR, secs

def sigma2(pR, secs, a, b, c):
    return 2 * meet(pR, secs[a]) - meet(secs[a], secs[b]) - meet(secs[a], secs[c]) + meet(secs[b], secs[c])

UNCOV_PAIRS = 21 - 5
def k_of(T):
    # mu_T: pairs all 2 except the three in T (1 each); one contact with D_T.
    # covered pair divisors: a2b1,a2b2,a2c2,b1c2,b2c2 ; covered triples: coaxial + events
    covered_pairs = {frozenset(p) for p in [(2, 3), (2, 4), (2, 6), (3, 6), (4, 6)]}
    covered_trip = {frozenset(t) for t in [(0, 1, 2), (0, 3, 4), (0, 5, 6)] + CHOSEN}
    k = 0
    for pr in itertools.combinations(range(7), 2):
        deg = 2 - (1 if T and set(pr) <= set(T) else 0)
        if frozenset(pr) not in covered_pairs: k += deg
    if T and frozenset(T) not in covered_trip: k += 1
    return k

def good_mob(T):
    while True:
        c = [[rnd() for _ in range(4)] for a in range(7)]
        if T:
            a0 = T[0]; z0 = c[a0][1] / c[a0][3]
            for b in T[1:]: c[b][1] = z0 * c[b][3]
        if any(ci[0] * ci[3] - ci[1] * ci[2] == 0 for ci in c): continue
        return [[md(v) for v in ci] for ci in c]

if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
    Ts = [()] + list(itertools.combinations(range(7), 3))
    if len(sys.argv) > 2: Ts = [()] + list(itertools.combinations(range(7), 3))[:int(sys.argv[2])]
    res = {}
    for T in Ts:
        pR, secs = lam_data(good_mob(T))
        s2 = [sigma2(pR, secs, *abc) for abc in [(0, 1, 2), (3, 4, 5), (1, 3, 5), (2, 0, 4)]]
        k = k_of(T)
        res[T] = (s2, k)
        print('T', [LAB[a] for a in T] if T else 'mu', 'Sigma_U^2 (4 event triples):', s2, 'k =', k,
              'lambda =', [s - k for s in s2], flush=True)
    if len(Ts) == 36:
        lam = {T: res[T][0][0] - res[T][1] for T in Ts}
        print('lambda_A(beta) = sum_T lambda(mu_T) - 32 lambda(mu) =',
              sum(lam[T] for T in Ts if T) - 32 * lam[()])
