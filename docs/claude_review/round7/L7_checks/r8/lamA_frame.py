# lambda_A(C) = deg F_X|_C for Mobius test curves C, computed in local frames of T_X (mod p = 2^61-1).
#
# Method. X = Mbar_{0,7} embeds in prod_{35 four-subsets} Mbar_{0,4} = (P^1)^35 (Keel), so T_X is a
# subbundle of the pulled-back tangent bundle and ord_t of a section of C^*T_X is the minimum over the
# 35 factors of the order of its component in a local frame of T_{P^1}.
# v(t) = d/du of the leaf through C(t) (leaf = ruling of Q_t containing R_t, in the G(2,4) model);
# v is a rational section of F_X along C (up to a rational factor, which does not change the degree).
# Component k: cross-ratio cr_k = N/D (quadratic in Plucker coords); W = N_u D - N D_u at u = 0;
# g = gcd(N0, D0). The component in the T_{P^1} frame is the rational function W/g^2 (orders at finite
# t), and its order at t = infinity is 2 deg h_k - deg(W/g^2), deg h_k = degree of t -> cr_k(C(t)).
# lambda(C) = sum_{t in A^1} min_k ord_t + min_k ord_inf   (valid if C avoids Sing(F_X)).
import sys, random, itertools
import flint
from star import point_of, quadric_through, LAB
from arr import PLANES, CHOSEN
from elamp import P, md, kernel_at, newton_interp, good_mob

def C(c): return flint.nmod_poly([c % P], P)
ZERO = C(0)
pl = [[md(c) for c in H] for H in PLANES]
qb = [[[C(md(Q[i, j])) for j in range(4)] for i in range(4)] for Q in quadric_through(
      [point_of([PLANES[i], PLANES[j], PLANES[k]]) for (i, j, k) in CHOSEN])]
PAIRS = [(i, j) for i in range(4) for j in range(i + 1, 4)]
# pi_ab(Plucker) = (l_a ^ l_b)(A ^ B) = sum_{i<j} (l_a,i l_b,j - l_a,j l_b,i) P_ij
PI = {}
for a in range(7):
    for b in range(7):
        if a != b:
            PI[(a, b)] = [(pl[a][i] * pl[b][j] - pl[a][j] * pl[b][i]) % P for (i, j) in PAIRS]

def gcd_list(L):
    g = None
    for k in L:
        if k.is_zero(): continue
        g = k if g is None else g.gcd(k)
    return g
def content_div(v):
    g = gcd_list(v)
    return [c // g for c in v]
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
def kernel_poly(mob):
    xs = list(range(3, 3 + 12))
    kers = [kernel_at(mob, x) for x in xs]
    return content_div([newton_interp(xs, [k[j] for k in kers]) for j in range(8)])

def leaf_data(mob):
    K = kernel_poly(mob)
    p, pp = K[:4], K[4:]
    rows = [[bil(Qk, u, v) for Qk in qb] for (u, v) in [(p, p), (p, pp), (pp, pp)]]
    cs = []
    for j in range(4):
        sub = [[r[k] for k in range(4) if k != j] for r in rows]
        d = det3(sub)
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
    P0 = wedge(p, q20)
    P1 = [x + y for x, y in zip(wedge(p, q21), wedge(w, q20))]
    P0 = content_div(P0); P1 = content_div(P1)
    return P0, P1

def lin(coef, Pv):
    s = ZERO
    for c, x in zip(coef, Pv):
        if c: s = s + x * c
    return s

def lam_of(mob, verbose=False):
    P0, P1 = leaf_data(mob)
    br0 = {k: lin(PI[k], P0) for k in PI}
    br1 = {k: lin(PI[k], P1) for k in PI}
    return lam_from_br(br0, br1, verbose)

def lam_from_br(br0, br1, verbose=False):
    # br0[(a,b)] = bracket [p_a p_b] along C (u = 0), br1 = its u-derivative along the leaf
    comps = []
    for (a, b, c, d) in itertools.combinations(range(7), 4):
        N0 = br0[(a, c)] * br0[(b, d)]; D0 = br0[(a, d)] * br0[(b, c)]
        Nu = br1[(a, c)] * br0[(b, d)] + br0[(a, c)] * br1[(b, d)]
        Du = br1[(a, d)] * br0[(b, c)] + br0[(a, d)] * br1[(b, c)]
        W = Nu * D0 - N0 * Du
        assert not N0.is_zero() and not D0.is_zero()
        if W.is_zero(): continue
        g = N0.gcd(D0)
        degh = max((N0 // g).degree(), (D0 // g).degree())
        num, den = W, g * g
        h = num.gcd(den); num = num // h; den = den // h
        comps.append((num, den, degh))
    L = C(1)
    for (num, den, _) in comps:
        L = L * den // L.gcd(den)
    G = gcd_list([num * (L // den) for (num, den, _) in comps])
    fin = G.degree() - L.degree()
    inf = min(2 * dh - (num.degree() - den.degree()) for (num, den, dh) in comps)
    if verbose:
        print('   fin', fin, 'inf', inf, 'ncomp', len(comps), 'degL', L.degree(), flush=True)
    return fin + inf

if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    lm = [lam_of(good_mob(()), verbose=True) for _ in range(reps)]
    print('mu', lm, flush=True)
    tot = 0
    out = {}
    for T in itertools.combinations(range(7), 3):
        ls = [lam_of(good_mob(T)) for _ in range(reps)]
        out[T] = ls
        print([LAB[a] for a in T], ls, ' diff from mu:', ls[0] - lm[0], flush=True)
        tot += ls[0]
    print('lambda_A(beta) = sum_T lambda(mu_T) - 32 lambda(mu) =', tot - 32 * lm[0], flush=True)
