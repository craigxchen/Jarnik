"""Star dictionary (ref4_bridge) in characteristic 2 for the Part V cluster, with explicit columns.
P_i = gcd(z_i, U), U U^bar = N, gcd(U,U^bar)=1; column (X_i,Y_i) in F_2[t]^2 with P_i = X_i + Y_i theta;
n_i = X^2+XY+Y^2; D_ij = X_iY_j + X_jY_i (= P_i Pbar_j + Pbar_i P_j in char 2).
Checks every hypothesis of (IV.2) and computes both sides, raw and after dividing by the global gcd.
Usage: python3 dict_check.py m L [npairs_sample]
"""
import sys, pickle, random, itertools, json
from f24 import *

def parity(x):
    return bin(x).count('1') & 1

m = int(sys.argv[1]); L = int(sys.argv[2])
nsample = int(sys.argv[3]) if len(sys.argv) > 3 else None
S = pickle.load(open(f'cluster_m{m}_L{L}.pkl', 'rb'))
d, Q, D, W, forms, H, sp, orient, val, bigdeg, vz, vN, z, N = (S[k] for k in
    ['d', 'Q', 'D', 'W', 'forms', 'H', 'sp', 'orient', 'val', 'bigdeg', 'vz', 'vN', 'z', 'N'])
Yv = {v: ((1 << L) ^ v, 1) for v in range(Q)}

def fpow(f, k):
    P = (1, 0)
    for _ in range(k):
        P = mul4(P, f)
    return P

big = {}
for v in range(Q):
    P = Yv[v]
    for f in sp:
        k = val.get((f, v), 0)
        if k:
            q, rr = divmod4(P, fpow(f, k))
            assert iszero4(rr)
            P = q
    assert deg4(P) == bigdeg[v]
    big[v] = P
# conj(big[v]) == big[v^1] (up to unit)
assert all(monic4(conj4(big[v])) == monic4(big[v ^ 1]) for v in range(2, Q))
reps = [v for v in range(2, Q) if v & 1 == 0]
U = (1, 0)
for v in reps:
    U = mul4(U, big[v])
for f in orient:
    U = mul4(U, fpow(f, vN[f]))
assert mul4(U, conj4(U))[1] == 0
UUb = mul4(U, conj4(U))[0]
lcU = 1
res = {}
res['U_Ubar_equals_N'] = (UUb == N)
res['gcd_U_Ubar_is_1'] = deg4(gcd4(U, conj4(U))) == 0 if m <= 32 else 'skipped (structural)'
Ub = conj4(U)
P = {}
for a in forms:
    Pa = (1, 0)
    for v in H[a]:
        if v and v in reps:
            Pa = mul4(Pa, big[v])
    for f in orient:
        Pa = mul4(Pa, fpow(f, vz[(f, a)]))
    P[a] = Pa
ok = dict(z_eq_P_Ub_over_Pb=True, primitive=True, n_div_N=True, Dnonzero=True, E0=True,
          deg_identity=True, degree_bound_delta0=True)
alpha = {}
n = {}
for a in forms:
    # z_a * Pbar_a == P_a * Ubar  (as F_4[t] polys; z monic, P,U monic -> exact)
    lhs = mul4(z[a], conj4(P[a])); rhs = mul4(P[a], Ub)
    if lhs != rhs:
        ok['z_eq_P_Ub_over_Pb'] = False
    X, Yc = P[a]
    if deg2(gcd2(X, Yc)) != 0:
        ok['primitive'] = False
    n[a] = norm4(P[a])
    alpha[a] = deg2(n[a])
    _, rr = divmod2(N, n[a])
    if rr:
        ok['n_div_N'] = False
pairs = list(itertools.combinations(forms, 2))
if nsample:
    random.seed(1)
    pairs = random.sample(pairs, nsample)
maxexcess = -10**9
for a, b in pairs:
    Xa, Ya = P[a]; Xb, Yb = P[b]
    Dab = mul2(Xa, Yb) ^ mul2(Xb, Ya)
    if Dab == 0:
        ok['Dnonzero'] = False
        continue
    g = gcd2(n[a], n[b])
    _, rr = divmod2(Dab, g)
    if rr:
        ok['E0'] = False
    dd = deg4(add4(z[a], z[b]))
    if 2 * dd != 2 * D + 2 * deg2(Dab) - (alpha[a] + alpha[b]):
        ok['deg_identity'] = False
    exc = deg2(Dab) - ((alpha[a] + alpha[b]) / 2 - W / 4)
    maxexcess = max(maxexcess, exc)
    if exc > 0:
        ok['degree_bound_delta0'] = False
res.update(ok)
res['pairs_checked'] = len(pairs)
res['max_excess_of_deg_D_over_budget'] = maxexcess
# layers: s_h over roots of N (over F_2bar); weight deg P
# raw
def iv2(mm, Lw, Qs, delta):
    return 0.75 * (mm - 24) * Lw + Qs, mm * (mm + 3) * delta
Qs = 0
glob = 0  # degree of global gcd from small primes (layers with s=0 or s=m)
for f in orient:
    e = 2 * deg4(f)
    for l in range(1, vN[f] + 1):
        s = sum(1 for a in forms if vz[(f, a)] >= l)
        Qs += e * (s - m / 2) ** 2
        if s == 0 or s == m:
            glob += e // 2   # common factor f-bar^.. or f^.. of degree e/2 in every z
# the big parts never give s in {0,m}
res['raw'] = dict(L=W, Q=Qs, delta=0, LHS=iv2(m, W, Qs, 0)[0], RHS=0)
# normalised: divide every z by g (deg glob); W' = W - 2 glob, layers with s in {0,m} removed
Qn = Qs - sum((m / 2) ** 2 * 2 * deg4(f) for f in orient for l in range(1, vN[f] + 1)
              if sum(1 for a in forms if vz[(f, a)] >= l) in (0, m))
Dn = D - glob
Kn = max(S['ddiff'].values()) - glob
deltan = Kn - Dn / 2
res['global_gcd_degree'] = glob
res['normalised'] = dict(D=Dn, L=2 * Dn, Q=Qn, delta=deltan, LHS=iv2(m, 2 * Dn, Qn, deltan)[0],
                         RHS=iv2(m, 2 * Dn, Qn, deltan)[1])
print(json.dumps(res, indent=1))
json.dump(res, open(f'dict_m{m}_L{L}.json', 'w'), indent=1)
