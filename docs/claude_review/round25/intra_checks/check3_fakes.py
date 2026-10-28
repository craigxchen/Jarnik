"""Check 3 (evidence + exact pieces): actual cluster versus same-profile angle models ("profile fakes"),
and a supercriticality scan of intra-cluster/route-5 supports at M = 10.

Part A (exact/numerical): for the record 10-point cluster, the angle system theta_x = sum_p (2a_xp - e_p) phi_p.
  Build profile fakes phi' = phi + t psi (psi random, t small so that the diameter stays <= 1.05 x actual).
  Every input of L_intra except the VALUE F(y) -- widths RHS_F, heights h_B(y_S), gcd data, the support of the
  multiplicative relation lattice -- is a function of the profile, hence identical for actual and fake (printed);
  the subtuple depths lambda(y_S) agree to the printed tolerance. Monomial relations Z^k = 1 (k in the relation
  lattice {k : <k,a_p> = 0 for all p}) hold for actual AND fake (exact for the actual point).
Part B (evidence): for supports S(P0, b) = {k : k_x in {-1,0,1} (x != o), <k,a_p> = 0 (p in P0), |<k,a_p>| <= b_p}
  compute reg(S) (upper bound via rank mod a 61-bit prime: rank_p <= rank_Q, so reg_p >= reg_Q), the depth
  L(S) = -log max_{k in S} |<k,delta>| and the comparison
     E(S) = reg(S) L(S) + log reg(S)! - RHS(S),   RHS(S) = (1/2) sum_p width_S(a_p) log p.
  For every G with supp G = S and ord_e G = m <= reg(S), the Taylor bound is log|G(Z)| <= log|c|_1 - m L(S) - log m!
  with |c|_1 >= 1, so a contradiction (or a forced exact relation) from such G with full width needs E(S) > 0.
  (Sub-supports are themselves scanned only through the family; this is evidence, not a proof.)"""
import math, random, sys, itertools
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from glib import *

P61 = (1 << 61) - 1

def rank_mod(rows, p=P61):
    rows = [[x % p for x in r] for r in rows]; rk = 0
    ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c]), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        inv = pow(rows[rk][c], p-2, p); rows[rk] = [x*inv % p for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c]:
                f = rows[i][c]; rows[i] = [(x-f*y) % p for x, y in zip(rows[i], rows[rk])]
        rk += 1
        if rk == len(rows): break
    return rk

def monos(nv, t):
    out = []
    for deg in range(t):
        for c in itertools.combinations_with_replacement(range(nv), deg):
            a = [0]*nv
            for i in c: a[i] += 1
            out.append(a)
    return out

def reg_upper(S, o):
    idx = [x for x in range(len(S[0])) if x != o]; t = 0
    while True:
        E = [[math.prod(k[x]**a for x, a in zip(idx, al)) for al in monos(len(idx), t+1)] for k in S]
        if rank_mod(E) == len(S): return t
        t += 1
        if t > 12: return t

def main():
    rng = random.Random(7)
    d, Z, Nn, g0 = normalised_cluster(1176852625, 10); M = 10
    prof = profile(Z, Nn)
    # actual angles
    th = [math.atan2(z[1], z[0]) for z in Z]
    o = 0
    delta = [math.remainder(th[x]-th[o], 2*math.pi) for x in range(M)]
    phis = [math.atan2(pi[1], pi[0]) for p, e, pi, a in prof]
    B = [[2*(a[x]-a[o]) for p, e, pi, a in prof] for x in range(M)]
    # rank of B (exact via Fractions)
    from fractions import Fraction as Fr
    def rk(m):
        m = [[Fr(v) for v in r] for r in m]; r = 0
        for c in range(len(m[0])):
            pv = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
            if pv is None: continue
            m[r], m[pv] = m[pv], m[r]
            for i in range(len(m)):
                if i != r and m[i][c] != 0:
                    f = m[i][c]/m[r][c]; m[i] = [x-f*y for x, y in zip(m[i], m[r])]
            r += 1
        return r
    rB = rk(B)
    print(f"Part A: 10-pt cluster, r = {len(prof)} primes, rank B = {rB}, free dimension m = r - rank B = {len(prof)-rB}")
    # check angle model reproduces actual relative angles (mod pi/2 units absorbed)
    def rel_angles(ph):
        return [sum(B[x][j]*ph[j] for j in range(len(prof))) for x in range(M)]
    ra = rel_angles(phis)
    dev = max(abs(math.remainder(ra[x]-delta[x], math.pi/2)) for x in range(M))
    print(f"  angle model at phi = arg(pi_p): max |rel angle - actual| mod pi/2 = {dev:.2e}")
    diam = max(abs(delta[x]-delta[y]) for x in range(M) for y in range(M))
    # profile fakes
    for trial in range(3):
        psi = [rng.gauss(0, 1) for _ in prof]
        dpsi = rel_angles(psi); spread = max(dpsi)-min(dpsi)
        t = 0.05*diam/spread
        php = [ph+t*ps for ph, ps in zip(phis, psi)]
        raf = rel_angles(php)
        dlt_f = [delta[x] + (raf[x]-ra[x]) for x in range(M)]
        diam_f = max(abs(dlt_f[x]-dlt_f[y]) for x in range(M) for y in range(M))
        lam_dev = 0.0
        for S in itertools.combinations(range(M), 5):
            la = -math.log(max(abs(2*math.sin((delta[x]-delta[y])/2)) for x in S for y in S if x != y))
            lf = -math.log(max(abs(2*math.sin((dlt_f[x]-dlt_f[y])/2)) for x in S for y in S if x != y))
            lam_dev = max(lam_dev, abs(la-lf))
        zdist = max(abs(dlt_f[x]-delta[x]) for x in range(M))
        print(f"  fake {trial}: diam ratio {diam_f/diam:.4f}; max |arg Z_f - arg Z| = {zdist:.3e} (Z_f != Z);"
              f" max_S |lambda_f(y_S) - lambda(y_S)| = {lam_dev:.4f}; widths/heights/gcds: profile-level, identical")
    # relation lattice: k with <k,a_p> = 0 for all p, k in {-1,0,1}^9 (anchor o)
    rel = 0; ok_rel = True
    for kp in itertools.product((-1, 0, 1), repeat=M-1):
        k = [0]*M; k[1:] = list(kp); k[0] = -sum(kp)
        if all(sum(kk*ax for kk, ax in zip(k, a)) == 0 for p, e, pi, a in prof) and any(kp):
            rel += 1
            v = char_value(Z, k)
            if v != (1, 0): ok_rel = False
    print(f"  monomial relations Z^k = 1 for k in relation lattice (k' in {{-1,0,1}}^9): {rel} found, all exactly 1: {ok_rel}"
          f" (fakes satisfy them identically: arg Z_f^k = sum_p 2<k,a_p> phi'_p = 0)")
    # Part B
    print("Part B: supercriticality scan (evidence)")
    lamZ = -math.log(max(abs(2*math.sin((delta[x]-delta[y])/2)) for x in range(M) for y in range(M)))
    print(f"  lambda(Z) = {lamZ:.4f}")
    primes = [p for p, e, pi, a in prof]
    best = None; rows = []
    allk = []
    for kp in itertools.product((-1, 0, 1), repeat=M-1):
        k = [0]*M; k[1:] = list(kp); k[0] = -sum(kp); allk.append(k)
    for r0 in range(0, len(prof)+1):
        for P0 in itertools.combinations(range(len(prof)), r0):
            for b in (0, 1, 2, 3):
                S = []
                for k in allk:
                    good = True
                    for j, (p, e, pi, a) in enumerate(prof):
                        v = sum(kk*ax for kk, ax in zip(k, a))
                        if (j in P0 and v != 0) or abs(v) > b: good = False; break
                    if good: S.append(k)
                if len(S) < 2 or len(S) > 260: continue
                m = reg_upper(S, o)
                L = -math.log(max(abs(sum(kk*dd for kk, dd in zip(k, delta))) for k in S if any(k)))
                rhs = 0.5*sum(width(S, a)*math.log(p) for p, e, pi, a in prof)
                E = m*L + math.lgamma(m+1) - rhs
                rows.append((E, [primes[j] for j in P0], b, len(S), m, L, rhs))
    rows.sort(reverse=True)
    for r in rows[:8]:
        print(f"  E={r[0]:8.3f}  P0={r[1]} b={r[2]} |S|={r[3]} reg<={r[4]} L={r[5]:.3f} RHS={r[6]:.3f}")
    print(f"  supports scanned: {len(rows)}; max E = {rows[0][0]:.3f}; any E > 0: {any(r[0] > 0 for r in rows)}")

main()
