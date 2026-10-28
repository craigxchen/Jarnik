"""Check 1 (exact): the intra-cluster n-point data are one point of T_M (Lemmas 2.1, 2.2, Thm 3.1).
For actual normalised clusters:
 (a) y_S computed from the gcd-reduced tuple zeta_S = z_S/g_S agrees with pr_S(Z): chi_s(y_S) = Z^{iota_S s};
 (b) v_pi(zeta_j) = a_j - min_S a  and  N/N(g_S) = prod p^{max_S a - min_S a}  (heights of subtuples = toric heights
     of projections; gcd reduction is a function of v_pi(Z));
 (c) for random F on T_5^n (with forced cancellations): F(y_S1..y_Sn) == G(Z), G = Phi^*F;
     width_{Q_F}((v_pi(y_i))_i) == width_{Phi^* supp F}(a_p) >= width_{supp G}(a_p);
     Norm(G(Z)) * prod_p p^{width_{supp G}(a_p)} is a positive INTEGER (toric Liouville on T_M, exact);
 (d) sharing: sum of the one-point right sides for the projections pr_i(supp F) versus the n-point right side.
Also prints subtuple heights in [W/2, W] (Lemma 1.2) and approximation exponents (evidence only)."""
import random, math, sys
from fractions import Fraction as Fr
from itertools import combinations
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from glib import *

def rand_s(d, rng):
    s = [0]*5
    for _ in range(d):
        if rng.random() < 0.8:
            a, b = rng.sample(range(5), 2); s[a] += 1; s[b] -= 1
    return tuple(s)

def pull(blocks, svec, M):
    k = [0]*M
    for S, s in zip(blocks, svec):
        for j, sj in zip(S, s): k[j] += sj
    return tuple(k)

def run(N, M, label, rng, trials=40):
    d, Z, Nn, g0 = normalised_cluster(N, M)
    prof = profile(Z, Nn)
    W = math.log(Nn); R = math.sqrt(math.sqrt(Nn))
    print(f"== {label}: N'={Nn} M={M} primes={[(p,e) for p,e,_,_ in prof]} W={W:.4f}")
    # (a),(b) all 5-subsets
    okA = okB = okH = True; hs = []; gam = []
    lamZ = -math.log(max(abs(complex(*z1)/complex(*z2)-1) for z1 in Z for z2 in Z if z1 != z2))
    for S in combinations(range(M), 5):
        gS = Z[S[0]]
        for j in S[1:]: gS = ggcd(gS, Z[j])
        zeta = [gexdiv(Z[j], gS) for j in S]
        for p, e, pi, a in prof:
            mn = min(a[j] for j in S); mx = max(a[j] for j in S)
            if [vpi(zz, pi) for zz in zeta] != [a[j]-mn for j in S]: okB = False
        HS = Nn // gnorm(gS)
        if HS != math.prod(p**(max(a[j] for j in S)-min(a[j] for j in S)) for p, e, pi, a in prof): okB = False
        h = math.log(HS); hs.append(h/W)
        if not (W/2 - 1e-12 <= h <= W + 1e-12): okH = False
        lam = -math.log(max(abs(complex(*zeta[x])/complex(*zeta[y])-1) for x in range(5) for y in range(5) if x != y))
        gam.append(h/lam)
        s = tuple(rng.choice([-2, -1, 0, 1, 2]) for _ in range(4)); s = s + (-sum(s),)
        k = [0]*M
        for j, sj in zip(S, s): k[j] = sj
        if char_value(zeta, s) != char_value(Z, k): okA = False
    hs.sort(); gam.sort()
    print(f"  (a) chi_s(y_S) via reduced tuple == Z^(iota_S s), all {len(hs)} 5-subsets: {okA}")
    print(f"  (b) v_pi(zeta_S)=a-min and N/N(g_S)=prod p^(max-min): {okB};  h_B(y_S)/W in [1/2,1] all: {okH}"
          f"  (min {hs[0]:.4f}, max {hs[-1]:.4f})")
    print(f"      lambda(Z)={lamZ:.4f}  W/lambda(Z)={W/lamZ:.4f} ; subtuple exponents h_B/lambda min/med/max ="
          f" {gam[0]:.4f}/{gam[len(gam)//2]:.4f}/{gam[-1]:.4f}  (evidence)")
    # (c),(d) random F
    okC = okW = okI = True; strict = 0; gains = []
    for t in range(trials):
        n = rng.choice([2, 3, 4]); blocks = [tuple(sorted(rng.sample(range(M), 5))) for _ in range(n)]
        degs = [rng.choice([1, 2, 3]) for _ in range(n)]
        supp = {}
        for _ in range(rng.choice([3, 5, 8])):
            sv = tuple(rand_s(dd, rng) for dd in degs)
            supp[sv] = (rng.randint(-3, 3), rng.randint(-3, 3))
        # forced cancellation: find a pair of block-tuples with the same pullback
        sh = [(i, j) for i in range(n) for j in range(n) if i < j and len(set(blocks[i]) & set(blocks[j])) >= 2]
        if sh:
            i, j = rng.choice(sh); com = sorted(set(blocks[i]) & set(blocks[j])); a, b = com[0], com[1]
            def unit(S, a, b):
                s = [0]*5; s[S.index(a)] = 1; s[S.index(b)] = -1; return tuple(s)
            base = list(next(iter(supp))); c0 = (rng.randint(1, 3), 0)
            v1 = base[:]; v1[i] = tuple(x+y for x, y in zip(v1[i], unit(blocks[i], a, b)))
            v2 = base[:]; v2[j] = tuple(x+y for x, y in zip(v2[j], unit(blocks[j], a, b)))
            supp[tuple(v1)] = c0; supp[tuple(v2)] = (-c0[0], -c0[1])
        supp = {s: c for s, c in supp.items() if c != (0, 0)}
        if not supp: continue
        # reduced tuples for each block
        zetas = []
        for S in blocks:
            gS = Z[S[0]]
            for j in S[1:]: gS = ggcd(gS, Z[j])
            zetas.append([gexdiv(Z[j], gS) for j in S])
        Fy = (Fr(0), Fr(0))
        for sv, c in supp.items():
            val = q_of(c)
            for zs, s in zip(zetas, sv): val = qmul(val, char_value(zs, s))
            Fy = qadd(Fy, val)
        G = {}
        for sv, c in supp.items():
            k = pull(blocks, sv, M); G[k] = (G.get(k, (0, 0))[0]+c[0], G.get(k, (0, 0))[1]+c[1])
        G = {k: c for k, c in G.items() if c != (0, 0)}
        GZ = (Fr(0), Fr(0))
        for k, c in G.items(): GZ = qadd(GZ, qmul(q_of(c), char_value(Z, k)))
        if Fy != GZ: okC = False
        if GZ == (0, 0): continue
        rhsF = rhsG = rhs1 = 0.0; cert = qnorm(GZ)
        for p, e, pi, a in prof:
            vy = []
            for S in blocks:
                mn = min(a[j] for j in S); vy.append([a[j]-mn for j in S])
            wF = max(sum(sum(x*y for x, y in zip(s, v)) for s, v in zip(sv, vy)) for sv in supp) - \
                 min(sum(sum(x*y for x, y in zip(s, v)) for s, v in zip(sv, vy)) for sv in supp)
            wPhi = width([pull(blocks, sv, M) for sv in supp], a)
            wG = width(list(G), a)
            w1 = sum(width([sv[i] for sv in supp], vy[i]) for i in range(n))
            if wF != wPhi or wG > wF: okW = False
            if wG < wF: strict += 1
            rhsF += 0.5*wF*math.log(p); rhsG += 0.5*wG*math.log(p); rhs1 += 0.5*w1*math.log(p)
            cert *= Fr(p)**wG
        if cert.denominator != 1 or cert < 1: okI = False
        gains.append(rhs1 - rhsF)
    print(f"  (c) F(y)==G(Z): {okC}; width_F == width_Phi*suppF >= width_suppG: {okW} (strict at {strict} prime-instances);"
          f" Norm(G(Z))*prod p^width integer>=1: {okI}")
    gains.sort()
    print(f"  (d) sharing gain sum_i RHS_i - RHS_F over {len(gains)} F: min {gains[0]:.3f} med {gains[len(gains)//2]:.3f}"
          f" max {gains[-1]:.3f} (>=0 always: {gains[0] >= -1e-12})")

if __name__ == "__main__":
    rng = random.Random(25)
    run(1176852625, 10, "record 10-pt cluster", rng)
    run(1176852625, 9, "9-pt subcluster", rng)
    run(5**2*13**2*17*29*37*41, 8, "8-pt window, squarefree-ish norm", rng)
