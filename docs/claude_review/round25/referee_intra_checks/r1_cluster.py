# Referee check r1: Lemma 2.1, Lemma 3.3, Theorem 3.1 on actual clusters, independent code.
import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round25/referee_intra_checks')
from rlib import *
from math import sqrt
random.seed(25)

def best_window(N, M):
    pts = lattice_points(N); R = sqrt(N); best = None
    n = len(pts)
    for i in range(n):
        w = [pts[(i+j) % n] for j in range(M)]
        L = max(sqrt(gnorm((a[0]-b[0], a[1]-b[1]))) for a in w for b in w)
        if best is None or L < best[0]: best = (L, w)
    return best[0]/sqrt(R), best[1]

def analyse(N, M, nF=60):
    C, z = best_window(N, M)
    g = z[0]
    for x in z[1:]: g = ggcd(g, x)
    print(f"== N={N} M={M} C={C:.4f} gcd norm={gnorm(g)}")
    if gnorm(g) != 1:
        z = [gdiv_exact(x, g) for x in z]; N = N//gnorm(g); print("   normalised N'=", N)
    R = sqrt(N); W = log(N)
    fac = factor(N); primes = sorted(fac)
    pis = {p: split_prime(p) for p in primes}
    for p in primes: assert p % 4 == 1
    a = {p: [vpi(x, pis[p]) for x in z] for p in primes}
    for p in primes:
        for x, ax in zip(z, a[p]): assert vpi(x, gconj(pis[p])) == fac[p]-ax
    lamZ = -log(max(sqrt(gnorm((u[0]-v[0], u[1]-v[1]))) for u in z for v in z)/R)
    ok1 = ok2 = ok3 = ok4 = True; ratios = []
    ys = {}
    for S in combinations(range(M), 5):
        gS = z[S[0]]
        for j in S[1:]: gS = ggcd(gS, z[j])
        zeta = [gdiv_exact(z[j], gS) for j in S]
        # (2) valuations
        for p in primes:
            mn = min(a[p][j] for j in S)
            for jj, j in enumerate(S):
                if vpi(zeta[jj], pis[p]) != a[p][j]-mn: ok2 = False
        # (3) heights
        NgS = 1
        for p in primes: NgS *= p**(min(a[p][j] for j in S)+fac[p]-max(a[p][j] for j in S))
        if gnorm(gS) != NgS: ok3 = False
        hB = sum((max(a[p][j] for j in S)-min(a[p][j] for j in S))*log(p) for p in primes)
        if abs(hB-log(N/gnorm(gS))) > 1e-9: ok3 = False
        ratios.append(hB/W)
        # (1) characters: s random in M_5
        for _ in range(3):
            s = [random.randint(-3, 3) for _ in range(4)]; s.append(-sum(s))
            k = [0]*M
            for jj, j in enumerate(S): k[j] = s[jj]
            if char_val(zeta, s) != char_val(z, k): ok1 = False
        lamS = -log(max(sqrt(gnorm((z[u][0]-z[v][0], z[u][1]-z[v][1]))) for u in S for v in S)/R)
        if lamS < lamZ-1e-12: ok4 = False
        ys[S] = zeta
    print(f"  L2.1(1) chars {ok1}  (2) vals {ok2}  (3) heights {ok3} hB/W in [{min(ratios):.4f},{max(ratios):.4f}]  (4) lam {ok4}")
    # Theorem 3.1 / Lemma 3.3: random F on n blocks, compare F(y) vs G(Z); widths; integrality
    blocks_all = list(ys)
    okv = okw = oki = True; strict = 0; gains = []
    for t in range(nF):
        n = random.choice([2, 3, 4]); Bs = random.sample(blocks_all, n)
        terms = {}
        for _ in range(random.randint(2, 6)):
            ss = []
            for i in range(n):
                s = [random.randint(-2, 2) for _ in range(4)]; s.append(-sum(s)); ss.append(tuple(s))
            c = (random.randint(-3, 3), random.randint(-3, 3))
            if c == (0, 0): c = (1, 0)
            terms[tuple(ss)] = c
        # forced cancellation: add a term and its negative image under a different block tuple with same pullback
        Fval = (Fr(0), Fr(0)); G = {}
        for ss, c in terms.items():
            v = (Fr(1), Fr(0)); k = [0]*M
            for i, S in enumerate(Bs):
                v = qmul(v, char_val(ys[S], list(ss[i])))
                for jj, j in enumerate(S): k[j] += ss[i][jj]
            Fval = qadd(Fval, qmul(tofr(c), v)); k = tuple(k)
            G[k] = (G.get(k, (0, 0))[0]+c[0], G.get(k, (0, 0))[1]+c[1])
        G = {k: c for k, c in G.items() if c != (0, 0)}
        Gval = (Fr(0), Fr(0))
        for k, c in G.items(): Gval = qadd(Gval, qmul(tofr(c), char_val(z, list(k))))
        if Gval != Fval: okv = False
        if Gval == (Fr(0), Fr(0)) or not G: continue
        RHS_F = RHS_G = sumproj = 0.0
        for p in primes:
            vF = [sum(sum(s[jj]*a[p][j] for jj, j in enumerate(S)) for s, S in zip(ss, Bs)) for ss in terms]
            wF = max(vF)-min(vF)
            wG = width(list(G), a[p])
            if wG > wF: okw = False
            if wG < wF: strict += 1
            RHS_F += 0.5*wF*log(p)*2/2; RHS_G += 0.5*wG*log(p)
            for i, S in enumerate(Bs):
                vv = [sum(ss[i][jj]*a[p][j] for jj, j in enumerate(S)) for ss in terms]
                sumproj += 0.5*(max(vv)-min(vv))*log(p)
        RHS_F *= 1.0
        gains.append(sumproj-0.5*RHS_F*2+0*RHS_F)
        # Lemma 3.3: Norm(G(Z)) * prod p^{width_suppG} is a positive integer
        val = qnorm(Gval)
        for p in primes: val *= p**width(list(G), a[p])
        if not (val.denominator == 1 and val >= 1): oki = False
    print(f"  Thm3.1(1) F(y)==G(Z): {okv}; width_suppG<=width_QF: {okw} (strict {strict}); L3.3 integral>=1: {oki}")
    print(f"  sharing gain sum_i RHS(pr_i) - RHS_F: min {min(gains):.3f} max {max(gains):.3f}  W={W:.3f}  lamZ={lamZ:.4f}")
    return z, N, a, primes, fac, lamZ

if __name__ == '__main__':
    analyse(1176852625, 10)
    analyse(1176852625, 9)
    analyse(5*13*17*29*37*41*53, 8)
