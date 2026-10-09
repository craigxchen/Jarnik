# E4: sanity checks of the Liouville step on actual short-arc lattice clusters.
#  * F(z) = sum_k c_k N^{(D-|k|_1)/2} prod z_j^{k_j^+} conj(z_j)^{k_j^-} is a Gaussian integer (exact);
#  * Ramana identity for Vandermonde forms (M = 2s'+1): Norm(F) * N^{s'^2} = prod_{i<j} Norm(z_i - z_j)
#    (up to the constant M!^2 when the measure is the full alternating orbit sum? -- checked exactly);
#  * explicit trigonometric Taylor bound |F(z)| <= R^D (D^m/m!) ||c||_1 delta^m, delta = max|theta_j-theta_1|
#    (floating point, illustration only);
#  * F vanishes at tuples with repeated points (so the grid lemma needs no distinctness).
import itertools, math, cmath
from fractions import Fraction

def perm_sign(p):
    s = 1; p = list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s
def vandermonde_measure(v):
    M = len(v); out = {}
    for p in itertools.permutations(range(M)):
        k = tuple(v[p[i]] for i in range(M)); out[k] = out.get(k, 0) + perm_sign(p)
    return {k: x for k, x in out.items() if x}

def gmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def gpow(a, k):
    out = (1, 0)
    for _ in range(k): out = gmul(out, a)
    return out
def conj(a): return (a[0], -a[1])
def gnorm(a): return a[0] ** 2 + a[1] ** 2

def F_value(meas, D, zs):
    N = gnorm(zs[0])
    tot = (0, 0)
    for k, c in meas.items():
        l1 = sum(abs(x) for x in k)
        assert (D - l1) % 2 == 0 and l1 <= D
        term = (c * N ** ((D - l1) // 2), 0)
        for z, kj in zip(zs, k):
            term = gmul(term, gpow(z, kj) if kj > 0 else gpow(conj(z), -kj))
        tot = (tot[0] + term[0], tot[1] + term[1])
    return tot

# Cilleruelo-Granville six-point family: prod_{j=1}^4 (a + j + i sigma_j), sum sigma = 0.
def cg_points(a):
    pts = []
    for sig in itertools.product((1, -1), repeat=4):
        if sum(sig) != 0: continue
        z = (1, 0)
        for j, s in enumerate(sig, 1): z = gmul(z, (a + j, s))
        pts.append(z)
    pts.sort(key=lambda z: cmath.phase(complex(*z)))
    return pts

F3 = vandermonde_measure((1, 0, -1)); D3, m3 = 2, 3
F5 = vandermonde_measure((2, 1, 0, -1, -2)); D5, m5 = 6, 10
worst = {3: -math.inf, 5: -math.inf}
nchecks = 0
for a in (10, 31, 100, 317, 1000, 3162):
    pts = cg_points(a)
    N = gnorm(pts[0]); R = math.sqrt(N)
    assert all(gnorm(z) == N for z in pts)
    th = [cmath.phase(complex(*z)) for z in pts]
    arcC = (max(th) - min(th)) * R / math.sqrt(R)
    for M, meas, D, m, sp in ((3, F3, D3, m3, 1), (5, F5, D5, m5, 2)):
        for sub in itertools.combinations(range(6), M):
            zs = [pts[i] for i in sub]
            Fv = F_value(meas, D, zs)
            prodnorm = 1
            for i, j in itertools.combinations(range(M), 2):
                prodnorm *= gnorm((zs[i][0] - zs[j][0], zs[i][1] - zs[j][1]))
            # Ramana identity (alternating orbit sum of v = (s',...,-s')): Norm(F) N^{s'^2} = prod Norm(z_i - z_j)
            assert gnorm(Fv) * N ** (sp * sp) == prodnorm, (a, sub)
            # Taylor bound (floats)
            ths = [th[i] for i in sub]
            delta = max(abs(t - ths[0]) for t in ths)
            c1 = sum(abs(c) for c in meas.values())
            bound_log = D * math.log(R) + m * math.log(D) - math.lgamma(m + 1) + math.log(c1) + m * math.log(delta)
            lhs = 0.5 * math.log(gnorm(Fv)) if Fv != (0, 0) else -math.inf
            worst[M] = max(worst[M], lhs - bound_log)
            nchecks += 1
    print(f"a={a}: N={N}, six points on an arc of normalized length C={arcC:.3f}; "
          f"log|F_3| ~ {0.5*math.log(gnorm(F_value(F3, 2, pts[:3]))):.2f} vs (D-m/2) log R = {(2-1.5)*math.log(R):.2f}")
print(f"Ramana identity exact on {nchecks} subtuples; max (log|F| - log Taylor bound): M=3 {worst[3]:.3f}, M=5 {worst[5]:.3f} (must be <= 0)")
# repeated points: F vanishes
pts = cg_points(31)
print("F_3 at (z,z,w):", F_value(F3, 2, [pts[0], pts[0], pts[1]]), " F_5 at (z,z,w,x,y):", F_value(F5, 6, [pts[0], pts[0], pts[1], pts[2], pts[3]]))
