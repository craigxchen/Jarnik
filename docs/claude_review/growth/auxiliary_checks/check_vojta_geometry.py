"""Exact checks supporting conditional.md, Sections 2-4.

(A) Singular locus of the affine cone {u_1 v_1 = ... = u_M v_M}: Jacobian rank
    at points with a prescribed number of vanishing pairs (exact rational rank).
(B) Eight-point translated Pell family (primitive_eight_point_translated_pell_family.md):
    for each 7-subtuple compute h = log max|coord| (after rational content removal),
    lambda_Y = log(max|coord| / max_j |coord(z_j - z_1)|), and the Vojta excess
    Phi = 5 lambda_Y - 2 h  (the lower bound for h_K in Lemma 4.1);  check
    lambda_Y >= (1/2) log R' - log C - (1/2) log 2 with C = max chord / sqrt(R).
(C) The Pell family is the image of a rational QUARTIC curve Gamma: P^1 -> X_8 and
    Gamma(+-sqrt5 : 1) lies on the diagonal Y (exact arithmetic in Q(sqrt5)(i)).
(D) Cilleruelo-Granville six-point polynomial family: degree-4 curve with contact
    order 2 with Y at a = infinity; lambda_Y/h -> 1/2.
"""
from fractions import Fraction
from itertools import combinations, product
from math import gcd, log, isqrt
import random, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from gauss import mul, conj, norm, sub, ggcd

# ---------------- (A) ----------------
def rank_Q(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = None
        for r in range(rk, len(rows)):
            if rows[r][c] != 0:
                piv = r; break
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for r in range(len(rows)):
            if r != rk and rows[r][c] != 0:
                f = rows[r][c] / rows[rk][c]
                rows[r] = [a - f*b for a, b in zip(rows[r], rows[rk])]
        rk += 1
    return rk

def jac(M, pt):  # pt = [(u_j, v_j)]; equations f_j = u_1 v_1 - u_j v_j, j = 2..M
    rows = []
    u1, v1 = pt[0]
    for j in range(1, M):
        r = [0]*(2*M)
        r[0], r[1] = v1, u1
        r[2*j], r[2*j+1] = -pt[j][1], -pt[j][0]
        rows.append(r)
    return rows

def random_point(M, zero_pairs, rng):
    pt = []
    t = 0 if zero_pairs else None
    if t is None:
        t = rng.randint(1, 50)
        for j in range(M):
            # u_j v_j = t with u_j | t
            divs = [d for d in range(1, t+1) if t % d == 0]
            u = rng.choice(divs); pt.append((u, t//u))
        return pt
    for j in range(M):
        if j < zero_pairs:
            pt.append((0, 0))
        else:
            pt.append((rng.randint(1, 9), 0) if rng.random() < 0.5 else (0, rng.randint(1, 9)))
    return pt

def check_A():
    rng = random.Random(1)
    print("(A) Jacobian rank of the cone {u_1v_1=...=u_Mv_M} (full rank = M-1)")
    for M in range(3, 8):
        res = {}
        for zp in range(0, 4):
            ranks = set()
            for _ in range(20):
                pt = random_point(M, zp, rng)
                ranks.add(rank_Q(jac(M, pt)))
            res[zp] = sorted(ranks)
        print(f"   M={M}: vanishing pairs -> ranks {res}")
    # points of Y: u_j = u, v_j = v, (u,v) != 0
    for M in range(3, 8):
        rk = {rank_Q(jac(M, [(u, v)]*M)) for (u, v) in [(1, 0), (0, 1), (2, 3), (5, -7)]}
        print(f"   M={M}: rank on Y = {sorted(rk)}  (Y in smooth locus iff all = {M-1})")

# ---------------- (B),(C) ----------------
LAM = (9, 4)  # 9 + 4 sqrt5
def pell(r):
    x, y = 1, 0
    for _ in range(abs(r)):
        x, y = 9*x + 20*y, 4*x + 9*y
    return (x, -y) if r < 0 else (x, y)
F = (1, 2)
def H(r):
    x, y = pell(r)
    return (x + y, 2*y)

def pell_points(n):
    D, A, B, C = H(n-2), H(n), H(n+2), H(n+4)
    pts = []
    for word in product((0, 1), repeat=4):   # 1 = unbarred
        p = sum(word)
        if p % 2 == 0: continue
        z = F if p == 1 else conj(F)
        for f, w in zip((D, A, B, C), word):
            z = mul(z, f if w else conj(f))
        pts.append(z)
    return pts

def tuple_stats(pts):
    coords = [c for z in pts for c in z]
    g = 0
    for c in coords: g = gcd(g, c)
    pts_r = [(z[0]//g, z[1]//g) for z in pts]
    N = norm(pts_r[0]); assert all(norm(z) == N for z in pts_r)
    mx = max(abs(c) for z in pts_r for c in z)
    z1 = pts_r[0]
    md = max(max(abs(z[0]-z1[0]), abs(z[1]-z1[1])) for z in pts_r)
    maxchord2 = max(norm(sub(a, b)) for a, b in combinations(pts_r, 2))
    h = log(mx); lam = log(mx) - log(md)
    logR = 0.5*log(N)
    return dict(g=g, h=h, lam=lam, logR=logR, maxchord2=maxchord2, N=N)

def check_B():
    print("(B) Pell eight-point family: 7-subtuples, Vojta excess Phi = 5 lambda_Y - 2h")
    for n in (61, 121, 181, 241, 301, 361):
        pts = pell_points(n)
        N = norm(pts[0]); assert all(norm(z) == N for z in pts) and len(set(pts)) == 8
        logR = 0.5*log(N)
        logC = 0.5*log(max(norm(sub(a, b)) for a, b in combinations(pts, 2))) - 0.25*log(N)
        worst_ineq = 1e9; ratios = []
        for S in combinations(range(8), 7):
            st = tuple_stats([pts[i] for i in S])
            # Lemma 4.1 inequality lambda_Y >= (1/2) log R' - log C - (1/2) log 2
            slack = st['lam'] - (0.5*st['logR'] - logC - 0.5*log(2))
            worst_ineq = min(worst_ineq, slack)
            ratios.append((5*st['lam'] - 2*st['h'])/st['h'])
        print(f"   n={n}: log R={logR:9.2f}  log C={logC:.3f}  min slack of Lemma-4.1 ineq={worst_ineq:.3f}"
              f"  Phi/h in [{min(ratios):.4f},{max(ratios):.4f}]")

# polynomial arithmetic: dict {(i,j): Gaussian coeff} for X^i Y^j, homogeneous
def pmul(P, Q):
    out = {}
    for (a, b), c in P.items():
        for (d, e), f in Q.items():
            k = (a+d, b+e); v = mul(c, f)
            o = out.get(k, (0, 0)); out[k] = (o[0]+v[0], o[1]+v[1])
    return {k: v for k, v in out.items() if v != (0, 0)}
def pconj(P):
    return {k: conj(v) for k, v in P.items()}
def peval(P, x, y):
    s = (0, 0)
    for (a, b), c in P.items():
        m = x**a * y**b
        s = (s[0] + c[0]*m, s[1] + c[1]*m)
    return s

def H_lin(k):
    # H_{n+k} as a Z[i]-linear form in (X,Y) = (x_n, y_n):  (x_{n+k}, y_{n+k}) = Mat^k (X,Y)
    a, b, c, d = 1, 0, 0, 1   # matrix [[a,b],[c,d]]
    step = (9, 20, 4, 9) if k >= 0 else (9, -20, -4, 9)
    for _ in range(abs(k)):
        a, b, c, d = step[0]*a + step[1]*c, step[0]*b + step[1]*d, step[2]*a + step[3]*c, step[2]*b + step[3]*d
    # x = aX + bY, y = cX + dY ; H = (x + y) + 2 i y
    return {(1, 0): (a + c, 2*c), (0, 1): (b + d, 2*d)}

def check_C():
    print("(C) Pell family as a rational quartic Gamma(X:Y) and Gamma(+-sqrt5:1) in Y")
    D, A, B, C = H_lin(-2), H_lin(0), H_lin(2), H_lin(4)
    words = []
    for word in product((0, 1), repeat=4):
        p = sum(word)
        if p % 2 == 0: continue
        P = {(0, 0): F if p == 1 else conj(F)}
        for f, w in zip((D, A, B, C), word):
            P = pmul(P, f if w else pconj(f))
        words.append(P)
    degs = {a+b for P in words for (a, b) in P}
    ok = all(peval(P, *pell(61)) == z for P, z in zip(words, pell_points(61)))
    print(f"   coordinate polynomials homogeneous of degree {sorted(degs)}; reproduce n=61 points: {ok}")
    # exact evaluation at (X,Y) = (s*sqrt5, 1): numbers a + b sqrt5 with Gaussian a,b
    for sgn in (1, -1):
        vals = []
        for P in words:
            re = [Fraction(0), Fraction(0)]; im = [Fraction(0), Fraction(0)]  # (rational, sqrt5 part)
            for (a, b), c in P.items():
                # (s sqrt5)^a = s^a 5^(a//2) sqrt5^(a%2)
                coef = sgn**a * 5**(a//2)
                part = a % 2
                re[part] += c[0]*coef; im[part] += c[1]*coef
            vals.append((tuple(re), tuple(im)))
        same = all(v == vals[0] for v in vals)
        print(f"   Gamma({'+' if sgn>0 else '-'}sqrt5:1): all eight words equal (point of Y): {same}")

# ---------------- (D) ----------------
def check_D():
    print("(D) Cilleruelo-Granville six-point family v_sigma(a)=prod_{j=1..4}(a+j+i sigma_j), sum sigma=0")
    sigmas = [s for s in product((1, -1), repeat=4) if sum(s) == 0]
    for a in (10**3, 10**5, 10**7, 10**9):
        pts = []
        for s in sigmas:
            z = (1, 0)
            for j, sj in enumerate(s, start=1):
                z = mul(z, (a + j, sj))
            pts.append(z)
        st = tuple_stats(pts)
        print(f"   a={a:>10}: points={len(set(pts))}, lambda_Y/h={st['lam']/st['h']:.4f}, log(chord/sqrt R)={0.5*log(st['maxchord2'])-0.25*log(st['N']):.4f}")
    # contact order at a = infinity: compare coefficient lists of the six degree-4 polynomials
    polys = []
    for s in sigmas:
        P = [(1, 0)]  # coefficients (Gaussian) of a^0..a^k
        for j, sj in enumerate(s, start=1):
            c = (j, sj)  # factor (a + c)
            Q = [(0, 0)] * (len(P) + 1)
            for k, pk in enumerate(P):
                Q[k + 1] = (Q[k + 1][0] + pk[0], Q[k + 1][1] + pk[1])
                m = mul(pk, c); Q[k] = (Q[k][0] + m[0], Q[k][1] + m[1])
            P = Q
        polys.append(P)
    agree = [len({P[k] for P in polys}) == 1 for k in range(5)]
    print(f"   coefficient of a^k identical for all six sigma, k=4,3,2,1,0: {agree[::-1]}  -> contact order 2 with Y at a=oo")

if __name__ == '__main__':
    check_A(); check_B(); check_C(); check_D()
