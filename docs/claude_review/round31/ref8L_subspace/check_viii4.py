"""Referee check of Astra Part VIII.4 (independent code; exact integer / Gaussian-integer arithmetic).

System (round8/weakest.md (5.1) = Astra (VIII.6)) over K=Q(i):
  |x|_inf = |x|^2, |x|_w = p^(-v_w(x)) for a Gaussian prime w above split p.
  at inf: forms X_1, X_i - X_1 ; at finite w in S: X_1..X_k.
  LHS = prod_v prod_i |L_iv(x)|_v / (max_i |x_i|_v)^k,  H = prod_v max_i |x_i|_v.
All quantities are computed PLACE BY PLACE (no product-formula shortcut), as exact Fractions.
The test LHS <= H^(-k-delta) with delta = a/q is checked exactly as LHS^q * H^(kq+a) <= 1.
"""
import sys, random, math, itertools
from fractions import Fraction
from sympy import primerange

random.seed(1)

# ---------- Gaussian integers as (a,b) ----------
def gmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gconj(x): return (x[0], -x[1])
def gnorm(x): return x[0]*x[0]+x[1]*x[1]
def gsub(x, y): return (x[0]-y[0], x[1]-y[1])
def gdivmod(x, y):
    n = gnorm(y); num = gmul(x, gconj(y))
    # rounded quotient
    qa = (2*num[0] + n) // (2*n); qb = (2*num[1] + n) // (2*n)
    q = (qa, qb); r = gsub(x, gmul(q, y)); return q, r
def gexact(x, y):
    q, r = gdivmod(x, y); assert r == (0, 0); return q
def ggcd(x, y):
    while y != (0, 0):
        _, r = gdivmod(x, y); x, y = y, r
    return x
def gval(x, pi):
    v = 0
    while True:
        q, r = gdivmod(x, pi)
        if r != (0, 0): return v
        x = q; v += 1

def gauss_prime_above(p):
    for a in range(1, int(math.isqrt(p))+1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2 and a > b: return (a, b)
    raise ValueError

def lhs_and_H(xs, primes):
    """xs: list of Gaussian integers (tuples), assumed supported above `primes` (checked).
    returns (LHS, H) exactly, place by place."""
    k = len(xs)
    # check support: every rational prime dividing the norms is in primes
    for x in xs:
        n = gnorm(x)
        for p in primes:
            while n % p == 0: n //= p
        assert n == 1, "coordinate not an S-unit"
    LHS = Fraction(1); H = Fraction(1)
    # archimedean
    mx = max(gnorm(x) for x in xs)
    forms = [xs[0]] + [gsub(x, xs[0]) for x in xs[1:]]
    num = 1
    for f in forms: num *= gnorm(f)
    LHS *= Fraction(num, mx**k); H *= mx
    # finite places above primes
    for p in primes:
        a, b = gauss_prime_above(p)
        for pi in [(a, b), (a, -b)]:
            vs = [gval(x, pi) for x in xs]
            vmin = min(vs)
            # |x_i|_w / max = p^-(v_i - vmin) ; forms are coordinates
            LHS *= Fraction(1, p**(sum(vs) - k*vmin))
            H *= Fraction(1, p**vmin)
    return LHS, H

def passes(LHS, H, k, a, q):
    # LHS <= H^(-k - a/q)  <=>  LHS^q * H^(kq + a) <= 1
    return LHS**q * H**(k*q + a) <= 1

def delta_eff(LHS, H, k):
    return -math.log(LHS.numerator/LHS.denominator) / math.log(H.numerator/H.denominator) - k if False else \
        -(math.log(LHS.numerator) - math.log(LHS.denominator)) / (math.log(H.numerator) - math.log(H.denominator)) - k

out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s); sys.stdout.flush()

# ---------- Part 1: rational squarefree products (VIII.4 main proof) ----------
def rational_instance(X, j, k, delta_num, delta_den, nsample):
    P = [p for p in primerange(X, 2*X+1) if p % 4 == 1]
    s = len(P)
    T = (2*X)**j
    b = Fraction(delta_num, delta_den) / (k-1)          # b = delta/(k-1)
    prods = sorted(math.prod(c) for c in itertools.combinations(P, j))
    assert len(set(prods)) == len(prods)                 # distinct (unique factorisation)
    # largest window with (max-min)/min <= T^(-b): exact test ((max-min)/min)^(1/b) <= 1/T
    # 1/b = bq/bp  -> ((max-min)^bq_inv ...) use: (max-min)^n * T^m <= min^n where b = m/n
    m, n = b.numerator, b.denominator
    def ok(lo, hi): return (hi-lo)**n * T**m <= lo**n
    best = (0, 0); i0 = 0
    for i1 in range(len(prods)):
        while not ok(prods[i0], prods[i1]): i0 += 1
        if i1 - i0 + 1 > best[1] - best[0]: best = (i0, i1+1)
    E = prods[best[0]:best[1]]
    # also the fixed-partition bin count (Astra's pigeonhole): bins of ratio 1+T^-b on [X^j,T]
    nbins = math.ceil(j*math.log(2) / math.log(1 + float(T)**(-float(b))))
    say(f"[rational] X={X} s={s} j={j} k={k} delta={delta_num}/{delta_den} b={b} T={T}: "
        f"C(s,j)={len(prods)}, #bins(fixed partition)={nbins}, C(s,j)/#bins={len(prods)/nbins:.1f}, "
        f"largest window |E|={len(E)}")
    # rigorous all-tuple sufficient condition (Astra's argument): ratio^(2(k-1)) * maxE^(2 delta) <= 1
    #   with ratio = (maxE-minE)/minE ; check exactly: ((maxE-minE)/minE)^(2(k-1) q) * maxE^(2 a) <= 1
    lo, hi = E[0], E[-1]
    a, q = delta_num, delta_den
    cond = (hi-lo)**(2*(k-1)*q) * hi**(2*a) <= lo**(2*(k-1)*q)
    say(f"   all-tuple sufficient condition (exact) holds: {cond}")
    # exact place-by-place on all or sampled ordered distinct k-tuples
    tuples = list(itertools.permutations(range(len(E)), k)) if math.perm(len(E), k) <= nsample else \
        [tuple(random.sample(range(len(E)), k)) for _ in range(nsample)]
    worst = 1e9; Hmin = None; Hmax = None; allpass = True; formula_ok = True
    for t in tuples:
        us = [E[i] for i in t]
        g = math.gcd(*us)
        xs = [(u//g, 0) for u in us]
        L, H = lhs_and_H(xs, [p for p in P if g and any(u % p == 0 for u in us)])
        # independent check of the closed form H^-k prod ((u_i-u_1)/u_i)^2 and of H=(max/g)^2
        cf = Fraction(1, 1) / H**k
        for u in us[1:]: cf *= Fraction(u - us[0], u)**2
        formula_ok &= (cf == L) and (H == Fraction(max(us)//g)**2)
        allpass &= passes(L, H, k, a, q)
        de = delta_eff(L, H, k); worst = min(worst, de)
        Hmin = H if Hmin is None or H < Hmin else Hmin
        Hmax = H if Hmax is None or H > Hmax else Hmax
    T2b = Fraction(T)**2  # compare H>=T^(2b): H^n >= T^(2m)
    lowok = Hmin**n >= Fraction(T)**(2*m)
    say(f"   checked {len(tuples)} ordered distinct tuples place by place: all solve (5.1): {allpass}; "
        f"closed form exact: {formula_ok}; min delta_eff={worst:.4f}; "
        f"H in [{float(Hmin):.3e},{float(Hmax):.3e}], H>=T^(2b): {lowok}, H<=T^2: {Hmax <= T2b}")
    NE = len(E)
    say(f"   subspace lower bound (|E|)_k/|E|^(k-1) = {math.perm(NE,k)/NE**(k-1):.2f}  (s = {s})")
    return NE, s

# ---------- Part 2: equal-norm pair monomials ----------
def equal_norm_instance(X, j, delta_num, delta_den, bnum, bden, mode, nmax_E=None):
    """mode 'pair5': five of the six pair monomials of 4 inputs (Astra, k=5);
       mode 'pair6': all six pair monomials (k=6, referee remark);
       mode 'lin4' : the 4 inputs themselves, coordinates G_i prod_{nu!=i} conj G_nu (k=4, referee remark)."""
    P = [p for p in primerange(X, 2*X+1) if p % 4 == 1]
    s = len(P); T = (2*X)**j
    pi = {p: gauss_prime_above(p) for p in P}
    b = Fraction(bnum, bden); delta = Fraction(delta_num, delta_den)
    # angular bins of width T^-b (float only used to SELECT E; all checks below are exact)
    w = float(T)**(-float(b))
    items = []
    for A in itertools.combinations(P, j):
        G = (1, 0)
        for p in A: G = gmul(G, pi[p])
        items.append((math.atan2(G[1], G[0]) * 2 % (2*math.pi), A, G))   # arg(G/conj G) = 2 arg G
    items.sort()
    # best window of angular width w (no wraparound needed for a lower bound)
    angs = [t[0] for t in items]; best = (0, 0); i0 = 0
    for i1 in range(len(items)):
        while angs[i1] - angs[i0] > w: i0 += 1
        if i1 - i0 + 1 > best[1]-best[0]: best = (i0, i1+1)
    E = items[best[0]:best[1]]
    if nmax_E: E = E[:nmax_E]
    NE = len(E)
    k = {'pair5': 5, 'pair6': 6, 'lin4': 4}[mode]
    say(f"[equal-norm {mode}] X={X} s={s} j={j} k={k} delta={delta} b={b} T={T}: C(s,j)={len(items)}, "
        f"angular window {w:.3e} holds |E|={NE}")
    pairs6 = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    pairs = pairs6[:5] if mode == 'pair5' else pairs6
    a, q = delta.numerator, delta.denominator
    good = 0; bad_eq = 0; allpass = True; worst = 1e9; lowok = True; highok = True; eqnorm = True
    for t in itertools.product(range(NE), repeat=4):
        Gs = [E[i][2] for i in t]; As = [E[i][1] for i in t]
        Gb = [gconj(G) for G in Gs]
        if mode == 'lin4':
            ys = []
            for i in range(4):
                y = Gs[i]
                for nu in range(4):
                    if nu != i: y = gmul(y, Gb[nu])
                ys.append(y)
        else:
            ys = []
            for (i, jj) in pairs:
                y = gmul(Gs[i], Gs[jj])
                for nu in range(4):
                    if nu not in (i, jj): y = gmul(y, Gb[nu])
                ys.append(y)
        if len(set(ys)) < len(ys): bad_eq += 1; continue
        g = ys[0]
        for y in ys[1:]: g = ggcd(g, y)
        ys = [gexact(y, g) for y in ys]
        nset = {gnorm(y) for y in ys}; eqnorm &= (len(nset) == 1)
        prs = sorted({p for A in As for p in A})
        L, H = lhs_and_H(ys, prs)
        allpass &= passes(L, H, k, a, q)
        worst = min(worst, delta_eff(L, H, k))
        # H >= T^(2b)/4  <=> (4H)^bden >= T^(2 bnum);  H <= T^4
        lowok &= (4*H)**bden >= Fraction(T)**(2*bnum)
        highok &= H <= Fraction(T)**4
        good += 1
    say(f"   parameter 4-tuples: {NE**4}, with a coordinate coincidence: {bad_eq}, good: {good}")
    say(f"   all good tuples solve (5.1) exactly: {allpass}; equal norms after gcd removal: {eqnorm}; "
        f"min delta_eff={worst:.4f}; H>=T^(2b)/4: {lowok}; H<=T^4: {highok}")
    if NE:
        say(f"   hyperplane lower bound good/(2|E|^3) = {good/(2*NE**3):.2f}")

if __name__ == "__main__":
    rational_instance(1000, 2, 3, 1, 2, 10**6)        # k=3, delta=1/2, b=1/4
    rational_instance(1000, 2, 5, 1, 1, 20000)        # k=5, delta=1,   b=1/4
    rational_instance(1000, 3, 3, 1, 2, 20000)        # j=3
    rational_instance(1000, 3, 5, 3, 2, 20000)        # k=5, delta=3/2 (arc range <2), b=3/8
    equal_norm_instance(1000, 3, 1, 2, 3, 10, 'pair5')   # k=5, delta=1/2, b=0.3
    equal_norm_instance(1000, 3, 1, 2, 3, 10, 'pair6')   # k=6 remark
    equal_norm_instance(1000, 3, 1, 3, 3, 10, 'lin4')    # k=4 remark, delta<6b/4=0.45
    with open(__file__.replace('.py', '_out.txt'), 'w') as f: f.write("\n".join(out) + "\n")
