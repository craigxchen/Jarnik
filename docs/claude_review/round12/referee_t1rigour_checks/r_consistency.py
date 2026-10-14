# Referee: Prop 5.5 scale consistency, author's o-bound vs the bound actually proved.
# EVIDENCE (double precision). r = n, C = 1, m' = ceil(4w).
import math
def params(n):
    d = n-1; w = math.sqrt(8*d*math.log(4*n)); j = int(w//2); m = math.ceil(4*w); b = (m-1)//(2*j)
    return d, w, j, m, b
def G(n, W, mode):
    d, w, j, m, b = params(n)
    abar = W/(2*(n-1))
    worst = math.inf; coef = math.inf
    for s in range(0, b+1):
        Emax = int(math.floor((m-1)/2 - j*s))
        for E in range(0 if s >= 1 else 1, Emax+1):
            q = 2**s
            if n - E - q <= 0: continue
            Q = q*n*(n-1)/((n-E)*(n-E-q))
            if mode == 'author':
                o = m - 2*j*s - 2*E
            else:  # proved: rho_A(F) >= 2 min(j, n - largest class); >= 2j if s>=1, >= 2 min(j,E) if s=0
                o = m - 2*j if s >= 1 else m - 2*min(j, E)
            if o <= 0: continue
            val = m*math.log(3) - o*(math.log(2) + Q*abar/2) - math.log(d/(m+1)) + W/4
            worst = min(worst, val)
            coef = min(coef, 0.25 - o*Q/(4*(n-1)))
    return worst, coef
for mode in ('author', 'proved'):
    last_inf = None; last_W = None
    n = 1000
    step = 1
    while n <= 600000:
        Wmin = 3*n*math.log(n) - 14*n
        g, c = G(n, Wmin, mode)
        if c <= 0: last_inf = n
        if g <= 0: last_W = n
        n += 1 if n < 20000 else 50
    print(mode, ": last n with W-coefficient <= 0:", last_inf, "; last n with G(3n ln n - 14n) <= 0:", last_W)
for n in (5185, 10000, 50000, 100000, 200000, 300000):
    print(n, "proved-bound G_coef = %.4f" % G(n, 1e9, 'proved')[1], " author G_coef = %.4f" % G(n, 1e9, 'author')[1])
