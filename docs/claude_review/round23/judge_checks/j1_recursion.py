# Judge J1 (independent code): exact recursive peeling for reg T_n(a,b) bound, n=2..5,
# checks (a) rho_n == F_n on [0,N]^2, (b) rho_5(m,m) == floor(10m/3),
# (c) peeling of slice bounds g(c) on S_m sliced by x_1 for two mean-width slice bounds.
from fractions import Fraction as Fr
import sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 70
def H(n,x,y):
    lo, hi = min(x,y), max(x,y)
    if n == 2: return 0
    if n == 3: return lo
    if n == 4: return min((x+y)//2, 2*lo)
    if n == 5: return min(lo + hi//3, 3*lo)
def F(n,x,y): return x+y+H(n,x,y)
def peel(vals):
    v = sorted(vals, reverse=True)
    return max(v[i] + i for i in range(len(v)))
rho = {2: {(a,b): a+b for a in range(N+1) for b in range(N+1)}}
for n in (3,4,5):
    r = {}
    for a in range(N+1):
        for b in range(N+1):
            V = [rho[n-1][(ap,b)] for ap in range(a+1)] + [rho[n-1][(a,bp)] for bp in range(b)]
            r[(a,b)] = peel(V)
    rho[n] = r
bad = [(n,a,b) for n in (3,4,5) for a in range(N+1) for b in range(N+1) if rho[n][(a,b)] != F(n,a,b)]
print("N=%d  mismatches rho_n vs F_n (n=3,4,5): %d" % (N, len(bad)), bad[:5])
print("rho_5(m,m)==floor(10m/3) for m<=N:", all(rho[5][(m,m)] == (10*m)//3 for m in range(1,N+1)))
print("rho_4(m,m)==3m:", all(rho[4][(m,m)] == 3*m for m in range(1,N+1)))
# (c) slice S_m=T_5(m,m) by x_1=c, c in [-m,m]; slice c is T_4(m-c,m) (c>=0) or T_4(m,m-|c|) (c<0);
# width-sum (pair functionals of 4 coords) of slice: each width = a'+b' = 2m-|c|.
def peel_frac(vals):
    v = sorted(vals, reverse=True)
    return max(v[i] + i for i in range(len(v)))
for m in (1,2,3,4,6,8,12,24,48):
    cs = range(-m, m+1)
    cbal = peel_frac([ (3*(2*m-abs(c)))//2 for c in cs])        # C_bal(4) bound, floored
    thmB = peel_frac([ (7*(2*m-abs(c)))//4 for c in cs])        # Theorem B(4) bound, floored
    true_ = peel_frac([ F(4, m-c, m) if c >= 0 else F(4, m, m-abs(c)) for c in cs])
    print("m=%2d  peel[C_bal(4)]=%4d (7m/2=%s)  peel[ThmB(4)]=%4d (15m/4=%s)  peel[true F_4]=%4d (10m/3=%s)" %
          (m, cbal, Fr(7*m,2), thmB, Fr(15*m,4), true_, Fr(10*m,3)))
