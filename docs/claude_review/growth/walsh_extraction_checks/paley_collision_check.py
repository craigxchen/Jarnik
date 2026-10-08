"""Literal fully-flipped capacity-two Paley profiles on actual nearby split primes:
compare, pair by pair, the small-prime part of the reduced chord norm with the available pair slack.
  slack_xy  = d_xy - (W+D)/2 + 2 log C      (C=1, D=0 here: primitive rows)
  coll_xy   = sum_{q<=Qmax} v_q(Norm(u_x-u_y)) log q,   u = z/gcd(z_x,z_y)
Necessary for an arc of constant C: coll_xy <= log Norm(u_x-u_y) <= slack_xy.  Also the aggregate form."""
import math, random, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gauss_common import *

def paley(q):
    sq = set((a * a) % q for a in range(1, q))
    ch = lambda a: 0 if a % q == 0 else (1 if (a % q) in sq else -1)
    M = q + 1
    H = [[1] * M for _ in range(M)]
    for x in range(1, M):
        for j in range(1, M):
            H[x][j] = -1 if x == j else -ch(x - j)
    return H

def nearby_split_primes(r, X):
    # densest window of relative log-width 1/r among split primes in [X, 2X)
    ps = [p for p in range(X | 1, 2 * X, 2) if p % 4 == 1 and is_prime(p)]
    best = None; j = 0
    for i in range(len(ps)):
        while j < len(ps) and math.log(ps[j]) - math.log(ps[i]) < 1.0 / r: j += 1
        if j - i >= r and (best is None or ps[i] < best[0]):
            best = (ps[i], ps[i:i + r]); break
    return best[1] if best else None

def smallprime_part(n, Qmax, plist):
    s = 0.0
    for q in plist:
        if q > Qmax: break
        while n % q == 0:
            n //= q; s += math.log(q)
    return s

random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
q = int(sys.argv[1]); b = int(sys.argv[3]) if len(sys.argv) > 3 else 5
H = paley(q); M = q + 1
labels = list(range(1, M))
r = b * (M - 1)
X = max(10 ** 5, M ** 4)
primes = None
while primes is None:
    primes = nearby_split_primes(r, X); X *= 2
cols = []; lab = []
for a in labels:
    for c in range(b): cols.append([H[x][a] for x in range(M)]); lab.append(a)
rows = list(range(M)); random.shuffle(rows)
for idx, x in enumerate(rows):          # capacity two: label idx mod (M-1), copy idx//(M-1)
    a = labels[idx % (M - 1)]
    j = (a - 1) * b + idx // (M - 1)
    cols[j][x] *= -1
pis = [gauss_prime(p) for p in primes]
orient = [random.choice((1, -1)) for _ in range(r)]
z = []
for x in range(M):
    v = random.choice(UNITS)
    for j in range(r):
        s = cols[j][x] * orient[j]
        v = gmul(v, pis[j] if s == 1 else gconj(pis[j]))
    z.append(v)
w = [math.log(p) for p in primes]; W = sum(w)
plist = [p for p in range(2, M + 1) if is_prime(p)]
worst = -1e9; aggL = 0.0; aggR = 0.0; viol = 0; npairs = 0; maxcoll = 0
for x in range(M):
    for y in range(x + 1, M):
        d = sum(w[j] for j in range(r) if cols[j][x] != cols[j][y])
        slack = d - W / 2
        g = ggcd(z[x], z[y])
        ux = gdivexact(z[x], g); uy = gdivexact(z[y], g)
        diff = (ux[0] - uy[0], ux[1] - uy[1])
        nrm = gnorm(diff)
        # consistency: log|g|^2 = W - d
        assert abs(math.log(gnorm(g)) - (W - d)) < 1e-6 * W
        coll = smallprime_part(nrm, M, plist)
        maxcoll = max(maxcoll, coll)
        worst = max(worst, coll - slack)
        viol += coll > slack
        aggL += coll; aggR += slack; npairs += 1
print("Paley order M=%d, b=%d, r=%d primes in [%d,%d], W=%.1f (W/(M log M)=%.2f)" % (M, b, r, primes[0], primes[-1], W, W / (M * math.log(M))))
print("  pairs=%d  max small-prime collision cost (q<=M)=%.2f  max(coll-slack)=%.2f  violations=%d" % (npairs, maxcoll, worst, viol))
print("  aggregate: sum coll=%.1f  sum slack=%.1f  (M^2/2) log M=%.1f" % (aggL, aggR, M * M / 2 * math.log(M)))
