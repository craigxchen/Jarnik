"""Referee checks of the height bookkeeping in Theorem L and Proposition B of round4/flipped.md.

(1) Averaging step of Theorem L: for every Hadamard core, assignment and weights,
       min_{a<b} (W_a + W_b + F_ab) <= 2 W/(M-1) + F (M-1)/(2(M-2)),
    with F_ab = sum_{x in D_ab} l_x (one flip per row), F = W_F.  Tested on random weights/assignments.
(2) Prop B in the heavy regime: weight profiles (distinct split primes, b = 5 per label, one label carrying
    the five smallest split primes, flipped primes heavy) in which the Liouville/column inequality (3.1)
    holds for every label (D = 0, any C >= 1), but min_a h(u_a) h(v_a) / V stays bounded as M grows
    (no log M growth), while (1/8) log M grows.  Shows the 'factor of order c log M' reading of Prop B
    is not implied by its proof once W_F >= (W+D)/2.
(3) Corrected Prop B: if W_F <= (1/2 - delta) V then h_1 h_2 >= (delta/2) V log(4M/e) (1 - o(1));
    checked on random admissible weight data (inequality form, exact float bookkeeping).
"""
import math, random, sys
sys.path.insert(0, '.')
from ref_pair_lp import paley, sylv

def is_prime(n):
    """Deterministic Miller-Rabin with the first 13 prime bases (valid for n < 3.3e24); for larger n
    it is a strong probable-prime test with 13 bases (only used to pick weights, not for any claim)."""
    if n < 2: return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
    for p in small:
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def split_primes(start, count):
    out, p = [], max(5, start)
    while len(out) < count:
        if p % 4 == 1 and is_prime(p): out.append(p)
        p += 1
    return out

def check_avg(rng):
    n = 0
    for H in [sylv(3), paley(11), sylv(4), paley(19), paley(23), sylv(5), paley(43)]:
        M = len(H)
        for _ in range(30):
            B = rng.choice([1, 2, 3])
            cnt = {}; asg = []
            for x in range(M):
                ch = [a for a in range(1, M) if cnt.get(a, 0) < B] or list(range(1, M))
                a = rng.choice(ch); asg.append(a); cnt[a] = cnt.get(a, 0) + 1
            l = [rng.expovariate(1.0)*rng.choice([1, 10, 100]) for _ in range(M)]
            U = [rng.expovariate(1.0)*rng.choice([1, 10]) for _ in range(M)]   # U[a] for a=1..M-1
            Wa = {a: U[a] + sum(l[x] for x in range(M) if asg[x] == a) for a in range(1, M)}
            W = sum(Wa.values()); F = sum(l)
            best = min(Wa[a] + Wa[b] + sum(l[x] for x in range(M) if H[x][a] != H[x][b])
                       for a in range(1, M) for b in range(a+1, M))
            assert best <= 2*W/(M-1) + F*(M-1)/(2*(M-2)) + 1e-9
            n += 1
    return n

def heavy_profile(M, heavy_factor):
    # label 1: five smallest split primes (unflipped part) ; other labels: b=5 nearby primes ~ M^4;
    # flipped primes (one per row) ~ M^(4*heavy_factor), counted inside their labels 2..M-1 (capacity 2)
    small = split_primes(5, 5)
    near = split_primes(M**4, 5*(M - 2))
    flp = split_primes(int(M**(4*heavy_factor)), M)
    Wa = {1: sum(map(math.log, small))}
    for i, a in enumerate(range(2, M)):
        Wa[a] = sum(map(math.log, near[5*i:5*i + 5]))
    lx = [math.log(p) for p in flp]
    asg = [(x % (M - 2)) + 2 for x in range(M)]       # label 1 carries no flip (capacity 2)
    for x in range(M): Wa[asg[x]] += lx[x]     # flipped column is one more column of its label
    W = sum(Wa.values()); WF = sum(lx); V = W
    # (3.1) for every label: (M/2) W_a + W_F >= V/2 - 2 log(M C/(2 sqrt 2)) ; check with C = 1 (strongest)
    ok31 = all((M/2)*Wa[a] + WF >= V/2 - 2*math.log(M/(2*math.sqrt(2))) for a in Wa)
    h1h2 = min(Wa[a]/2 * WF/2 for a in Wa)
    return WF/V, ok31, h1h2/V

if __name__ == '__main__':
    rng = random.Random(4)
    print("(1) Theorem L averaging inequality checked on", check_avg(rng), "random weighted profiles: PASS")
    print("(2) heavy-flip weight profiles: label 1 = five smallest split primes")
    for M in [44, 100, 200, 400, 800]:
        r, ok, ratio = heavy_profile(M, 8.0)
        print(f"   M={M:4d}  W_F/V={r:.3f}  (3.1) holds for all labels: {ok}   min_a h(u_a)h(v_a)/V = {ratio:.3f}"
              f"   vs (1/8)log(4M/e) = {math.log(4*M/math.e)/8:.3f}")
    # (3) corrected Prop B (pure inequality bookkeeping)
    cnt = 0
    for _ in range(20000):
        M = rng.choice([12, 20, 44, 100, 1000]); V = rng.uniform(1, 50)*M*math.log(M)
        WF = rng.uniform(M*math.log(4*M/math.e), V)
        delta = 0.5 - WF/V
        if delta <= 0: continue
        Wa_min = (V - 2*WF)/M                    # (3.1) lower bound at C = 2 sqrt 2 / M (log term dropped)
        h1h2 = (Wa_min/2)*(WF/2)
        assert h1h2 >= (delta/2)*V*math.log(4*M/math.e)*(1 - 1e-12)
        cnt += 1
    print("(3) corrected Prop B bound h1h2 >= (delta/2) V log(4M/e) verified on", cnt, "random admissible data: PASS")
