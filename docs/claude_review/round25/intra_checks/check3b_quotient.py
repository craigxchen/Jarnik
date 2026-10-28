"""Check 3b (evidence): supercriticality scan modulo the multiplicative relations, record 10-point cluster.
The cluster point Z lies on the subtorus T' = {x^k = 1 : k in Lambda_rel}, Lambda_rel = {k : <k,a_p> = 0 for all p}
(check3: Z^k = 1 exactly for the 50 short relation vectors). Characters of T' are l(k) = (<k,a_p>)_p in Z^r.
For boxes S = {l : |l_p| <= b_p} (l in the image of k' in {-1,0,1}^9): reg(S) (mod-p rank, an upper bound for reg_Q),
theta(l) = <k,delta> mod 2pi (checked consistent across k with equal l), L(S) = -log max_{l in S}|theta(l)|,
RHS(S) = (1/2) sum_p width_p(S) log p, and E(S) = reg(S) L(S) + log reg(S)! - RHS(S).
For G supported on S with ord m <= reg(S) and |c|_1 >= 1, the Taylor bound log|G(Z)| <= -m L - log m! + log|c|_1 is
compared with -RHS; a forced non-monomial relation needs E > 0 (necessary, not sufficient; evidence only)."""
import math, itertools, sys
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


def reg5(S):
    nv = len(S[0]); t = 0
    while True:
        E = [[math.prod(l[x]**a for x, a in enumerate(al)) for al in monos(nv, t+1)] for l in S]
        if rank_mod(E) == len(S): return t
        t += 1

def main():
    d, Z, Nn, g0 = normalised_cluster(1176852625, 10); M = 10; prof = profile(Z, Nn)
    th = [math.atan2(z[1], z[0]) for z in Z]; o = 0
    delta = [th[x]-th[o] for x in range(M)]
    tab = {}; incons = 0
    for kp in itertools.product((-1, 0, 1), repeat=M-1):
        k = [0]*M; k[1:] = list(kp); k[0] = -sum(kp)
        l = tuple(sum(kk*ax for kk, ax in zip(k, a)) for p, e, pi, a in prof)
        t = math.remainder(sum(kk*dd for kk, dd in zip(k, delta)), 2*math.pi)
        if l in tab:
            if abs(math.remainder(tab[l]-t, 2*math.pi)) > 1e-9: incons += 1
            if abs(t) < abs(tab[l]): tab[l] = t
        else: tab[l] = t
    print(f"characters of T' reached: {len(tab)}; theta(l) consistent mod 2pi across k: {incons == 0}")
    primes = [p for p, e, pi, a in prof]; r = len(primes)
    rows = []; S_of = {}
    for b in itertools.product(range(3), repeat=r):
        S = [l for l in tab if all(abs(l[j]) <= b[j] for j in range(r))]
        if len(S) < 2 or len(S) > 150: continue
        m = reg5(S)
        L = -math.log(max(abs(tab[l]) for l in S if any(l)))
        rhs = 0.5*sum((max(l[j] for l in S)-min(l[j] for l in S))*math.log(primes[j]) for j in range(r))
        E = m*L + math.lgamma(m+1) - rhs
        rows.append((E, b, len(S), m, L, rhs)); S_of[b] = S
    rows.sort(reverse=True)
    print(f"primes {primes}; boxes scanned {len(rows)}")
    for x in rows[:6]:
        print(f"  S={sorted(S_of[x[1]])} E={x[0]:.6e} b={x[1]} |S|={x[2]} reg<={x[3]} L={x[4]:.3f} RHS={x[5]:.3f}")
    print(f"  max E = {rows[0][0]:.3f}; any E > 0: {any(x[0] > 0 for x in rows)}")

main()
