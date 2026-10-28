"""Independent recomputation of the place-by-place slack table (section 4.2) and of Lobatto diameters."""
import math, itertools
from fractions import Fraction
from sympy import factorint, primerange
import numpy as np
from numpy.polynomial import legendre as Lg
from kgauss import *
from kclusters import configs, normalise

def E(n, Q):
    # min number of colliding pairs placing n items in Q classes
    if Q >= n:
        return 0
    b, r = divmod(n, Q)
    return r*math.comb(b+1, 2) + (Q-r)*math.comb(b, 2)

def lobatto_delta(M):
    if M == 2:
        return 2.0
    # Fekete points of [-1,1]: +-1 and zeros of P'_{M-1}
    c = [0]*(M-1) + [1]
    d = Lg.legder(c)
    x = np.sort(np.concatenate(([-1.0, 1.0], Lg.legroots(d))))
    s = sum(math.log(abs(x[i]-x[j])) for i in range(M) for j in range(i+1, M))
    return math.exp(s/math.comb(M, 2))

def pigeon_eps(q, M, inert):
    tot = 0; a = 1
    while True:
        m = (q+1)*q**(a-1) if inert else (q-1)*q**(a-1)
        e = E(M, m)
        if e == 0:
            break
        tot += e; a += 1
    return 2*math.log(q)*tot  # both places (split) or Norm q^2 (inert)

def psi_free(p, e, M):
    # min over level vectors (n_0..n_e), sum M, of Q_p(n) + 2 sum_t sum_a E(n_t,(p-1)p^{a-1}); exhaustive for small M,e
    best = None
    for comp in itertools.combinations(range(M+e), e):
        # stars and bars
        prev = -1; n = []
        for c in comp:
            n.append(c-prev-1); prev = c
        n.append(M+e-1-prev)
        S = 0; Qv = Fraction(0)
        for t in range(e):
            S += n[t]
            Qv += (Fraction(S) - Fraction(M, 2))**2
        coll = 0
        for nt in n:
            a = 1
            while True:
                ee = E(nt, (p-1)*p**(a-1))
                if ee == 0:
                    break
                coll += ee; a += 1
        val = Qv + 2*coll
        if best is None or val < best:
            best = val
    return best

def run():
    lines = ['%-22s %2s %6s %8s %8s %8s %8s %8s %8s %8s' % ('cluster', 'M', 'W', 'prodf', 'slack', 'ex_inf', 'fekgain', 'ex_q<=M', 'ex_p|N', 'ex_q>M')]
    for name, n0, pts, C, isarc in configs():
        if not isarc:
            continue
        q, g = normalise(pts)
        M = len(q); N = gnorm(q[0]); fac = factorint(N); K = math.comb(M, 2)
        diffs = {(i, j): gsub(q[i], q[j]) for i, j in itertools.combinations(range(M), 2)}
        act_inf = -sum(math.log(gnorm(d)) for d in diffs.values())
        L = max(math.sqrt(gnorm(d)) for d in diffs.values())  # chord between extreme points
        dM = lobatto_delta(M)
        eps_inf = -2*K*math.log(L/2*dM)
        fek = 2*K*math.log(2/dM)
        # finite places: all primes dividing some Norm(z_i-z_j)
        primes = set([2])
        for d in diffs.values():
            primes |= set(factorint(gnorm(d)).keys())
        primes |= set(r for r in primerange(3, M+3))
        act = {}; eps = {}
        for r in primes:
            a_r = 0.0
            for d in diffs.values():
                n = gnorm(d); v = 0
                while n % r == 0:
                    n //= r; v += 1
                a_r += v*math.log(r)
            act[r] = a_r
            if r == 2:
                eps[r] = K*math.log(2)
            elif r in fac:
                eps[r] = math.log(r)*(fac[r]*(M*M/4 - M/2) + float(psi_free(r, fac[r], M)))
            else:
                eps[r] = pigeon_eps(r, M, r % 4 == 3)
        prodf = act_inf + sum(act.values())
        exq_small = sum(act[r]-eps[r] for r in primes if r not in fac and (r <= M or r == 2))
        exq_big = sum(act[r]-eps[r] for r in primes if r not in fac and r > M and r != 2)
        exN = sum(act[r]-eps[r] for r in primes if r in fac)
        exinf = act_inf - eps_inf
        slack = exinf + exq_small + exq_big + exN
        lines.append('%-22s %2d %6.2f %8.1e %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f' % (name, M, math.log(N), prodf, slack, exinf, fek, exq_small, exN, exq_big))
    lines.append('lobatto delta_M: ' + ' '.join('%d:%.5f' % (m, lobatto_delta(m)) for m in (8, 10, 50, 200)))
    return lines

if __name__ == '__main__':
    out = run()
    open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round20/kill_arakelov_checks/kadelic_out.txt', 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))
