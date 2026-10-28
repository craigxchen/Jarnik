"""Independent recount of section 4.3 / 6.1: large-prime pair and multi-point contacts on the arc clusters."""
import math, itertools
from sympy import factorint, primerange
from kgauss import *
from kclusters import configs, normalise

def run():
    lines = ['%-22s %2s %8s %6s %6s %8s %8s %7s %s' % ('cluster', 'M', 'gam>M', 'pairs', 'trip', 'pexp', 'texp', 'depth2', 'details')]
    tot_p = tot_t = 0; tot_pe = tot_te = 0.0
    for name, n0, pts, C, isarc in configs():
        if not isarc:
            continue
        q, g = normalise(pts)
        M = len(q); N = gnorm(q[0]); fac = factorint(N)
        pairs = list(itertools.combinations(range(M), 2))
        cn = {}
        primes = set()
        for i, j in pairs:
            d = gsub(q[i], q[j]); gij = ggcd(q[i], q[j]); c = gdiv_exact(d, gij)
            cn[(i, j)] = gnorm(c)
            for r in factorint(gnorm(c)):
                primes.add(r)
        big = sorted(r for r in primes if r > M and r != 2 and r not in fac)
        npairs = ntrip = depth2 = 0
        gam = 0.0
        det = []
        for r in big:
            # classes mod r (for split r with equal norm prime to r, z_i = z_j mod varpi iff mod r)
            cls = {}
            for j in range(M):
                cls.setdefault((q[j][0] % r, q[j][1] % r), []).append(j)
            for k, v in cls.items():
                if len(v) >= 2:
                    npairs += len(v)*(len(v)-1)//2
                    det.append('%d:%s' % (r, v))
                if len(v) >= 3:
                    ntrip += math.comb(len(v), 3)
            for i, j in pairs:
                d = gsub(q[i], q[j])
                if d[0] % (r*r) == 0 and d[1] % (r*r) == 0:
                    depth2 += 1
            # gamma_{>M} contribution: sum over pairs of log of r-part of Norm c_ij
            for i, j in pairs:
                n = cn[(i, j)]
                while n % r == 0:
                    n //= r; gam += math.log(r)
        # random-model expectations
        pe = te = 0.0
        absc = {k: math.sqrt(v) for k, v in cn.items()}
        cmax = max(absc.values())
        for r in primerange(M+1, int(cmax)+2):
            if r in fac:
                continue
            m = r - 1 if r % 4 == 1 else r + 1
            for k in pairs:
                if absc[k] >= r:
                    pe += 1.0/m
            for t in itertools.combinations(range(M), 3):
                if all(absc[(a, b)] >= r for a, b in itertools.combinations(t, 2)):
                    te += 1.0/m**2
        tot_p += npairs; tot_t += ntrip; tot_pe += pe; tot_te += te
        lines.append('%-22s %2d %8.2f %6d %6d %8.2f %8.3f %7d %s' % (name, M, gam, npairs, ntrip, pe, te, depth2, ' '.join(det)))
    lines.append('TOTAL pairs %d (exp %.2f), triples %d (exp %.3f)' % (tot_p, tot_pe, tot_t, tot_te))
    return lines

if __name__ == '__main__':
    out = run()
    open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round20/kill_arakelov_checks/kmulti_out.txt', 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))
