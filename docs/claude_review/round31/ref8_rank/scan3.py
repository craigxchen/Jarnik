"""Dump per-window records (arc constant, chord constant, M, omega, rank) and tabulate defects
by threshold.  Theorem assertions are those of scan2.analyse_window (run on every window)."""
import sys, pickle
from math import atan2, sqrt, pi, log
from sympy import factorint, primerange
import scan2
from scan2 import points_of, newstats, gen_split_N

def windows(fac, CMAX, out, N, minM=3):
    ps, es, pts = points_of(fac)
    q = N ** 0.25
    R = sqrt(N)
    ang = sorted((atan2(z[1], z[0]) % (2 * pi), z, a, zp) for z, a, zp in pts)
    n = len(ang)
    st = newstats()
    for i in range(n):
        if ang[i][0] >= pi / 2:
            break
        j = i + 1
        while j < i + n:
            t = ang[j % n][0] + (2 * pi if j >= n else 0)
            C = R * (t - ang[i][0]) / q
            if C > CMAX:
                break
            M = j - i + 1
            if M >= minM:
                win = [(ang[t % n][1], ang[t % n][2], ang[t % n][3]) for t in range(i, j + 1)]
                M_, om, rk, Cp = scan2.analyse_window(ps, es, win, st, Nval=N)
                out.append((N, M, om, rk, C, Cp))
            j += 1

if __name__ == '__main__':
    mode = sys.argv[1]
    out = []
    if mode == 'all':
        X = int(float(sys.argv[2])); CMAX = float(sys.argv[3])
        for N, fac in gen_split_N(X):
            windows(fac, CMAX, out, N)
    elif mode == 'struct':
        # N <= X built from the first k split primes with many points
        X = int(float(sys.argv[2])); CMAX = float(sys.argv[3]); kp = int(sys.argv[4]); dmin = int(sys.argv[5])
        P1 = [p for p in primerange(5, 10**4) if p % 4 == 1][:kp]
        cnt = 0
        def dfs(idx, cur, fac, d):
            global cnt
            if d >= dmin:
                cnt += 1
                windows(dict(fac), CMAX, out, cur)
            for t in range(idx, len(P1)):
                p = P1[t]
                if cur * p > X: break
                c = cur; e = 0
                while c * p <= X:
                    c *= p; e += 1; fac[p] = e
                    dfs(t + 1, c, fac, d * (e + 1))
                del fac[p]
        dfs(0, 1, {}, 1)
        print('circles scanned', cnt)
    pickle.dump(out, open(sys.argv[-1], 'wb'))
    print('windows', len(out))
