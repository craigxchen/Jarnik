# Independent scans for Cor 3.7 (sep vs radius), Prop 5.3 (n_A, n_B), at r = n, n ln n, n^2.
import math, numpy as np, sys
def params(n, rkind):
    d = n-1
    r = {'n': n, 'nlnn': n*math.log(n), 'n2': n*n}[rkind]
    w = math.sqrt(8*d*math.log(4*r))
    j = math.floor(w/2); mp_ = math.ceil(4*w)
    return d, w, j, mp_
def worst(n, rkind, which):
    d, w, j, mp_ = params(n, rkind)
    N2 = n*(n-1)/2
    b = (mp_-1)//(2*j)
    W = 3*n*math.log(n) - 14*n
    abar = W/(2*(n-1))
    rhs = W/4 + mp_*math.log(4) - math.log(d/(mp_+1))
    worstA = -1e300; worstB = -1e300; worstC = -1e300
    for s in range(0, b+1):
        emax = (mp_-1)//2 - j*s
        if emax < 0: continue
        e = np.arange(0, emax+1, dtype=float)
        den = (n-e)*(n-e-2**s)
        if np.any(den <= 0): return None
        Q = 2**(s+1)*N2/den
        o = mp_ - 2*j*s - 2*e
        worstA = max(worstA, float(np.max(mp_*Q)) - d)
        worstB = max(worstB, float(np.max(o*Q)) - d)
        worstC = max(worstC, float(np.max(o*(Q*abar/2 + math.log(2)))) - rhs)
    return worstA, worstB, worstC, b, w
if __name__ == "__main__":
    for rk in ['n','nlnn','n2']:
        lastB = lastC = None
        for n in range(1000, 20001):
            res = worst(n, rk, None)
            if res is None: lastB = lastC = n; continue
            A,B,C,b,w = res
            if B >= 0: lastB = n      # n_B uses <= d ; record strict failure B>0 separately
            if not (C < 0 and B < 0): lastC = n
        print(f"r={rk}: last n<=20000 failing (B>=0): {lastB}; last failing Cor3.7 (C<0 and B<0): {lastC}")
        sys.stdout.flush()
