# Froeberg-type heuristic top degree of Q[k]/(k_x^3 (d of them), r generic forms of degree 2L+1)
# (1+t+t^2)^d (1-t^{2L+1})^r truncated before first nonpositive coefficient. Heuristic only.
import sys
def series(d, r, L, N):
    a = [0]*(N+1); a[0] = 1
    for _ in range(d):
        b = [0]*(N+1)
        for i in range(N+1):
            if a[i]:
                for j in range(3):
                    if i+j <= N: b[i+j] += a[i]
        a = b
    D = 2*L+1
    for _ in range(r):
        b = a[:]
        for i in range(D, N+1): b[i] -= a[i-D]
        a = b
    return a
def topdeg(d, r, L):
    a = series(d, r, L, 2*d)
    for i, x in enumerate(a):
        if x <= 0: return i-1
    return 2*d
for L in (1, 2, 3):
    print('L=%d: need rho >= 8L = %d' % (L, 8*L))
    for d in range(7, 64):
        r = d  # r = n-1 columns (Hadamard-like), d = n-1
        t = topdeg(d, r, L)
        t2 = topdeg(d, (d+1)//2, L)
        print('  d=%d  r=d: %d   r=d/2: %d' % (d, t, t2))
