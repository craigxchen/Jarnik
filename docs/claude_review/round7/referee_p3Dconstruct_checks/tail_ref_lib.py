import math
def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    for p in [3,5,7,11,13,17,19,23,29,31,37]:
        if n % p == 0: return n == p
    d = n-1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in [2,3,5,7,11,13,17,19,23,29,31,37]:
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True
def leg(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q-1)//2, q) == 1 else -1)
def mul(x, y, Q): return ((x[0]*y[0]-x[1]*y[1]) % Q, (x[0]*y[1]+x[1]*y[0]) % Q)
def gpow(x, e, Q):
    r = (1, 0)
    while e:
        if e & 1: r = mul(r, x, Q)
        x = mul(x, x, Q); e >>= 1
    return r
def ratio(z, Q):  # z / conj z mod Q
    n = (z[0]*z[0]+z[1]*z[1]) % Q; ni = pow(n, -1, Q)
    zz = mul(z, z, Q)  # z/conj z = z^2 / Norm z
    return (zz[0]*ni % Q, zz[1]*ni % Q)
