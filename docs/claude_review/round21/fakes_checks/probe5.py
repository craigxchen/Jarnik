from mpmath import mp, mpf, asin, pi, sqrt, atan2, sin, floor
mp.dps = 80
t0 = atan2(1, 2)
q = pi/2
def viol(t, p=5, nmax=200):
    out = []
    for n in range(1, nmax+1):
        x = n*t
        xm = x - q*floor(x/q)
        dax = min(xm, q - xm); ddg = abs(xm - q/2)
        ea = asin(mpf(p)**(-mpf(n)/2)); ed = asin(mpf(p)**(-mpf(n)/2)/sqrt(2))
        if dax < ea: out.append((n,'ax', mp.nstr((ea-dax)/ea,5)))
        if ddg < ed: out.append((n,'dg', mp.nstr((ed-ddg)/ed,5)))
    return out
for eta in ['1e-3','1e-6','1e-10','1e-20','1e-40']:
    for s in [1,-1]:
        v = viol(t0 + s*mpf(eta))
        print(eta, s, v[:4])
# tight n for genuine
for n in range(1, 40):
    A = (2+1j)**n
    X, Y = int(round(A.real)), int(round(A.imag))
    m = min(abs(X), abs(Y), abs(X+Y), abs(X-Y))
    if m <= 2: print('n',n,'X',X,'Y',Y,'min(|X|,|Y|,|X+-Y|)=',m)
