from common import *
import sys
rows = []
print("n r d w j m' R1(m'<=m_max) W_min abar ln(1/eps0) m'ln(2d) CW_log10(bound) inv_err_lower MNV_range")
for n, rname in [(2900,'n'),(4096,'n'),(5809,'n^2'),(10**4,'n'),(10**4,'n^2'),(65536,'n'),(10**5,'n'),(10**6,'n'),(10**9,'n'),(10**12,'n')]:
    r = {'n': n - 1 if n in (4096, 65536) else n, 'n^2': n * n}[rname]
    d, w, j, m = wjm(n, r)
    ok = mmax_ok(n, m) if n <= 10**6 else 'n/a'
    W = Wmin(n)
    abar = W / (2 * d)
    le = ln_inv_eps0(n, m, W)
    # CW (optimistic constant 1): Pr <= m * (d^m eps0 / sqrt(m!))^{1/m}
    lcw = (log(m) + (m * log(d) - le - mp.loggamma(m + 1) / 2) / m) / log(10)
    # MOO invariance error lower bound m * tau^{1/(8m)}, tau >= 1/d
    inv = m * exp(-log(d) / (8 * m))
    mnv = log(r) / log(log(r))
    print(n, rname, d, mp.nstr(w, 6), j, m, ok, mp.nstr(W, 8), mp.nstr(abar, 5), mp.nstr(le, 8),
          mp.nstr(m * log(2 * d), 8), mp.nstr(lcw, 4), mp.nstr(inv, 5), mp.nstr(mnv, 4))
