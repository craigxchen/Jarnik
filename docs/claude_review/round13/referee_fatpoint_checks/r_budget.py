import itertools, numpy as np
d = 7
Pm = np.array([[1 if bin(a & (x + 1)).count('1') % 2 else 0 for x in range(d)] for a in range(1, 8)])
pts = np.array(list(itertools.product((-1, 0, 1), repeat=d)))
for L in (2, 3):
    A = pts[(np.abs(pts @ Pm.T) <= L).all(1)]
    tau0 = len(np.unique(A.sum(1))) - 1
    Tc = sum(len(np.unique(A @ Pm[i])) - 1 for i in range(7))
    mx = max(m for m in range(1, 40) if m * m <= (2 * d + min(tau0, m)) * Tc / 8)
    print('Syl(8,%d): tau0=%d T_core=%d  budget allows m <= %d' % (L, tau0, Tc, mx))
