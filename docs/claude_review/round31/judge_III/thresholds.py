from fractions import Fraction as F
from math import comb
def FM(M): return F((M-1)*(M-4),8*(M-2)*(M-3))
def astra_max(D):  # largest M with (III.7) failing => max marks = smallest passing M - 1
    D=F(D)
    for M in range(5,2000):
        if (D-1)*FM(M) > F(3)*D/(2*(M-1)):
            # check monotone afterwards
            assert all((D-1)*FM(m) > F(3)*D/(2*(m-1)) for m in range(M,M+400))
            return M-1
def rho(M): return F(M*M//4, comb(M,2))
def pair_max(D):
    D=F(D)
    for M in range(2,5000):
        if D > 4*rho(M): return M-1
    return None
for D in ['11/10','3/2','2','21/10','11/5','5/2','3','4','5','10','100','1000000']:
    print(D, 'quartet(III.7) max', astra_max(F(D)), ' pair-count max', pair_max(F(D)))
print(float(3*FM(18)), float(F(3*4,2*17)), float(3*FM(17)), float(F(12,32)))
