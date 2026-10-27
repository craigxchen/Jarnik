# Certified (interval-partition) lower bound for G(p) on (0,1) with lambda=min(p/2,1-p); and the author's lambda=p/2 version near p=1.
from mpmath import mp, mpf, log, binomial
mp.dps=30
q=3251; j=257; m=2057; P=(q+1)//2
lb=log(binomial(P,j+1))/2
def Glow(p1,p2):
    # lower bound of G on [p1,p2]: ln(1/lambda) >= ln(1/max lambda); tail term >= (j+1)/2 ln(p1/2) + (P-j-1)/2 ln(1-p2)
    lam_max = max(min(p1/2,1-p1), min(p2/2,1-p2)) if not (p1<=mpf(2)/3<=p2) else mpf(1)/3
    return m*log(3)+(m/mpf(2))*log(1/lam_max)+lb+(j+1)/mpf(2)*log(p1/2)+(P-j-1)/mpf(2)*log(1-p2)
# partition: (0,1e-30] handled analytically; [1e-30, 1-1e-30] by geometric grids
worst=None
pts=[mpf(10)**(-30+ k*mpf(30)/20000) for k in range(20001)]  # 1e-30 .. 1
pts=[p for p in pts if p<=mpf('0.5')]+[mpf('0.5')]
pts2=[1-mpf(10)**(-mpf(k)*30/20000) for k in range(1,20001)]  # 0.5.. 1-1e-30 approx (reversed)
pts2=sorted(p for p in pts2 if p>mpf('0.5'))
grid=pts+pts2
for a,b in zip(grid[:-1],grid[1:]):
    g=Glow(a,b)
    if worst is None or g<worst[0]: worst=(g,a,b)
print('certified lower bound of G on [1e-30, 1-1e-30]:',worst[0],'on',worst[1],worst[2])
# ends: p<=1e-30: G >= m ln3 + (m/2)ln(2/p) + lb + (j+1)/2 ln(p/2) + (P-j-1)/2 ln(1-p); since m > j+1 the p-terms give ((m-j-1)/2) ln(2/p) >0
print('p<=1e-30: G >= m ln3 + lb + ((m-j-1)/2) ln(2e30) + (P-j-1)/2 ln(1-1e-30) =', m*log(3)+lb+((m-j-1)/mpf(2))*log(2*mpf(10)**30)+(P-j-1)/mpf(2)*log(1-mpf(10)**-30))
# p>=1-1e-30: lambda=1-p: G >= m ln3 + lb + (j+1)/2 ln(1/4)... + ((m - (P-j-1))/2) ln(1/(1-p)) with m > P-j-1
print('m=',m,' P-j-1=',P-j-1,' so (m/2)ln(1/(1-p)) dominates (P-j-1)/2 ln(1-p) near 1; value at 1-1e-30 lower bound:', m*log(3)+lb+(j+1)/mpf(2)*log(mpf(1)/4)+((m-(P-j-1))/mpf(2))*log(mpf(10)**30))
# author's formula lambda=p/2 at p near 1:
for p in [mpf('0.891'),mpf('0.95'),mpf('0.99'),mpf('0.999'),1-mpf(10)**-6]:
    G=m*log(3)+(m/mpf(2))*log(2/p)+log(binomial(P,j+1)*(p/2)**(j+1)*(1-p)**(P-j-1))/2
    print('author formula p=',p,'G=',G)
