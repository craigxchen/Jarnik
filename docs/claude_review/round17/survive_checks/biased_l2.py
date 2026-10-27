# Method check: L2 restriction with a p-biased product measure nu_p (P(k_x=+-1)=p/2) conditioned on A.
# The chain E[h^2 1_A] >= E h^2 - ||h||_4^2 nu(A^c)^{1/2} with the standard (2,4)-hypercontractive bound
# ||h||_4 <= (sqrt3 * lambda^{-1/4})^{m'} ||h||_2, lambda = p/2 (min atom), closes only if
#   G(p) := m' ln 3 + (m'/2) ln(2/p) + (1/2) ln nu_p(A^c) < -ln 2.
# We use an exact LOWER bound for nu_p(A^c): one pattern, exactly j+1 nonzeros in P, all +1:
#   nu_p(A^c) >= binom(|P|, j+1) (p/2)^{j+1} (1-p)^{|P|-j-1}.  If G(p) > 0 with this lower bound, the chain cannot close.
from mpmath import mp, mpf, log, binomial
mp.dps=40
q=3251; j=257; m=2057; P=(q+1)//2
best=None
for e in range(1,400):
    p=mpf(10)**(-mpf(e)/20)
    tail_low=binomial(P,j+1)*(p/2)**(j+1)*(1-p)**(P-j-1)
    G=m*log(3)+(m/mpf(2))*log(2/p)+log(tail_low)/2
    if best is None or G<best[0]: best=(G,p)
print(f"q={q} j={j} m'={m}: min over p in [1e-20,1] of G(p) (with exact tail LOWER bound) = {float(best[0]):.1f} at p={float(best[1]):.3g}; chain needs G < -ln2")
