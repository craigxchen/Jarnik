# Sylvester k=12 (n=4096, r=4095) at m'=ceil(4w): (R1) check, and the level g* above which Prop 5.1(ii)
# (with the proved order o >= m' - 6j of the 3-column cheap span) forces MIN*(A) to fail.
import math
from params import j_mprime
from r1_check import m_max
k = 12; n = 2**k; d = n-1; r = n-1
w, j, mp_ = j_mprime(n, r); w = float(w)
mm = m_max(n)
print(f"n={n} r={r} w={w:.4f} j={j} m'={mp_} m_max={mm} (R1) ceil(4w)<=m_max: {mp_ <= mm}  w>=10: {w>=10}")
o = mp_ - 6*j
N2 = n*(n-1)//2; P = 8*((n//8)*(n//8-1)//2)
slope = o*N2/(2*P)
print(f"3-col span: P={P}, N2/P={N2/P:.5f}, proved order o>={o}, o*N2/P={o*N2/P:.1f} vs d={d}")
delta_y = 1/7; a_y = 2*math.log(3/delta_y); Gam_L = a_y + 2*math.log(6*n/delta_y)
# MIN* needs ln(1/Lambda(g)) < (d/2) g + ln((m'+1)/d); Prop 5.1(ii): ln(1/Lambda(g)) >= slope*(g-Gam_L) - o ln(4n) - m' ln 4
const = slope*Gam_L + o*math.log(4*n) + mp_*math.log(4) + math.log((mp_+1)/d)
gstar = const/(slope - d/2)
print(f"Gamma_L={Gam_L:.3f}; MIN* fails for every g > g* = {gstar:.2f}  (i.e. W = 2 d g > {2*d*gstar:.4g});")
print(f"near-extremal level abar ~ 1.5 ln n = {1.5*math.log(n):.2f}; W_min = 3 n ln n - 14 n = {3*n*math.log(n)-14*n:.4g}")
