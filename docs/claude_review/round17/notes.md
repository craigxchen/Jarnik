# round17 cap notes (working; not a report)

## sdc_q definition (referee_barrier_checks/sdc_sym.py)
rows: base inf + F_q. Relative patterns (base omitted):
 P_a  = {a} u {x : chi(a-x)=1}  = {a} u (a - QR) = {a} u (a + NQR)   (q=3 mod 4, -1 NQR)
 P'_a = {a} u {x : chi(a-x)=-1} = {a} u (a + QR)
|P|=(q+1)/2, 2q patterns, each x in q+1 patterns, each pair x!=y in exactly (q+1)/2 common... check: #sep = q+1.
A = {k in {-1,0,1}^q : |k(P)| <= j all P}.

## slab coordinates
S_a := k(P'_a) = k_a + k(a+QR);  S^-_a := k(P_a) = k_a + k(a+NQR); S_a + S^-_a = K + k_a, K = sum k.
Fourier hat f(b) = sum_x f_x zeta^{bx}.  hat S(b) = hat k(b) * (1 + sum_{y in QR} zeta^{-by}).
G(b) := sum_{y in QR} zeta^{by} = (chi(b) g - 1)/2 with g = sum_x chi(x) zeta^x = i sqrt q (q=3 mod 4).
sum_{y in QR} zeta^{-by} = G(-b) = (chi(-b) g - 1)/2 = (-chi(b) i sqrt q - 1)/2.
So hat S(1) = hat k(1) (1 - i sqrt q)/2   => hat k(1) = 2/(1 - i sqrt q) * sum_a zeta^a S_a.
|2/(1-i sqrt q)| = 2/sqrt(q+1).   [verify numerically]
Similarly hat S^-(1) = hat k(1)(1 + G(1)) = hat k(1) (1 + i sqrt q)/2.

## scales at q=3251
d=q=3251, r=6502, w=514.2, j=257, m'=2057, 2j+1=515, |P|=1626, sd of k(P) on cube = sqrt(2*1626/3)=32.9
threshold ln T* = -q tau/2 + ln q - 1.5 ln(m'+1) - ln 2 - pi^2/(4q) (C=1)
tau_D (c0=0.184) ~ 296 -> ln T* ~ -4.8e5 ; tau>= 2 ln(1000 q n^2) ~ 62 anyway -> <= -1.0e5.

## key structural facts found
* products of linear forms (zonotopes) in inverse system: each integer (even real) factor v hits >= q+1 of 2q
  pattern directions -> total directional degree >= m(q+1) vs budget 4jq -> m <= 4j < m'. ratio 2 vs needed 4.
* H_A(m') > 0 by counting (|A| >= 3^q/2 > #{alpha in {0,1,2}^q, |alpha| < m'}), so non-product G exist abstractly.
* Isotypic: l_{u_theta}^{m'} = 2^{-m'} sum_i binom(m',i) e^{i theta(2i-m')} M_i, M_i = l_psi^i l_psibar^{m'-i};
  translation eigenvalue zeta^{(m'-2i)} (or inverse), distinct for i=0..m' since m' < q.
  => Psi_A(u_theta) = 4^{-m'} sum_i binom(m',i)^2 N_i, constant in theta. (Hermitian L2 on A, translation-invariant)
* N_i rational? Galois sigma_s maps psi -> psi_s; box AGL-invariant -> N_i(psi_s) = N_i(psi) (multiplication by s is a box
  automorphism mapping l_psi to l_{psi_s} up to relabel). careful: x -> s x maps P_a -> P_{sa} for s QR, swaps families for NQR.
  Need: N_i(psi_1) defined via rational q_A, sigma_s(N_i(psi_1)) = N_i(psi_s). So N_i in Q.
* slab-only reduction acts on negligible multinomial mass.
* cube-only and S-grid-only lower bounds >> threshold.
