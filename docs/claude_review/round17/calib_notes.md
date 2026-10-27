# calib notes

## N1. Setup (from flib.py / referee sdc_sym.py)
sdc_q: coordinates x in F_q (base inf omitted), patterns P_a = {a} u {x: chi(a-x)=1}, P'_a = {a} u {x: chi(a-x)=-1},
2q patterns of size (q+1)/2. Box A_L = {k in {-1,0,1}^q : |k(P)| <= L all P}. AGL(1,q)-invariant.
khat_b(k) = sum_x zeta^{bx} k_x, translation (gk)_t = k_{t-1}: khat_b(gk) = zeta^b khat_b(k).
u_theta = Re(e^{i theta} psi_1): <k,u_theta> = (e^{i th} khat_1 + e^{-i th} khat_{-1})/2.
l^m = 2^{-m} sum_i binom(m,i) e^{i th (2i-m)} F_i,  F_i = khat_1^i khat_{-1}^{m-i}, F_i in component c_i = 2i-m mod q.
m < q, q odd => c_i distinct. Proj onto Poly_{<m}|_A commutes with translation, so residuals R_i orthogonal:
  Psi_A(u_theta) = 4^{-m} sum_i binom(m,i)^2 ||R_i||^2   (mu_A-norm), CONSTANT in theta (unnormalised u_theta).
  ||R_{m-i}|| = ||R_i|| (conjugation). ||R_i||^2 rational (Galois sigma_s = multiplication automorphism of A).
Normalised hat u = u/|u|_inf, |u_theta|_inf in [cos(pi/q),1] => Psi(hat u_theta) = Psi/|u_theta|_inf^{2m}.
Cheb bracket with c = residual r (in V): ||r||^2/||r||_1 <= Cheb_A(u) <= max_A |r|  (counting norms).
Generators of Poly_{<=e} component c: Poly_{<=e-1,c} + sum_b khat_b * top_{e-1, c-b}.

## N2. Exact runs (exact_iso.py, multimodular, cross-checked by full_check.py on sdc7 L=1 m=2,3: agree to 12 digits)
L=1: A_1 = {0, +-e_x, e_x - e_y}, |A| = 1+q(q+1), rho(A_1) = 4 (dim V(A,4) = 0) for all q in 7..103.
 m=1: Psi = (U-1)/U; m=2: Psi ~ 1.06 (q=7) .. 1.2356 (q=103), rel 0.52..0.553;
 m=3: Psi ~ 1.209 (q=7) .. 1.1355 (q=103), rel 0.222 .. 0.1833. dim V(A,3) = U - binom-ish (q=43: 903).
L=2 q=7: rel = 1, .677, .373, .155, .0401, .00578 for m=1..6; Psi(C_1) = 2.19, 10.07, 60.4, 348, 1427, 3468.
 dim V = 756,749,721,644,483,232 (m=1..6); h_A = 1,7,28,77,161,251.
full_check sdc7 L=1: Cheb(LP) m=2: 1.807 (theta=0), 1.901 (theta=pi/7); m=3: 1.713, 1.853. random u Psi 0.3-0.8 (|u|inf=1).
## N3. sdc11 L=2 exact: rel(C_1) m=1..8: 1, .6756, .3954, .2182, .1072, .03658, .006483, 6.0e-6; dim V = ..., 2162 (m=7), 5 (m=8, all trivial char) => rho=9.
 Psi(C_1): 2.77, 15.96, 131.8, 1416, 1.70e4, 1.67e5, 9.50e5, 3.01e4.
 cube rel (q=11, C_1): 1, .642, .356, .178, .0822, .0349, .0137, .0049, .0016, .00047
sdc7 L=2 m=8 (m>q): Psi const 23.127 on circle, Cheb(LP)=7.65625 const, LLL c height 42, |c|_1=1008, Phi_c = 7717.5 = barrier value. OK.
sdc7 L=2 m=4 LP: kappa_V(hat u) in [2.0910e-2, 2.0956e-2] over theta: nearly constant (ratio .998); Cheb/sqrt(Psi) ~ 2.4-2.7
## N4. LP kappa_V(hat u) on circle nearly constant: inf/sup >= 0.994 (sdc7 L2 m=2..6), >= 0.856 (L=1). Cheb/sqrtPsi ~ 2-3.
kappa_H: H_red = 1 for all L=1 (q<=23,m=2,3); sdc7 L2 m=5: H_red=3 (profile 261/472/483), m=6: H_red 4; m=8: 42.
 kappa_H(hat u, H=1) >= 6.3e-4 (m=5) vs sqrtPsi/d^m 2.25e-3 -> H=1 already within factor 4.
 L=1 limits q->inf: rel(m=2) -> .553, rel(m=3) -> .183. Psi(C_1) ~1.1-1.24.
 regime_need: required decay -ln rel/m' >= 466 (q=3251), 613, 841, 2806 (q=1e5), 9321 (1e6); cube ~0.69.
 tau* (toy) <= 2.2 everywhere vs tau_D 164-221.
## N5. q=11 L=2 float == exact (10 digits). Circle rel > max random rel at EVERY m (q=7, q=11): circle is the least tangent among samples.
 A rel(C_1) >= cube rel(C_1) for m <= rho-2 (q=7: m<=4; q=11: m<=6).
 sdc19 L=2 c=0 dims Poly_<m,0 (m=1..7) = 1,2,12,81,447,2024,6858 / 13911 (2024 matches record 13911-11887).
## N6 (resume run, after container restart)
 q19 L2 m<=6 complete (float_iso_sym, out/float_q19_L2.txt); deficiencies match record exactly (13464; 11887/11886).
 q19 L2 rel/rel_cube on C_1: m=2..6: 1.01, 1.02, 1.04, 1.03, 0.88 -> at fixed L, growing q, toy approaches cube.
 relaunched: float_iso_sym 19 2 7 6 -> out/float_q19_L2_m7run.txt ; float_iso_sym 23 2 5 6 -> out/float_q23_L2.txt (timeout 3600 each).
 cheb q11 L2 m5: only theta=0 finished: kappa_V(hat u_0)=3.516e-3, Cheb/sqrtPsi=4.34 (not rerun: >40 min/theta).
 exact_iso_big.py (p<2^24) q19 L2 m<=4 EXACT = float to 10 digits. exact_iso_big23 (p<2^23) q23 m<=3 exact.
 q23 L2 float m<=5: rel .658,.381,.2075,.1035 ; cube .655,.380,.205,.106 -> ratios 1.00,1.00,1.01,0.98.
 q23 m=5 circle rel inside random range [.0882,.106] (below max) -> §5.2 nuance.
 -ln rel/m at m=5 L=2: q=11 .447, q=19 .453, q=23 .454 (converging in q).
 running: float q19 m7 (m7run), float q23 m6 (m6run), exact q19 m5 (m5run), exact q23 m4.
