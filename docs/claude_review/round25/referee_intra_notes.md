# referee notes (intra)
## pass 1 (reading intra.md fully)
- L2.1 (1)-(4): re-derived, OK. h_B >= W - 2log L ok; W/2 needs L<=sqrt R (C<=1); for C>1, kappa<=2+O(logC/W).
- L2.2 ok (toric bound for F(zeta) == toric bound for G(Z) after common factor).
- L3.2, L3.3, Thm 3.1 re-derived OK. Reduction near-tautological: same algebraic number, same valuation data.
- Thm 4.1: Phi^* basis map ok; s_j in 4T R_5 ok (|s|_1<=8T). (4): padding => Phi^* surjective, Phi closed immersion (writeup says "surjection onto subtorus": wording).
- Thm 5.1 algebra ok; Cor 5.2 arithmetic ok: 2^{M-3}(1-eta)-(2^{M-1}-1)(1+eta)/4=(1-(2^M-1)eta)/4.
- Prop 6.1 algebra ok (E=(8/3)sum, 40 = (8/3)*15).
CONCERNS
 C1 Cor 5.2(2)/5.3: w threshold depends on |C(F)| (template dependent); Cor 5.3 "every inference consistent with ONE datum" has quantifier order for-all-F-exists-w. Fix: only negative part of C matters; Taylor bound C(G) >= -c*ord (ord<=reg<2 K1 by Thm B) -> uniform.
 C2 "settles open point 1 for intra-cluster configs": only for Z cube-balanced on ALL 2^{M-1}-1 cuts (r>=2^{M-1}-1 primes). Subtuple-balanced sparse designs (random cuts, r=O(log M/eta^2)) have all y_S balanced at mu~4 but Z not cube; Thm 5.1 silent; route 5 supercritical there (Cor 7.1). So open pt 1 settled only for cube class.
 C3 Def 3.0 allows "any U"; Thm 5.1/Cor 5.3 only for Dyson/Taylor-type U. Cor 4.2 "= route 5's framework": route 5 bound isotropic (L4.0) while L_intra anisotropic. Precision.
 C4 Def 1.5 "beats cap = contradiction or exact relation": monomial relations (sparse profiles) trivially count.
## pass 2 (context)
- round24 open pt 1: "cap for (*) on configurations with shared primes"; assessment there "suggests same cap".
- round10 Thm 2.1: slab box at T=n, ARBITRARY profile beta_p, omega_p; supercritical when n-1>=300ln4r. Thm 4.1 rigidity for any cluster with (4.1).
- => C2 sharpened: for 5-wise-uniform designs (linear cuts u.c_x, any 5 c_x independent) every subtuple EXACTLY balanced,
  r=2^m-1 << 2^{M-1}-1; in route-5 regime the intra method (star blocks) is supercritical -> exact relations
  at data where all y_S are balanced with mu<4(1-eta). So NO data-level cap for shared intra-cluster configs in general;
  cap only for cube Z. Headline "settles open point 1 ... cap 15/4 even with maximal sharing" overclaims; "maximal sharing" undefined.
- C(G) fix: width_S(1_A)<=K1:=max|k|_1 (|k(A)|<=|k|_1/2) => reg<2K1 ; Taylor bound bracket >= m(1-log4)-O(log m) using |<k,t>|<=(|k|_1/2)diam.
## pass 3 (computations)
- r1: L2.1(1)-(4), Thm3.1(1), L3.3 exact on record 10-pt (C=8.2386, lamZ=3.1127 matches), 9-pt, 8-pt N=5*13*17*29*37*41*53. all True.
  sharing gain is F-distribution dependent (mine 0..134); author's 3.5-64 is evidence only.
- r2/r2b: (5.1) 0/95 random + structured (R_5:3, 2R_5:6, cube slices) ok; reg<2K1 ok (max reg/K1 = 1 on cube6 slice? printed 0.5 for random; cube6: reg 6 = K1).
- r3: 5-wise linear designs: all subtuples EXACTLY balanced, cut count 32/63, 64/255, 64/127, 256/32767. ext-BCH t>=13 (M=8192, r=2^27-1) meets n-1>=300ln4r -> route-5 LP supercritical (margin m/2-w>=m/4 in Norm units) => intra (*) violated at data level with all subtuples balanced, mu<4. => NO data-level cap for shared intra-cluster configs in general.
- C3: Cor 5.3 "every inference of L_intra" but Def 3.0 allows any U; Thm 5.1 only for Dyson/Taylor U depending on (idx/ord, lambdas). Direction-sensitive bounds (route 5's Phi_c(u)) not covered.
- check3b: reg mod-p is an upper bound -> E upper bound -> "max E<0" conservative; but family restricted (slab boxes, k' in {-1,0,1}^9) -> "nothing supercritical at M=10" only within scanned family.
