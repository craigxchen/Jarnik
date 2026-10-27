# pass2 notes
- prior referee report saved as round17/referee_hlight_pass1.md (not trusted; independent re-check)
- Hand re-derivation OK: L1.1, L1.2(1-4), L2.1, L2.2, L2.3 (Euler+Vandermonde; binom(m,j)!=0 since m<q), P2.4 (phi concave,
  n=4 nu>=q-3), P2.5 (M bound, g(n) algebra, (n+2)(g-(q-3)/2) identity), P2.6 (c = (q/2)(d1-d0) - f + e), L2.7 Walsh
  Phi = 3/2 + e2(x)/2 - x1x2x3x4/2, constants -6-9+4 = -11.
- L_T for sdc_q: lambda in L_T => h(a) = lambda(P_a)+lambda(P'_a) = lambda_a - lambda_inf in Z => lambda = mu + t*1, mu integral,
  t in (2/(q+1))Z (from lambda(P_a) in Z, |P_a| = (q+1)/2). relevant (k = e_1) iff t in Z iff integral. => L_T/Z^n_0 = Z/((q+1)/2). OK.
- q=7: LB(7)=2 < 7/3 (n>=6 F=0 bound (q-3)/2=2); "c0=1/6 every q" needs n-1>=5 at q=7.
- Table proved LB/r uses ceil(sqrt q): conservative (q=199: true LB/r=.1692, table .1667). harmless.
- MILP "exact global minima" q<=43: HiGHS is floating point; but q>=19 min follows from exact 4-term C enumeration + Thm;
  q=7,11 need exact exhaustive (do it).
- c1: plib built from exact.md Thm4.1 (Paley-I, DKD normalisation, base inf): |S| = (q+1)/2, 2q distinct, formula == definition
  (300 random lambda each q in 7..59), sigma_xy = 1 all pairs, L_T/Z^n_0 = Z/((q+1)/2) by SNF (q=7,11,19,23).
- c2: exact global minima q=7 (|lam|<=2 all 8 coords) = 3, q=11 (ternary) = 5; with Lemma 2.1 => exact global min. q=7 F=0 n>=6 min 5.
  q=11 class-by-class minima all >= proved bounds; F=0 n=10 min 9 = n-1 (Prop 2.5(1) tight there).
- r4b.c (own C, normal form d0+d1-dc-de): 4-term minima identical to author's at all checked q incl 3251:1571, 3259:1575,
  3271:1581, 5479:2667. 0 violations of L2.7 bound; max|S_T| <= 2sqrt q + 1 (Hasse) everywhere (5479: 149 vs 149.04).
- c3: normal form == brute force (q=19: 7, q=23: 9); Lemma 2.3 0 violations (6884 random incl multiples of q; exhaustive ternary q=11);
  dual pairs sigma=(q-1)/2 by definition q=7..83, formula 3251,3259,5479; base 4-terms take {(q-1)/2,(q+3)/2} (q=7..103).
- c4: first regime prime 3251 (3203 fails by 4.8), m'=2057<q, 2m'<q false; true Weil LB at 3251 = 1528.97 = 0.23515 r; o>=3 / o>=4.
- [L] LN Thm 5.41 located via citing papers arXiv 1902.04736, 2012.01483 (d = distinct roots, f monic not m-th power).
