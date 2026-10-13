# Referee notes (compute.md)

## Reading pass 1 (lines 1-200)
- Prop 2.1: zonotope certificate. Check: c in V(A,m) via forward diff -- OK. |c|_1 <= 2^m OK. H>=2^m for T>=2 OK; T=1 "check directly" -- need to check H for T=1: H=(U*1^{m-1})^2 = U^2. 2^m <= U^2? For M8: m=8, 2^8=256, U = |A| ... check.
  BUT: Prop 2.1 conclusion "every point of Bs satisfies w^{v_i}=1 for some i" - only for i with alpha_i>0. Fine.
  Also note: degenerate in the sense of complex torus: w^{v}=1 means <v,u> in 2pi Z, not = 0. "near e" OK.
- Prop 2.2: uses Prop 5.3(iii) corrected constant. |v_i|_1: v_i in Z^d has |v_i|_inf <= 2T/L_i actually; v as element of Z^n_0: last coord -sum v, |.|_1 <= 2*sum|v_x| <= 2*(n-1)*2T = 4n(n-1) = 2D. OK. Check final inequality W <= 8 log(C D/sqrt2) arithmetic.
  Character inequality from Prop 5.3(iii): X_v >= W/2 - 2 log(C|v|_1/(2 sqrt2)). With |v|_1<=2D: X >= W/2 - 2 log(C D/sqrt2). Combined with X <= W/4: W/4 <= 2 log(CD/sqrt2), W <= 8 log(CD/sqrt 2) = 8 log C + 8 log(2n(n-1)) - 4 log 2 < 8logC + 16 ln n + 8 ln2 - 4 ln 2 = 8 log C + 16 ln n + 2.77. OK (<3). But if C<1 need log^+... fine.
  Issue: "some X_{v_i} <= W/4" -- from sum L_i X_{v_i} <= (sum L_i) W/4, yes a weighted average. But the inequality chain: sum_p omega_p h_Z(beta_p) <= sum omega_p 2w <= 2w <= m/2. And omega_p |<v,beta_p>| = (e_p log p/W) |<v,beta_p>| = (e_p log p / W)(2/e_p)|<v,a_p>| = (2 log p/W)|<v,a_p>|. So sum_i L_i (2/W) X_{v_i} <= m/2 -> sum L_i X <= mW/4. OK.
  Subtle: Prop 5.3(iii) needs <v,delta> != 0 mod 2pi; and also uses the realisation (phi, kappa) of actual cluster. Fine.
  Also subtle: the cluster positions in the regime are nondegenerate -- Prop 2.2 assumes nondegenerate positions (for the full cluster). OK as hypothesis.
- Prop 2.3 section principle: fine.  Claim "rho(A') <= 2Tj" in regime - section in box, coordinate functionals... t coordinates? A' in t-space; the functionals k_x restricted give linear forms in t; need j independent ones each with <=2T+1 values -> complete intersection, rho <= 2Tj. OK-ish.
- Prop 2.4: n-1=2. Need to check carefully. "(ii) h_A(beta_p)<m for every column" -- uses supercritical in every column. Two cutting patterns case: rho(A) <= 4L. Claim: l_i^{2L+1} in J. Then (l1^{2L+1}, l2^{2L+1}) contains all forms of degree >= 4L+1. Yes. And m > h_A(beta_p) = 2 width(l_i) = 2*2L = 4L. So m >= 4L+1 > rho. Need: width of l_i on A = 2L when slab cuts. "box symmetric and contains k with l_i(k)=L" -- if slab cuts, there's a box point with |l_i| > L, then ... does A contain a point with l_i = L exactly? Box is lattice points of [-T,T]^2; l_i = k_1 or k_2 or k_1+k_2 (up to sign). Values of l_i on box are contiguous integers; the points of the box with l_i(k)=L: exist. But do they satisfy the other slab |l_j| <= L? Hmm! With two cutting patterns, the point with l_1 = L must also satisfy |l_2| <= L and |l_3|<=L if cutting. E.g. l1=k1, l2=k2: point (L,0) fine. l1 = k1, l3 = k1+k2: (L, 0): l3 = L fine. OK so with only coordinates and sum, (L,0) or (0,L) or (L,-? )... l1=k1,l2=k2,l3=k1+k2 all three: point (L,0): k2=0, sum=L fine. Point with k1+k2=L: (L,0) fine. So width = 2L for all cutting ones. But what if width of some cutting pattern is less than... fine. But what if only some patterns cut and "every column" includes non-cutting ones with h_A < m: non-cutting means |l|<=L automatically on box -> width <= 2L anyway.
  Hmm, but wait: profile columns at n=3: beta in {+-1}^3 non-constant -> patterns {1},{2},{3} up to complement. Constant columns? beta=all +1 gives <k,beta>=0. Fine.
  Also "one L for all columns" -- restriction from author's model class. Is that faithful to Theorem 2.1? Thm 2.1 uses one w for all. OK.
  But the README asks also for e_p=... squarefree only. OK.
  One cutting pattern {3}: hexagon. Claim zeta >= 2T+L and rho <= 2T+L. Check zonotope: (-T,T-L) + t1(1,-1) + t2(1,0) + t3(0,1), t1 in [0,2T-L], t2,t3 in [0,L]. k1 = -T + t1 + t2 in [-T, -T+2T-L+L] = [-T,T]. k2 = T-L - t1 + t3 in [T-L-(2T-L), T-L+L] = [-T, T]. k1+k2 = -L + t2 + t3 in [-L, L]. OK. Sum = 2T - L + 2L = 2T+L. Good.
  rho <= 2T+L: line-by-line Newton interpolation. Term j: prod_{i<j}(k1+k2-s_i) * q_j, q_j has degree <= (#points on line s_j) - 1 = 2T - |s_j|. j-th term degree j + 2T - |s_j|, |s_j| = ceil(j/2). Max over j <= 2L: j - ceil(j/2) = floor(j/2) <= L. So degree <= 2T + L. OK. But wait, does the Newton form interpolation work? We need: on line s_j, previous terms are given and we need q_j to fit the residual divided by prod_{i<j}(s_j - s_i) which is nonzero. Line s_j: points with k1+k2=s_j, |k1|,|k2|<=T: number of points 2T+1-|s_j|; parametrize by k1; univariate interpolation of degree 2T-|s_j|. OK.
  Missing case: zero cutting patterns -> slabs don't cut, excluded by (i). Also "two or more cutting" includes {1},{3} etc. OK.
  Also: what about L >= 2T for pattern {3}? then doesn't cut. L<2T fine.
  One cutting {1}: L < T; A = [-L,L]x[-T,T]; zeta >= 2L+2T, rho <= 2L+2T. (ii) needs m > h_A(beta_{1}) = 2*2L... but also non-cutting columns {2},{3} need h_A < m: h_A for {2} = 2*2T=4T, for {3}: width of k1+k2 = 2(L+T), h = 4(L+T)? Hmm wait: h_A(beta) = 2 width(k(S)). So (ii) requires m > 4(L+T) if {3} in profile. Then m > rho. Anyway conclusion fine if profile includes them; if the profile is just {{1}} then m> 4L and rho = 2L+2T... m <= rho needed, so zeta=rho>=m. Fine, conclusion "m <= zeta = rho" holds.

## Reading pass 2 (lines 200-613): candidate issues
1. M8, M8' satisfy only (ii') (range < m), not README (ii) (h_A < m), and far from the regime LP condition
   sum omega h_A <= m/2 (Thm 2.1(d)). M8: h=8=m; M8': h=12>11. "faithful" = only necessary condition zeta<m.
   -> major interpretive: "first faithful models" overstated.
2. Prop 2.4 (n-1=2) proved under (ii) h_A<m, but n-1=3,4 scans and n=8 models use (ii'). Apples/oranges.
   Need: does n-1=2 have zeta<rho under (ii')? (2D centrally symmetric hexagons are zonotopes -> probably zeta=rho anyway; check.)
3. ZE-braid: likely FALSE in regime. Prop 5.2 strata delta = sum_{p in S} t_p a_p with sum_S h_A(a_p) < m lie in Bs;
   if rows restricted to S are distinct (|S| >= log2 n, e.g. Sylvester basis columns) generic points have no collisions.
   In regime m/w ~ 0.24 sqrt(n/ln4r) >> log2 n. => ZE-braid false for large n. Check at model size: H8, T=1, L=1:
   range 2, 3 columns -> 6 < m <= rho(A_1)? compute rho(A_1).
4. Thm 5.2 indexing: "x,y in {1..8}" vs 0-based 0..7 with base 7. minor.
5. Lemma 5.1 requires all two-term factors of form w^a - w^{a'} (coef -1); check w^a + w^{a'} absent. 
6. Thm 5.2(b) uses 24 of 42 polys: fine logically (superset).
7. "no gap-admissible fake position lies on Z_P" OK given nondegenerate.
8. Prop 6.2 proof checked logically OK; verify exact inputs.
9. "n-1=3: all 355" restricted T<=4, L<=3; "n-1=4": 12 random profiles -> overgeneralised "model sizes n-1=2,3,4 are all unfaithful".

## Computations done (referee's own code)
- r_had8.py: H8 T=1: L=1 |A|=57 rho=3; L=2 |A|=1037 rho=8 dimV(8)=42; L=3 |A|=1877 rho=11 dimV(11)=7. (mod p ranks; exact below)
- r_kernel.py: exact kernels, M8' dim 7 (sup<=4), M8 dim 42 (sup<=6). LLL. Verified M c = 0 exactly.
- r_factor.py: M8' each basis poly = monomial * linear-in-x binomials * ONE x-form (3 cubics, 4 quartics). Matches author.
  M8: 24 polys with only x-form non-binomial factors, 18 with non-form factor; 27 distinct two-term factors, all w^a - w^b;
  my basis has 7 non-linear two-term factors (w0w3-w4, w0w5-w2, w0w6-w1, w1w3-w5, w1w6-w0, w2w3-w6, w2w5-w0) vs author's 6. Basis-dependent; harmless.
- r_saturate.py (flint Groebner + Bayer, own code; G_j = product of factors): M8' unit ideal after x0..x6, x0-x1..x2-x3 (30-element GB of forms, same as author). M8: same list -> unit. THEOREM 5.2 CONFIRMED.
  (Author's cone_reduce.py also uses products G_j: checked.)
- r_sdc7.py / r_sdc7_c.py: sdc q=7 built from Paley; 14 cuts Plotkin tight; |A|=757 for every base point; exact dim V(A,8)=1, rho=8;
  c even sup 42 supp 687; Phi = 8! g, g 196 terms coeffs {-2,-1,1,3}, irreducible over Q; g<0 at all 65604 nonzero-valued pts of [-2,2]^7 (12520 zeros);
  positive witness (25,6,-12,1,-11,-14,-20) -> 76047300 (base=7, my coordinates). zeta<=7 rigorous via min_v sum|v(S_p)| = 8 and 14*4/8.
  PROP 6.2 INPUTS CONFIRMED.
- ZE-braid falsity: Sylvester n=2^k, k>=16: k*w < 0.7 n(n-1) <= m, so the k-column subtorus (Prop 5.2, kappa=0) lies in Bs; generic points collision-free.
- r_scan.py n=3 (T<=8,L<=7) under (ii'): 55 distinct admissible A, all zeta=rho (witness); none under (ii).
  2D general proof: P={|k1|<=a,|k2|<=b,|k1+k2|<=s}: rectangle/parallelogram (CI) or hexagon = zonotope with edge lengths
  s-a+b, a+b-s, s+a-b (sum a+b+s); line interpolation gives rho<=a+b+s. So zeta=rho for EVERY 2D slab-box, no (ii) needed.
- r_scan.py n=4 (T<=4,L<=3): 355 distinct admissible A (matches author), all zeta=rho with my own witnesses. (instances count 3664 incl. m; author 980 = different convention)
- r_tc_m8.py: TC claims of §4 confirmed (exp<=2, >=4 vars, 7+21+35 subspaces dim<=3, exactly 7 4-planes = cuts avoiding base, none dim>=5). max coeff 161280 in my basis.
- r_had8_T2: T=2,L=1 |A|=57 rho=3; T=2,L=2 |A|=2507 rho_p=8 (h=8): README (ii) infeasible at n=8 Hadamard T<=2.
- r_braid_threshold: Sylvester k=16 (n=65536): regime holds, k*2floor(w/2)=2.682e9 < m=3160214892. => ZE-braid false. k=12..15 not by this argument.
## Minor issues list
- Prop 6.1: A' arbitrary finite subset: |v_i|_1 <= 2D bound needs A' in box B_T; else log diam enters.
- "every model with m<=2T is unfaithful": coordinate segment of length 2T lies in A only if not cut (L>=T or x,base not separated).
- Thm 5.2(a) labels {1..8} vs 0..7.
- "six four-term relations": basis dependent (mine: 7).
- "Each single P_c has non-toric zeros through e with nondegenerate real points": not proved (P_c(e^{iu}) complex-valued: 2 real eqs).
- §6 "zonotopal whenever |J|<=5": overgeneralised from samples.
