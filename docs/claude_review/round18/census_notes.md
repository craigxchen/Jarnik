# Notes (parking derivations) for census.md

## N1 Formulation
- (L_8) ptolemy §7 / L8 §1: anchored P_0=(1,0), det(P_i,P_j)=t_ij prod_{T∋i,j} n_T, n_T two nested chains per prime,
  |log n_T - w| <= eta_0 w; H = sum log|Y_j| + sum log|t_ij| >= c w - O(1) off Z.
- Visible T: 2<=|T|<=m-1 (119 at k=8). Singletons/top invisible and free (L8 Thm B proof: top by diag(n,1), same point of M_0k).
- Lemma 5.2 anchoring: every Phi has exact form with Y=±1, t=±1, n_T = d_T (contact integers), top = delta.
  (round4 L7 Thm 1 = k=7 case.)  So in the integral chart (P_i=(X_i,1)) H=0 always; census = contact profile.
- Lipschitz (Thm B(iv)(2)): |l_T(p) - v_p(n_T)| <= 2 max_ij v_p(t_ij) => |log d_T - log n_T| <= 2H.
- Sandwich: (L_8') "off Z, all |log d_T - w|<=eta w and w large impossible" implies (L_8) (eta'=eta/2, c=eta/4);
  (L_8) implies (L_8') restricted to admissible-tree Phi (two chains at every prime), which includes all
  configurations with pairwise coprime d_T.

## N2 Heuristic (lattice) counts, units of w
k: visible blocks V, index I = sum(|T|-1), pair-blocks b = 2^{m-2}-1, dim m-1.
 count exponent  V + (m-1) b - I:  k5: 10+9-14=5; k6: 25+28-45=8; k7: 56+75-124=7; k8: 119+186-315=-10.
 (ptolemy table 1-K.beta: 6,9,8,-9: differ by 1 = shift class / top block.)
 per random assignment success prob exp((m-1)b - I) : k5 -5, k6 -17, k7 -49, k8 -129.
 => random-assignment CRT search cannot find bounded-H points even at k=5,6 (calibration failure).

## N3 Families threshold
7K = -2 D_2 + D_3 + 2 D_4 (28,56,35). If gamma.D_S in [(1-eta)a,(1+eta)a]: 7K.gamma >= a(70-182 eta).
 K.gamma <= 0 => eta >= 5/13. Covering genus<=1 families have K.gamma<=0 (rational: <= -2).
 Points version via Lemma 6.2/Step 3: h_K >= 10w - 282 eta w - 60 H - c; on fixed curve with K.C<=0: h_K <= o(h).
 => eta >= 5/141 - o(1) if H = o(w).
 Product charts (round4 L7 Thm 2(b)): beta in cone{gamma_k}, K.gamma_k<=-2, K.beta=10: impossible.

## N4 Exclusion: Y cap D_T = empty => d_T bounded on Y(Q) (Z-model, finitely many primes, bounded mult.)
 - symmetric loci Y_sigma (sigma != id): exists T with T, sigma T incompatible -> D_T cap Y = empty. (check all S_8)
 - forgetful fibres pi_Q^{-1}(c), c in M_{0,Q}(Q): D_T with T|Q nontrivial split.
 - Y_sigma example (12)(34)(56)(78): orbits of splits 67; T={1,2,3} type (pair+one) 24 splits & {1,2,3,5}-type 24 splits untouched.

## N5 pull-back of E_t under forgetting 2 points: each 6-contact must split 4 ways balanced; 19 new splits.
 heuristic: x must satisfy congruences mod ~e^{50w}, |x|~e^{31w}: count ~ e^{-19w} per new point.

## Progress log
- census.md §1 (formulation, Lemmas 1.1-1.2, Prop 1.3) and §2 (Theorem E, emptiness table) written.
- scripts done: census_lib.py, test_lib.py (284752 cmp ok), check_symmetric_loci.py, emptiness.py(.out)
- NEXT: k5 exhaustive calibration; k6 E_t calibration (repo scripts threshold/L8/k6fam); lattice census;
  families Prop F + exclusion props S, Phi; coverage census small heights (C?).
- k6 E_t census (k6_family_census.py): N=2..14: w* 5.6,15,32.4,53.8,82.2,112.8,152.4,195.2,243.7,298.5,358.1,423.1,493.9;
  eps 1,.893,.388,.353,.211,.236,.149,.152,.119,.118,.095,.103,.078 ; Sigma-free eps .460(N4) ... .080(N14);
  Sigma-part log 46-79; maxlog|X| ~ 7 w* (3431 at N=14). H=0 representation (anchored at x1).
- census_small k=5 B=100: 3,646,069 configs, full coverage 775,434 (21%), best eps .1333 at w*=2.77 ([0,22,35,46,85]); frontier
  eps<=.1345 for w*>=3, .3235 for w*>=4.
- lattice k=5 A=50 (w=3.91) 20 samples: minH/w 7.34, mean 8.81 (H=0 balanced points exist at this scale: calibration failure).
- census_small done: k6 B60 max cov 18/25; k7 B40 max 23/56; k8 B34 max 29/119 (5,359,838 configs). Tables filled.
- lattice census done: H/w ~ 8, 42, 150, 464 (k=5..8). verify_lattice_example ok.
- census.md complete (sections 0-9).
