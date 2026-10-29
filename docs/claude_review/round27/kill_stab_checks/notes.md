# kill_stab notes (scratch)

## Reading pass 1 (stab.md)
- Thm 3.3 proof steps 1-3 re-derived by hand: Step1 peeling for dims OK (kernel = A_t(slice), image in A_{t-1}(rest)).
  Step2 OK incl. t-i+1<=0 edge case (forces singleton slice, MW=0). Step3 (*) OK: P'=empty term weight 2^{1-n},
  exact value c_max-c_min vs bound c'-c. Two-critical contradiction j-1<j-1 OK; j=i+1 case contradiction is immediate.
- Cor 3.4: 8 MW_4 = sum over 15 cuts of width (Cor B1 with M=5: 2^{-4} * 2 * sum_cuts). width_{P_D}(1_A) <= d_A + d_{A^c}. OK.
  NOTE: for D=-K, D.beta=30, 8t>=30 => t>=15/4 > 10/3=tau: Cor 3.4 is vacuous for multiples of -K; content only on face.
- Thm 3.2 = h(phi_L(y)) <= h_L(y)+O(1) with h_L(y) ~ (L.gamma_hat /2) w at balanced near-tight points; i.e. Northcott
  for the Iitaka map of L = pi^*D - tE. s=0 case needs kappa(L)>=1 on gamma_hat^perp face.
- Lemma 3.1 re-derived: D_S uses Newton-polygon min_S<m,a(p)>, stronger than toric-divisor psi_D (additive on chains).
  chain deficit = sum_p [psi_{D_S}(a(p)) - width_S(a(p))] log p = extra p-adic proximity to Bs of the section.
- Root binomials x_a/x_b - 1: width additive along chains -> chain deficit identically 0.
- Prop 1.3 trivial reformulation, correct.
- Sign slip §1.5(3)/Prop 3.8: R|iota+i => x = -it = iota*t mod R, not -iota*t. Harmless (criterion symmetric in sign).

## Record cross-check (novelty)
- round25/intra.md:523-524 open item 3: "Non-simple balanced profiles: for chains sum_t 1_{A_t} the width is not additive
  and Theorem 5.1 is not proved". So the chain non-additivity is ALREADY in the record; Prop 3.5 makes it explicit (gain).
- threshold/construct.md §3 item 3: "Constant-level effects in the equality case of Theorem B. These are closed under
  Conjecture E or C_bal." -> the pencil principle (two tight sections, bounded c's) is a constant-level equality-case
  effect = a mapped loophole. Thm 3.3 closes its pencil part unconditionally (not the single-tight-non-root-form part).
- threshold/referee_construct.md §2: Theorem B == Cayley transform (w=(i-t)/(i+t)) of multiplicity Schwartz-Zippel on
  {i,-i}^n (DKSS 2009), recorded as hermitian_polynomial_barrier.md ineq (1). So Thm 3.3 = equality-case uniqueness
  for multiplicity SZ on the 2-point grid (Cayley picture). stab.md does not mention this.
- referee_construct.md §4: finite reduction + exhaustive: 2097 tight kernels (cut closures in small boxes), ALL dim 1.
  So "dim A_MW = 1" was already verified as evidence in the record; Thm 3.3 is the general proof.

## Computations
- cox.py -> cox.out: 400 random genuine one-class profiles, pi=a+bi with b>a>0 (opposite of author's): param, primitivity,
  norms, conj-primitivity of Q_ab, zeta=G'Q, |z_a-z_b|^2 N(Q)=4t^2N', triangle Q_ab Q_bc = G^(b) Q_ac (all 60 ordered
  triples), Pi_a conj(Pi_b)=G_ab Q_ab, x_b t_a - x_a t_b = G t_ab, quadratic identity: 0 failures. 181 coincident pairs (t=0).
  (Ch) <=> primitive&canonical 3000/3000. Prop 1.1(iii) example reproduced (canonical blocks pibar^2 in {3}, pibar in {1,3,4}).
- thm33.py (running): exhaustive boxes (4,1),(1,1,1),(2,1,1),(5,1),(1,1,1,1) = 75003 sets, 386 tight, 0 fails of B / 3.3.
- thm33.out FINAL: 82703 sets (exhaustive (4,1),(1,1,1),(2,1,1),(5,1),(1,1,1,1); random n=3,5,6; zonotope unions n=3,4,5;
  cut closures n=3,4), 1638 tight, 0 failures of Thm B, 0 of Thm 3.3; mod-p certified (dim_p >= dim_Q).
- cor34.out: (A) 120 root-product classes D.beta=8t: h^0 = 1 exactly, 0 fails. (C) defect search 3000 supports: min s=4
  (triangle {0,e1-e2,e1-e3} type), found s in {4,6,7,8,...}; (D) 1736 classes D.beta-8t in {1,2,3}: no pencils.
- cor34b: sums of 2-3 root-product divisors: P_D = zonotope sum exactly (gen. permutohedra), h^0=1.
- cor34c: 8794 random (mostly non-nef) classes on the face D.beta=8t: h^0 = 0 always.
  => effective part of face gamma_hat^perp in samples = cone of the 10 diagonal divisors (Conj E structure half, open).
- split_check.out: all 4 Pell orbits of z^2-56x^2=65, 25 points x<=4e10. Author's orbit reproduced exactly (alpha 6.5143,
  4.6651, 4.3679, 4.1969, 4.0970, 4.0071, 3.9913; residues [8,5,9,5],[1,5,1,1],[8,1,1,1],[1,1,9,5]). Twin rep holds.
  Max chain deficit over ALL orbits 11.15 (x=98382962, w=9.0), not <=6.8 (author's orbit only). Never >= s w = 4w. Bounded.
  2-6 nonsimple (bounded) primes per point. Balance deviation ~0.25-0.48 at x~1e10 (slow).
- k4_check.out: Im P_x = -420,-180,-240 (sympy). 38 values of s (random up to 1e9): gcd norm <= 10368000; max|c_ab|^2 =
  230400 (author: "<=14400" was a sample max over s=10^k+7, all but one giving 162); true bound 4*420^2. balance -> 0.09;
  chain deficit max 7.2 vs s_S w ~ 80: Prop 3.5 k=4 analogue never fires. Bounded residues confirmed.
- lemma31.out: 300 non-simple cases: D_S|F, N(F)/N(D_S) = |f(y)|^2 prod p^width, subadditivity: 0 fails.
  Example S={0,e1-e2+e3-e4}: sum widths 12, s=4. By hand: +-2 cuts {13},{24} never cancel (subsets of {1,3} have <m,.> >=0,
  supersets <=... ), +-1 cuts pair into 4 nested opposite-sign chains {1}<{124},{3}<{234},{2}<{123},{4}<{134} -> max Delta ~ 8w > 4w.
  So Prop 3.5 hypothesis is satisfiable by balanced profiles (needs chain primes carrying >= (2+c/2)w log weight).
- record: round25/README.md:73 "Direction-sensitive bounds (route 5's leading form Phi_c(u)) are not covered" -> Prop 3.7
  placement right, citation "intra item 3" wrong (that item = non-simple profiles).
- Literature (WebSearch): no characterization of equality in multiplicity Schwartz-Zippel found; no mean-width regularity
  uniqueness found. Castravet-Tevelev (arXiv 1311.7673; CLTU arXiv 2009.14298): rational contractions
  Bl_e LM_{n+1} -> M_0n -> Bl_e LM_n; so Cor 3.4 is about a face of Eff of a birational model of M_{0,7} (unchecked pointer).
