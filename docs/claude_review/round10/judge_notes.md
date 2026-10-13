# Judge notes (round 10), working file

## Context recap
- Best bound M <= (2/3+eps) log R/loglog R. Growth = W/(M log M) -> infinity.
- Round 9: route 3 = L4 (inverse) + Thm A/B (Linnik/Heath-Brown) + L7 (Poisson tails); L5 twin channel; L1/L2.
- P3-all* open; bracket 3 M log M - O(M) <= W_min^all <= (66A+o(1)) M log^2 M/loglog M.
- R1: residue caps only at fixed A; regime (1.4) p_j <= M^{h_0}, A chosen after c_0.

## L4K
- Thm 1 (proved, refereed): norm<=K budgets force collisions only for n_q < |S_h| ~ (2r)^h/h!, h=floor(K/2);
  free above A_K. Sharpens Prop H M^{K} -> M^{floor(K/2)}.
- Thm 2 (proved, refereed w/ minor repairs): L4(K) as posed in round9 §6 is false for every K>=2.
  Not fakes: violate norm-(2h+2) budgets at every level on positive proportion of free moduli.
- Lemma P (level 1): equal collision lattices => c' = u c.
- M1: threshold 2L(h0+h) wrong; via Thm B (Hamming) K ~ 2(h0+h) suffices; rigid-core form needs neither A nor B.
- M2: displayed L4^forced uses F_K={n_q<=A_K} mostly free -> fails first test (coherent share 0.13-0.25).
  Norm-5 budgets shrink K=4 flexible top: 107,113,127,151 rigid one-modulus; 167,199 keep 96, 40 non-genuine;
  no joint non-genuine top found. 211 (first norm-5-forced above norm-4 range) flexible (832).
- Referee replacement: L4*(K) on exact forced set, budgets all norms K' > K.
- FAKE TEST (referee §7, proved): every fake with genuine labels of its own primes satisfies hypotheses AND conclusions
  of L4, L4^forced, L4*. So the surviving lemma uses no fake-violated input. Separation must come from label-angle step (L7/slabs).
- Score thoughts: correctness high (refutation proved, Thm 1 proved); fake-violation: none (0); growth: none,
  but kills route-3 as formulated & relocates. ~4/10?


## analytic
- Prop A (fibre real-analytic statements = functions of angle vector -> capped). Referee major: arithmetic fibre inputs
  (N^{2k}S_k(N) in Z, residues) separate (dist to Z 0.095-0.417 on 4 explicit fakes); per-prime form p^{2k}e^{4ik phi} in Z[i]
  pins phi to finite set -> original problem; product/residue form = P3-all (open).
- Thm B graft barrier (proved): robust prime-supported ensemble inputs (error >= c X^gamma; LS const >= c(X+K X^{1/2}))
  satisfied by grafted fakes; A>4 needed for lw fakes (2<A<=4 only sharp.md fakes, const bA). gamma, A fixed before c0.
  All-ideal large sieve / Type I-II not covered.
- Prop C: absolute-value Fourier counting rotation-invariant -> never certifies emptiness.
- Thm D (sharpened): LCG_n(phi;n) iff phi = arg beta, beta primitive, Norm beta <= n. Verified fakes violate at 142/142 primes.
  = original problem.
- D.1/D.2 conditional on GPA (strip, not known to follow from GRH).
- Sharp (X+K) large sieve = WCG_c (max-scale cross-gap) up to constants; does not collapse; status open (between).
- Open: P_+(N) < M^{2-delta} regime (no fakes known; Thm B doesn't apply).
- Score thoughts: correct barrier (survives), fake-violated input: LCG yes, but = original problem; no growth. ~3-4.

## additive
- Collapse thms 2.1-2.3, 2.6, Prop 2.7 correct: lattice-ness in full = problem (norm-form for arbitrary imag quad order).
- AG(H) gap. Referee: (i) only a.e.-fake-violated (generic); universal = Lemma (AG), strictly stronger than growth
  (restriction to actual clusters IS growth), contains circle problem for all orders. (ii) heights wrong: l1 height
  M^{1+c0/2} (box M^{c0/2}); curvature rectangle -> Thm R2: relation-free violates AG at M^{1+c0/4+o(1)};
  rho <= 2 + c0 M/(4(A-1)); lattice coset at e^{W/4+O(M)}.
- Referee Prop R1: monomial gaps already give AG on all c with w(v(c)) <= W + log2; exact relations have w >= 3W/2 - O(log).
- Referee Prop R3: admissible positions (rank-M profile) empty interior, positive measure (fat Cantor).
- Prop 3.1 (capped, A_0) answers round9 L1L2 open pt 1 with delta=1/4. Thm 3.2 low-position (R1-arch) only, actual only.
- result_survives: FALSE (headline). Status partial.
- Score thoughts: correct parts but headline fails; generic-only; Lemma(AG) stronger than growth. ~2-3.

## determinant
- Thm 2.1 (proved): slab-box LP ratio >= 0.12 sqrt((n-1)/ln 4r); Hadamard <= sqrt n (Thm 2.3); LP-critical needs r >= e^{0.0144(n-1)}/4.
  => corrects round8 §2.4 / threshold "forced divisibility exactly critical": only for form fixed before cluster.
- Thm 4.1 rigidity (proved): F_c(z)=0 for all c in V_H; needs n >= 2137 (exact), maybe vacuous.
- Cor 4.3: a.e. fake violates; 'every fake' is special case of ZE. Referee: generic like "points in Z[i]";
  separating input = integrality + forced divisibility of multi-term forms (outside.md §6 item 1, outside barrier);
  new content = LP supercritical.
- Lemma 3.1 squarefree rank >= n/4-2; Lemma 3.2 nondegenerate full-rank subcluster; Thm 5.7: ZE => n = O(loglog R) squarefree.
- Referee Prop R1: ZE holds when slabs don't cut (Z_P = sparse strata); all difficulty in cuts; ZE must exclude translates gH.
- Score: 5 overall (correct 8, fake 4, growth 3).

## Judge's own analysis for next step (README remark, unrefereed)
- Z_P ⊇ Bs|A - m e| (real torus near e) since V_H ⊂ V(A,m). Toric components: 1-param subgroups u rational with width_A(u) < m
  (e.g. u = beta_p: width <= 2w < m). These are degenerate. Prop 5.2 low rank = subtori of dim <~ sqrt(n/log r).
- Non-toric components: analogue of negative curves on blowups of toric surfaces at identity (Gonzalez Anaya-Gonzalez-Karu;
  Castravet-Tevelev/ Goto-Nishida-Watanabe) -- quoted from memory; known source of non-toric base loci at a fat identity point.
  Risk for real ZE; then need algebraic-point ZE or intersection over all admissible supports (open pt 2).
- First test: n-1 = 2..5, random ±1 columns (squarefree), r ~ n..2n, choose T,w,m with m > 2w (slabs cut) and V(A,m) != 0;
  compute Bs near e (exact linear algebra for V(A,m); Groebner / numerical AG for components); test tangent cone Phi_c;
  test intersection over several (T,w); test gap-admissible fake positions on Z_P (ZE => none).
