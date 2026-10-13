# Referee working notes (L1L2)

## Read so far: sections 0-2.1
- Thm 1.1: algebra of four cases re-derived by hand: OK. Normalisation (G absorbs inert primes, (1+i)) OK.
- Cor 1.2: (1.1) OK by norm identity |z_x-z_y|^2 = Norm g * nu_s B^2.
- Lemma 1.3: re-derived. Argument: Re(z conj v) integer, arc contains bisector direction, so
  Re(z conj v) in [R|v|cos(L/R), R|v|], <= |v|L^2/(2R)+1 integers, <=2 pts each. OK.
  Need to check outside.md Prop C.2 statement, and arc-contains-bisector claim (short arc).
- Thm 2.1: essentially a one-line rewrite: L1 on near-extremal <=> max mu = O(log M) (since mubar=Theta(log M) by (F))
  <=> max |B|/t = M^O(1) by (1.1). Correct; (b) constant O_C(1/log M) actually depends on a,b also.
- Thm 3.2: correct (one line: product of c window primes <= |B| <= (C/sqrt2)M^b). It is the contrapositive of
  Thm C's proof applied to near-extremal clusters: ANY two hypotheses that jointly contradict near-extremality are
  incompatible on every existing near-extremal cluster. Not a barrier.
- HEADLINE ISSUE: "L1 and L2 are each other's negation" is FALSE: they are contraries, not contradictories.
  not-L1 does not imply L2; and if growth holds (no near-extremal clusters for large M) BOTH L1 and L2 hold
  vacuously (write-up itself says this in Cor 3.4.2). Same for R1-arch vs R1-arith (R1-arith is far stronger
  than not-R1-arch: forall b forall cluster vs forall b exists cluster).
- Thm 3.3: re-derived: W >= 2(M-1) mubar holds (column sum_{x<y}|a_x-a_y| <= e M^2/4). OK. With L2:
  leading coefficient 4(1-eta)(2-eps/4)=8(1-eta)-eps(1-eta) >= 8(1-eta)-eps. OK. It is p3lower B.1 + identity (1.1).
- Cor 3.4: "route 1 not weaker than growth": tautology for ANY pair of lemmas universally quantified over
  near-extremal ACTUAL clusters (vacuous truth). Not specific to route 1. Round-7 L1 was meant as "a valid input
  forcing..." i.e. derived in a larger class; the actual-cluster rewrite loses that.
- R1-arith alone does not give growth (needs uniformity). Write-up says "together". OK.
- Lemma 1.3 re-derived independently (worst-case arc position: Re(z conj v)/|v| in [R cos(L/R), R]). Correct.
  Analogous chord-direction argument gives only LOWER bounds on |beta| (consistent with write-up claim).
- Sum of three margins: mu_xy+mu_yw+mu_xw = W/2 - 2 log Norm g_xyw <= W/2 (trivial upper bound on triples).
- Prop 3.5(d) heuristic: bold conclusion "L2 cannot follow from genuineness + profile structure" contradicts its own
  caveat that near-Plotkin-tightness of adapted rows is open (profile structure is exactly that).
- Prop 3.5 status para: "B_xy of size M^Theta(log M/loglog M)" only a lower size for L2 pair; upper is C e^{W/4}.
- 5.1 item 4: "cannot be pursued as independent lemmas ... each is the other's negation" -- logical error.
  In a proof by contradiction both lemmas are proved for a hypothetical cluster; incompatibility is the point.
- Prop 5.1(b) derivation re-done: correct (circumradius, triple ranges). Prop 5.2 covolume/Minkowski re-done: correct.
- Lemma 5.3 and 4-point remark re-done: correct (sorted values, triple ranges a4-a2, a4-a1, a4-a1, a3-a1).
- MAJOR: Prop 2.2 / Open point 6.3 / the "window omega(M log M) <= lambda_max(D) <= O(M^{3/2} log M)":
  p3lower Thm A (refereed, valid for A_X hence actual) gives W >= 2 lambda_max(D) - 2(M-1)c_C^+, so on
  near-extremal clusters lambda_max(D) <= (c_0/2) M log M + (M-1)c_C^+. The window is EMPTY; open point 3 is
  already answered (yes, O(M log M)). Prop 2.2 survives only as a bound on max_x D_x (point totals), where
  it is a genuine but different statement (lambda_max <= max_x D_x is the wrong direction to transfer).
- Thm 1.1 novelty: p3lower Def 1.1 bullets + sharp Prop 1.2 already have: reduced chord, q^a | reduced chord <=> collision,
  and for monomials q^a | D_s(A_v) nonzero rational integer. Thm 1.1 = pair specialisation with EQUALITY beta = i^k gamma_s D_s.
  Notation clash: p3lower's B_xy is the Gaussian zhat_x/g (|B|=e^{d/2}); L1L2's B_xy is the rational D_s(A_v).
- Edge cases of Thm 1.1 checked by hand: v=0 (rotation by i^s, B=+-1); v!=0 => A_v not real/not on diagonal => B!=0.
- Data reading 2 ("every margin grows linearly in W") contradicted by author's own known_out: min_mu = -2.73 constant
  along golden family, -1.76 along Pell.
- stats.py "K quartiles" are percentiles [25,75,95]; matches write-up wording. M=4 clusters with mubar<=0 give K=inf;
  "without complementary pairs" uses n=22784 of 22802 (18 dropped, unexplained).
- Cor 3.4.1 drops "M >= M_0(C)" of Thm 2.1(b): small near-extremal clusters with mubar<=0 violate L1 for every K.
## Computations done
- records.py (own toolkit refgauss.py): 9pt, 10pt (C=8.2393, NOT in write-up), 7pt, golden-8 records: Thm 1.1 formula,
  exact (1.1) in integers, gcd formula, bisector parallelism u_x+u_y || g or g(1+-i), congruence mod q^a (q<=200,a<=3)
  for normalised AND original points, Prop 3.1 dlog code (q<=60): all pass. Table values reproduced (A, mubar, max, K).
  Natural "K without complementary pairs" (exclude from max and mean): 30.7 (9pt), 25.0 (golden) vs author's K2 9.0, 5.17.
- scan_ref.py: independent re-scan N<=1e9 (own method: numpy generation, midpoint-octant dedup): 2,621,071 N,
  22,907 clusters (22,802 M=4; 105 M=5), IDENTICAL cluster set to author's s9.pkl (pickle checked: pure data opcodes);
  max |C diff| 1.6e-13.
- stats_ref.py: 137,862 pairs: Thm1.1 failures 0, bisector failures 0, Lemma 1.3 violations 0, min slack 4.16.
  minC/minA/median K/percentiles/complementary fractions all reproduced. 19 M=4 clusters have mubar<=0 (K undefined).
  Author-def K2 median 1.17 reproduced (n=22783 vs author 22784). Natural K2 median 5.00 (M=4; 8372 undefined since
  non-comp mean <=0), 4.17 (M=5). Clusters WITHOUT any complementary pair: median K 2.67 (M=4, n=7194), 3.54 (M=5, n=72).
  => "non-uniformity driven by symmetry / drops to 1.17" is an artefact of the K2 definition.
- Prop 2.3(c) exponent at log M=100, eps=.4, eta=.5: 0.692 (lam^k/k! form) -- matches "about 0.69".
- records2.py: Thm A ceiling vs Prop 2.2 on records (9pt: lam_max(D)=14.2, ThmA ceiling 37.9, Prop2.2 first line 57.4);
  area identity exact on all triples; relation lattice rank M-3, LLL shortest |c|=2.83..6.93 << Minkowski bound;
  Lemma 5.3 no violation.
- golden7.py: golden n=7 member: 24 of 28 margins IDENTICAL to n=4 member (-2.73..2.73); only 4 grow (12.61, 3x W/2).
  => Reading 2 "every margin grows linearly in W" FALSE. K->7 = 28/4 is a counting artefact.
- neg_mubar.py: 26 M=4 clusters with mubar = 0 (to rounding), A down to 2.27; L1 fails for every K there (small M).
