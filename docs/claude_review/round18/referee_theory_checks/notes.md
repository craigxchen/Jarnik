# Referee notes (theory.md)

## Read theory.md lines 1-400
- Lemma 2.2 index: per-level |C|-1 conditions; T class at l_T(p) levels; index prod d_T^{|T|-1}. Looks right. sums 196, 294 ok; H.beta=31 (31 visible T contain {i,j}), K.beta=-186+196=10 ok.
- Kapranov K = -6H + sum(|T|-2)E_T: codim L_T = |T|-1 in P^5 -> discrepancy |T|-2. ok.
- Prop 2.3: bounds kappa (finite part) only, NOT h_K (archimedean lambda_{E_T,inf} terms omitted). Summary item (4) says "h_K <= -sum ell_T": overstated. Minor.
- Trivial bound 30 = 31-1 (exclude pair). ok.
- Lemma 2.4: needs check of def of matching coeffs (3.1) and of nested-chain hypothesis (is it on blocks divisible by p? or p-parts?).
- Thm 2.5 Ridout: derived myself via Subspace n=2: S-part xyz <= C H^eps |L1 L2|, so S-free >= H^{1-eps}/C. ok. 16 blocks per matching ok. combination ok.
- TODO: Prop 2.6 abc gain; identity; §3; check "Bal^pc nonempty?" and whether (L_8)=>(L_8^pc) correct (census Prop 1.3(2)).
- TODO: check §1.1 verbatim vs ptolemy §7.

## Read theory.md 400-594; derivations checked by hand
- Prop 2.6 proof OK (Y in (u)^nu over Z by unimodular change; factors of transl.-inv. poly are transl.-inv.; accidental p|d_T dividing b accounted in Y cofactor). gcd(A,B,C)=gcd(A,C).
- identity gain=10-K^log.beta: K^log = K+F_red-L. ok.
- morphisms: L.beta=16 deg h, n>=deg h+2 by RH -> gain<=-32. ok (B-M cited from abstract).
- Mason line bound: re-derived, ok.
- multiplicity form 3L=F: ok; cross ratio 48*(-2/3)=-32 ok.
- F-curves K.F with c(s)=s(8-s)/7-2: (2222)=2,(3221)=1,(3311)=0,(4211)=0,(5111)=-1. ok.
- expected dims: 186d+2-196d; genus1 186d; k=7: 75d+1-68d. ok.
- OBSERVATION (strengthening): random specialization can only OVERestimate nu_T; overestimating nu_T(Y) never lowers computed gain (net 0 or +1). So computed max -32 => true gain <= -32 rigorously, given enumeration complete & factorization exact.
- ISSUE A (rate dichotomy, item 3): f-version => (L8^c) for CATERPILLAR Phi only; (L8^c)^all => (L8). Gap: (L8^c)^cat => (L8)? not shown (data with residues can have non-caterpillar Phi). So "no easier than (L_8)" overstated: correct = f-version sandwiched (L8^c)^all => L8 => f-version => (L8^c)^cat. CHECK census Prop 1.3.
- ISSUE B (§3.1 iff and §3.4 "exactly a refutation of V_8"): Bal dense => not(L8) needs caterpillar. "exactly" overstated; only not(L8) => not V_8 direction.
- ISSUE C (§3.3(1)): "covering families excluded (proved)": genus>=2 covering families each with finitely many points (Faltings) are NOT excluded as a source of Zariski-dense balanced points. Wording.
- ISSUE D: Prop 2.3 bounds kappa (finite part) not h_K; summary claims h_K.
- ISSUE E: "what remains" dichotomy (rigid curves / sporadic) is a description, not a theorem; every point lies on lines.
- TODO: check (L8) verbatim in ptolemy §7, item-532 definition, census Prop1.3, referee_census/lift; independent code.

## Formulation check (ptolemy §1, §3, §6.2, §7 read)
- item-532 chain condition (ptolemy §1): "blocks with p | n_T form two nested chains on disjoint row sets" -> theory (a) matches.
- (L_8) §7 balance "n_T satisfy |log n_T - w|<=eta_0 w" (all T ambiguous); visible-only version equivalent (top block via P_j->(nX_j,Y_j) fresh prime, singletons fresh primes) - census §1.1. OK.
- census Prop 1.3(2) is only for ADMISSIBLE Phi -> confirms ISSUE A/B.
- L8.md §7 already has genus-0 expected dim -10d+2 (novelty: theory table partly not new). L8 Cor A already: Liouville can't prove f(w)->inf either.
- abc barrier scope (round8 weakest Prop 4.2 'Precise scope'; conditional Prop 8.3 needs squarefree pairwise coprime n_i and bounded cofactor): referral to aux-poly barrier = Liouville only. Gap claim (7) is fair.
- Prop C (L_5) uses INTEGRAL chart (inf unmarked), Prop 2.3 is Kapranov chart only. Hand count for k=8 integral chart box: positive density needs sum C(|C_S|,2)/4 >= sum(|C_S|-1): 101.5+3.5a+1.75b >= 245+4a+2b impossible. So conclusion OK for balanced points in that chart too, but Prop 2.3 as stated is chart-specific.
## Independent computations
- ref_caterpillar.py: 26/236/2752/39208 trees; caterpillars 26/221/2227/25698; (a)anchor0=(b)all anchors=(c)=(d internal degree<=2): 0 mismatches. CONFIRMED.
- ref_lattice.py: 80 configs, index via HNF = prod d_T^{|T|-1}, (1.1) holds: 80/80. CONFIRMED.
- ref_gain.py (own code, exact compose nu_T, flint factor): D=1: 40 pencils (-32:20,-48:20); D=2: 1180, max -32 (24); D=3: 19800, max -32 (20), next -40(4),-44(2),-48(208): EXACT MATCH with author. D=4: 22 orbits confirmed; full run in progress (first attempt died: my cache leak, not flint).
- ISSUE F (new): enumeration covers only lam = s = +-1 in f = lam A/C. Unit pencils with other lam in Q^* (e.g. 4z(1-z): B=(x0+x1-2x2)^2, lam=-1/4 in my param) are NOT enumerated, so "all unit pencils of degree <=4 excluded" is overstated. ref_gain_lambda.py: special rational lam (non-reduced B / extra pair factor / higher order on L_T): D<=2: 6 special lam, all gain -32. D=3 to run.
- ISSUE G: §5(2) "equivalent in substance to (L_8)" for (L_8^pc) - overclaim (pc is a subcase); lattice set slightly larger than Bal^pc but inside Bal^cat, so still a consequence of (L_8). minor.
- §3.3(4) 'rigid' curves: should be 'non-covering'; and 'sporadic' dichotomy = census-referee's major issue repeated in weaker form (theory restricts to curves carrying infinitely many points? It says "lying on no such curve" where such = rigid genus<=1 with eta->0). Not exhaustive: curves moving in non-covering families omitted.
- Thm 2.5 checked vs (4.1): Y_i t_jl N_jl - Y_j t_il N_il + Y_l t_ij N_ij = 0, each residue once -> |abc|<=e^{H}. ok.
