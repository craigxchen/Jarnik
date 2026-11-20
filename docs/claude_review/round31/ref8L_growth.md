# Referee report on Astra Part VIII, sections VIII.5 and VIII.6 (growth block)

Scope. This report covers sections VIII.5 and VIII.6 of `astra_lemma_VIII/answer.md` (lines 1870 to 1914) and items 9 to 11 of her summary. To read them I also used her VIII.4, for the definition of (VIII.6), together with `round8/weakest.md` §5, `lit18/README.md` §2.2–2.3 and `overview/overview.md` §6.7. All scripts are my own and live in `round31/ref8L_growth/`. I did not import Astra's script, which in any case contains no check of VIII.5 or VIII.6.

## 0. Verdicts

| claim | verdict | label |
|---|---|---|
| (1) Incidence bound. A ratio u ≠ 1 is realised by at most 1+C² ordered pairs on one critical arc. A hyperplane contains at most max(2,1+C²)M^{k-2} distinct tuples. The covering number lies between c_{k,C}M² and M(M−1). | CORRECT | Proved; sanity-checked by computation |
| (2) Conditional growth theorem with relative cutoff H ≥ H_0(∏p)^γ, 0<γ<1/2 | CORRECT | Proved |
| (2') "Is γ < 1/2 exactly what is needed?" | The proof needs exactly γ < 1/2. Sharpness is not shown (see §2.3). | Proved / Gap |
| (3) VIII.6 reduction: m^{2ℓ} ≪ F_δ(P), hence M ≪ 1+(log R/loglog R)^{A/(2ℓ)} | CORRECT | Proved |
| (4) Integer form of F_δ | CORRECT | Proved; computed two independent ways |
| (5) Plausibility of F_δ(P) ≪ (1+\|P\|)^A | Not refuted and heuristically true, even with A = 0. It is not implied by abc, and no known conjecture implies it. | Heuristic + Computed |
| (6) Novelty | The VIII.5 hypothesis is the record's QS^prod. VIII.5 relaxes the needed bound from o(s) to o(s²). VIII.6 is new. | — |

Main result. Every proof in VIII.5–VIII.6 is correct as written. The sharpened incidence bound is new: it lowers the requirement of the record's open QS^prod hypothesis (lit18 §2.3) from f(s)=o(s) to f(s)=o(s²). The reduction in VIII.6 is also new. It replaces QS^prod with a one-variable count F_δ(P) of norm-one P-units near 1. Through ℓ-fold products, any polynomial bound with exponent independent of δ gives growth at every power of log R. This hypothesis is a weaker statement than QS^prod(2,δ,s^A), which implies it.

Status of the hypothesis. A random-angle model predicts F_δ(P) = O_{δ,γ}(1) for every P I can model. A = 0 would already give the uniform conjecture through VIII.6, so the hypothesis is not a "weaker standard conjecture". It is open, non-standard, and heuristically as safe as the uniform conjecture itself. The model also predicts that the implicit constant is astronomically large. In that model F_δ grows exponentially in |P| for P equal to the first s split primes, up to s ≈ 10^4–10^6 when δ, γ are near 1/2, and much further for small γδ. So no finite computation can support the hypothesis. My exact computations agree with the model to within a few percent and show no anomaly. Every structured P I tried behaves like the model or produces only O_{δ,γ}(1) extra solutions. This covers P built from norms t²+1 with y = 1, from unions of unrelated clusters, and from PTE-type cancellation. I found no refutation of either the VIII.5 or the VIII.6 hypothesis.

## 1. The incidence bound (VIII.5, first PROVED block): Proved

Setting. X is a set of M distinct Gaussian integers of norm N lying in an arc of angular width θ ≤ C N^{-1/4}, which is a critical arc of length C√R with R = √N. Fix k ≥ 2 and M ≥ 2k.

Ratio multiplicity. Let u = z/z' ≠ 1 with z, z' ∈ X, and write u = A/B with A, B coprime in Z[i]. Then the following hold.
- |A| = |B|, and |u−1| = |A−B|/|B| ≥ 1/|B|. Also |u−1| ≤ θ ≤ C N^{-1/4}, because the chord is at most the arc. Hence |B| ≥ N^{1/4}/C.
- From zB = z'A and gcd(A,B) = 1 it follows that A | z. So z = gA and z' = gB.
- The multipliers g lie on the circle of radius √N/|B|, in an arc of length at most θ√N/|B| ≤ C N^{1/4}·C N^{-1/4} = C².
- Distinct lattice points are at distance at least 1, and an arc is at least as long as its chord. So there are at most 1 + C² multipliers g, hence at most 1 + C² pairs.

Hyperplanes. Take a hyperplane Σc_iX_i = 0. If only one coefficient is nonzero, it contains no tuple, because marks are nonzero. Otherwise pick two indices i_0, i_1 with nonzero coefficients and freeze the other k−2 coordinates (at most M^{k-2} choices). This leaves az + bz' = S.
- If S ≠ 0, then az lies on two circles, |w| = |a|√N and |w−S| = |b|√N. Their centres are distinct, so there are at most 2 solutions.
- If S = 0, the ratio z/z' = −b/a is fixed. Distinct pairs need ratio ≠ 1, so there are at most 1+C² pairs.

So each hyperplane contains at most max(2,1+C²)M^{k-2} distinct tuples. Every proper subspace lies in a hyperplane. Since (M)_k ≥ (M/2)^k for M ≥ 2k, at least c_{k,C}M² subspaces are needed. For the upper bound, the hyperplanes X_1 = uX_2, one for each realised ratio u = z_1/z_2 with z_1 ≠ z_2, cover every distinct tuple. There are at most M(M−1) of them.

Computed (sanity checks only, since the statement is elementary).
- `g1_arcs.py`: 400 random circles with N a product of 3 to 7 split primes and exponents in {1,2,3}, at C = 1, 2, 3. Pairs per ratio were counted per arc, over every window that starts at a mark. There were no violations.
- `g1b_ratio.py`: targeted circles N = n(a²+1) with n ∈ {25, 65, 325, 1105, 5525}. Here ratio multiplicities above 1 occur, and the maximum of (multiplicity − 1)/C² was 0.125 ≤ 1.
- Caution. A first version of my check pooled ratios over the whole circle and reported multiplicities up to 192. That is not a counterexample: the pairs (z,z') and (εz, εz') share a ratio but lie on different arcs. The claim is per arc and is correct.

Comparison with the record. Round 8 (Theorem 5.2) used the generic bound of M^{m-1} tuples per hyperplane and so needed f(s) = o(s). The arc geometry gives M^{k-2} per hyperplane, so f(s) = o(s²) suffices. This improvement is genuine and is new.

## 2. The conditional growth theorem with relative cutoff (VIII.5, second PROVED block): Proved

### 2.1 Statement (as I verified it)

Fix k ≥ 2, 0 < δ < (k−1)/2, 0 < γ < 1/2, a number H_0 > 0 and a function f. Assume the following: for every finite set P of split rational primes, the projective solutions x ∈ P^{k-1}(Q(i)) of (VIII.6) with H(x) ≥ H_0(∏_{p∈P}p)^γ lie in at most f(|P|) proper linear subspaces. Here (VIII.6) uses the system of round 8 §5: at infinity the forms X_1, X_i − X_1; at the places above P the coordinate forms; normalised absolute values; and H(x) = ∏_v max_i|x_i|_v. Then every critical arc of every circle of radius R satisfies M ≤ max(2k, c_{k,C}^{-1/2}f(s)^{1/2}) once N' ≥ N_0(C,k,δ,γ,H_0). Here N' is the normalised norm and s = ω(N') ≤ (2+o(1)) log R/loglog R. Consequently f(s) = o(s²) gives M = o(log R/loglog R), and f(s) ≪ (1+s)^{2−ξ} gives M ≪ (1+log R/loglog R)^{1−ξ/2}.

### 2.2 Step-by-step check

1. Normalisation by the Gaussian gcd G of all marks. Put N' = N/|G|². If an inert prime q divided N', then q would divide z z̄ and hence z, for every normalised mark z. The same holds for 1+i, since 1+i and its conjugate are associate. Either contradicts gcd 1. So every prime factor of N' splits. The angular width is unchanged, θ ≤ C N^{-1/4} ≤ C N'^{-1/4}. No Gaussian gcd cost is discarded: the whole cluster is divided, and the tuple gcd is accounted for in step 2.
2. Height window. For a distinct k-tuple, the tuple gcd g divides a nonzero chord of length < θ√N' ≤ C N'^{1/4}. So H = N'/|g|² satisfies N'^{1/2}/C² ≤ H ≤ N'.
3. Relative cutoff. P is the set of primes dividing N', so ∏_{p∈P}p ≤ N'. Then H ≥ H_0(∏p)^γ holds as soon as N'^{1/2−γ} ≥ C²H_0. This is uniform in the tuple and the arc, depends only on (C, H_0, γ), and needs γ < 1/2.
4. (VIII.6) by the product formula. x = z/g is primitive, so max_i|x_i|_w = 1 at every finite place. Each x_i is supported above P with norm H, so the finite part equals H^{-k} exactly. The archimedean part is ∏_{i≥2}|z_i − z_1|²/N' ≤ (C²N'^{-1/2})^{k-1}. The total is at most C^{2(k−1)}H^{-k-(k-1)/2}, which is ≤ H^{-k-δ} once H^{(k−1)/2−δ} ≥ C^{2(k−1)}. That follows from step 2 for large N', since δ < (k−1)/2.
5. Homogeneity. Membership of z in a linear subspace is equivalent to membership of z/g, so the cover passes to the raw grid. §1 then gives c_{k,C}M² ≤ f(s) when M ≥ 2k.
6. Counting primes. The s primes are distinct and the j-th is at least j+1, so s! ≤ N' ≤ R². Hence s ≤ (2+o(1)) log R/loglog R. If N' < N_0, the circle has at most 4d(N') = O(1) points.

The rates follow. One point on quantifiers: in "f(s) = o(s²) gives M = o(log R/loglog R)", write √f(s) ≤ εs + K_ε, so the o(·) holds as R → ∞ for fixed C.

### 2.3 Is γ < 1/2 exactly what is needed?

Proved: the argument works for every γ < 1/2 and for no γ > 1/2 without new input. Its only lower bound for H is H ≥ N'^{1/2}/C², while ∏p can be as large as N' when N' is squarefree. At γ = 1/2 the argument works only if H_0 ≤ C^{-2}, which is not uniform in C.

The pair-level bound is attained before normalisation (Computed, `g5_gamma.py`). Take z_1 = g(a+i) and z_2 = g(a−i) with a even, |g|² = n a prime ≡ 1 (mod 4) just below a², and N = n(a²+1). Then the angular distance times N^{1/4} equals 2.0000, the pair gcd is g, H = a²+1, N is squarefree, and log H/log rad N = 0.500000 for a = 10³, 10⁵, 10⁷, 10⁹.

Gap: this does not prove that γ < 1/2 is necessary for the theorem. On a two-mark arc the arc-gcd normalisation removes g, so after normalisation H = N'. Sharpness would need an arc of arc-gcd 1 that contains a pair of height O(N'^{1/2}) with N' squarefree and N' → ∞. My search (`g5b_gamma.py`, 3000 random squarefree N with 3 to 6 split primes below 3000) found arc-gcd-1 arcs at C = 3 with min log H/log N = 0.38. That is consistent with N^{1/2}/C² at these small sizes but is not an asymptotic example. For k ≥ 3 the tuple gcd is usually much smaller, so a larger γ may be admissible for most tuples. The answer to the question is therefore: γ < 1/2 is exactly what this proof uses; whether γ ≥ 1/2 is excluded is undetermined.

### 2.4 Scope paragraph ("PROVED (scope of the test)")

This is correct. In VIII.4 log H ≤ 4j log(2X) = o(s log X), while log∏P ≥ s log X. So the VIII.4 constructions lie below the relative cutoff and leave the relative hypothesis open. The record already said this for Props A and B (lit18 §2.3: Prop A heights exp(O(√(s log s))), Prop B heights (∏p)^{1/s}).

## 3. VIII.6, the reduction to the norm-one torus: Proved

Definitions (as verified). F_δ(P) is the number of u ∈ Q(i) with |u| = 1, u ≠ 1, and v(u) = 0 at every finite place not above P, satisfying H(u:1) ≥ H_0(∏_{p∈P}p)^γ and |u−1|² ≤ H(u:1)^{-δ}. For u = A/B reduced, H(u:1) = |B|². Fix 0 < γ < 1/2, ℓ ≥ 1 and 0 < δ < 1/(2ℓ).

Steps.
1. Subdivision. Cut the arc into ⌈C⌉+1 pieces. One piece of angular width < N'^{-1/4} retains m ≥ M/(C+2) marks, and every chord in it is < N'^{1/4}.
2. Gcd formula. For p | N' with p = ππ̄ and e_p = v_p(N'), a normalised mark has v_π = a_{ip} and v_{π̄} = e_p − a_{ip}. The gcd exponents are therefore min(a_{ip}, a_{jp}) at π and e_p − max(a_{ip}, a_{jp}) at π̄, and |gcd(z_i,z_j)|² = exp(W − D_ij). Since the gcd divides a chord of length < N'^{1/4}, D_ij > W/2 for all i ≠ j.
3. Affine independence (barycenter argument). Suppose Σλ_i a_i = 0 and Σλ_i = 0 with λ ≠ 0. Then there are probability measures μ⁺, μ⁻ with disjoint supports and equal barycenters. Draw i ~ μ⁺ and j ~ μ⁻ independently, so i ≠ j and E D_ij > W/2. At each p put X = a_{ip}/e_p and Y = a_{jp}/e_p in [0,1]. They have the same mean μ and are independent. The inequality |X−Y| ≤ X+Y−2XY holds on the unit square: for X ≥ Y it reduces to 2Y(1−X) ≥ 0. Hence E|X−Y| ≤ 2μ − 2μ² ≤ 1/2, and summing over p gives E D ≤ W/2, a contradiction. Computed: all 756 clusters of width < N^{-1/4} in `g1_arcs.py` have full affine rank and D > W/2. These clusters have m ≤ 2 on random circles, so the check is weak; the proof is what matters. A corollary is m ≤ s+1.
4. Distinct ratios. Take I, J disjoint with |I| = |J| = ℓ. Equal ratios give equal π_p-valuations, so (1_I − 1_J) − (1_{I'} − 1_{J'}) is an affine relation, because its coefficients sum to 0. It therefore vanishes, and disjointness recovers I and J. u = 1 is excluded in the same way. There are C(m,ℓ)C(m−ℓ,ℓ) ≫_ℓ m^{2ℓ} ratios when m ≥ 4ℓ.
5. Closeness. Pair numerator and denominator marks. Each quotient has real argument difference < N'^{-1/4}, so the real arguments sum to less than ℓN'^{-1/4} and |u−1| ≤ ℓN'^{-1/4}. The common arc centre cancels exactly, so no drifting phase is replaced by zero.
6. Height bounds. H = |B|² ≤ N'^ℓ. Also |B| ≥ 1/|u−1|, so H ≥ N'^{1/2}/ℓ². The relative cutoff is passed for large N', because ∏p ≤ N' and γ < 1/2. Moreover |u−1|² ≤ ℓ²N'^{-1/2} ≤ ℓ²H^{-1/(2ℓ)} ≤ H^{-δ} for large H, because δ < 1/(2ℓ). Each u is a P-unit of modulus 1.
7. Conclusion. m^{2ℓ} ≪ F_δ(P) ≪ (1+s)^A, so M ≪ 1 + (log R/loglog R)^{A/(2ℓ)}. For the final claim choose ℓ > A/(2τ) first and then δ < 1/(2ℓ). This uses that A is independent of δ, as stated. The constants depend on (C, ℓ, δ, γ, H_0, A) and on the implicit constant of the hypothesis; there is no hidden dependence on R.

Remark (Proved, immediate). Dimension two in (VIII.6), with x = (b̄, εb), has left side |u−1|²H^{-2}. So QS^prod(2, δ, f) implies F_δ(P) ≤ f(|P|). The VIII.6 hypothesis is therefore implied by the record's QS^prod at m = 2. It restricts to equal norms, which lit18 §2.3 identified as the input QS^prod over-asks.

## 4. Integer form: Proved and Computed

A norm-one u = A/B in lowest terms with all primes in P has the following structure. For each p, if π^a ∥ A then π ∤ B, so π̄^a ∥ B, and π̄ ∤ A. Hence A = ε'b and B = ε''b̄ with gcd(b, b̄) = 1. The prime 1+i cannot occur, because equal valuation would contradict coprimality. So u = εb/b̄, H = h = |b|², and |u−1|² = |εb − b̄|²/h. Replacing b by ηb changes ε by η̄², which is ±1. Each u therefore corresponds to exactly 4 pairs (b, ε), and u ≠ 1 is equivalent to εb ≠ b̄. For b = x+iy the four values of |εb − b̄|² are:
- ε = 1: (2y)², that is 4y²;
- ε = −1: 4x²;
- ε = i: |−(x+y) + i(x+y)|² = 2(x+y)²;
- ε = −i: 2(x−y)².

These are verified by hand. Computed: the meet-in-the-middle exact counter `g2_F.py` enumerates signed exponent vectors, filters with floats under an explicit angle margin, and verifies every candidate exactly with integers. The brute-force counter `g2_brute.py` loops over x, y directly in the integer form. They agree: P = {5,13,17}, h ≤ 3·10⁶, γ = 1/4 gives 70 for δ = 2/5 and 102 for δ = 3/10; P = {5,13,17,29}, h ≤ 2·10⁶, δ = 49/100 gives 92.

## 5. Plausibility of the hypothesis F_δ(P) ≪ (1+|P|)^A

### 5.1 What it would give, and what implies it

- Proved. A = 0, meaning F_δ(P) ≤ K_δ, gives the uniform Cilleruelo–Granville bound through VIII.6, because m^{2ℓ} ≪ F. So the hypothesis is not "weaker standard": with A = 0 it is at least as strong as the target.
- Proved (abc is insufficient). Truncated Vojta, which is abc on P¹ minus {0, 1, ∞}, applied to u = εb/b̄ gives only rad(h) ≥ h^{δ/2−ε}. The relative cutoff gives only rad(h) ≤ ∏P ≤ h^{1/γ}. These are compatible for every γ < 1/2 and δ < 1/2, and abc in any case bounds heights, not counts. This agrees with Astra's statement that she could not derive the count from abc.
- Proved (CG does not imply it). The points b̄ and εb lie on the circle of radius √h at distance ≤ h^{(1−δ)/2}. That is an arc much longer than √R, and different u lie on different circles.
- I know of no conjecture in the literature that implies a polynomial-in-|P| count above a relative height cutoff. I found none in a short search.

Literature, read versus snippet. I did not read any paper for this section. From search snippets only: Konyagin–Soundararajan, "Two S-unit equations with many solutions", arXiv:math/0604453, construct S with superpolynomially many solutions of a+b = c (≥ exp(s^{2−√2−ε}) as quoted elsewhere), improving Erdős–Stewart–Tijdeman. Those solutions have heights far below (∏S)^γ, so, like VIII.4, they do not touch the relative cutoff. That last statement is Heuristic, because I did not read their construction.

### 5.2 Random-angle heuristic

Model. arg b mod π/4 is uniform and independent over the signed exponent vectors e (b = ∏π_p^{e_p}, h = ∏p^{|e_p|}). Then

E F_δ(P) = Σ_{e : h ≥ (∏P)^γ} (8/π) arcsin(h^{-δ/2}/2) ≈ (4/π) Σ_{h ≥ Y_0} 2^{ω(h)} h^{-δ/2}.

Rankin's trick with σ = δ/2 gives

log E ≤ log(4/π) + min_{0<η<σ} [ −ηγ log∏P + Σ_{p∈P} log((1+p^{−(σ−η)})/(1−p^{−(σ−η)})) ].

Dense P (all split primes up to X, s ≍ X/log X; Heuristic). Take η = σ − 1/log X. The prime sum is O(s loglog X), while ηγ log∏P ≍ σγ s log s. So log E → −∞ like −(σγ/2+o(1)) s log s. The model predicts F_δ → 0 for large s, that is F_δ(P) = O_{δ,γ}(1), consistent even with A = 0. The decay sets in only for log s ≫ loglog s/(σγ). Computed, `g3_heur.py large`, values of (log E)/s:
- δ = 0.49, γ = 0.49: +0.98 at s = 10, +0.17 at s = 10³, +0.012 at s = 3·10³, −0.15 at s = 10⁴, −0.69 at s = 10⁶;
- δ = 0.45, γ = 0.25: turns negative near s ≈ 10⁵;
- δ = 0.30, γ = 0.10: still +0.17 at s = 10⁶.

So in the model the admissible constant K_δ behaves like exp(c·s_*(δ,γ)) with s_* ≈ exp(C/(γδ)). This matches the lit18 judge's finding for QS^prod (exp(Ω(s)) up to s ≈ 5·10³) and is more extreme for small γδ.

Sparse P (s primes near Q, Heuristic). Solutions above the cutoff need Σ|e_p| ≥ γs, so E ≲ C(s,γs)(2Q^{−σ})^{γs} → 0, both as Q → ∞ at fixed s and as s → ∞. Computed (Rankin bound, `g3_heur.py sparse`, available in the script): the bound is negative already for Q ≥ 10⁵.

### 5.3 Exact counts, P = the first s split primes (Computed, `g2_F.py`, h ≤ e^{80})

The model predicts that more than 99.5% of the mass has log h < 80 at these δ (checked with `expected()` up to log h = 2500). Every solution was verified exactly with integers.

| s | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| δ=0.49, γ=0.10 | 68 | 170 | 394 | 900 | 2082 | 4570 | 9196 |
| δ=0.49, γ=0.25 | 66 | 166 | 386 | 888 | 2062 | 4538 | 9138 |
| δ=0.49, γ=0.49 | 60 | 154 | 358 | 832 | 1946 | 4292 | 8646 |
| model δ=0.49, γ=0.25 | 61.5 | 161 | 389 | 920 | 2040 | 4410 | 9150 |
| δ=0.45, γ=0.25 | 76 | 212 | 560 | 1432 | 3516 | 8220 | 18340 |
| δ=0.45, γ=0.49 | 70 | 200 | 532 | 1374 | 3386 | 7934 | 17724 |
| model δ=0.45, γ=0.25 | 79.1 | 223 | 581 | 1480 | 3530 | 8200 | 18300 |

Growth in s. F roughly doubles with each added prime, which is exponential in this range. The model reproduces it to within about 5%. This is the pre-asymptotic hump of §5.2, not evidence against the hypothesis: the hypothesis is asymptotic with a δ-dependent constant, and the model predicts eventual decay. The relative cutoff barely matters at these s, because (∏P)^{0.49} ≤ e^{15} for s = 9.

### 5.4 Adversarial and structured P (Computed, `g4_struct.py`)

Each line gives the exact count up to h ≤ e^{log∏P+12}, then the model expectation over the same range, then the model tail beyond.

(a) P = the primes of t²+1 for L even t starting at t_0, where b_t = t+i has y = 1.
- t_0 = 10³, L = 6, s = 10: δ = 0.49 gives 152 / 54 for γ = 0.25 / 0.49, against model 107 / 22. δ = 0.3 gives 5244, against model 4817.
- t_0 = 3·10⁴, L = 4, s = 7: δ = 0.49 gives 22, against model 3.2. δ = 0.3 gives 154, against 72.
- t_0 = 10⁶, L = 3, s = 10: δ = 0.49 gives 480, against 454. δ = 0.3 gives 21786, against 21510.

The structure adds a bounded excess on top of the model.

(c) P = {t²+1 prime}, with L even values t ≈ 10⁶. All primes are ≈ 10¹², so the model predicts 0 and every solution comes from structure (products ∏b_t^{±1}).

| L = s | 2 | 4 | 6 | 8 | 10 |
|---|---|---|---|---|---|
| δ=0.49, γ=0.10 | 6 | 20 | 44 | 82 | 110 |
| δ=0.49, γ=0.25 | 6 | 16 | 32 | 34 | 14 |
| δ=0.49, γ=0.40 | 6 | 12 | 2 | 10 | 4 |
| δ=0.30, γ=0.25 | 12 | 126 | 546 | 1528 | 3540 |

This matches a structural bound (Heuristic, with an explicit mechanism). Write u_c = ∏(b_{t_j}/b̄_{t_j})^{c_j} with c ∈ Z^L and t_j = t_0 + d_j, d_j ≪ t_0^{o(1)}. If the first V moments Σc_j d_j^r (r < V) vanish, then arg u_c ≍ t_0^{−(V+1)}. The approximation condition then forces |c|_1 ≤ (V+1)/δ, and the cutoff forces |c|_1 ≥ γL. A ±1 vector with V vanishing moments needs |c|_1 ≥ 2V, so the count is O_{δ,γ}(1) and vanishes for L > (1/δ+1)/γ·O(1). At δ = 0.49, γ = 0.25 the observed peak is at L = 8 (34 solutions), followed by 14 at L = 10.

(b) P = the union of the primes of r unrelated two-mark critical clusters (N_j = n_j(a_j²+1)), r = 1..4. The counts follow the model: for example r = 3, s = 10, δ = 0.3, γ = 0.25 gives 6886 against model 6847, and δ = 0.49 gives 152 against 152.1. r = 5 exceeded memory and was not run.

### 5.5 The subquadratic-subspace hypothesis of VIII.5 with several clusters (Heuristic)

Take P = P_1 ∪ … ∪ P_r for r clusters of norms N_j ≈ N_0. Single-cluster tuples have H ≤ N_j < (∏N_j)^γ once r > 1/γ, so they fall below the cutoff. Products of marks across clusters have archimedean factor ≈ (Σθ_j)^{2(k−1)} ≈ N_0^{−(k−1)/2}. They need this to be ≤ H^{−δ} with H ≥ N_0^{γr}, which forces r ≤ (k−1)/(2γδ). So for fixed (k, δ, γ) only boundedly many clusters combine, and the covering number stays bounded by a function of the largest known cluster sizes. This does not refute f(s) = o(s²). The k = 2 computation in (b) shows nothing beyond the random model.

### 5.6 Conclusion on plausibility

The hypothesis is not trivially false: no construction I know of, nor the record's (VIII.4, Props A and B, PTE chains), exceeds O_{δ,γ}(1) above the relative cutoff. It is heuristically true even with A = 0. It is not implied by abc or by the uniform CG conjecture, and no known conjecture implies it. It is a non-standard uniformity in S, and at A = 0 it is as strong as the target. The random model predicts constants so large that no computation can test it. I could not rule out adversarial P: a union bound over all subsets P fails, so the heuristic does not control sup_P F_δ(P).

## 6. Novelty relative to the record

- The relative-height hypothesis is not new. It is the record's QS^prod (`lit18/README.md` §2.3; `overview/overview.md` §6.7: "for some m ≥ 2, 0<δ<(m−1)/2, 0<c<1/2 and f(s)=o(s) … H(x) ≥ max(H_0,(∏p)^c)", which says "open, non-standard … Do not spend proof effort on it"). Astra's γ < 1/2 is the record's c < 1/2, for the same reason.
- New: the per-hyperplane bound M^{k−2}, which weakens the requirement from f = o(s) to f = o(s²).
- New: VIII.6. The reduction to a count of norm-one P-units near 1 is implied by QS^prod(2), restricts to equal norms, and gives an arbitrarily small log-exponent from any polynomial count with A independent of δ. The affine-independence lemma is credited to the parallel report. My searches of the record for "norm-one", "S-unit", "relative height", round 8 and round 10 found no prior statement of the torus count or of the ℓ-fold product amplification.
- Route 5 (the determinant method, rounds 10–17) is unrelated.
- The record's statement that "no intermediate growth bound is known from a weaker standard conjecture" stands. The new hypothesis is non-standard, and at A = 0 it implies the uniform bound.

## 7. Files (`round31/ref8L_growth/`)

- `g1_arcs.py`, `g1b_ratio.py`: incidence bound, affine independence, ratio distinctness and height window on genuine circles.
- `g2_F.py`: exact meet-in-the-middle counter for F_δ(P). `g2_brute.py`: independent brute force in the integer form. Outputs: `out_F_049.txt`, `out_F_045.txt`.
- `g3_heur.py`: random-angle model (exact DP histogram and Rankin bound). Output: `out_heur_large.txt`.
- `g4_struct.py`: structured P. Outputs: `out_struct_a.txt`, `out_struct_b.txt`, `out_struct_c.txt`.
- `g5_gamma.py`, `g5b_gamma.py`: the γ = 1/2 boundary.
