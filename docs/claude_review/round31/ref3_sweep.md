# Referee report, round 31, Part III, block "sweep"

Claims under review: astra_lemma_III/answer.md III.6 (first theorem, "bounded block degree excludes density", together with its corollary on leading cofactor shapes), and the last two results of III.10 ("growing block-degree cutoff" for ten-pair templates; "the same degree-cutoff principle reaches fifteen cuts"). These are summary items 5 and 6.

Scripts are in `ref3_sweep/`:
* `k4det.py` and `k4det.out`. These test the log-determinant mechanism on the record's k = 4 family (round20 gauss.md §4.1), using exact SymPy arithmetic.
* `degcount.py` and `degcount.out`. These check the contact and degree arithmetic.

## Verdicts in brief

| Claim | Verdict |
|---|---|
| III.6. For fixed d, the union of all strict degree-d ten-pair template ratio curves has Zariski closure of dimension ≤ 2 in G_m^4 (char 0) | **CORRECT** |
| III.6. Leading cofactor shapes [h_1:…:h_4] lie in a closed set of dimension ≤ 1 | **CORRECT** (an immediate tangent-cone corollary) |
| III.10. Ten-pair templates of degree ≤ d_0 lie on a hypersurface of degree O(6^{24 d_0}) | **CORRECT** |
| III.10. Fifteen-cut templates with contact ≥ ⌈D/2⌉ and max d_A ≤ d_0 lie on a hypersurface of degree O((d_0+1)^{15} 8^{30 d_0+4}) = m^{o(1)}, so they make up o(m^5) tuples | **CORRECT**, with a scope correction (below) |
| Relation to hypothesis (5) of astra_free_III Part I §5 | (5) is neither proved nor reduced. The result only shows that one structured subclass of tuples cannot carry a counterexample |

## 1. The parameter space and why the union is not dense in a 3-fold

Fix d. Take the split coordinates: roots r_{e,ν} and "conjugate" roots r̄_{e,ν}, which are independent (e ranges over the ten pairs, ν = 1..d), so 20d coordinates in all. Monic blocks β_e and β̄_e give stars P_i = ∏_{j≠i} β_{ij}. The template conditions say that Z_{ij} − Z̄_{ij} is a constant 2i t_{ij} for the four anchored pairs, where Z_{ij} = P̄_i P_j / (β_{ij} β̄_{ij}) has degree 6d. These are closed polynomial conditions. The strictness conditions (twenty distinct roots, t ≠ 0) are open. The result is a finite-type scheme T_d. Ordering the roots makes it a finite cover of coefficient space, so T_d has finitely many irreducible components.

The source x ranges over A^1 minus the roots. The map is Φ: T_d × A^1 → G_m^4, u = ∏_e ∏_ν ((x−r_{e,ν})/(x−r̄_{e,ν}))^{b_e}, where b_e ∈ Z^4 is the exponent column of edge e in u_i = q_i/q_5. Explicitly, b_{ij} = e_i + e_j for i,j ≤ 4, and b_{i5} = e_i − (1,1,1,1). These columns span Q^4. The "union of all template ratio curves" is exactly Im Φ.

By Chevalley, Im Φ is constructible. Its closure is the finite union of the closures Φ(component × A^1). In characteristic 0, each of those closures has dimension equal to the generic rank of dΦ on that component. So the whole question is a pointwise rank bound, rank ≤ 2. There is no room for a union that is dense in a 3- or 4-dimensional subvariety, because there are only finitely many components and each has an image of dimension ≤ 2. For unbounded d, a countable union of surfaces can be dense, and Astra says so.

## 2. The proof of rank ≤ 2, checked line by line

1. **Residue structure.** Each column ∂ log u is a vector rational function in x. This holds for the source direction ∂_x and for each parameter direction ∂_p. The poles are simple, at the r_{e,ν} and r̄_{e,ν}, with residues in ⟨b_e⟩. For example, ∂_p log(x−r) = −∂_p r/(x−r). Correct.
2. **Order at infinity.** Monicity and Z − Z̄ = 2it give q_j/q_i = Z/Z̄ = 1 + 2it/Z̄ = 1 + O(x^{-6d}). So u = 1 + O(x^{-6d}), and this holds identically on a component. Hence ∂_p log u = O(x^{-6d}) and ∂_x log u = O(x^{-6d-1}). Correct.
3. **Projection along an edge e.** Let π_e: K^4 → K^4/⟨b_e⟩ ≅ K^3, where K is the function field of the component. The residues at the 2d poles of edge e are killed, so the projected columns are regular there. At each of the other 18d pole locations, every column has a residue on one line π_e⟨b_f⟩. By multilinearity, the 3×3 determinant Δ_e of (∂_x, ∂_{p_1}, ∂_{p_2}) therefore has at most a simple pole there. This also holds under coincidences of roots between blocks: the pole order at a point is at most the rank of the residue span, which is at most the number of root-sets meeting there. So the total is still ≤ 18d. Since Δ_e = O(x^{-18d-1}), writing Δ_e = N/∏(x−ρ) gives deg N ≤ 18d − 18d − 1 < 0, so Δ_e ≡ 0. Correct.
4. **Conclusion.** Suppose the rank is ≥ 3. Then ∂_x log u ≠ 0, since otherwise u ≡ 1, which is impossible for distinct stars. So ∂_x can be completed by two parameter directions to an independent triple spanning a 3-dimensional V. Because the b_e span K^4, some b_e lies outside V, and π_e is injective on V. Then Δ_e ≠ 0, a contradiction. Correct.

**The margin is exactly one, and it comes from the source derivative.** With three parameter directions the order would be only 18d, against 18d poles, and no contradiction would follow. The extra x^{-1} from ∂_x is essential.

I tested this bookkeeping independently on the record's k = 4 family (gauss.md §4.1). That family has 7 linear blocks, 14 poles and contact 4, and it covers T_4, so no determinant there should vanish. `k4det.out` shows:
* The 3×3 determinant (∂_s, ∂_X, ∂_m) has exactly 14 simple poles and order 13 = 3·4+1 at infinity, with numerator of degree 1 = 14 − 13.
* Each projected determinant (modulo 1_T) has 12 poles (the two of edge T are removed) and order exactly 9 = 2·4+1.

So the predicted pole counts and the +1 are reproduced exactly, and the inequality is sharp there.

**Comparison with the record.**
* The covering-family inequality (gauss.md Thm 4.3(ii), toric.md Prop 4.3) says covering ⇒ (k−1)μ_0 ≤ B·C. Here B·C = 20d and μ_0 = 6d, so 24d > 20d and dimension 4 is excluded. Astra states this.
* The record bounds the dimension at 3 only qualitatively: the union of bounded-degree curves *with infinitely many A-points* lies in a proper closed set.
* III.6 lowers the bound to 2 and does not need A-points. The projection-along-a-residue-line device does not occur in the record (grep finds no log-determinant argument of this kind). The result is **new**, but it is an incremental sharpening.
* The bound is far from the expected dimension. The parameter count 20d − 2 − 4(6d−1) is negative, so templates are expected to be empty. III.6 is a statement conditional on existence, not evidence of existence.

**Leading-shape corollary.** Along a template curve, u − 1 ~ 2i(t_{5i})_i x^{-6d}, and this vector is nonzero. So its projective class lies in the projectivised tangent cone at e of the closure of the image. That cone has dimension ≤ 1. Correct and immediate.

## 3. III.10, ten-pair degree cutoff

* There are 4(6d−1) contact equations, each of coefficient degree ≤ 6 (a coefficient of a product of six monic blocks). Heintz–Bézout gives a component-degree sum ≤ 6^{4(6d−1)}.
* There are four graph equations w·Z̄(0) − Z(0) of degree 7.
* The degree of the closure of a linear image is at most the source degree (Heintz 1983, Lemma 2).
* A variety of dimension ≤ 3 and degree δ lies on a hypersurface of degree ≤ δ.
* Translation x → x + c preserves the template class, so evaluation at 0 covers every finite source value.

All of these steps are correct. The cutoff c < 1/(24 log 6) ≈ 0.0233 is correct.

## 4. III.10, fifteen cuts

**The dimension step.** Here there are 2D poles in total, each simple (one residue line per block), and the contact is k at infinity. The 4×4 determinant is O(x^{-4k-1}). With 4k ≥ 2D, i.e. k ≥ D/2, it vanishes, so the image has dimension ≤ 3. This is the record's Prop 4.3 inequality, again sharpened by the +1: Prop 4.3 allows covering when 4μ_0 = B·C, and this argument excludes it.

That sharpening is *needed*. For equal degrees with d even, k = 15d/2 = D/2 exactly (`degcount.out`). The +1 mechanism verified in §2 covers it.

**The degree count.**
* Each pair crosses 8 of the 15 cuts: the two singletons and 2·3 doubletons. So each reduced word has degree 8d.
* A cofactor of degree ≤ ⌊d/2⌋ gives contact 8d − ⌊d/2⌋ = ⌈15d/2⌉, which I checked for d ≤ 6.
* There are at most 4(k−1) contact equations, each of degree ≤ 8, and four graph equations of degree 9.
* With k ≤ ⌈15d_0/2⌉, so 4k ≤ 30d_0 + 2, and at most (d_0+1)^{15} degree vectors, the bound O((d_0+1)^{15} 8^{30 d_0+4}) follows. The constant 9^4·8^2 is absorbed.
* This is m^{o(1)} for d_0 = o(log m). It is o(m) already for d_0 ≤ c log m with c < 1/(30 log 8) ≈ 0.016. Astra does not claim this, but it is true.
* Homogenising in the anchor and applying Schwartz–Zippel on S^5 gives ≤ m^{4+o(1)} tuples.

**CORRECT.**

**Scope correction.** "Evaluations of templates" means u(x_0) for split templates and x_0 ∈ C. A Gaussian template evaluated at a *non-real* Gaussian integer x_0 = a + bi does not give u(x_0). It gives the real evaluation at a of the family with oppositely shifted roots (r − bi, r̄ + bi), and that family need not have the contact. Such evaluations are covered only when the shifted family is itself in the class. If it is not, then generically |u − 1| is of order b/|x|^2 there, so they are not short-arc for large |x|. The statement is therefore correct as stated, about the image of the split family. Applying it to lattice tuples needs real (or effectively real) source values.

## 5. Relation to hypothesis (5) of astra_free_III Part I §5

Hypothesis (5) asks, at each parent, for one polynomial F_{N,B} of degree o(W/log W) = o(m) that vanishes on all but o(m^5) classes of B_B(N). The fifteen-cut result gives such a polynomial, of degree m^{o(1)}, *for the subclass of classes that are evaluations of fifteen-cut monic templates with max block degree ≤ d_0 = o(log m) and contact ≥ D/2*. The covering is uniform in the template coefficients, and that uniformity is the real content: a single parent needs about m^5 different templates. Schwartz–Zippel then shows directly, without (5), that on a failure sequence at most o(m^5) of the ≫ m^5 extracted tuples can be of this kind. So template-generated growth counterexamples *of exactly this type* are excluded.

The contact hypothesis is the natural one. For a fixed template with |x| → ∞, N ≍ |x|^{2D} up to bounded gcd factors, and the displacement is ≍ |x|^{-k}. So the short-arc condition ≤ ½N^{-1/4} forces k ≥ D/2, i.e. α = B·C/μ_0 ≤ 4. The cofactor condition H ≤ B log W forces constant cofactors in this regime, which is stronger than required.

The result does **not** exclude the following:
1. Templates of block degree ≳ log m. A parent with bounded x and large d is of this kind.
2. Templates whose short-arc property occurs at special finite x_0 rather than through contact at infinity. For example, large coefficients, or high-degree cofactors that happen to be small at the evaluation point. The asymptotic argument forcing contact needs |x_0| large relative to the coefficients, and the extraction gives no such control.
3. Multi-variable families.
4. The generic case: isolated tuples with no polynomial deformation. Astra says this explicitly.

The extraction in astra_free Part I §5 produces arbitrary genuine tuples and gives no degree bound or polynomial representation. So nothing here proves (5) or reduces it to templates.

## 6. Significance

* Growth rate: unchanged.
* Uniform bound M(C): unchanged.
* Question A: unchanged.
* Route 5: unchanged.
* Question B: the refutation-by-density route through bounded-degree ten-pair templates was already blocked in dimension 4 by the record. It is now blocked in dimension 3 as well. Whether any template exists remains open, so this changes no status.

**Peripheral.** The one genuinely useful technical item is the elementary log-determinant inequality (k−1)μ_0 + 1 ≤ B·C for covering families, with its projected version. It is cleaner than the record's sheaf argument, slightly sharper, and effective in the degree.

None of the standing cautions is violated: no fixed-S Subspace theorem is used, no gcd is discarded, no phase is set to zero, and no threshold is claimed to be uniform beyond what is stated.
