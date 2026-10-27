# OpenAI "local smoothing in all dimensions" (Oct 26 2026): what exactly was claimed

Lens: facts and statement only. Written 2026-10-27. Sources were read directly where possible. Where only a
search-engine snippet was available, this is flagged. Nothing here has been peer reviewed. I did not compile the
Lean proof.

## 0. What I could and could not access

| Source | Status |
| --- | --- |
| GitHub repo `openai/math`, shallow clone at commit `7ef5c0cc5ca4b0e072b7a050fca09ab24bc6fa3a` (2026-10-26 14:18 -0700) | **Read directly.** Local copy: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/local_smoothing/dl/math` |
| Paper: `preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/` (`paper.pdf` plus LaTeX `build/paper.tex`, `build/sections/*.tex`, about 4,200 lines) | **Read directly**: intro, history, the transfer/sharpness/consequences/Kakeya statements, and the planar appendix |
| `history.md` (Oct 26 entry), `README.md`, `overview.tex` (family 079) | Read directly (WebFetch of github.com also worked) |
| Lean: `lean/ComparatorChallenges/LocalSmoothing.lean` (statements), `lean/OAI/Analysis/Fourier/LocalSmoothing/*` (proof), `lean/formalization.yaml` | Read directly. I scanned the import closure textually. **Not compiled**: there is no Lean toolchain here, and Mathlib plus about 190k lines would need to be built |
| openai.com blog post "Working with academic communities to build scientific understanding" | **Not fetched** (DNS failure). Known only from search snippets |
| X posts (Dan Roberts, OpenAI), news (Scientific American, the-decoder, Quanta "Transformation", mlq.ai, datastudios, tech-insider, qz, kingy.ai) | **Not fetched** (DNS failure for WebFetch and curl). Search-engine snippets only, so treat their content as unverified |

## 1. The exact claimed theorems

Notation (paper): `U(T) = exp(-i T |D|)` is the half-wave propagator on R^n, `I = [1,2]`, `J = (1-Delta)^(1/2)`,
`p_c = 2n/(n-1)`, `s_+(p) = max{0, (n-1)/2 - n/p}`.

**Theorem 1.1 (Euclidean local smoothing).** For `n >= 2`, `2 <= p < infinity` and `eps > 0`,
`||U(T) f||_{L^p(R^n x [1,2])} <= C_{n,p,eps} ||J^{s_+(p)+eps} f||_{L^p(R^n)}`, with `eps = 0` allowed at `p = 2`.

* This is Sogge's 1991 conjecture in its standard form: every gain below `1/p` over the fixed-time order
  `(n-1)(1/2-1/p)` for `p >= p_c`, and arbitrarily small loss for `2 < p < p_c`.
* The operator is the Euclidean half-wave propagator (constant coefficients). The claim is **not** for variable-coefficient
  FIOs or manifolds. The paper notes that some 3-dimensional metrics need `p >= 10/3` (Minicozzi-Sogge).
* It is **sharp up to eps**. It is **not eps-free** above `p_c`. The paper's own refinement is logarithmic (below).

**Theorem 1.2 (critical annular estimate; the real core).** For `n >= 2`, `p = p_c`, and `P_N = phi(D/N)`,
`||U(T) P_N f||_{L^{p_c}(R^n x I)} <= C [log(2+N)]^(2/p_c) ||f||_{L^{p_c}}`.
* Lower bound (Prop. "Radial focusing"): `(log N)^(1/p_c)` is necessary, so the true power lies in `[1/p_c, 2/p_c]`.
  A weak-type obstruction `(log N)^((n-1)(1/2-1/p))` holds for `2 < p <= p_c`.
* Interpolation gives `N^(s_+(p)) (log N)^(beta(p))` on each annulus, with `beta(p) = min{(n-1)(1-2/p), 2/p}`.
  Corollary "Logarithmic Sobolev estimate": the loss `J^(s_+(p)+eps)` can be replaced by
  `J^(s_+(p)) log(2+J)^B (loglog)^C` for `B > beta(p) + 1/2 - 1/p`. So above `p_c` the paper claims **log loss, not
  zero loss**. The eps-free endpoint conjecture for `p > p_c` (known for `n >= 4`, `p > 2(n-1)/(n-3)` by
  Heo-Nazarov-Seeger 2011) is **not claimed** in general.

**Consequences claimed (paper abstract, Sections 9-12):**
1. Bochner-Riesz conjecture: `sup_r ||B_r^alpha f||_p <= C||f||_p` for all `1 <= p <= infinity` and
   `alpha > max{0, n|1/p-1/2| - 1/2}`. Also the maximal Bochner-Riesz bound for `2 <= p < infinity`, `alpha > s_+(p)`.
2. Stein's square function: bounded for `2 <= p < infinity`, `alpha > 1/2 + s_+(p)` (Lean form: `alpha > max{1/2, n(1/2-1/p)}`).
3. Fourier restriction (extension) for the sphere, in the full **open** range `L^a(S^(n-1)) -> L^b(R^n)` with
   `b > 2n/(n-1)` and `b > (n+1)/(n-1) a'`. The same holds for the truncated cone in R^(n+1) by slicing and the compact
   paraboloid by rescaling, and gives Schrodinger local smoothing in `n-1` space dimensions. Endpoints are excluded
   ("interior of the conjectured range"). There is a log-loss endpoint on balls:
   `||Eg||_{L^{p_c}(B_R)} <= C (log R)^(2/p_c) ||g||_{L^{p_c}(S^(n-1))}`.
4. Kakeya and Nikodym maximal conjectures: `||K_delta f||_{L^n(S^(n-1))} <= C (log 1/delta)^(4/p_c) ||f||_{L^n}`, and every
   Kakeya set in R^n has Hausdorff dimension n.
5. **Not claimed:** decoupling (no new decoupling theorem), Mizohata-Takeuchi, endpoint (eps-free) restriction, Minkowski
   or Kakeya results over finite fields, and anything arithmetic.

The deductions from local smoothing to Bochner-Riesz, to restriction and to Kakeya follow known implication chains
(Sogge 1991, Beltran-Hickman-Sogge 2021 Props 39-40, Tao's "Bochner-Riesz implies restriction", Bourgain's
"restriction implies Kakeya"). The paper re-derives them with explicit log powers. **All new mathematics is in
Theorem 1.2.**

## 2. The method ("novel auxiliary PDE")

From the introduction and Sections 2-7, as summarised by the paper itself:
* **Imaginary translation along light rays.** `X -> X - i l (1, omega)` damps the frequencies of a cap Gaussianly:
  `exp(-l|xi| |omega - xi/|xi||^2 / 2)`, with angular width `(N l)^(-1/2)`.
* **Conserved current.** The fluxes of the symmetric first-order wave system, continued to complex spacetime with two
  auxiliary time variables `tau in R^2`, give a positive semidefinite tensor. Each of its real columns is
  divergence-free in (physical plus auxiliary) variables, and its time component detects the wave on the cap. Precursors
  are Friedrichs (bilinear flux), Kaiser (analytic-signal currents) and Serre (compensated integrability: a
  divergence-free positive tensor paired with a Monge-Ampere potential).
* **The auxiliary PDE.** On each spatial fiber a **quaternionic Monge-Ampere Dirichlet problem** is solved:
  `det H_fib(lambda|Phi|^2 + u) = lambda mu`, `u = 0` on the boundary. Here `mu` is the regularised current density and
  `lambda|Phi|^2` is a rank-one background. The energy `F = int(-u) mu`, with an Aubin-Mabuchi-type correction,
  detects `mu^(n/(n-1))`. Since `mu ~ |wave|^2`, this produces the exponent `2n/(n-1) = p_c`.
* **Scalar inequality.** `Delta_{T,tau} F <= C s^(-2) e^(-cs/h) F` on `2a <= |tau| <= 3`. Two "circle comparisons" in
  the planar `tau` variable lose one log each, which gives `log^2` and hence the `log^(2/p_c)` power.
* **n = 2 (Appendix "planar").** The determinant equation becomes a **linear** Dirichlet problem on a disk. So the paper
  contains a self-contained PDE proof of 2-D local smoothing at `p = 4` with `(log N)^(1/2)` loss.
* Also: a generalisation to a simple extremal eigenvalue of a real symmetric pencil (Appendix "pencil").

## 3. The Lean formalization

* Comparator challenge `lean/ComparatorChallenges/LocalSmoothing.lean` uses Mathlib-only definitions and states six
  theorems: `euclidean_local_smoothing` (exactly Theorem 1.1 for Schwartz `f`, `n >= 2`, `2 <= p < infinity`, eps > 0 or
  p = 2 with eps = 0), `bochner_riesz` (uniform, full strict range including `p = infinity`), `spherical_restriction`
  (strict mixed range), `kakeya_maximal` (`C delta^(-eps)`), `kakeya` (`dimH F = n` for every Kakeya set, measurability
  not assumed) and `stein_square`. Permitted axioms are `propext`, `Quot.sound` and `Classical.choice`.
* Solution: `lean/OAI/Analysis/Fourier/LocalSmoothing/Comparator.lean`, which goes through `Main.lean` and
  `WavePackets/PositivePacketSolvability.lean`, with quaternionic-Hessian, Dirichlet and current modules. The import
  closure is **2,261 modules, about 190k lines**. My textual scan found no `sorry`, `axiom`, `admit`, `native_decide`,
  `implemented_by` or `extern`. Toolchain: `leanprover/lean4:v4.34.1`. **I did not build it.**
* **The formalized log power is weaker than the paper's.** In Lean, `LogarithmicRefinement` is the annular estimate
  with `(1 + log(1/h))^(2n+5)` instead of `log^(2/p_c)`, and `LogarithmicSobolev` needs `B > 2n+6`. So the eps-level
  statements are formalized, but the sharp log powers of the paper (and the `(log)^(4/p_c)` Kakeya bound) are only
  stated in the paper.
* Same-repo context: on Oct 7, three unformalized manuscripts were withdrawn for a sign error. This one does have a
  formal proof, so the remaining risk is mainly that it fails to compile or that the statements were transcribed
  wrongly. The challenge statements look faithful to the standard conjectures as I read them.
* OpenAI had already released, on Sept 23-24, "Critical local smoothing for the three-dimensional wave equation"
  (the n = 3 case, by a different-looking method), "The Kakeya maximal conjecture in three dimensions",
  "Every four-dimensional Kakeya set has full Hausdorff dimension", "Bochner-Riesz in three dimensions" and
  "Elliptic capacity propagation and Fourier restriction to the sphere". The Oct 25 paper is the all-dimensions version.

## 4. Verification status and reactions (as found)

* Not peer reviewed. `history.md` (Oct 26) says: "We are releasing the model output directly." The repo README says
  the manuscripts are "raw research data" and that "some of the unformalized results could have issues".
* Dates in news snippets (unverified): natural-language proof Oct 18, Lean formalization Oct 19, announced Oct 26.
* No named verification group for the local smoothing result was found. The Impagliazzo-Paturi-Williams group is for the
  ETH disproof (confirmed in `history.md`). An IAS advisory group statement concerns Yang-Mills/QCD (snippet only).
* Expert reactions found are about OpenAI's math releases in general, and all come from **snippets only**:
  * Tao is reported as critical of mass AI "harvesting" of open problems and as chairing a group that calls for a
    boycott.
  * Andrew Sutherland: treat one-shot claims as unverified until replicable.
  * Ellenberg and Litt were cautiously positive about the earlier ~350/700-paper release.
  * Guth is quoted (tech-insider) on the 4-D Kakeya difficulty: tubes concentrating along quadrics.
* **No expert assessment of this specific proof was found, positive or negative.** No identified gap was found.

## 5. Dimension 2 (the plane) before the announcement

Here `n = 2` (circles in R^2), so `p_c = 4`.

| Statement at n = 2 | Status before Oct 2026 |
| --- | --- |
| Local smoothing, all `2 <= p < infinity`, eps loss | **Theorem.** Guth-Wang-Zhang, Ann. Math. 192 (2020), via the sharp square function for the cone in R^3 (arXiv 2019) |
| Restriction for the circle (full open range) | **Theorem**, classical: Fefferman 1970, Zygmund 1974, Carleson-Sjolin 1972 |
| Bochner-Riesz in R^2 | **Theorem**, Carleson-Sjolin 1972 (also Fefferman, Cordoba) |
| Kakeya in R^2 | **Theorem**: Davies 1971 (dimension 2), Cordoba 1977 (maximal, with `(log 1/delta)^(1/2)`). The new paper's 2-D bound `(log)^(4/p_c) = log^1` is **weaker** than Cordoba's |
| Cone restriction in R^3 (n = 2 cone) | **Theorem**, Barcelo 1985 |
| Parabola restriction / 1-D Schrodinger | Classical (I did not pin down a specific reference for sharp-order 1-D Schrodinger local smoothing) |
| Maximal Bochner-Riesz and Stein square function, n = 2 | Carbery 1983 for maximal BR at `p >= 4`. The full p >= 2 range follows from GWZ via the known reductions (BHS 2021). I believe the square-function range was also known (S. Lee and others) but **did not verify** this |
| Log power at `p = 4` for 2-D local smoothing | Not in GWZ (`N^eps`). The paper cites Cushman-Demeter-Wu 2026 (arXiv:2609.13105) for the sharp `log^(1/4)` in the square-function inequality for canonical sectors. Combined with a log-loss data-side square function bound (Cordoba-type), this may already give a log power. **Not verified** |

**Net for n = 2:** every eps-level consequence was already a theorem. At most, the plane gains:
* a `(log N)^(1/2)` bound at `p = 4`, with lower bound `(log N)^(1/4)`;
* a new, short-ish linear-PDE proof (Appendix "planar").

## 6. What is genuinely new (if correct)

* **n >= 3, everything critical.** Before this, local smoothing in all dimensions held only for
  `p >= 2(n+1)/(n-1)` (Bourgain-Demeter 2015). For n = 3: `p >= 10/3` (Gan-He-Li-Wu 2026 per the paper's
  bibliography), or a partial gain at p = 3 (Gao-Liu-Miao-Xi 2023). For n >= 4 there was also the eps-free large-p
  result (Heo-Nazarov-Seeger 2011). The paper claims Bochner-Riesz and restriction for n >= 3 (open since the 1970s),
  Kakeya sets for n >= 4 (Wang-Zahl 2025 settled n = 3), and Kakeya maximal for n >= 3 (open).
* **Quantitative.** The eps loss is replaced by an explicit `(log N)^(2/p_c)` at `p_c` (the paper shows a log loss is
  necessary there), and by log-Sobolev losses above `p_c`. The eps-free endpoint above `p_c` is still not claimed.
* **Method.** A positivity/energy argument (divergence-free positive tensor from complexified waves, plus a quaternionic
  Monge-Ampere potential on fibers) replaces induction on scales, decoupling and polynomial partitioning. If it holds
  up, this is a new tool in harmonic analysis.

## 7. Facts that bear directly on our problem (no speculation beyond what is stated)

1. **Every statement is Euclidean, translation-invariant, and holds for all `f` in L^p.** Applied to exponential sums,
   it holds for **any** finite set of frequencies on a real circle with any coefficients. So it is satisfied by our fake
   circles (real angles on a circle). Any argument whose only inputs are these inequalities is capped at
   `log R/loglog R` by the fake-circle barrier, in the same way as round 10 Thm B.
2. **The paper's torus is a fixed-size device** (`R^n/(100 Z^n)`, time `<= 2`, no wraparound). It carries no arithmetic.
   Long-time or periodic (Strichartz-on-torus) statements, where lattice structure enters, are **not** claimed.
3. **Our arcs live inside one canonical sector.** At frequency `N ~ R` the cone and circle decompositions use sectors of
   angular width `N^(-1/2)`, which is arc length `~ sqrt R`. An arc of length `C sqrt R` therefore meets `O(1 + C)`
   sectors, and its points lie in a `C sqrt R x C^2/8` rectangle. Square-function, decoupling and local-smoothing
   estimates have no content inside a single sector: the sectors are their atoms. The new theorem changes inter-sector
   bounds only.
4. **The losses are far too large.** The losses are `N^eps` or polylog in `N = R`. Our proved bound is already
   `O(log R/loglog R)`, and the target is `o(log R/loglog R)` or `O(1)`.
5. **In dimension 2, nothing at eps level is new.** GWZ has been available since 2019, and the record contains no route
   that was waiting on 2-D local smoothing.
6. At `p = 4` on a single circle, the `L^4` (additive energy) count of points on a circle is exact and elementary
   (`xi_1 + xi_2 = xi_3 + xi_4` on a circle forces trivial or antipodal solutions). This is pair-chord information, which
   is already capped.

Verdict from this lens: **not useful** for the uniform bound or for a growth improvement. The claimed result is a
major development in Euclidean harmonic analysis, mainly for n >= 3. Its 2-D content is pre-existing up to log factors,
its inequalities cannot tell genuine clusters from fakes, and it acts only between sectors, while our problem sits
inside one sector.

## Sources

* https://github.com/openai/math (commit 7ef5c0c, cloned 2026-10-27)
* https://github.com/openai/math/blob/main/history.md
* https://github.com/openai/math/blob/main/preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/paper.pdf
* https://github.com/openai/math/blob/main/lean/ComparatorChallenges/LocalSmoothing.lean
* https://github.com/openai/math/blob/main/preprints/Critical-local-smoothing-for-the-three-dimensional-wave-equation-September-24-2026/paper.pdf
* https://openai.com/index/building-scientific-understanding-with-academic-communities/ (snippet only)
* https://x.com/danintheory/status/2114833685417689321 (snippet only)
* https://www.scientificamerican.com/article/openai-unleashes-hundreds-more-math-results-upon-a-field-already-in-shock/ (snippet only)
* https://the-decoder.com/some-mathematicians-call-for-openai-boycott-after-ai-generated-proofs-flood-their-field/ (snippet only)
* https://tech-insider.org/openai-722-math-claims-mathematician-scrutiny-2026/ (snippet only)
* https://www.datastudios.org/post/openai-models-exponential-time-hypothesis-disproof-kakeya-yang-mills-proofs (snippet only)
* Guth-Wang-Zhang, Ann. of Math. 192 (2020) 551-581; Bourgain-Rudnick arXiv:1308.4247 (states the Cilleruelo-Granville
  conjecture and the `B_lambda << log lambda` bound for individual lambda)
