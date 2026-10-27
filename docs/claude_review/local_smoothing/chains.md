# From Euclidean harmonic analysis to lattice points on short arcs: what the claimed local smoothing proof adds

Lens: every known route from restriction, decoupling, Kakeya, square-function or local-smoothing estimates to
lattice points on circles and arcs. Written 2026-10-27, one day after the announcement.

Labels used below:
* **[verified]**: I read it in a primary source.
* **[snippet]**: seen only in search-engine summaries of the source.
* **[memory]**: my recollection, not re-checked this session.
* **[proved here]**: an elementary argument given in full in this note.
* **[speculation]**: a guess, not an argument.

## 0. Verdict

**Not useful for us.** By our strict criterion it gives no plausible route to `M = o(log R/loglog R)` or to the
uniform bound, and it unblocks none of route 5, P3-all\* or H^gen/L7. Four reasons, each argued below:

1. **Scale.** An arc of length `C sqrt(R)` is exactly O(C) *canonical caps* of the circle of radius `R` at lattice
   resolution (sagitta `C^2/8`). Restriction, decoupling, square-function, local-smoothing and Kakeya estimates
   treat such a cap as flat and cannot count points inside it (§2).
2. **Dimension and loss.** The 2D Euclidean theory is already complete:
   * Fourier restriction for curves in `R^2` is classical;
   * decoupling for the circle and cone is known (Bourgain-Demeter 2015);
   * local smoothing in 2+1 dimensions is known (Guth-Wang-Zhang 2020).

   The new content is `n >= 3`, with `eps`-losses (log-losses at the critical exponent). In 2D any `R^eps`-loss
   estimate for lattice points on a circle is vacuous, because `#E_R = R^(O(1/loglog R))` (§3). The paper itself
   shows that an `eps`-free version fails at the critical exponent (§1).
3. **Cap theorems.**
   * *(Prop. 1)* One explicit real configuration satisfies every Euclidean and exact-additive input. It has
     `floor(C sqrt R)` points, exactly on the circle and in the arc.
   * *(Prop. 2)* Any argument that bounds `M` through `L^p(T^2)` norms, `Lambda(p)` constants or k-fold additive
     energies can reach a contradiction only when `M >= c log R`. So it is capped at `Theta(log R)`, worse than the
     proven `(2/3 + eps) log R/loglog R`.
4. **Direction.** In every chain I found, arc counts are an *input* to harmonic analysis, never an output:
   * Bourgain-Rudnick (restriction of toral eigenfunctions);
   * decoupling-based discrete restriction for spheres;
   * the additive-energy papers;
   * Cilleruelo-Granville's link to squares in arithmetic progressions.

The fake-circle barrier is not even needed for most chains: they are already capped at the trivial bound
`C sqrt(R)` (Prop. 1) or at `log R` (Prop. 2). Every gap inequality valid for actual circles is in the fake class
by definition, including unit separation. So fakes satisfy every Euclidean input as well (§4, column c).

## 1. What was announced, and what I could verify

**Primary source, located and read [verified].** The GitHub repository `openai/math`:
* `CONTENTS.md`, entry 079, "Local smoothing for the Euclidean wave equation";
* `history.md`, entry dated October 26, 2026;
* the preprint `preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/paper.pdf` ("OpenAI",
  dated October 25, 2026, 54+ pages);
* `lean/docs/079.md`;
* the formal statement file `lean/ComparatorChallenges/LocalSmoothing.lean`.

Local copies are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/local_smoothing/dl/`
(`CONTENTS.md`, `ls_paper.pdf`, `ls_paper.txt`, `lean079.md`, `LocalSmoothing.lean`).

### Exact statements [verified, from the paper and the Lean statement file]

* **Theorem 1.1.** For `n >= 2`, `2 <= p < infinity` and `eps > 0`,

  ```text
  || e^{-iT|D|} f ||_{L^p(R^n x [1,2])} <= C_{n,p,eps} || (1-Delta)^{(s_+(p)+eps)/2} f ||_{L^p},
  s_+(p) = max(0, (n-1)/2 - n/p).
  ```

  Here `eps = 0` is allowed only at `p = 2`. So this is Sogge's conjecture **with arbitrarily small loss**: sharp up
  to `eps`, over the full range `2 <= p < infinity`. The Lean comparator statement `euclidean_local_smoothing` has
  exactly this form.
* **Theorem 1.2 (critical annular estimate).** At `p = p_c = 2n/(n-1)`, frequency-`N` data satisfy
  `||U(T) P_N f||_{L^p(R^n x I)} <= C (log(2+N))^{2/p} ||f||_p`.
* **Endpoint obstructions (Props. 8.1-8.2).** The optimal log power lies between `1/p_c` and `2/p_c`. Radial
  focusing forces at least `(log N)^{1/p_c}`. So an **`eps`-free (endpoint) estimate is false at `p_c`** for
  unrestricted data.
* **Stated consequences.** All in open ranges, with `eps`-losses where applicable:
  * Bochner-Riesz;
  * maximal Bochner-Riesz for `p >= 2`;
  * Stein's square function;
  * restriction for the sphere, compact paraboloid and truncated cone "throughout the interior of the conjectured
    range" (`b > 2n/(n-1)`, `b > (n+1)a'/(n-1)`; **not** the endpoints);
  * Schrodinger local smoothing;
  * Kakeya and Nikodym maximal bounds with loss `delta^(-eps)`;
  * full Hausdorff dimension of Kakeya sets.
* **Method.** Section 1.4 of the paper:
  * an imaginary translation along light rays gives each mode a Gaussian directional weight;
  * a positive semidefinite, divergence-free tensor is built from the continued wave in two auxiliary time
    variables;
  * on each spatial fibre one solves a **quaternionic Monge-Ampere equation** `det H_fib(lambda|Phi|^2 + u) =
    lambda mu`, with an Aubin-Mabuchi-type energy correction;
  * a scalar differential inequality in the auxiliary times, compared across circles of radius `2a` and 1, costs
    two logarithms.

  The core estimate is proved on a **fixed spatial torus** `R^n/(100 Z)^n` for times `T in [1,2]`, then
  transferred to `R^n` (§7). The torus is used for compactness, not arithmetic: this is the short-time regime,
  which is equivalent to the Euclidean one.

### Not verified

* **Correctness.** I did not check the proof. I did not read or build the Lean development. In the comparator file
  the statements end in `sorry`; the proofs live elsewhere in the repository.
* **Expert commentary.** Searches on 2026-10-27 found only press coverage and OpenAI's own post
  (`openai.com/index/building-scientific-understanding-with-academic-communities/`, seen as a [snippet]). The post
  and the press describe the material as raw model output, not peer reviewed, and handed to outside experts. I
  found no comment from Tao, Guth, Wang, Zahl or Demeter.
* **Access.** `openai.com`, arXiv, `annals.math.princeton.edu`, the news sites and the authors' home pages were
  unreachable from this sandbox (DNS or proxy refusal). Only GitHub was reachable.

### Context

* The same repository lists 374 numbered results claiming many famous problems. Examples: entry 003, a
  "quasi-Riemann hypothesis" with zero-free region `Re s > 7/8`; entry 074, Kakeya maximal in `R^3`; entry 196, a
  counterexample to Kaplansky's zero-divisor conjecture; entry 028, the Gaussian moat conjecture.
* This is a reason for caution about correctness. It does not matter here, because nothing below depends on whether
  the theorem is true.
* No entry of the repository concerns lattice points on circles, short arcs or sums of two squares (I searched
  `CONTENTS.md` for these terms).

### Previously known [verified from the paper's introduction; snippets agree]

* Guth-Wang-Zhang (Annals 2020): `n = 2`, via the sharp square function for the cone in `R^3`.
* Bourgain-Demeter cone decoupling: full gain for `p >= 2(n+1)/(n-1)`.
* Gan-He-Li-Wu: `n = 3`, `p >= 10/3`.
* Cushman-Demeter-Wu: sharp `log^(1/4)` power for the cone square function in `R^3`.
* Wang-Zahl (2025): Kakeya sets in `R^3`.

## 2. The scale fact that decides most chains

Take an arc of length `L = C sqrt(R)` on `|xi| = R`. Its sagitta is `R(1 - cos(L/2R)) ~ L^2/(8R) = C^2/8`. For
`T^2`-periodic functions, frequency resolution is 1, i.e. relative thickness `delta = 1/R` after scaling the circle
to radius 1.

The canonical caps of the `delta`-neighbourhood have relative length `delta^(1/2) = R^(-1/2)`, i.e. length
`sqrt(R)` before scaling. These caps are the atoms of every Euclidean estimate here:
* the caps of `l^2` decoupling;
* wave-packet directions;
* the sectors of the Guth-Wang-Zhang and Cordoba square functions;
* the plates or tubes dual to Kakeya.

**So an arc of length `C sqrt(R)` is `O(C)` canonical caps, each flat at the uncertainty scale.** An estimate stated
for functions Fourier-supported in the `O(1)`-neighbourhood of the circle cannot distinguish two cases:
* `M` lattice points on the arc;
* `M` separated points in a `C sqrt(R) x O(1 + C^2)` rectangle, where `M` can be `~ C sqrt(R)`.

This is why, for example in decoupling-based discrete restriction for spheres, the number of lattice points per cap
must be supplied as an *input* (§4, row B).

### Proposition 1 (Euclidean and exact-additive inputs are capped at the trivial bound) [proved here]

Let `omega = e^(i/R)` with `R >= 1` an integer, and let

```text
Z = { R omega^k : 0 <= k < C sqrt(R) } ⊂ C = R^2.
```

Then:
* `Z` has `floor(C sqrt R)` (or one more) points, exactly on the circle `|z| = R`;
* the points lie in an arc of length `<= C sqrt(R)`;
* consecutive points are at distance `2R sin(1/2R) in (1 - 1/(24R^2), 1)`.

By Lindemann-Weierstrass, `e^a` is transcendental for algebraic `a != 0`, so `omega = e^(i/R)` is transcendental.
Hence `sum_k c_k omega^k = 0` with algebraic `c_k` forces `c = 0`. So:

1. **Euclidean inputs hold for `Z`.** Every estimate valid for (essentially) 1-separated frequency sets in the
   `O(1)`-neighbourhood of a circle holds for `Z`. These estimates are perturbation-stable by the uncertainty
   principle. They include:
   * Fourier restriction and extension;
   * `l^2` and `l^p` decoupling, with canonical or small caps (Demeter-Guth-Wang), including the log-loss variants;
   * the cone square function;
   * local smoothing on the 3D frequency set `{(xi, |xi|)}`;
   * any Kakeya-derived bound;
   * all with any loss (`eps` or log).
2. **Exact-additive inputs hold for `Z`.**
   * Every exact additive relation among points of `Z` is trivial, so every k-fold energy `E_k(Z)` equals its
     diagonal count.
   * Kronecker's theorem then makes the Besicovitch (almost-periodic) `L^p` norm of `sum a_k e(z_k . x)` equal to
     `||sum a_k e(theta_k)||_{L^p(T^M)}`. That is the Steinhaus value, `<= C sqrt(p) ||a||_2`, the best possible
     `Lambda(p)` behaviour for every `p`.
   * Cushman-Demeter-Wu's `J_3 <= C_eps |X|^(3+eps)` for points on convex curves holds trivially for `Z`.

Therefore no conjunction of such inputs implies a bound better than `M <= C sqrt(R) + 1`. The integer lattice can
enter only in two ways:
* through periodicity on `T^2`, which is the box pigeonhole of Prop. 2;
* through Gaussian-integer arithmetic, i.e. the input classes our record already caps with fakes.

## 3. L^p, Lambda(p) and energy chains are capped at `Theta(log R)`

### Proposition 2 [proved here; the Khintchine lower bound is standard]

**Setting.** Let `S ⊂ Z^2` be `M` lattice points in an arc of length `L = C sqrt(R)` of `|xi| = R`. Let
`f_a = sum_{xi in S} a_xi e(xi . x)` on `T^2`. Let `K_p(S) = sup_a ||f_a||_p/||a||_2` be the `Lambda(p)` constant,
and `E_k(S)` the k-fold additive energy. Put `V = 36 C sqrt(R) max(1, C^2/2)`.

**(a) Box lower bound.** Take `xi_0` an endpoint of the arc. Then `S - xi_0` lies in a rectangle with sides `L`
(tangential) and `C^2/2` (normal). On the rectangle `|x_t| <= 1/(12L)`, `|x_n| <= 1/(12 max(1, C^2/2))`, every phase
satisfies `|<xi - xi_0, x>| <= 1/6`, so `|f_1(x)| >= M/2`. That rectangle has measure `1/V`. Hence for every
`p >= 2`:

```text
||f_1||_p >= (M/2) V^(-1/p),   so   K_p(S) >= (sqrt(M)/2) V^(-1/p).
```

**(b) Energy lower bound.** By Cauchy-Schwarz over the k-fold sumset,
`E_k(S) >= M^(2k)/#(kS) >= M^(2k)/((kL+1)(k C^2/2+1))`.

**(c) What any `Lambda(p)` hypothesis yields.** Suppose some hypothesis gives `K_p(S) <= kappa`. Then (a) gives
`M <= 4 kappa^2 V^(2/p)`. With the strongest conceivable form, `kappa = A sqrt(p)` for all
`2 <= p <= 2 ln V`, the choice `p = 2 ln V` gives

```text
M <= 8 e A^2 ln V = (4 e A^2 + o(1)) ln R.
```

The same computation with (b) and the best possible (`B_k`) energy `E_k <= k! M^k` gives
`M <= (k!)^(1/k) (k^2 V)^(1/k)`. That is `M <= (1/2 + o(1)) ln R`, at `k ~ ln V ~ (1/2) ln R`. This is the order of
the Cilleruelo-Cordoba bound `k_R = ceil(log R/2)`.

**(d) The chain cannot do better.** To rule out `M`-point arc sets, the chain needs a theorem
`K_p(S) <= kappa < (sqrt(M)/2) V^(-1/p)` at some `p`. Two cases:
* **`2 <= p <= c_1 M`.** Flat coefficients with random signs give `K_p(S) >= c_0 sqrt(p)`, by the central-limit /
  large-deviation lower bound for Rademacher sums. So such a `kappa` can exist only if
  `M > 4 c_0^2 min_{2 <= p <= c_1 M} p V^(2/p)`, which forces `M >= c_2 log R` for an absolute `c_2 > 0`.
* **`p > c_1 M`** (interpretation, not a theorem). No harmonic-analysis lower bound obstructs `kappa`. But the required inequality says that the
  `L^p` mass of the dual-box peak is smaller than it must be, which is equivalent to the absence of `M` points in
  the arc. This is the conclusion restated, with no leverage.

**Conclusion.** **No upper bound on `L^p(T^2)` norms, `Lambda(p)` constants or energies of arc sets, however sharp,
can prove `M = o(log R)` through this chain.** In particular none can improve our
`(2/3 + eps) log R/loglog R`.

Two further remarks:
* **Large `k`.** To beat `log R/loglog R` by energies one would need near-diagonal control of `E_k` at
  `k >= (log R)^(1+delta)`. At that level, near-diagonality is equivalent to the pigeonhole statement itself, i.e.
  to the arc bound. It is arithmetic, not harmonic analysis.
* **Vacuity of `R^eps` losses (2D).** `||sum_{E_R} a e||_p <= |E_R|^(1/2) ||a||_2` and
  `|E_R| = R^(O(1/loglog R))`. So every `R^eps`-loss estimate for lattice points on one circle in `Z^2` is implied
  by the divisor bound. Decoupling, discrete restriction and every `eps`-loss local smoothing statement therefore
  carry no information on our problem.

## 4. Chain-by-chain survey

Columns:
* **(a)** Does the chain need the new content (`n >= 3`, or an `eps`-free estimate)?
* **(b)** The best bound it could give for arcs of length `C sqrt(R)`.
* **(c)** Are its inputs satisfied by fake circles (hence capped)?

| # | Chain (literature) | Direction in the literature | (a) new content? | (b) best for `C sqrt R` arcs | (c) fakes / cap |
|---|---|---|---|---|---|
| A | 2D restriction, Stein-Tomas, Zygmund `L^4` for toral eigenfunctions; Bourgain's eigenfunction `L^p` bounds; Bourgain-Rudnick restriction to curves and surfaces (AHP 2011; GAFA 22 (2012) 878-937, arXiv 1105.0018) [snippet, memory] | arc and cap counts are **inputs** (Jarnik, Cilleruelo-Cordoba, cluster structure). The output is restriction or nodal-set bounds. | No. 2D restriction is classical; Bourgain-Rudnick work in `T^2`, `T^3`. | `L^4` (Zygmund constant `2^(1/4)`): `M <~ R^(1/4)`. Any `L^p`: `<= Theta(log R)` (Prop. 2). | Zygmund `L^4` holds for any real points exactly on a circle (Besicovitch norm), including `Z` and fakes. Capped, trivially. |
| B | Decoupling for circle, parabola and cone (Bourgain-Demeter, Annals 182 (2015) 351-389 [snippet]); small-cap decoupling (Demeter-Guth-Wang 2020 [memory]); log-loss decoupling for the parabola (Guth-Maldague-Wang [memory]); discrete restriction for spheres (Bourgain-Demeter; ranges `p >= 2n/(n-3)`, then improved [snippet]) | caps of length `~ sqrt(R)` are the decoupling atoms. Per-cap lattice counts are **inputs** to discrete restriction [memory]. | No. Needed decoupling is known in every dimension; the new theorem adds nothing to decoupling. | Output is `L^p` bounds, so `<= Theta(log R)` (Prop. 2). With `R^eps` loss: vacuous (divisor bound). | Perturbation-stable, so `Z` and fakes satisfy it. Capped at trivial. |
| C | Square function and local smoothing: Guth-Wang-Zhang 2D (Annals 192 (2020)) [verified via paper and search]; claimed `n >= 3` (openai/math #079) [verified statement]; Cushman-Demeter-Wu `log^(1/4)` [verified citation in paper] | No chain to lattice points on circles found. Data on one circle `|xi| = R`: `e^(-iT|D|)` acts by the phase `e^(-iTR)`, so (1.1) is an identity. The paper's torus is fixed-size and short-time, i.e. Euclidean. | `n >= 3` is irrelevant to `Z^2`. `eps`-free is not offered and is false at `p_c` (Props. 8.1-8.2). | Nothing. Single circle gives an identity; for several circles see Prop. 2 / row B. | Satisfied by every configuration. Capped at trivial. |
| D | Kakeya: Wang-Zahl `R^3` (2025) [snippet]; claimed maximal Kakeya in all `n` [verified statement]; arithmetic Kakeya (Bourgain; Katz-Tao) [memory] | No chain to lattice points on arcs found. | Not applicable. | Nothing. A cluster has `O(M^2) = O(log^2 R)` chord directions in a sector of angle `~ C R^(-1/2)`. Kakeya concerns `delta^(-(n-1))` directions. | Chord-direction information is the "pair chords" class, capped in the record. |
| E | Additive energy on circles: Bombieri-Bourgain, "A problem on sums of two squares", IMRN 2015, 3343-3407 [snippet] (sixth correlation of `x^2+y^2 = m`, motivated by Krishnapur-Kurlberg-Wigman nodal length; `E_3 <~ \|E\|^(7/2)` also in arXiv 1401.1873 [snippet]); Bourgain-Demeter via decoupling [snippet]; Mudgal (JLMS 2022) [snippet]; **Cushman-Demeter-Wu, arXiv 2608.12316 (Aug 2026)** [snippet] | Energy bounds are outputs; arcs do not appear. CDW: `J_3(gamma(X)) <=_eps \|X\|^(3+eps)` for **every** finite `X` and strictly convex `gamma`, settling a conjecture of KKW13/BB15 on lattice points on dilates of the circle. | No. CDW is already 2D and does not need the new theorem. | With pigeonhole (b): `M^(3-eps) <~ C^3 sqrt R`, so `M <~ R^(1/6+eps)`. Energies in general: `<= Theta(log R)` (Prop. 2). | CDW holds for all real points on a circle, so for fakes and `Z`. Upper energy bounds valid for actual circles also hold for fakes with algebraically independent prime angles: those satisfy only the formal Laurent-polynomial relations, a subset of the actual ones. (Genericity of the record's fakes after perturbation is not checked.) Capped. |
| F | Rudin `Lambda(p)` for squares (open for `2 < p < 4` [snippet: Chang]); squares in AP: Rudin's `sigma(k) << sqrt k`; Szemeredi `o(k)`; Bombieri-Granville-Pintz `k^(2/3+o(1))` (Duke 1992); Bombieri-Zannier `k^(3/5+o(1))` (2002); Cilleruelo-Granville 2007 (CRM Proc. 43, arXiv math/0608109) [content via our paper's §related, snippet]; Gabdullin (arXiv 2205.13611: `L^4` for squares in short ranges, Cilleruelo-Cordoba conjecture, Ruzsa's divisor conjecture) [snippet] | CG: arc / short-interval boundedness conjectures **imply** Rudin-type conjectures. The direction is arcs ⇒ squares. BGP and BZ use rational points on curves (algebraic geometry), not Fourier analysis. | No. No source links local smoothing or restriction to `Lambda(p)` for squares; the 1D `Lambda(p)` problem is not a consequence of any known Euclidean estimate [speculation]. | Nothing. The reverse implication (squares ⇒ arcs) is not known. | Curve-counting methods are in our record's Vojta / uniformity family, the only conditional route that works. Harmonic analysis plays no part. |
| G | Almost-all radii: Bourgain-Rudnick separation `R^(1-eps)` outside a sparse set (AHP 2011); Kurlberg-Lester Poisson spacing on a density-one subsequence (Invent. 2025) [via our paper] | Arithmetic and statistical, not Euclidean. | No. | Nothing uniform: density-one statements miss worst-case `R` (compare our round 10 Theorem B). | Statistical / robust inputs are capped by round 10 Thm B. |
| H | Cilleruelo-Cordoba, PAMS 115 (1992) 899-905: arcs of length `sqrt2 R^(1/2 - 1/(4[k/2]+2))` contain `<= k` points [via our paper] | The integrality gap of Gaussian-integer half-angles (pair chords) gives the arc bound. I believe "Trigonometric polynomials" in the title refers to applications from arcs to trigonometric polynomials, but have not checked it [memory, unverified]. | No. | `<= (1/2) log R + 1` per sub-arc; the same order as Prop. 2(c), which is no coincidence. | Pair chords: capped by fakes (record). |
| I | Gauss circle problem through decoupling / Bombieri-Iwaniec (Bourgain-Watt) [memory] | Counts in disks, not on circles. | No. | `r(n) <= P(n^+) - P(n^-) <~ n^theta`: useless. | N/A. |

**A reverse benefit (for others, not for us).** Rows A and B take arc and cap counts as input. So our proven
`(2/3 + eps) log R/loglog R` could slightly sharpen toral-eigenfunction restriction or nodal estimates whose loss is
controlled by arc clustering. The commented-out passage near the end of §related in `paper/jarnik_sublog.tex`
already notes this direction.

## 5. Could anything in the new work touch our open routes? [speculation, clearly marked]

* **Route 5 / MIN^pair.** MIN^pair is a *discrete uniform-approximation* (Chebyshev or Remez type) statement.
  `kappa_H(u)` is the best uniform approximation error of `<k',u>^m` by lower-degree polynomials on the slab-cut
  ternary box (round 12 §2.3, "Duality"). The relevant tools would be Remez inequalities, sublevel-set and
  Lojasiewicz bounds, and polynomial partitioning. The paper uses a conserved current for a continuous wave and a
  fibrewise quaternionic Monge-Ampere equation. I see no transfer. Its "determinant" is a Hessian determinant,
  unrelated to the Vandermonde or Ramana determinants of route 5.
* **P3-all\*.** A residue-design and rigidity question, combinatorial and arithmetic. No harmonic-analysis entry
  point.
* **H^gen / L7** (Poisson tails of collision statistics of Gaussian-prime residues). This is analytic number theory
  over `Q(i)`; round 10 Theorem B caps robust ensemble inputs.
  * The same repository's entry 003 claims a zero-free half-plane `Re s > 7/8` for Dirichlet `L`-functions and for
    finite-order Hecke `L`-functions over `Q(sqrt -3)`.
  * Even if correct, it is a robust input with error `X^(7/8)`, covered by Theorem B.
  * It does not concern the infinite-order Grossencharaktere `(alpha/|alpha|)^(4k)` of `Q(i)` that control angles.
  * Outside this lens; noted only so it is not chased separately.
* **Long-time torus dynamics.** Schrodinger or wave on `T^2` at times `>> 1` does see arithmetic (Bourgain's `T^2`
  Strichartz uses divisor bounds). But the arithmetic is the input, losses are `eps`-type, and Prop. 2 caps any
  resulting `L^p` route at `log R`. Local smoothing is a short-time statement and says nothing here.

## 6. Sources

**Read directly (GitHub, reachable):**
* `https://github.com/openai/math`: `CONTENTS.md` (entries 003, 028, 074, 077, 079), `history.md`,
  `preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/paper.pdf` (abstract, §§1-2, §8),
  `lean/docs/079.md`, `lean/ComparatorChallenges/LocalSmoothing.lean`.
* Our own `paper/jarnik_sublog.tex`, §related, for the Cilleruelo-Cordoba, Cilleruelo-Granville,
  Bombieri-Granville-Pintz, Bombieri-Zannier, Bourgain-Rudnick and Kurlberg-Lester statements and references.

**Search snippets only** (hosts unreachable for fetching):
* OpenAI post `openai.com/index/building-scientific-understanding-with-academic-communities/`.
* Press: cryptobriefing.com, unite.ai, tech-insider.org, qz.com, mlq.ai, kucoin.com.
* arXiv 2608.12316 (Cushman-Demeter-Wu).
* IMRN 2015 (Bombieri-Bourgain), with the abstract via OUP.
* arXiv 1403.5335 / Annals 182 (Bourgain-Demeter).
* arXiv 1401.1873.
* arXiv 1105.0018 (Bourgain-Rudnick).
* arXiv math/0608109 (Cilleruelo-Granville).
* arXiv 2205.13611 (Gabdullin).
* M.-C. Chang, "On problems of Erdos and Rudin" (`Lambda(p)` for squares open).
* Guth-Wang-Zhang (Annals 2020, IAS and Annals listings).

**Memory, flagged in the text:**
* the exact cap-count input in Bourgain-Demeter's discrete spherical restriction;
* Demeter-Guth-Wang small caps;
* Guth-Maldague-Wang log-loss decoupling;
* Bourgain-Watt;
* arithmetic Kakeya.
