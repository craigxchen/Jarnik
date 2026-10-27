# Local smoothing / Kakeya / restriction in all dimensions: high-dimensional lifts of our open statements

Lens: can the claimed theorem (or its eps-free cores, or its auxiliary-PDE method) supply the missing step of
route 5 (MIN^pair), of P3-all\* rigidity, or of H^gen / L7, possibly after lifting the planar arc problem to
higher dimension?

Labels: **[V]** verified from the primary source by me; **[P]** proved here (elementary, all quantifiers);
**[L]** published literature I located in this session; **[M]** literature cited from memory, not re-checked
this session; **[S]** speculation / heuristic.

## 0. Verdict

**Not useful by our strict criterion.** No implication from the new results to MIN^pair, to P3-all\*, to
H^gen/L7, or to a growth improvement is visible, and four structural mismatches (§2) make one very unlikely.
These mismatches hold even if the whole local smoothing / Bochner-Riesz / restriction / Kakeya family were known
at the endpoint with no loss. They are:

1. Single-shell vacuity **[P]**. Our lattice points lie on one circle, i.e. in one Laplace eigenspace. On any
   function whose frequencies lie on one sphere, every local smoothing and Bochner-Riesz inequality is a scalar
   identity and says nothing.
2. Resolution **[P]**. A fixed-`p` L^p inequality can see an arc cluster only if the cluster has size at least
   `R^(c/p)`, which is a power of `R`. To see clusters of size `O(log R)` one needs `p ~ log R/loglog R`, with
   arithmetic Lambda(p) constants on a single circle. That is an open planar arithmetic problem, and the new
   theorem does not touch it.
3. Wrong direction and wrong kind of object. Route 5 needs a pointwise lower bound: no near-zeros of explicit
   polynomial forms on a region, at exponential precision. The new theorems give average (L^p) upper bounds,
   and these are blind to measure-zero zero sets.
4. Wrong regime and wrong geometry. Our lifts live in dimension `d ~ n -> infinity`, at bounded radius, with
   flat (ternary-cube) frequency sets. Restriction theory is about fixed dimension, scale tending to infinity
   and curvature. Its constants `C_n` (or `C_(n,p,eps)`) carry no stated dependence on `n`.

In addition, every conclusion of measure-theoretic type ("sublevel sets are small", "most configurations",
Hausdorff dimension) is already capped by `round10/analytic.md` Cor. B.1, which covers Haar-generic and robust
ensemble inputs. The auxiliary-PDE method (a quaternionic Monge-Ampère potential with an Aubin-Mabuchi energy)
is a compensated-integrability device for L^p upper bounds. It shares only the word "determinant" with route
5's interpolation determinant.

## 1. What I could verify about the news

**Primary source [V].** I cloned the public repository `github.com/openai/math` (commit `7ef5c0cc`,
2026-10-26 14:18 -0700). WebFetch failed with DNS errors on openai.com, arxiv.org and every news host, but git
over the proxy worked. The repository contains:

* the preprint `preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/` (author "OpenAI",
  dated 2026-10-25), with `paper.pdf` and its LaTeX sources;
* the Lean challenge file `lean/ComparatorChallenges/LocalSmoothing.lean` and its `.json`;
* the solution tree `lean/OAI/Analysis/Fourier/LocalSmoothing/` (74 files, about 5.5k lines, no `sorry` in that
  directory).

I did **not** build the Lean project, check its transitive imports, or check the mathematics.

**Exact statements (preprint, Theorem 1.1 and Theorem 1.2 / Prop. "periodic cap", Cor. "Spherical extension
and Kakeya") [V].** Notation: `U(T) = e^(-iT|D|)`, `J = (1-Delta)^(1/2)`,
`s_+(p) = max{0, (n-1)/2 - n/p}`, `p_c = 2n/(n-1)`.

* Local smoothing: for `n >= 2`, `2 <= p < infinity` and `eps > 0`,
  `||U(T) f||_(L^p(R^n x [1,2])) <= C_(n,p,eps) ||J^(s_+(p)+eps) f||_(L^p(R^n))`, with `eps = 0` allowed at
  `p = 2`. This is **sharp up to eps** (focusing example), not eps-free.
* The core estimate has a logarithmic loss at the critical exponent: for `p = p_c`,
  `||U(T) P_N f||_(L^p(R^n x [1,2])) <= C_(n,phi) [log(2+N)]^(2/p) ||f||_p`. The optimal log power lies in
  `[1/p_c, 2/p_c]`.
* The proof runs through a periodic version on `T^n = R^n/(100Z)^n`, for frequencies in a fixed-aperture cap of
  the annulus `|xi| ~ 1/h`: `(int_1^2 ||U(T) f||_(L^(p_c)(T^n))^(p_c) dT)^(1/p_c) <= C [log(2/h)]^(2/p_c) ||f||_(p_c)`.
* Consequences:
  * Bochner-Riesz in the full open range;
  * Stein's square function for `2 <= p < infinity`, `alpha > 1/2 + s_+(p)`;
  * global Fourier extension from `S^(n-1)`, `L^a -> L^b`, in the open range `b > 2n/(n-1)`, `b > (n+1)a'/(n-1)`;
  * local extension with a log loss, `||Eg||_(L^(p_c)(B(z,R))) <= C_n [log(2+R)]^(2/p_c) ||g||_(L^(p_c)(S^(n-1)))`;
  * the Kakeya maximal bound `||K_delta f||_(L^n(S^(n-1))) <= C_n [log(2/delta)]^(4/p_c) ||f||_(L^n)`, hence
    `delta^(-eps)` and full Hausdorff dimension of Kakeya sets;
  * Nikodym, paraboloid extension and Schrödinger local smoothing.
* The Lean comparator certifies six **statements**:
  * local smoothing with `eps > 0`;
  * Bochner-Riesz (strict range, including `p = infinity`);
  * spherical restriction (strict mixed range);
  * Kakeya maximal with `delta^(-eps)`;
  * Kakeya dimension `= n`;
  * Stein square function (strict range).

  The log-loss forms are **not** among the comparator targets.
* Method (paper §1.4):
  * a conserved, positive semidefinite, divergence-free tensor built from the wave continued to imaginary
    light-ray displacements, with two auxiliary time variables `tau in R^2`;
  * on each spatial fiber, a quaternionic Hessian equation `det H_fib(lambda|Phi|^2 + u) = lambda mu`;
  * an Aubin-Mabuchi-type corrected energy `F(T, tau)`;
  * the scalar inequality `Delta_(T,tau) F <= C s^(-2) e^(-cs/h) F`;
  * circle-mean comparisons in the `tau`-plane, each of which loses one logarithm.

  It cites Serre's compensated integrability and Kaiser's positive Klein-Gordon current as antecedents.
* Prior work cited in the preprint's history section: Sogge 1991; Mockenhaupt-Seeger-Sogge; Wolff 2000;
  Łaba-Wolff; Bourgain-Demeter 2015 (cone decoupling, `p >= 2(n+1)/(n-1)`); Heo-Nazarov-Seeger 2011;
  Guth-Wang-Zhang 2020; Gao-Liu-Miao-Xi 2023; Gan-He-Li-Wu 2026 (`p >= 10/3` in `R^3`); Wang-Zahl 2025
  (Kakeya in `R^3`). The same repository also holds earlier OpenAI preprints, dated 2026-09-23/24:
  * "The Kakeya maximal conjecture in three dimensions";
  * "Every four-dimensional Kakeya set has full Hausdorff dimension";
  * "Critical local smoothing for the three-dimensional wave equation";
  * "Elliptic capacity propagation and Fourier restriction to the sphere".

**Secondary sources [L] (search snippets only; pages not fetchable).**
* Announcement date: 2026-10-26, on OpenAI's page "Working with academic communities to build scientific
  understanding".
* Timeline: AI output on Oct 18, Lean formalization on Oct 19.
* Expert reaction:
  * a reported "historic development" remark attributed to Alex Iosevich (Rochester);
  * a statement by an Advisory Group on Mathematics and AI, reportedly dated Oct 25, criticising the release
    process.
* Not found: no comment from Guth, Wang, Zahl, Demeter or Tao, and no independent verification.

**Status for our purposes: unrefereed and unverified. None of the conclusions below depends on whether it is
correct.** They would hold even if the endpoint forms were true.

## 2. Four structural facts that decide the lens

### F1. Single-shell vacuity [P]

Let `f` on `R^n` or on a torus `R^n/(LZ)^n` have Fourier support on one sphere `{|xi| = rho}`. This is the case
for every function built from the lattice points of `x^2 + y^2 = N`: they form one Laplace eigenspace of `T^2`.
Then:

* `U(T) f = e^(-iT rho) f`, so `int_1^2 ||U(T)f||_p^p dT = ||f||_p^p`;
* `J^s f = (1 + rho^2)^(s/2) f`.

Every local smoothing inequality therefore reads `||f||_p <= K (1+rho^2)^(s/2) ||f||_p`, which holds for any
`K >= 1`. Likewise `B_r^alpha f = (1 - rho^2/r^2)_+^alpha f` and the Stein square function act as scalars.

The same holds approximately for frequencies in a shell `rho <= |xi| <= rho + eta` with `eta = o(1)`: the time
window `[1,2]` cannot separate shells closer than `O(1)`.

The new theorem's periodic-cap estimate (Prop. "periodic cap", on `T^n`, fixed aperture) is therefore true but
empty on every object our problem produces. Local smoothing is a statement about dispersion **across many
shells**. A single circle has none.

### F2. Resolution of L^p [P]

Take `M` lattice points on `|z| = R` (`R^2 = N`) inside an arc of length `C sqrt R`, and put
`F(x) = sum_x e(<z_x, x>)` on `T^2`.

* Chord differences are at most `C sqrt R` tangentially and at most `C^2/2` normally (the sagitta).
* So `|F| >= M/2` on a rectangle of area `>= c C^(-3) R^(-1/2)`.
* Hence `||F||_p^p >= c (M/2)^p C^(-3) R^(-1/2)`.
* `||F||_2 = M^(1/2)`.

Any bound `||F||_p <= K_p ||F||_2` therefore yields only `M <= (c' K_p^p C^3 R^(1/2))^(2/p)`, i.e.
`M <~ R^(1/p)`. A cluster becomes visible to L^p only when `M^(p/2) > R^(1/2)`, i.e. `p > log R/log M`. At
the scale that matters, `M ~ log R/loglog R` (so `log M ~ loglog R`), this means
`p > (1 + o(1)) log R/loglog R`, with constants `K_p` uniform in `R`. For bounded `M` it means
`p >> log R`.

That is an arithmetic Lambda(p) problem for one circle with `p` growing, Bourgain's eigenfunction problem on
`T^2` in an extreme range. The implication is known to run the other way: arc-cluster bounds feed
eigenfunction and nodal estimates **[L]**. Bourgain-Rudnick state their results with `B_lambda`, the maximal
number of lattice points on arcs of length `sqrt lambda`. Cilleruelo-Granville predicts `B_lambda = O(1)`;
Bourgain-Rudnick prove it for most `lambda` and `B_lambda << log lambda` in general.

By F1, the new theorem says nothing about such single-circle norms. **Even an eps-free, log-free form
would not change this**, because the obstruction is the volume factor `R^(-1/2)`, not the loss.

### F3. Direction and kind of conclusion [P]

* **The conclusions are upper bounds on averages.** The new theorems bound L^p norms (or Hausdorff
  dimension). An L^p bound is unchanged when the function is modified on a null set. So it can never
  exclude a zero, or a near-zero at precision `e^(-W/4)`, at a particular point, or on a lower-dimensional
  real-algebraic set. The real tangent-cone points of route 5 are such sets.
* **Measure-type consequences are capped.** What an L^p upper bound can give is a measure-type
  statement: sublevel sets are small, most configurations behave. `round10/analytic.md` Cor. B.1 caps exactly
  these inputs. Haar-generic properties and robust ensemble statements are satisfied by fakes with
  `W = O(M log M)`.

### F4. Regime and geometry [P/S]

**The dimension grows.** In every lift of our problem the dimension is a growing parameter:
* `d = n - 1` for the slab box;
* `r = Theta(M)` for the exponent lattice;
* `2M` or `M^2` for tensor lifts.

The new theorems are stated for fixed `n`, with `C_n` or `C_(n,p,eps)` of unstated (in the proofs, at best
super-exponential **[S]**) growth in `n`. Route 5 works at precision `e^(-W/4)`, with `W ~ n log n` in the
regime. A dimension-free or explicitly `e^(o(n log n))` constant would be needed, and nothing of the kind is
claimed.

**No scale separation.** Restriction, decoupling and Kakeya gain from curvature at scale `R -> infinity`
in fixed dimension. Our frequency sets are:
* finite subsets of the ternary cube `{-1,0,1}^d`;
* cut by linear slabs `|<k, beta_p>| <= w`;
* of radius `<= sqrt d`, at mutual separation `>= 1`.

So there is no scale separation and no curvature. For flat frequency sets the sharp L^p bounds factor
(tensor products), and harmonic analysis has nothing beyond Plancherel to add.

## 3. The three open statements, one by one

### 3.1 Route 5: MIN^pair (round 12 §2.3 item 6; round 13 §0, §7)

**Needed (record, refereed formulation).** Take the cluster's own `T = 1` slab box
`A = {k in Z^n_0 : |k_x| <= 1, |<k, beta_p>| <= w}` with `d = n - 1`, `w = sqrt(8 d ln 4r)`, and order
`m' = ceil(4w)`. The statement is

```text
min { kappa_H(u) : |u|_inf = 1, |u_x - u_y| >= e^(-eps_xy/2) for all pairs }  >  eps_0 = C d e^(-W/4)/(m'+1),
kappa_c(u) = |Phi_c(u)|/(|c|_1 d^(m')),   Phi_c(u) = sum_(k in A) c_k <k',u>^(m'),   c in V_H.
```

Here `V_H` consists of the integer measures on `A` with height `<= H` that annihilate all polynomials of
degree `< m'`.

Equivalently (nearstrata duality): for every `u` in the gap region, the function `k -> <k',u>^(m')` on the
finite set `A` lies at `l^inf(A)`-distance `> d^(m') eps_0` from the polynomials of degree `< m'`. The
trigonometric form `sum_k c_k e^(i<k, delta>)` enters only through its order-`m'` Taylor coefficient at
`delta = 0`.

**Given (new theorem).** L^p upper bounds for `e^(-iT|D|) f`, for the extension operator from `S^(n-1)` and
for Kakeya maximal functions. They hold in fixed dimension, with curved frequency support, and with losses
`eps`, `delta^(-eps)` or `log^(O(1))`.

**Match check.**

| aspect | needed | given | match |
| --- | --- | --- | --- |
| direction | lower bound (no near-zero) | upper bound | **no** |
| locality | pointwise, on a region of directions, near `delta = 0` | spacetime or spatial L^p averages | **no** (F3; and by the uncertainty principle, a function with frequencies in a set of radius `<= sqrt d` is polynomial-like below scale `d^(-1/2)`, while MIN^pair is an order-`m'` Taylor statement at `0`) |
| objects | polynomial forms / exponential polynomials on `(C^*)^d` with integer coefficients, frequencies in a flat finite set | Schwartz functions, frequencies on cone / sphere | **no** (F4) |
| dimension | `d -> infinity`, precision `e^(-W/4)`, `W ~ n log n` | fixed `n`, unspecified `C_n` | **no** |
| losses | `e^(o(n log n))` at most | `eps`, `log` | moot, given the rows above |

**What could conceivably bridge it, and why it does not.**
1. **Lower-bound tools for exponential sums [M].** Turán's lemma, Nazarov's form of it, and multivariate
   Remez inequalities bound `sup_I |f| <= (C|I|/|E|)^(K-1) sup_E |f|` for `K`-term sums. Here
   `K = |A| ~ 3^d`, so the loss is `exp(3^d)`. Remez / Bernstein-Walsh bounds have exponent of order the
   degree times a capacity ratio, and they bound sublevel-set **measure**, not emptiness.
2. **These are classical, not new.** Nearstrata Prop. 2.7 already shows that no bound depending only on
   dimension, degree, height and the qualitative zero set reaches the needed exponent `d/2`. A sublevel-set
   measure bound is also capped by F3: fakes sit exactly on the small set, as fake Thm 2.1 shows.
3. **The auxiliary PDE does not transfer.** Its mechanism is Alexandrov-type: mass concentration of a
   density forces a Monge-Ampère potential to be deep. That is anti-concentration, which yields upper bounds
   for L^p norms of the density. Pluripotential theory, the classical home of Monge-Ampère ↔ polynomial
   inequalities [M], is old and gives the measure-type statements of item 1.
4. **No algebraic content.** Nothing in the method produces auxiliary polynomials, interpolation determinants
   or information about the leading-form ideal `J(A)`. Round 13 identified `J(A)` as the only structural
   handle on the regime tangent cone; Loc(b0) is open there and refuted for column projections.

**Verdict: no.** The method does not decide the fake test (MIN^pair on (Mar) fake pattern boxes) either: a
universal real-analytic inequality applies to fake boxes and actual boxes alike.

### 3.2 P3-all\* rigidity (round 14; round 15 Rig^SL)

**Needed.** A combinatorial-arithmetic statement about label systems: residues `ell_j(q) in Z/m_q` at moduli
up to `3M^A`, parity constraints, monomial budgets `m_(v,s)^2 <= Norm Z_v 2^[s odd]`, and one angle vector
`phi in R^r` (slab system 1.2). The goal is either a forced lower bound `W >> M log M` (or rigidity plus
H^gen + L1), or a fake (C1) in regime (1.4).

**Given.** Euclidean tube-overlap and oscillatory-integral bounds.

**Match.** There is no Euclidean or continuous structure to apply them to. The only continuous datum is the
single angle vector `phi`. The constraints on it are linear slab inequalities (Lemma J0 box test), which are
flat, finite and per-monomial.

The finite-ring analogues of Kakeya are older and independent of the Euclidean theorem [M]: Dvir 2008 over
`F_q`, and Dhar-Dvir style results over `Z/N`. They concern sets containing lines, not label systems. The
known bridge between arithmetic and Euclidean Kakeya runs arithmetic → Euclidean (Bourgain 1999;
Katz-Tao, sum-difference estimates) [M]. As far as I know, a proof of Euclidean Kakeya does not yield the
arithmetic (sum-difference) statements.

**Verdict: no.**

### 3.3 H^gen / L7 (round 9; round 10 Thm B)

**Needed.** Poisson-type tails for collision statistics of the genuine Gaussian-prime residues of the cluster
primes, at moduli up to `M^A` with `A` chosen after `c_0`. These are pointwise or selection-sensitive
statements about `O(M)` primes among roughly `M^A/log M`.

**Given.** Nothing about primes. Harmonic-analysis inputs that do concern primes are robust ensemble
statements, with error `>= c X^gamma` or large-sieve shape: Green's restriction theorem for primes, Bourgain's
Lambda(p) for primes, Hecke large sieve [M]. Theorem B (graft barrier) caps every such input, because grafting
`O(M)` angles or residues is invisible to them.

**Verdict: no.**

## 4. High-dimensional lifts of the planar arc problem

1. **Tensor / product lifts** (`(z_i, z_j) in Z^4` on `|w|^2 = 2N`) **[P]**.
   * The `M^2` product points lie in a cap of radius `~ C N^(1/4)` on `S^3(sqrt(2N))`.
   * That sphere carries `~ N` lattice points (`r_4(m) >> m`).
   * The expected number in such a cap is `~ (C N^(1/4))^3/(2N)^(1/2) ~ C^3 N^(1/4) >> M^2`. The lifted
     cluster is **not rare**.
   * It is also a product configuration on a 2-torus. Configurations concentrated on lower-dimensional sets
     are exactly the extremisers (Knapp-type) for which restriction and Kakeya bounds are sharp, so the lift
     gains nothing.
   * Higher tensor powers behave the same way. Multiplying cluster points in `Z[i]` doubles the angle and
     moves the arc length to `R'^(3/4)`.
2. **Exponent-lattice lift** **[P]**. `arg z_x = <2a_x - e, theta> (mod 2pi)`. A cluster is `n` points of the
   box `prod [0, e_j] subset Z^r` inside a thin slab `{v : ||<2v, theta> - c|| <= Delta}`, with
   `Delta = C N^(-1/4)`.
   * The slab is flat; its normal is the Gaussian-prime angle vector.
   * The content is Diophantine (small values of `<v, theta>`): Baker / Lang-Waldschmidt / GAP-Fourier, all
     capped in the record.
   * No curvature appears, so restriction and decoupling are inapplicable.
3. **Spheres in dimension `>= 3`** **[L/M]**.
   * Lattice points in caps on spheres (Bourgain-Rudnick, cap bound via Jarník's theorem) are bounded by
     reducing to planar sections. A plane meets the sphere in a conic of a 2-dimensional sublattice, so the
     reduction runs **planar → spherical**, and the planar binary-form problem is a sub-case.
   * Discrete restriction to lattice points on spheres needs decoupling (Bourgain-Demeter 2015, `N^eps`
     loss) [M]. Euclidean restriction does not imply it. In dimension 2 everything relevant was already
     classical: Zygmund's L^4, and Bombieri-Bourgain on L^6 [M].
4. **The slab box itself as a spherical set** **[P]**.
   * The weight-`s` layer of `A subset {-1,0,1}^d` lies on the sphere `|k|^2 = s`.
   * So `A` is a union of `<= d` concentric lattice-point spheres of radius `<= sqrt d` in dimension `d`.
   * That is the opposite corner from restriction theory (fixed dimension, radius tending to infinity), and
     F1 shows the wave / Bochner-Riesz side is scalar on each layer.
   * The paraboloid lift `k -> (k, |k|^2)` used in decoupling has no scale separation here: the radius is
     `~ sqrt d` and the separation is `1`.
5. **Many shells.** The only place in the record where several shells enter is the pointwise lattice
   separation `|sin(theta_j - arg alpha)| >= (p_j Norm alpha)^(-1/2)` (round 10 Theorem D, external lattice
   points `alpha`). It is an exact integrality statement. It already collapses fakes onto actual clusters, and
   it returns the original problem. An averaged (L^p) version of it is weaker (F3).

## 5. What a harmonic-analysis input would need, to help

The following necessary conditions are [P] where stated as consequences of F1-F4, and [S] as a checklist. An
input can help only if it:

* (a) is applied to an object that sees integrality and spans **many** frequency shells. Single eigenspaces
  are vacuous (F1);
* (b) yields **pointwise or exact-counting** information, or a measure statement that fakes violate. Plain
  measure statements are capped (F3, Cor. B.1);
* (c) after lifting, has dimension dependence `e^(o(n log n))` and precision `e^(-W/4)` (F4);
* (d) for single-circle L^p, works at `p ~ log R/loglog R` with uniform constants (F2).

The new results fail (a), (b) and (c), and do not address (d). Their sharper features fail the same tests:
the log-loss critical estimates, the periodic-cap estimate on `T^n`, and the PDE method.

## 6. Bottom line for the user's question

The claimed theorem, if correct, is a landmark in Euclidean harmonic analysis. It does not bear on the
Cilleruelo-Granville problem, its growth rate, or any of our three open routes. Every link we can name either:
* is vacuous on a single circle (F1);
* runs in the wrong direction: cluster bounds imply eigenfunction and nodal bounds, not conversely (F2);
* asks for a lower bound where the theorem gives an upper bound (F3);
* lives in fixed dimension with curvature, where our lifts are high-dimensional and flat (F4).

Measure-type corollaries are already capped by the fake-circle barrier (Cor. B.1).

No action recommended. If a referee confirms the result, the log-loss Kakeya and extension bounds are worth
recording as background only.

## Sources

* `https://github.com/openai/math` (commit `7ef5c0cc`):
  * `preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/` (`paper.pdf`, `build/sections/*.tex`);
  * `lean/ComparatorChallenges/LocalSmoothing.lean`.
* OpenAI announcement page (snippet only; fetch failed):
  `https://openai.com/index/building-scientific-understanding-with-academic-communities/`.
* Dan Roberts (OpenAI) on X, summary of the claim (snippet): `https://x.com/danintheory/status/2114833685417689321`.
* Press, snippets only: `https://www.unite.ai/openai-says-its-models-produced-three-proofs-in-math-cs-and-physics/`;
  `https://qz.com/openai-ai-discoveries-math-computer-science-physics-102726`;
  `https://www.engadget.com/2297941/openai-new-discoveries-mathematics-computer-science-physics/`;
  `https://www.mindstudio.ai/blog/openai-ai-discoveries-math-physics-cs`.
* Prior state of the art (search results): arXiv 2502.05973 (local smoothing for `p >= p(n)`); Guth-Wang-Zhang,
  Annals 192 (2020).
* Cilleruelo-Granville and eigenfunctions (search snippets): Bourgain-Rudnick arXiv 1012.3843 (nodal lines;
  `B_lambda`; CG for most `m`), arXiv 0907.4824 / 1105.0018 (restriction to hypersurfaces; caps via Jarník),
  arXiv 1308.4247, Granville-Wigman arXiv 1612.07819.
* Record: `docs/claude_review/README.md`; `round9/README.md`; `round10/analytic.md` (Thm B, Cor. B.1, Thm D);
  `round11/README.md` §1; `round12/README.md` §2.3, §4.2-4.3, §7; `round13/README.md` §0; `round14/README.md`
  §0-1; `round15/README.md`.
