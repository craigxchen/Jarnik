# Route 5, round 16: polynomial anti-concentration / inverse Littlewood-Offord theory versus AC(u)

## 0. Summary and verdict

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` stands. `AC(u)` is neither proved
nor refuted for any regime box in the form MIN^pair needs (`AC_A` on the gap region `G`). Setting throughout:
`T = 1`, `w = sqrt(8(n-1) ln 4r)`, `j = floor(w/2)`, `m = m' = ceil(4w)`, regime `n - 1 >= 316 ln 4r`, scale
`delta_0 = d^m eps_0`, `C = 1` where numbers are given.

What this round establishes (all proved unless marked):
1. **Exact bookkeeping (Prop 1.2, Remark 1.3).** The sup form `AC_A(u)` *is* `kappa_V(u) > eps_0` (MIN^pair at
   `H = infinity`); the steelman's density form `AC(u) = AC_T(u; 1/2)` implies it and is strictly stronger; the weakest
   density form that still implies it uses the threshold `mu(A)`. Dual form: two measures on `A` with equal moments of
   order `< m` and different `m`-th moments along `u`. The `kappa_H`/`kappa_V` gap cannot be closed by counting
   (Minkowski + Vaaler leave codimension `E^{1-o(1)}` even at the largest height rigidity tolerates).
2. **Literature (§2).** Every anti-concentration or inverse Littlewood-Offord theorem located (Carbery-Wright; MOO;
   DSTW regularity; Meka-Nguyen-Vu; Razborov-Viola; Fox-Kwan-Sauermann; Kwan-Sauermann; Jin-Kwan-Sauermann-Wang;
   Grebennikov-Kwan; Kane; Nguyen-Vu; Filmus-Mossel; dimension-free Remez) is bounded-degree (MNV: non-trivial only
   up to degree `c log r/loglog r`, by their own account), continuous, or degree 1; the two that allow growing degree
   (Filmus-Mossel slice invariance, dimension-free Remez) are distributional comparison or norming statements on
   product-type sets and give no small-ball bound at scale `e^{-Theta(n ln n)}` on `A`. Our degree is
   `m ~ 11.3 sqrt(d ln 4r)` (`>= 1864`). The structures in the inverse theorems (low rank, few polynomials, bounded Chow
   rank) are structures of `Q`; route 5 needs structure of the leading form `l_u^m` alone, against an arbitrary
   lower-degree adversary. No located theorem has this "prescribed leading form" shape.
3. **Degree barrier (Theorem 3.1).** At our degree there are maximally regular polynomials on the ternary cube,
   with all `L^2` mass on the top chaos and `b` disjoint maximal monomials, that vanish on all but a `10^{-328}`
   (`n = 2900`) to `10^{-2088}` (`n = 10^6`) fraction of the cube (exact). So no hypothesis of degree / influence /
   chaos-weight / rank type can give anti-concentration at degree `m`; only the power-of-a-linear-form structure can.
   The degree is not an artefact of the record's slab width: by a heuristic count it stays `Theta(sqrt d)` for every
   admissible width when `r >= n/4` (Remark 3.3, evidence).
4. **The density forms are false at regime sizes (Theorem 4.4, finite `n`).** Through an interpolation lemma
   (Lemma 4.1) the density form has strictly more exceptional directions than `A`. On the full cube at density
   `1/2` these include partition strata with up to `~50-90` classes (Prop 4.2) and **bit strata** spanned by `q` binary
   digits, with `2^q` classes: `q` up to 11 at `n = 2^12` (classes of size 2) and 14 at `n = 2^22` (Prop 4.3). For
   `(n, r) = (2^k, 2^k - 1)`, `k = 12..22`, at `W = 3n ln n`, the steelman's AC fails at
   explicit points of the (MIN) region `{Gamma <= abar}`; the best **profile-free** density form fails there for
   `n = 2^12..2^17`. At the Sylvester-specific density `mu(A) >= 1 - 2r e^{-3w^2/(8n)}` the exceptional bit strata are
   exactly the record's cheap range (`<= 4` columns) for all `k <= 60`. So an anti-concentration statement that is
   true at regime sizes must know `mu(Cube \ A)` to polynomial accuracy in `r`, i.e. it must see the slabs. Whether
   the constructed points lie in `G` depends on pair data (for genuine Sylvester data `G` is empty), so **MIN^pair is
   not touched**.
5. **The task's template (a)-(d) (§5).**
   * (a)+(c) work exactly: **Lemma 5.1 (fibre restriction)** moves any point near a rational stratum to a fibre
     problem at loss `eta^m`, with certificates extended by zero (height preserved, so compatible with `kappa_H`).
     The loss fits in the room for 16-class cheap spans from `n ~ 10^6` at the worst exponent `m`.
   * (d) as a heuristic reproduces the conjecture exactly: if conditional moments given the slab statistics were
     polynomial, `kappa_A` would reduce to the grid box `{-j..j}^r`, whose exceptional set is exactly the spans of
     `<= floor((m-1)/(2j)) = 4` columns, with exponent `m` (**Lemma 5.3**, proved for the grid). Quantitatively it fails:
     the lift needs precision `e^{-(m ln(2de/m) + m ln n)}` (`e^{-1.9*10^4}` at `n = 2900`), fixed-order LLT/Edgeworth
     give `n^{-O(1)}`, and in the main regime `r >= d` the fibres are points. Off the column span, `r < d` gives
     exact kernel-zonotope certificates, but the actual point lies in the column span ([D] Lemma 3.2), where they
     vanish.
   * (b) is the whole problem: the **regular core** needs a non-product dual certificate at degree `m` (two measures
     on a fibre of `A` matching moments `< m`). Product (zonotope) certificates do not exist in LP-supercritical boxes,
     collision-form ones do not exist for `W > 8 ln C'` (record), and every located theorem is proved by differencing
     (which produces products) or by continuous comparison (precision `n^{-O(1)}`).
6. **Ledger (§6).** The room is asymptotically huge (`(3/4) n ln n` against `O(m ln n)`, so polynomial separation
   `n^{-K}` with `K <~ 0.066 sqrt(n/ln 4r)` is affordable), but it opens late: even an ideal exponent-`m` inequality at
   distance `1/n` does not fit for `n <= 10^4` (`L0 - Lc - Lpoly < 0`) and fits from `n = 2^16`. Nothing in the chain is
   limited by the room; it is limited by the absence of any certificate at degree `m`.

**The missing theorem (§5.4), stated so it can be attacked:** an inverse Littlewood-Offord theorem at degree
`m ~ sqrt(d ln r)`, for polynomials with *prescribed leading form `l_u^m`*, on the slab box itself: smallness of
`l_u^m + Poly_{<m}` on `A` at level `kappa` forces `u` within `(kappa/kappa_*)^{1/m}` of a projection-cheap rational
subspace, with `ln(1/kappa_*) = O(m ln n)`. With the record's separation lemma this implies MIN^pair (`H = infinity`)
for width-cheap strata from `n ~ 10^6`; collision counts of the other projection-cheap strata remain open. By fake
Thm 2.1 it is equivalent in difficulty to the barrier question. No literature result comes close.

**Fake-test (rubric b).** Nothing here uses an input that fakes violate: AC is a statement about the box, and the only
arithmetic used (Step 0, separation, the column-span relation) is the record's. The column-span relation of §5.3(3) is
fake-violated only when the profile rank is `< d`; nothing here exploits it (it only shows that kernel certificates
vanish at the actual point).

| item | outcome | where |
| --- | --- | --- |
| survey of polynomial anti-concentration / inverse LO | all located results bounded-degree, continuous or linear | §2 |
| best applicable theorem | none applies; MNV's own range `c log r/loglog r` | §2, §3 |
| degree barrier | proved (Theorem 3.1) | §3.2 |
| density form of AC | false at regime sizes on the (MIN) region, `n = 2^12..2^22` (`theta = 1/2`), `2^12..2^17` (profile-free) | Theorem 4.4 |
| structured part (a)+(c) | proved reduction (Lemma 5.1), loss in the room from `n ~ 10^6` | §5.1 |
| conditioning (d) | grid model proved (Lemma 5.3); transfer needs `e^{-O(m ln n)}` precision: gap | §5.3 |
| regular core (b) | no method; obstruction explained | §5.2 |
| AC_A / MIN^pair | open | §5.4 |

## 1. The target, restated exactly

**Setting (fixed throughout unless a statement says otherwise).** `T = 1`; `d = n - 1`; `w = sqrt(8 d ln 4r)`;
`j = floor(w/2)`; `m := m' = ceil(4w)`; regime `d >= 316 ln 4r` (so `m <= m_max`); any profile `beta_p in [-1,1]^n`;
base point `n`, `u_n := 0`; `Cube := {-1,0,1}^d` (coordinates `k' = (k_1..k_d)`, `k_n = -sum k'`), `mu` = uniform
measure on `Cube`; `s_p(k') = <k, beta_p> = sum_{x<n} k_x (beta_{p,x} - beta_{p,n})`;
`A = {k' in Cube : |s_p(k')| <= w for all p}`; `l_u(k') = <k', u>`; `Poly_{<m}` = real polynomials of degree `< m`
in `k'`; `eps_0 = C d e^{-W/4}/(m+1)`; `delta_0 := d^m eps_0` (the *scale*); `G` = the gap region
`{u : |u|_inf = 1, |u_x - u_y| >= e^{-eps_xy/2} for all pairs}` of the cluster's pair data; `Gamma`, `abar`, `N_2`,
`P(L)` as in `round12/README.md` §1. *Proved* = proved here with all quantifiers, or refereed in the record and cited.
*Exact* = integer/rational computation. *Evidence* = anything else. **[L]** = literature located this session by
WebSearch (snippets/abstracts only; arxiv.org was not reachable by WebFetch), **[M]** = quoted from memory.

**Definition 1.1 (the anti-concentration statements).** For `u in R^d`, `|u|_inf = 1`, and `Q` ranging over
`l_u^m + Poly_{<m}` (polynomials with leading form `l_u^m`):
* `AC_A(u)` (sup form): every such `Q` has `max_A |Q| > delta_0`;
* `AC_T(u; theta)` (density form, threshold `theta in (0,1]`): every such `Q` has `mu{k' in Cube : |Q(k')| <= delta_0} < theta`.

The steelman's `AC(u)` (`local_smoothing/steelman.md` S3) is `AC_T(u; 1/2)`.

**Proposition 1.2 (what implies what; proved).**
1. `AC_A(u)` is equivalent to `kappa_V(u) > eps_0`, where `kappa_V` is `kappa` taken over all of `V(A, m)`.
2. `AC_T(u; theta)` implies `AC_A(u)` whenever `theta <= mu(A)`; in particular `AC_T(u; 1/2) => AC_A(u)` (since
   `mu(A) >= 1/2`, t1rigour §2).
3. `AC_T(u; theta) => AC_T(u; theta')` for `theta <= theta'` (smaller `theta` = stronger statement). The steelman
   form `theta = 1/2` is a strong choice; `theta_A := mu(A)` is the **weakest** density form that still implies `AC_A`.
4. (Moment-matching form.) `d^m kappa_V(u) = sup { (E_{nu+} l_u^m - E_{nu-} l_u^m)/2 }` over pairs of probability
   measures `nu+, nu-` on `A` whose moments `E k'^alpha` agree for all `|alpha| < m`.
5. `AC_A(u)` for all `u in G` is MIN^pair in the `H = infinity` form. MIN^pair itself (Theorem 4.1 of nearstrata as
   strengthened by its referee) needs `kappa_H(u) > eps_0`, and `kappa_H <= kappa_V`.

*Proof.* 1 is nearstrata Prop 2.1 (`d^m kappa_V(u) = dist_inf(l_u^m|_A, Poly_{<m}|_A)`), since
`Q = l_u^m - P` ranges over `l_u^m + Poly_{<m}`. 2: if `mu{|Q| <= delta_0} < theta <= mu(A)` then some `k' in A` has
`|Q(k')| > delta_0`. 3 is monotonicity. 4: the constant `1` lies in `Poly_{<m}`, so every `c in V(A,m)` has
`sum c_k = 0`; write `c = (|c|_1/2)(nu+ - nu-)` with `nu+-` the normalised positive and negative parts. Then
`c perp Poly_{<m}` iff `nu+` and `nu-` have equal moments of order `< m`, and
`|Phi_c(u)|/|c|_1 = |E_{nu+} l_u^m - E_{nu-} l_u^m|/2`; take the supremum (nearstrata (1.x) definition of `kappa`).
5 is the definition plus nearstrata Prop 2.1 (`kappa_H <= kappa_V`). QED

**Remark 1.3 (the `H`-gap is real; proved bound).** Let `E = dim Poly_{<m}|_A` (`<= binom(m+d-1, d)`) and `U = |A|`.
The lattice `Lambda = V(A,m) cap Z^A` (rank `N = U - E`) is primitive, so
`covol(Lambda) = covol(V^perp cap Z^A) <= prod (row norms of an E x U row basis of the {-1,0,1} moment matrix) <= U^{E/2}`.
Minkowski's second theorem for the body `V tensor R cap [-1,1]^A`, whose `N`-volume is `>= 2^N` by Vaaler's cube-slicing
theorem, gives successive minima `1 <= lambda_1 <= ... <= lambda_N` (sup norm) with `prod lambda_i <= covol(Lambda) <= U^{E/2}`
[M: Minkowski II, Vaaler 1979, `covol(Lambda) = covol(Lambda^perp)` for primitive lattices; Bombieri-Vaaler's Siegel
lemma gives the same shape]. Hence

```text
codim_{V(A,m)} V_H  <=  #{ i : lambda_i > H }  <=  (E/2) ln U / ln H .
```

Rigidity tolerates `ln|c|_1` up to about `wW/2 - m ln(d C') + ln m!` (nearstrata Lemma 1.1), so `H` may be taken as
large as `exp(Theta(n^{3/2} (ln n)^{3/2}))`; but `ln E = Theta(w ln n)`, so the bound leaves a codimension of size
`E^{1-o(1)}`. **No counting argument closes `kappa_H` versus `kappa_V`.** Everything below is about `AC_A`/`AC_T`, i.e.
the `H = infinity` form; any proof must additionally produce its certificate `c` with integer entries of height
`<= exp(wW/2 - O(m ln n))` (this is how a certificate-producing proof would close the gap; a Hahn-Banach proof does
not).

**Remark 1.4 (scale bookkeeping).** The three scales that matter, all at `T = 1`, `m = ceil(4w)`:
* *Step-0 scale* `ln(1/eps_0) = W/4 - ln(d/(m+1)) >= (3/4) n ln n - (7/2) n - O(ln n)` (`C <= 1`, `W >= W_min = 3n ln n - 14 n`);
* *certificate scale*: a certificate of "zonotope type" (nearstrata Lemma 2.3) gives `kappa >= (m!/(2d)^m) prod|<v_i,u>|^{alpha_i}`,
  i.e. `ln(1/kappa) = m ln(2de/m) + m ln(1/eta)` at distance `eta` from the strata;
* *room* := `ln(1/eps_0) - m ln(2de/m)`; it is `Theta(n ln n)` against a loss budget `O(m ln(1/eta)) = O(sqrt(n ln n) ln(1/eta))`.
Table 5.1 gives the numbers. Any argument whose loss is `e^{-O(m ln n)}` at polynomially separated `u` has room to
spare; the problem is that no argument of any loss is available at degree `m`.

## 2. Literature survey: what is actually proved, measured against our parameters

Our parameters (numbers in §3.1 and Table 6.1): degree `m = ceil(4w) ~ 11.3 sqrt(d ln 4r)` (1864 at `n = 2900`, 44112 at
`n = 10^6`); dimension `d -> infinity`; domain the ternary cube restricted to `A` (measure `>= 1 - (4r)^{-1/2}/4`
profile-free, `>= 1 - 2r(4r)^{-3d/n} ~ 1 - 2r(4r)^{-3}` for balanced sign columns); scale `delta_0/sigma` with
`ln(sigma/delta_0) ~ W/4 ~ (3/4) n ln n`; required conclusion: a structural description of the exceptional `u`
(collision/cheap strata) with **quantitative** distance control at scale `e^{-O(log n)}` per pair.

| result | setting | degree range where non-trivial | form of the structure | applicable? |
| --- | --- | --- | --- | --- |
| Carbery-Wright 2001 [L, as quoted in arXiv 2108.04268] `mu{|f| <= t} <= C min(d_f, n) (t/||f||_2)^{1/d_f}` | log-concave / convex bodies, continuous | all degrees | none needed (no exceptional set) | only the Gaussian proxy (steelman S3); see §3.1 for the room |
| Mossel-O'Donnell-Oleszkiewicz, Ann. Math. 171 (2010) 295-341 [L]: invariance for low-influence, bounded-degree multilinear polynomials | product spaces | error `O(d_f tau^{1/(8 d_f)})` [M for the exponent]: needs `d_f << log(1/tau) <= log n` | "regular" = small influences | **no**: error `>= m e^{-ln d/(8m)} ~ m` (`params_out.txt`) |
| Diakonikolas-Servedio-Tan-Wan, regularity lemma for degree-d PTFs, CCC 2010 / Theory Comput. 10 (2014) 27-53 [L] | Boolean cube | constant `d_f` (number of pieces grows like `(1/tau)^{O(d_f)}` [M]) | restriction to few variables + regular remainder | **no** (degree) |
| Meka-Nguyen-Vu, Theory Comput. 12 (2016) art. 11, arXiv 1507.00829 [L] | independent variables, multilinear, "rank" = number of disjoint max-degree monomials | **"non-trivial for degrees up to `c log r/loglog r`"**; they note no non-trivial bound is possible near degree `log_2 r` [L, quoted from the arXiv text] | none (forward bound only) | **no**: `log r/loglog r <= 8.3` at `n = 10^12`, against `m >= 1864` |
| Razborov-Viola, "Real advantage" (ECCC 2012/134) [L] | Boolean cube | degree `~ (1/2) log_2 log_2 n` | many disjoint maximal monomials (Costello-Tao-Vu type) | **no** |
| Fox-Kwan-Sauermann, Math. Proc. Camb. 171 (2021) 227-248, arXiv 1905.12142 [L] | Bernoulli, nonnegative coefficients / Poisson-type | fixed degree | — | **no** (sign condition, fixed degree) |
| Kwan-Sauermann, Discrete Analysis 2020:12, arXiv 1909.02089 [L]; Kwan-Sauermann, Compositio 2025 [L] | Rademacher, quadratic | degree 2 | point probability `>> n^{-1+o(1)}` forces closeness to a **low-rank quadratic form**; optimal `O(n^{-1/2})` | degree 2 only |
| Jin-Kwan-Sauermann-Wang, arXiv 2505.23335 (2025) [L] | Rademacher, degree `d_f` | fixed degree | MNV improvable unless `f` is close to a polynomial with special algebraic structure ("robust irreducibility") | fixed degree; point probabilities, not small balls at scale `e^{-n ln n}` |
| Grebennikov-Kwan, arXiv 2505.24699 [L] | bounded **Chow rank** (polynomial of boundedly many linear forms) | fixed degree | lattice point counting in low dimension | the closest *structural* match (our structured part is a function of few linear forms) but fixed degree and bounded rank |
| Grebennikov, arXiv 2610.08708 (Oct 2026) [L, seen only through a MathWorld entry; **not verified**] | polynomial LO | "for fixed degree" `C r(F)^{-1/2}` | — | fixed degree |
| Kane, Ann. Probab. 45 (2017) 1612-1679, arXiv 1204.0543 [L] | **Gaussian** chaos | fixed degree; decomposition size effective but huge | poorly anticoncentrated `p` is approximated by a function of a bounded number of polynomials ("diffuse decomposition") | **no**: Gaussian (our continuous proxy has no exceptional set at all, §3.1) |
| Nguyen-Vu, Adv. Math. 226 (2011) 5298-5319, arXiv 1004.3967 [L]; Halasz, Sarkozy-Szemeredi [M] | linear forms | degree 1 | GAPs | only for the degree-1 statistic `l_u` (cannot see `Poly_{<m}`) |
| Filmus-Mossel, PTRF 175 (2019) 721-782, arXiv 1507.02713 [L] | Boolean **slice** vs cube | harmonic multilinear, **degree `o(n)`** | — | the only high-degree transfer found; one linear constraint, Boolean, harmonic: see §5.4 |
| Slote-Volberg-Zhang (Discrete Analysis 2025:4, arXiv 2305.10828); Becker-Klein-Slote-Volberg-Zhang, arXiv 2310.07926 [L] | dimension-free Remez: sup over a product domain `<= C(d_f, K) sup` over a product test set ("norm designs") | constants depend on degree, not dimension | — | wrong direction and wrong sets: they bound the **continuous** sup by a discrete product sup; `A` is not a product and is not norming at degree `>= 3` (`k_x^3 - k_x` vanishes on it; steelman S4) |

**What the theory proves, in one sentence.** Every discrete anti-concentration or inverse theorem located is
either (i) a statement about bounded degree (or degree `<= c log n/loglog n`), or (ii) a continuous (Gaussian,
log-concave) statement, or (iii) a degree-1 statement, or (iv) a growing-degree comparison/norming statement
(Filmus-Mossel, dimension-free Remez) with polynomial precision on product-type sets. Our degree is `m = Theta(sqrt(d ln r))`, and §3 proves that
this is not a technicality: at this degree the hypotheses those theorems use (degree, influences, variance, weight
on the top chaos, number of disjoint maximal monomials) are compatible with a polynomial vanishing on all but a
`10^{-328}` fraction of the cube.

**Best applicable theorem, by component.** (1) For the continuous proxy: Carbery-Wright, with room only from
`n ~ 2 * 10^4` (§3.1). (2) For the degree-1 statistic `l_u`: Halasz / Nguyen-Vu (GAP structure); irrelevant because
`Poly_{<m}` is invisible to it. (3) For the structured part (`u` near a span of few integer directions, `Q` then a
function of few linear forms on fibres): the bounded-Chow-rank lattice-point method of Grebennikov-Kwan is the right
*shape*, but at fixed degree; in our setting this part is handled exactly anyway (Lemma 5.1). (4) For slab
conditioning: Filmus-Mossel (one Boolean slice, degree `o(n)`), the only growing-degree transfer, with CDF-type
precision. **For the regular core nothing applies** (§3, §5.2).

**The form of the structure in the inverse theorems, against what we need.**
* Quadratic (Kwan-Sauermann): low rank. Gaussian (Kane): a function of few lower-degree polynomials. Bounded Chow rank
  (Grebennikov-Kwan): a function of few linear forms.
* What route 5 needs: `u` (not `Q`) is close to a span of few integer directions whose projections of `A` have
  interpolation degree `< m` (the projection-cheap spans; t1rigour §5, round 13 loctheory). In our class `Q` is never
  structured: the adversary `P in Poly_{<m}` is arbitrary. The structure must be read off from the *leading form*
  `l_u^m` alone. No located theorem has a hypothesis or conclusion of this "prescribed leading form" type.

## 3. Why the located theorems cannot apply: three barriers at degree `m`

### 3.1 The continuous proxy has room only asymptotically

For the Gaussian proxy the steelman's computation (S3) gives `Pr(|Q| <= delta_0) <= C m (delta_0/sigma)^{1/m}` with
`sigma >= sqrt(m!)`. With the optimistic constant `C = 1` and `W = W_min` (`params.py`, 50-digit arithmetic):

| `n` | `r` | `m` | `log10` of the CW bound (needs `< log10(1/2) = -0.30`) |
| --- | --- | --- | --- |
| 2900 | `n` | 1864 | +3.64 |
| 4096 | `n-1` | 2256 | +3.35 |
| `10^4` | `n` | 3683 | +1.98 |
| `10^4` | `n^2` | 5035 | +3.13 |
| `2^16` | `n-1` | 10231 | -6.36 |
| `10^6` | `n` | 44112 | -59.0 |

(The skeptic's table, `local_smoothing/skeptic.md` §3.2, uses `W = 3 n ln n`; the crossover moves from `~6000` to
`~2 * 10^4` at `W_min`.) This is evidence about a proxy and is not used below.

### 3.2 The degree barrier (proved)

**Theorem 3.1 (regular, top-chaos, almost-everywhere-vanishing polynomials at degree `m`).** Let `1 <= m <= d`,
`b = floor(d/m)`, `B_1, ..., B_m` disjoint `b`-subsets of `[d]`, and

```text
Q(k') = prod_{i=1}^m  k(B_i) / sqrt(2b/3),        k(B) = sum_{x in B} k_x,   k' uniform on Cube.
```

Then:
1. `Q` is multilinear and homogeneous of degree `m`; in the orthonormal basis `phi_0 = 1`, `phi_1(k) = sqrt(3/2) k`,
   `phi_2(k) = (3k^2 - 2)/sqrt 2` of `L^2({-1,0,1})`, **all** of its `L^2` mass lies on the top chaos (degree exactly `m`);
   `E Q = 0`, `E Q^2 = 1`.
2. Every coordinate `x in B_i` has influence `Inf_x(Q) = E Var_{k_x}(Q) = 1/b`; the total influence is `m`, so the
   largest influence is a fraction `1/(mb) ~ 1/d` of the total: `Q` is maximally regular.
3. `mu{Q != 0} = (1 - p_b)^m` with `p_b = T_b/3^b` (`T_b` the central trinomial coefficient); `p_1 = 1/3` and
   `p_b >= 0.4/sqrt b` for `2 <= b <= 3000` (exact check). Hence `mu{Q != 0} <= exp(-0.4 m^{3/2}/sqrt d)` whenever
   `2 <= b <= 3000`.
4. At the regime parameters (`m = ceil(4w)`, exact `p_b`, `barrier_product.py`):

| `n` | `r` | `m` | `b` | `log10 mu{Q != 0}` |
| --- | --- | --- | --- | --- |
| 2900 | `n` | 1864 | 1 | -328.2 |
| 4096 | `n-1` | 2256 | 1 | -397.3 |
| `10^4` | `n` | 3683 | 2 | -648.5 |
| `2^16` | `n-1` | 10231 | 6 | -955.1 |
| `10^5` | `n` | 12850 | 7 | -1105.4 |
| `10^6` | `n` | 44112 | 22 | -2088.5 |

*Proof.* 1: the factors are independent, centred, with `E k(B_i)^2 = 2b/3`; `k(B_i) = sqrt(2/3) sum_{x in B_i} phi_1(k_x)`,
so `Q` is a combination of products of `m` distinct `phi_1`'s, which have degree exactly `m`. 2: for `x in B_i`,
`Var_{k_x} Q = (2/3) prod_{i' != i} k(B_{i'})^2/((2b/3)^m)`, whose mean is `(2/3)/(2b/3) = 1/b`. 3: `Q = 0` iff some
`k(B_i) = 0`; the events are independent with probability `p_b = T_b 3^{-b}`; `1 - p <= e^{-p}`. 4: exact `T_b`
(integer recursion), logs to 40 digits. QED

**Corollary 3.2 (what the barrier excludes; proved).** Any theorem of the form "if `f` has degree `<= m`, variance
`1`, every influence `<= tau`, top-chaos weight `>= 1 - epsilon`, and (after scaling the top coefficients to modulus
`>= 1`) at least `R` pairwise disjoint maximal-degree monomials, then `mu{|f| <= t} <= F(m, tau, epsilon, R, t)`" must
satisfy `F(m, 1/b, 0, b, 0) >= 1 - e^{-0.4 m^{3/2}/sqrt d}` (for `2 <= b <= 3000`). Indeed `Q` of Theorem 3.1 has
influences `1/b`, no lower chaos, and `b` disjoint transversal monomials of equal coefficient (rank `b` in the sense of
Meka-Nguyen-Vu). So such a theorem gives nothing once `m >> d^{1/3}`; our `m ~ 11.3 sqrt(d ln 4r)` is far beyond.
This matches MNV's own remark that their bound is non-trivial only up to degree `c log r/loglog r` [L].

**What survives.** `Q` is a product of `m` linear forms, not an element of `l_u^m + Poly_{<m}`. The only feature of
route 5's class that Theorem 3.1 does not defeat is the **prescribed leading form `l_u^m`** (a power of a single
linear form). Any proof of `AC` must use that feature; regularity, invariance and chaos-weight arguments do not see
it. §4 shows that the class itself, at density `1/2`, already contains polynomials that vanish on half of the cube
for collision-free `u`.

### 3.3 Remark: is the degree forced? (evidence)

The degree enters through two constraints of the `T = 1` method (nearstrata Lemma 1.1): rigidity needs
`m > 2w` (Norm `>= N^{2d} e^{-wW}` against Taylor `e^{-mW/2}`), and Siegel needs `|A| >= 4 binom(m+d-1, d)`. The record's
`w = sqrt(8 d ln 4r)` is the union-bound choice that makes `mu(A) >= 1/2`; it is not forced. A thinner slab keeps
Siegel alive as long as `ln(1/mu(A)) <= d ln 3 - m ln(ed/m) - ln 4`. If the `r` slab events behaved independently
(true for the Gaussian proxy by the Gaussian correlation inequality; **not proved** for the ternary cube), this allows
`w ~ sigma exp(-d ln 3/r)` with `sigma ~ sqrt(4d/3)`, i.e. `m ~ 2w ~ c(r/d) sqrt d` with `c(rho) ~ 2.3 e^{-1.1/rho}`.
For squarefree norms `r >= n/4 - 2` ([D] Lemma 3.1), so `m >= ~0.03 sqrt d`: **the degree stays `Theta(sqrt d) >> log d`
in the main regime**, and the barrier of §3.2 stands. Degrees `O(log d)` would need `r = O(d/log d)`, i.e. low-rank
profiles, which [D] Prop 5.2 already shows are invisible to Theorem 4.1. (Heuristic; not used.)

## 4. The density form `AC_T` has a strictly larger exceptional set, and is false at regime sizes

### 4.1 The interpolation lemma

**Lemma 4.1 (transversal smallness on interpolation sets; proved).** Let `v_1, ..., v_q in Z^d`, `pi(k') = (<k',v_i>)_i`,
`Y` a subset of `Cube`, and `rho` an integer such that every function on `pi(Y)` is the restriction of a polynomial of
degree `<= rho` on `R^q`. Let `y in span_R{v_i}` and `v in R^d` with `|y|_inf, |v|_inf <= 1`, `0 < eta <= 1/2`,
`u = y + eta v`, and `o = m - rho >= 1`. Then there is `Q in l_u^m + Poly_{<m}` with

```text
|Q(k')| <= 2^m eta^o d^m      for every k' in Y .
```

If moreover `|u|_inf >= 1` then `Q/|u|_inf^m in l_{hat u}^m + Poly_{<m}` and the same bound holds for it; so
`AC_T(hat u; theta)` fails as soon as `mu(Y) >= theta` and `2^m eta^o <= eps_0` (and `AC_A(hat u)` fails if `Y = A`).

*Proof.* `l_u^m = sum_{i=0}^m binom(m,i) eta^i l_y^{m-i} l_v^i`. Since `y = sum t_i v_i`, `l_y^{m-i} = F_i o pi` with
`F_i(z) = (sum t_i z_i)^{m-i}`; let `R_i` be a polynomial of degree `<= rho` with `R_i = F_i` on `pi(Y)`, so that
`R_i o pi` has degree `<= rho` in `k'`. Put `Q = l_u^m - sum_{i<o} binom(m,i) eta^i (R_i o pi) l_v^i`. Each subtracted
term has degree `<= rho + i < m`, so `Q in l_u^m + Poly_{<m}`. On `Y`, `Q = sum_{i>=o} binom(m,i) eta^i l_y^{m-i} l_v^i`,
and `|l_y|, |l_v| <= |k'|_1 <= d`, so `|Q| <= d^m sum_{i>=o} binom(m,i) eta^i <= 2^m eta^o d^m`. QED

(This is nearstrata Lemma 5.0 with `A` replaced by an arbitrary `Y`; exact toy check in `lemma41_check.py`.)
The useful `Y` are boxes in few statistics: if `pi(Y) subset prod_i [a_i, a_i + 2c_i] cap Z^q`, tensor Lagrange
interpolation gives `rho <= sum_i 2c_i`.

### 4.2 Two families of `AC_T`-exceptional subspaces

Tails below use the rigorous sub-Gaussian bound `mu(|S_N| >= t) <= 2 exp(-3t^2/(4N))` for a sum of `N` independent
uniform `{-1,0,1}` variables (from `(1 + 2 cosh s)/3 <= e^{s^2/3}`, which holds termwise:
`2/(3 (2i)!) <= 1/(3^i i!)` for all `i >= 1`).

**Proposition 4.2 (partition strata; proved).** Let `X_1, ..., X_q` be disjoint subsets of `[d]` of size `>= N`, and
`L_pi = span{1_{X_i}}` (the stratum of `u` constant on each `X_i` and `0` elsewhere, a collision subspace with `q + 1`
classes). If integers `c_i` satisfy `prod_i (1 - 2 e^{-3(c_i+1)^2/(4N)}) >= theta` and `sum_i 2c_i <= m - 1`, then
every `u in L_pi`, `|u|_inf = 1`, violates `AC_T(u; theta)` (with `Q` vanishing exactly on a set of measure
`>= theta`). For `theta = 1/2` and equal `c_i` the largest admissible `q` is (`act_strata.py`):

| `n` | `r` | `m` | `q_max >=` |
| --- | --- | --- | --- |
| 2900 | `n` | 1864 | 49 |
| 4096 | `n-1` | 2256 | 49 |
| `10^4` | `n` | 3683 | 52 |
| `10^4` | `n^2` | 5035 | 89 |
| `2^16` | `n-1` | 10231 | 58 |
| `10^6` | `n` | 44112 | 69 |
| `10^8` | `n` | 503517 | 86 |

*Proof.* The construction of Lemma 4.1 with `eta = 0` (`u = y`, `Q = l_y^m - R_0 o pi`, which vanishes on `Y`), `Y = {|k(X_i)| <= c_i for all i}` (independent coordinates, so
`mu(Y) = prod_i mu(|k(X_i)| <= c_i)`), and `rho <= sum 2c_i`. QED

So on the full cube at density `1/2` the exceptional set contains collision strata with up to `~50-90` classes,
against `<= 16` classes for the cheap spans of `A` (t1rigour Thm 6.3). These are collision-rich
(`P(L_pi) >= N_2/(q+1) - O(n)`), so they are harmless for the separation argument; the next family is not.

**Proposition 4.3 (bit strata; proved).** Let `n = 2^k`, identify `[n]` with `F_2^k` (base point `0`), and for
`q <= k` let `X_i = {x : x_i = 1}` (`|X_i| = n/2`, `i <= q`) and `L_q = span{1_{X_1}, ..., 1_{X_q}}` (a collision
subspace with `K = 2^q` classes of size `n/K`; `L_k` is collision-free). Let `c` be the least integer with
`q * 2 e^{-3(c+1)^2/(2n)} <= 1 - theta`, `rho_q = 2cq`, `o_q = m - rho_q`. If `rho_q <= m - 1`, then every `u in L_q`
violates `AC_T(u; theta)` exactly (some `Q` vanishes on a set of measure `>= theta`), and for every `u = y + eta v`
with `y in L_q`, `|y|_inf <= 1`, `|v|_inf <= 1`, `|u|_inf >= 1`, `0 < eta <= 1/2` and `2^m eta^{o_q} <= eps_0`, the
normalised `hat u` violates `AC_T(hat u; theta)`.

*Proof.* Lemma 4.1 with `Y = {|k(X_i)| <= c, i <= q}`, `mu(Y) >= 1 - q * 2e^{-3(c+1)^2/(4(n/2))} >= theta` (union bound). QED

Three thresholds are relevant: `theta = 1/2` (the steelman's form); `theta_0 = 1 - (1/4)(4r)^{-1/2}`, the
**profile-free** lower bound for `mu(A)` (sub-Gaussian bound with `sum_x (beta_{p,x} - beta_{p,n})^2 <= 4d`); and
`theta_S = 1 - 2r e^{-3w^2/(8n)}`, the lower bound for `mu(A)` for a Sylvester profile (`sum_x (...)^2 = 2n`). On `A`
itself the bit statistics of a Sylvester profile are profile columns, cut at `|s| <= w`, i.e. `|k(X_i)| <= j`.

### 4.3 `AC_T` fails on the averaged region at regime sizes

**Theorem 4.4 (proved; finite `n`).** Let `n = 2^k`, `r = n - 1`, `C = 1`, `m = ceil(4w)`. For each `theta` and `k`
listed below there is `q` and a point `hat u` with `|hat u|_inf = 1` and `Gamma(hat u) <= abar = W/(2d)` at which
`AC_T(hat u; theta)` fails, at the stated `W`:

| `theta` | `k` (all checked: `12 <= k <= 22`) | failing `q` (examples) | `W` |
| --- | --- | --- | --- |
| `1/2` (steelman) | every `k in [12, 22]` | e.g. `k = 12`: `q = 2..7`; `k = 16`: `q = 4..10`; `k = 22`: `q = 8..13` | `3 n ln n` |
| `theta_0` (profile-free `mu(A)`) | `k in [12, 16]` at `3 n ln n`; `k = 17` at `4 n ln n` | `k = 12`: `q = 2..6`; `k = 16`: `q = 6` | `3 n ln n` / `4 n ln n` |
| `theta_S` (Sylvester `mu(A)`) | none (`q_max = 4` for all `12 <= k <= 60`) | — | — |

At `W = W_min = 3n ln n - 14n` none of the constructed points lies in `{Gamma <= abar}` for `k <= 22`.

*Proof.* `act_bits_theta.py` (rigorous tail bounds; logs in 50-digit arithmetic, margins visible in
`act_bits_theta_out.txt` and `act_strata_out.txt`). For each `(k, theta, q)` with `rho_q <= m - 1`:
take `eta = (eps_0 4^{-m})^{1/o_q}`; check `eta <= 1/(2(K-1))`; take `y = sum_{i<=q} 2^{i-1} 1_{X_i}/(K-1)` (values
`t/(K-1)`, `t = 0..K-1`, one per class, base class `0`) and `v` = rank within its class divided by the class size
`s = n/K` (`v_0 = 0`). Then `|u|_inf in [1, 1 + eta]`; pairs in one class have `|u_x - u_y| >= eta |Delta rank|/s`;
pairs in different classes have `|u_x - u_y| >= 1/(K-1) - eta >= 1/(2(K-1))`. Hence

```text
(N_2/2) Gamma(hat u) <= K [ binom(s,2) ln(s(1+eta)/eta) - sum_{t=1}^{s-1} (s-t) ln t ] + (N_2 - K binom(s,2)) ln(2(1+eta)(K-1)),
```

which the script evaluates and compares with `abar`. Proposition 4.3 then gives the failure. QED

**Corollary 4.5.** (i) For `(n, r) = (2^k, 2^k - 1)`, `12 <= k <= 22` (e.g. Sylvester profiles), at `W = 3 n ln n`,
the steelman's `AC_T(.; 1/2)` is false at points of the (MIN) region `{|u|_inf = 1, Gamma(u) <= abar}`; so (MIN)
cannot be proved by establishing `AC_T(.; 1/2)` on that region. (ii) The same holds for the best profile-free density
form `AC_T(.; theta_0)` for `12 <= k <= 17` (at `W = 3n ln n` for `k <= 16`, at `4 n ln n` for `k = 17`). (iii) Whether these points lie in the
gap region `G` depends on the pair data: they need excess `eps_xy >= 2 ln(s/eta) ~ 2(W/4 + m ln 4)/o_q` on the
`N_2/K` within-class pairs, which the averaged budget allows exactly when `K o_q >~ d` (the necessary condition of
nearstrata Prop 5.1). For genuine Sylvester pair data (all `eps_xy` nearly equal) `G` is empty (referee_nearstrata §2.3), so
nothing is said about MIN^pair.

**Asymptotics of the bit-strata test (exact integer `c`, rigorous tails; `act_bits_theta_out.txt`).** The necessary
condition `K o_q > d` for some exceptional `q` holds for `k <= 36` at `theta = 1/2`, for `k <= 19` at `theta_0`, and
never at `theta_S`. So for these families the density forms are false at finite sizes only; asymptotically the bit
strata are harmless, as are the partition strata (`P >= N_2/(q+1)`, `q <= 90`).

### 4.4 Reading

* **The density threshold carries the arithmetic of the slabs.** At `theta = 1/2` the bit statistics may be cut at
  `c ~ 0.2-0.3 j`; at the profile-free `theta_0` at `~0.45-0.5 j`; only at the Sylvester-specific `theta_S` at
  `~0.84 j`, which reproduces the record's cheap range (`<= 4` long directions, t1rigour Thm 6.3). So an
  anti-concentration statement strong enough to be true at regime sizes must know `mu(T \ A)` up to a polynomial in
  `r` — i.e. it must see the profile. A profile-free statement of the steelman type is false at regime sizes
  (Theorem 4.4).
* **The sup form is MIN itself.** `AC_A` is literally `kappa_V > eps_0` (Prop 1.2). Between the two sits
  `AC_T(theta_A)`, which quantifies over every set of measure `>= mu(A)`; it is still strictly stronger than `AC_A`
  (it includes the slab boxes of every other profile with the same `mu`, fake profiles included).
* **Fake test.** `AC_T(theta)` for every set of measure `>= theta` contains `AC_A` for the (Mar) fake pattern boxes of
  that density. By fake Thm 2.1 any proof of it decides the barrier question (as the steelman noted for `AC_A`).

## 5. Proof attempt: regularity + invariance + exact algebra + slab conditioning

The template the task proposes is: (a) decompose; (b) invariance/Carbery-Wright on the regular part; (c) exact
algebra on the structured part (cheap families); (d) condition on the `r` slab statistics with a local limit theorem.
In route 5's class the structure lives in the **direction** `u`, not in `Q` (the adversary `P in Poly_{<m}` is
arbitrary), so the decomposition is `u = y + eta v` with `y` in a rational subspace `L` of low projection cost. Below:
what is proved for each step, and exactly where the chain breaks.

### 5.1 (a)+(c): the structured part costs nothing beyond `eta^m` (proved)

**Lemma 5.1 (fibre restriction).** Let `v_1, ..., v_q in Z^d`, `pi(k') = (<k', v_i>)_i`, `L = span_R{v_i}`, and for
`t in pi(A)` let `F_t = A cap pi^{-1}(t)`. Let `u = y + eta v` with `y in L`, `eta > 0`. Then

```text
dist_inf( l_u^m |_A , Poly_{<m}|_A )  >=  eta^m  dist_inf( l_v^m |_{F_t} , Poly_{<m}|_{F_t} )      for every t in pi(A),
```

and a certificate `c in V(F_t, m)` for the right side, extended by `0`, is a certificate in `V(A, m)` for the left side
with the same `|c|_1` and the same entries (so integer height is preserved: the reduction is compatible with
`kappa_H`).

*Proof.* Restriction to a subset can only decrease the distance (a polynomial approximating on `A` approximates on
`F_t`). On `F_t`, `l_y = <t, tau>` is a constant `c_t` (`y = sum tau_i v_i`), so
`l_u^m = sum_{i<m} binom(m,i) c_t^{m-i} eta^i l_v^i + eta^m l_v^m`, and the first sum lies in `Poly_{<m}|_{F_t}`.
Hence `dist(l_u^m|_{F_t}, Poly_{<m}) = eta^m dist(l_v^m|_{F_t}, Poly_{<m})`. A vector `c` supported on `F_t` that
annihilates `Poly_{<m}|_{F_t}` annihilates `Poly_{<m}|_A`. QED

**Corollary 5.2 (structured part; proved implication).** Suppose that for some fibre `F_t` of some cheap family
`F` (t1rigour Thm 6.3: `<= 4` long directions, `<= 16` classes at `m = ceil(4w)`), and for all `v` in a set `R_F`
("regular relative to `F`"), `dist(l_v^m|_{F_t}, Poly_{<m}) >= d^m kappa_core`. Then every `u = y + eta v` with
`y in L_F`, `v in R_F`, satisfies `kappa_V(u) >= eta^m kappa_core`. With `eta >= sigma_F = (1/2) e^{-Q_F abar/2}`
(separation, nearstrata Cor 3.5) the loss is `m (Q_F abar/2 + ln 2) <= m (8 abar + ln 2)(1 + O(w/n))`.

So (a)+(c) are not the difficulty: exact algebra (Lemma 5.1) moves every near-stratum point to a fibre problem with an
explicit, `u`-independent loss. Table 6.1 shows the loss fits in the room for `n >~ 10^6` with the worst exponent `m`
(and much earlier with the true transversal orders; nearstrata §5.3). **But the fibre problem is again `AC` at full
degree `m` for a direction `v` that is regular relative to `F` — the core problem in a box of dimension
`d - O(1)`.** The decomposition does not reduce the degree.

### 5.2 (b): the regular core — no available tool (proved obstruction + identified gap)

Theorem 3.1 shows that degree, influence, variance, chaos-weight and rank hypotheses cannot give anti-concentration
at degree `m`; Theorem 4.4 shows that even within the class `l_u^m + Poly_{<m}`, density-only hypotheses (any set of
measure `theta <= theta_0`) admit collision-poor exceptional directions at regime sizes. The record adds:
* product (zonotope) certificates do not exist at all in LP-supercritical regime boxes (t1rigour 2.1(7):
  `zeta <= 2w(n-1)/n < m`), and when they exist one generator is a relation at the actual point (nearstrata
  Remark 2.10);
* collision-form certificates do not exist for `W > 8 ln C'` (round 13 fatpoint Thm 3.1).

So the regular core needs a **non-product dual certificate** at degree `m`: by Prop 1.2(4), two probability measures
on `A` (or on a fibre) with equal moments of all orders `< m` whose `m`-th moments in direction `v` differ by
`>= 2 d^m kappa_core`. Every located anti-concentration theorem is proved by (i) decoupling / differencing, which
produces exactly product certificates (finite differences along parallelepipeds), or (ii) comparison with a continuous
measure, whose error is polynomial in `1/n` or `tau^{1/(8m)}`. Neither produces such a pair.

### 5.3 (d): slab conditioning — the heuristic is right, the precision is not available

**Lemma 5.3 (the grid model; proved).** Let `Gamma_j = {-j, ..., j}^r` and `z in R^r`. Then
`dist_inf((z . sigma)^m|_{Gamma_j}, Poly_{<m}|_{Gamma_j}) = 0` iff `|supp z| <= floor((m-1)/(2j))`, and if at least
`s := ceil(m/(2j))` coordinates of `z` have `|z_p| >= eta`, then
`dist_inf((z . sigma)^m, Poly_{<m}|_{Gamma_j}) >= (m!/2^m) eta^m`.

*Proof.* If `|supp z| <= (m-1)/(2j)`, `(z . sigma)^m` is a function of `sigma_{supp z}`, and every function on
`{-j..j}^{supp z}` is a polynomial of degree `<= 2j |supp z| < m` (tensor Lagrange). Otherwise choose `s` coordinates
with `|z_p| >= eta` and exponents `alpha_p <= 2j`, `sum alpha_p = m`; the lattice zonotope
`-j sum_p e_p + sum_p [0, 2j] e_p` lies in `Gamma_j`. The mixed difference `c = prod_p Delta_{e_p}^{alpha_p}` based
at its corner annihilates `Poly_{<m}`, has `|c|_1 = 2^m`, and `<c, (z . sigma)^m> = m! prod_p z_p^{alpha_p}` (the identity
of nearstrata Lemma 2.3). Duality (nearstrata Prop 2.1) gives `dist >= m! prod |z_p|^{alpha_p}/2^m >= (m!/2^m) eta^m`.
If `|supp z| >= s`, take `eta = min_{p in supp z} |z_p| > 0`. QED

*Heuristic content of (d).* If `k'` were Gaussian, `E[k'^alpha | s(k') = sigma]` would be a polynomial of degree
`<= |alpha|` in `sigma`, so a pair of measures on the grid `s(A)` with equal moments `< m` would lift to a pair on `A`
with equal moments `< m`, and `kappa_A(u) >~ kappa_{Gamma_j}(z)` with `z` the coordinates of `u` in the column basis.
Lemma 5.3 then says: the exceptional set is exactly the spans of `<= floor((m-1)/(2j)) = 4` columns (at `w >= 10`),
i.e. **the cheap spans**, with exponent `m` — which is precisely the conjectured shape of MIN^pair. This is a
derivation of the conjecture from a Gaussian proxy *with lattice-valued slab statistics*; it is not evidence about
regime boxes.

*Why it fails quantitatively (identified gap).*
1. **Precision.** The lift is exact only if conditional moments are exact polynomials. With approximate matching to
   precision `epsilon` (in the dual-sup sense: `|<c, p>| <= epsilon |c|_1 sup_A |p|` for `p in Poly_{<m}`), the bound
   survives iff `epsilon < kappa_core`, i.e. `ln(1/epsilon) >~ m ln(2de/m) + m ln(1/eta)`: `~ 1.9 * 10^4` nats at
   `n = 2900`, `~ 8 * 10^5` at `n = 10^6` (Table 6.1, column `Lc + Lpoly`). Local limit theorems and Edgeworth expansions
   of fixed order give precision `n^{-O(1)}`; an expansion of order `~ m` would be needed, and none is available.
2. **Degenerate fibres.** In the main regime `r >= d` (Hadamard-type profiles, `r = n - 1`), `s` is injective on
   `Cube` up to the kernel, so there is nothing to condition on: the fibres are points.
3. **Exact algebra off the column span.** If `r < d`, the kernel of the profile supplies exact certificates: if
   `v_1, ..., v_{ceil(m/2)} in {-1,0,1}^d` have disjoint supports and `s_p(v_i) = 0` for all `p, i`, then
   `k_0 = -sum v_i` and the zonotope `k_0 + sum [0,2] v_i` lie in `A`, and nearstrata Lemma 2.3 gives
   `kappa_H(u) >= (m!/(2d)^m) prod |<v_i, u>|^{alpha_i}`. But at an actual cluster every such `v_i` is a relation:
   `sum_x v_{i,x} delta_x in (pi/2) Z` ([D] Lemma 3.2 proof: `sum v_x delta_x = <2 sum v_x a_x, phi> - (pi/2)<v, kappa>`),
   hence `= 0`; so `hat u` lies in the column span and these certificates vanish there. (This is [D] Prop 5.2's
   "low-rank subtori are invisible" in certificate form.)
4. **The one high-degree transfer known** (Filmus-Mossel, slice vs cube, harmonic multilinear, degree `o(n)` [L])
   handles one linear constraint, compares distributions (CDF-type closeness), and gives no `e^{-m ln n}` precision;
   the one-cut case it would cover is already solved exactly (nearstrata Thm 2.5).

### 5.4 The statement that would be needed

**Conjecture ILO_m (inverse Littlewood-Offord at degree `m`, prescribed leading form).** There are `c_0, K_0 > 0` such
that for every regime profile (`n >= n_1`) and every `u` with `|u|_inf = 1`: if some `Q in l_u^m + Poly_{<m}` has
`max_A |Q| <= d^m kappa`, then there is a rational subspace `L` with `rho_A(L) < m` (projection-cheap) and

```text
dist_inf(u, L)  <=  ( kappa / kappa_* )^{1/m},        ln(1/kappa_*) <= K_0 m ln n .
```

*What it gives (proved implication).* With `kappa = eps_0` and the separation lemma: if every projection-cheap `L` has
`P(L) >= N_2/Q_L` with `m Q_L abar/2 + K_0 m ln n < ln(1/eps_0) - m ln 2`, then MIN^pair (`H = infinity` form) holds,
since `dist(hat u, L) >= (1/2) e^{-Q_L abar/2}` (nearstrata Lemma 3.3) forces `kappa(hat u) > eps_0`. For width-cheap
`L` (`Q_L <= 16(1 + O(w/n))`) and `ln(1/kappa_*) = m ln(2de/m)` (the value in the zonotope models) the inequality
holds for `n >~ 10^6` (Table 6.1, `room_16 > 0`). For projection-cheap `L`
the collision counts are open (round 13: 32-class spans `L_{2,1}` lie in `TC` for Sylvester `k >= 5`).

*Status.* Open. It is the "inverse theorem" form of MIN^pair; it is a statement about the box only, so by fake Thm 2.1
it fails on any (Mar) fake box with a nondegenerate real `TC_H` point satisfying (2.2). No located theorem implies
it at any `m >> log n`; Theorem 3.1 shows why the standard hypotheses cannot; Theorem 4.4 shows that its density
version is false at regime sizes unless the density is profile-specific.

## 6. Loss ledger against the room

**Table 6.1** (`ledger.py`; `C = 1`, `W = W_min = 3n ln n - 14n`, `m = ceil(4w)`; natural logs; exact integer
parameters, 50-digit arithmetic). `L0 = ln(1/eps_0)` is the Step-0 budget. `Lc = m ln(2de/m)` is the base cost of an
exponent-`m` certificate (`m!/(2d)^m`, as in the zonotope models). `Lpoly = m ln n` is the extra cost at distance
`1/n` from the strata. `L16 = m(8 abar + ln 2)` is the extra cost at the separation radius `(1/2)e^{-8 abar}` of a
16-class cheap span, with the worst exponent `m`.

| `n` | `m` | `L0` | `Lc` | `L0 - Lc - Lpoly` | `L0 - Lc - L16` |
| --- | --- | --- | --- | --- | --- |
| 2900 | 1864 | 7190 | 3979 | -11650 | -72051 |
| 4096 | 2256 | 11216 | 5165 | -12714 | -94380 |
| `10^4` | 3683 | 34077 | 9914 | -9759 | -179223 |
| `2^16` | 10231 | 315735 | 36323 | +165947 | -516341 |
| `10^5` | 12850 | 513467 | 48123 | +317404 | -599266 |
| `10^6` | 44112 | 6.86e6 | 2.12e5 | +6.04e6 | +1.78e6 |
| `10^9` | 1682268 | 1.20e10 | 1.36e7 | +1.20e10 | +1.17e10 |

**Losses of each step of §5, with their status.**

| step | loss (nats) | fits the room? | status |
| --- | --- | --- | --- |
| Step 0 (record) | budget `L0 = W/4 - O(ln n) ~ (3/4) n ln n` | — | proved |
| separation from a stratum with `P(L) = N_2/Q` | `m (Q abar/2 + ln 2)` at exponent `m`; `o_L (...)` at the true order `o_L` | `Q <= 16`: from `n ~ 10^6` at exponent `m`; from `n_B ~ 5 * 10^3` with proved orders (nearstrata Prop 5.3) | proved |
| fibre reduction (Lemma 5.1) | `m ln(1/eta)`, included above | same | proved |
| core certificate on a fibre (§5.2) | `ln(1/kappa_core)`; `m ln(2de/m)` in every model where it is known | yes if `ln(1/kappa_core) = O(m ln n)` and `n >~ 2 * 10^4` | **no method** |
| Gaussian / Carbery-Wright (proxy only) | `m ln(1/t) ` scale: needs `W/(4m) >> ln(Cm)` | from `n ~ 2 * 10^4` | proved for the proxy, irrelevant to `A` (§3, §4) |
| invariance transfer (MOO) | additive error `>= m d^{-1/(8m)} ~ m` | never | vacuous |
| MNV / polynomial LO | non-trivial only for `m <= c log r/loglog r` | never (`m >= 1864`) | vacuous |
| LLT/Edgeworth for slab conditioning (§5.3) | needs precision `e^{-(Lc + Lpoly)}` (`~e^{-1.9*10^4}` at `n = 2900`, `e^{-8*10^5}` at `10^6`) | fixed-order expansions give `n^{-O(1)}` | **gap** |
| `kappa_H` vs `kappa_V` (Remark 1.3) | certificate height `<= exp(wW/2 - O(m ln n))` | counting leaves codimension `E^{1-o(1)}` | **gap** |

**Reading.** The room `Theta(n ln n)` is enormous asymptotically: at polynomially separated `u` (`eta >= n^{-K}`)
with an exponent-`m` core bound, the budget allows `K <= (3/4)(n/m)(1 - o(1)) ~ 0.066 sqrt(n/ln 4r)`.
Nothing in the chain is limited by the room; the chain is limited by the absence of any core certificate at degree `m`.
At the record's thresholds (`n ~ 3 * 10^3 - 10^4`) even an ideal exponent-`m` Lojasiewicz inequality at distance
`1/n` does **not** fit (`L0 - Lc - Lpoly < 0`); the room opens between `10^4` and `2^16`.

## 7. Files

* Write-up: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/antitheory.md`
* Notes: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/antitheory_notes.md`
* Scripts (each writes `<name>_out.txt`), in
  `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/antitheory_checks/`:

| script | content | status |
| --- | --- | --- |
| `common.py` | regime parameters `w, j, m, (R1)`, `W_min`, `ln(1/eps_0)` | exact integers, 50-digit logs; `(R1)` exact big-integer binomials for `n <= 10^6` |
| `params.py` | Table §3.1 (CW proxy), MOO error, MNV range | 50-digit arithmetic; CW constant set to 1 (optimistic) |
| `barrier_product.py` | Theorem 3.1: exact central trinomial coefficients, `mu{Q != 0}`; check `p_b >= 0.4/sqrt b`, `2 <= b <= 3000` | exact integers; logs 40 digits |
| `act_strata.py` | Prop 4.2 (partition `q_max`), Prop 4.3 and Theorem 4.4 at `theta = 1/2` | rigorous sub-Gaussian tails, exact integer `c`; `Gamma` by 50-digit summation |
| `act_bits_theta.py` | Theorem 4.4 at `theta = 1/2, theta_0, theta_S`, `k = 12..60` (Gamma test `k <= 22`) | as above |
| `lemma41_check.py` | Lemma 4.1 on a toy (`d = 8`, `q = 2`, `c = 1`, `m = 7`): `deg(Q - l_u^m) < m` and the identity on `Y` | exact (sympy rationals) |
| `ledger.py` | Table 6.1 | 50-digit arithmetic |
