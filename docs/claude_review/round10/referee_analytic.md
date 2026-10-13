# Referee report on round10/analytic.md

Adversarial referee, round 10. Scripts: `referee_analytic_checks/`. Each runs as `python3 <script>` from that
directory, and its output is next to it as `*_out.txt`. Labels: **exact** means integer arithmetic, or mpmath with a
certified margin far above the working precision; **illustration** means floating point that no conclusion relies on.

## 0. Verdict summary

**The result survives: barrier, no growth improvement.** I found no error that changes the status. Theorem B,
Proposition C, Theorem D, Proposition D.1 and Lemma D.2 are correct as stated, apart from the quantifier slips
listed below. The claimed fake-violated input (LCG) really is violated by the fakes:
* I built explicit realisations of the round-5 `construct.md` fake (Paley-I flipped, `q = 11`, `b = 8`, `r = 100`)
  and of the `lw.md` 6.1 skew double core (`q = 3, 7, 11`, `P = (40M)^4`).
* Every one of their 142 primes violates LCG with an `alpha` of norm `<= p_j`.
* The genuine angles of the same primes satisfy LCG (exact check).

Two claims in the write-up's **classification** of analytic inputs are wrong. Both overstate how complete the
dichotomy is; neither gives growth.

1. **(major) Fibre inputs are not all capped.** Proposition A's consequence, the §2 slogan and the first row of the
   §6 table say that Hecke characters on the norm fibre "add nothing", and that "the input that separates `theta` from
   `phi` must involve primes outside `N`". This is false.
   * For actual clusters the fibre sums are *arithmetic*: `N^(2k) S_k(N) = sum_(|z|^2=N) z^(4k)` is a rational
     integer, and the same holds for its residues mod `q`.
   * This uses only the primes of `N`, and it separates genuine angles from fakes. On the four explicit fakes the
     distance of `N^2 S_1(N)/4` to `Z` is between 0.095 and 0.417.
   * These are prime-level residue and integrality inputs on every monomial, i.e. the P3-all class / "the actual
     primes" of `outside.md` §2.4. Their cap is **open**.
   * Proposition A is correct only for fibre statements that are identities or inequalities valid for every real
     angle vector.
2. **(major) The sharp large sieve is misclassified.** Remark 3.4.2 and the §6 table list "the sharp `(X+K)` large
   sieve" as a pointwise prime separation that "collapses under GPA (Prop. D.1)". It does not.
   * The sharp large sieve controls the spacing of the points `4 theta(varpi)` only at scale `1/X`. So it can enforce
     only `|sin(phi_p - arg varpi)| >~ 1/max(p, Norm varpi)`, not PCG's `(p Norm varpi)^(-1/2)`.
   * Composite lattice directions `beta` with `p < Norm beta <= 4p` satisfy the former exactly (12,546 cases), yet
     violate LCG and, in 9,360 cases, PCG.
   * The optimal large-sieve constant does not change when such an angle is grafted in (illustration).
   * Proposition D.1 therefore does not apply. Whether the sharp large sieve is capped or fake-killing is open.

Minor slips:
* Theorem B needs `A > 4`, not `A >= 4`, for the `lw.md` fakes. The summary's "`W <= (2A+o(1)) M log M` for any
  `A > max(2, 1/gamma)`" is false for `2 < A <= 4`: only the `sharp.md` fakes, with constant `bA` (`b >= 33`), exist
  there.
* Theorem B's caps bar an argument only when `gamma` and the residue exponent `A` are fixed before the target constant.
  This is the round-9 R1 caveat, and it is not stated.
* Theorem B is proved only for prime-supported inputs. The standard all-ideal Hecke large sieve and Type I/II sums are
  not covered by Lemma 3.2's argument.
* Lemma 5.3 needs only `Norm alpha <= n`. In that form LCG is *equivalent* to "`phi_j` is the argument of a Gaussian
  integer of norm `<= p_j`".
* The GPA region is a strip of width `2 h^(-1/2)`, not a narrow sector at one scale. "Beyond GRH" should read "not
  known to follow from GRH".
* The Kubilius-type "never empty" is a heuristic. In the summary it sits among the proved items.

## 1. Proposition A (fibre inputs)

**Re-derivation.**
* (2.1) holds. By unique factorisation, the elements of norm `N'` are `i^u prod pi_j^(a_j) conj(pi_j)^(e_j - a_j)`, and
  `lambda_k(i) = 1`. Also `sum_(a=0)^e e^(i x (2a - e)) = sin((e+1)x)/sin x` with `x = 4k theta_j`.
* The Parseval statement of Remark A.1 holds, including the factor `4 r(N')`. The `r(N')/4` associate classes each
  contribute `16` to `|S_k|^2` on average over `k` when `theta` is generic.
* As literally stated, Proposition A is a tautology and is correct. Any inequality derived for actual clusters from
  identities and inequalities that hold for every real angle vector also holds for every angle model.

**The overstatement (major).** The write-up then says three things:
* "So Hecke characters restricted to the norm fibre are the GAP Fourier analysis, and add nothing."
* "The input that separates `theta` from `phi` must involve primes outside `N`: the ensemble."
* §6, table row 1: fibre inputs do not separate actual from fake, and are capped.

All three are false once one allows **arithmetic** statements about the fibre values, which are just as much
"Hecke characters restricted to the norm fibre":
* `N'^(2k) S_k(N') = sum_(|z|^2 = N') z^(4k)` is a rational integer, being a sum of Gaussian integers invariant under
  conjugation. With `e_j = 1` it equals `4 prod_j 2 Re(pi_j^(4k)) = 4 prod_j 2 p_j^(2k) cos(4k theta_j)`. Each
  factor is itself an integer.
* Its residues modulo every `q` are determined by the residues of the `pi_j`. They are prime-level residue data
  imposed on a monomial-type sum.

These inputs use only the primes of `N`, are true for actual circles, and are false for generic fakes:

| fake (explicit realisation, `lcg_fakes.py`) | `r` | `|N^2 S_1(N)/4|` | distance to `Z` | min per-prime distance of `2 p_j^2 cos 4phi_j` to `Z` |
|---|---|---|---|---|
| `construct.md` Paley-I flipped, `q = 11`, `b = 8`, primes near `10^7` | 100 | `~10^1401` | 0.095 | 0.005 |
| skew double core `q = 3`, `P = (40M)^4` | 6 | `~10^106` | 0.236 | 0.011 |
| skew double core `q = 7` | 14 | `~10^278` | 0.196 | 0.020 |
| skew double core `q = 11` | 22 | `~10^469` | 0.417 | 0.027 |

The products were computed with mpmath at `sum log10(2p^2) + 60` digits, from a realisation computed at that precision.

These are exactly the inputs that `outside.md` §2.4 calls "residue information for all monomials at all moduli = the
actual primes". At bounded moduli they are the P3-all class, whose cap is **open** (README of round 9, bracket
`3 M log M - O(M) <= W_min^all <= (66A + o(1)) M log^2 M/loglog M`). So the "fibre" bin of the dichotomy splits in two:
* real-analytic uses of fibre Hecke sums (sizes, sector counts, Riesz products), which are capped by Proposition A;
* arithmetic uses (integrality, residues), which are **not capped**, since they are P3-all.

The write-up's parenthesis "(The trivial identity `theta_j = theta(pi_j)` is such an input; it is all of the problem.)"
files that identity under the ensemble. It belongs to the fibre (`pi_j` divides `N`).

*Correction.* Restrict the slogan to fibre statements valid for all real angle vectors. Add that arithmetic fibre
inputs are the P3-all class, which is open. This weakens the "dichotomy" in §0 and §6. It does not create a growth
argument: P3-all is the already-known open class.

## 2. Theorem B (graft barrier)

**Re-derivation of Lemma 3.2.** All three parts are correct.
* (1) `|D_k| <= 2 sum_(varpi above S) |a_varpi| <= 2 (2 s_[1,X])^(1/2) ||a||_2`, so
  `sum_k |D_k|^2 <= 8 s (2K+1) ||a||^2`. Then Minkowski in `l^2`.
* (2) Each changed term is at most `2 ||f|| log X_1`, and there are two elements per prime.
* (3) Each changed element enters or leaves `B` at most once.

**Re-derivation of the three estimates.**
* (R1) `8 s(2K+1) <= 24 s K` for `K >= 1`. Then `t = 24 s/(c P_0^(1/2)) <= 48 M^(1 - A/2)/c` and
  `(1 + t^(1/2))^2 - 1 <= 3 t^(1/2)`.
* (R2) The ratio `4 s log X_1/(c X_1^gamma)` decreases for `X_1 >= e^(1/gamma)`, so it is at most
  `8 A M^(1 - A gamma) log M/c`.
* (R3) The ratio is `2 s_J/(c x_B^gamma) <= 4 M^(1 - A gamma)/c`.
* For lower bounds, `count_graft >= L - 2 s_J >= (1 - eta) L`.

All correct. Corollary B.1 is valid barrier logic: the fakes exist for infinitely many `M` (Paley orders), and any fixed
relaxation `eta > 0` exceeds `eta_M` for large `M`.

**Slip 1 (minor, constants): `A >= 4` must be `A > 4`.**
* `lw.md` Thm 6.1 requires `P >= P_1 = max{(40M/C_0)^4, (64M)^(2/eps)}`. Its proof also uses
  `sqrt(P') >= 1600 M^2`. With `C_0 <= 1`, `(40M/C_0)^4 >= 2.56*10^6 M^4 > M^4`, so `P = M^A` is admissible only for
  `A > 4` (and `M` large).
* So the summary's "`W <= (2A + o(1)) M log M` for any `A > max(2, 1/gamma)`" is false for `2 < A <= 4`. In that range
  only the `sharp.md` P3-char fakes are available. They need `Z = max(3M^A, 21 r^2)(log M)^2` and
  `r = b(M-1)`, `b >= 33`, and give `W <= (bA + o(1)) M log M`.
* Corrected:
  * `lw.md` fakes: `A > max(4, 1/gamma)`, `W <= (2A + o(1)) M log M`;
  * `sharp.md` fakes: `A > max(2, 1/gamma)`, `W <= (bA + o(1)) M log M`.
* Only the constant changes. For the GRH quality (`gamma = 1/2`) and for Harman-Lewis (`gamma ~ 0.6`) this gives
  `2A > 8`.

**Slip 2 (minor, quantifiers): the barrier needs `(c, gamma)` and the residue exponent fixed before the target
constant.**
* The fake constant is `2 max(4, 1/gamma)` (or `b max(2, 1/gamma, A_res)`). So Corollary B.1 bars an argument that,
  for each target `c_0`, uses robust inputs with fixed `gamma`. It does not bar one that chooses `gamma = gamma(c_0) -> 0`
  as `c_0 -> infinity`.
* For the published inputs that §3 lists, `gamma` is bounded below by an absolute constant (about `1/2` for polar-box
  sieve bounds, whose lattice-count error is the perimeter), so the conclusion stands for them.
* The P3-char clause ("residue data at all moduli `<= 3M^A` ... still allow `W = O(M log M)`") inherits exactly the
  round-9 caveat R1. All these fakes have `p_j > X = 3M^A`, and the constant `bA` grows with `A`. The clause therefore
  bars only arguments that fix `A` before `c_0`. The write-up should say so, as the round-9 README does.

**Slip 3 (minor-moderate, scope): only prime-supported inputs are covered.**
* Definition 3.3 and Lemma 3.2 concern sums and counts over Gaussian **primes**. §3 says the theorem covers "the
  large sieve for Hecke Grossencharacters in the form `<< (X + K^2)||a||^2`". That inequality is usually stated, and
  used in sieve arguments, for `a` supported on **all** ideals. Type I/II sums are also composite-supported.
* For such inputs a fake must be extended to composites. The natural extension is multiplicative: shift the angle of
  every multiple of `varpi_p` by `v(n)(phi_p - theta_p)`.
* Lemma 3.2's Cauchy-Schwarz bound then gives a change of about `(2K+1) 4 #{n <= X : S | n} ||a||^2 ~ K X s/P_0`.
  This is not `o(K X^(1/2))` once `X >> (P_0/s)^2`, so the proof does not extend.
* `composite_graft.py` (illustration, all Gaussian integers of norm `<= 3000`, graft at `p = 101, 197`):
  * the Lemma-3.2-type bound exceeds the whole constant `Lambda` (e.g. `56280` vs `10990` at `K = 100`);
  * the true optimal constant changes by at most `7.5%` (`K = 30`) and by `< 0.4%` for `K >= 100`.
* So the extension is plausibly true for `P_0 >= M^A`, but it needs a dual-large-sieve argument (points per arc of
  length `1/K`), not Lemma 3.2.
* Correction: "covers every published Hecke-type statement" should read "covers the prime-supported forms of these
  statements".

**Selection remark (3.4.1) and Remark 3.4.3 (residues).** Correct.

## 3. Proposition C (rotation invariance) and §4

**Re-derivation.**
* The Fourier expansion is valid: `<k, Cmat theta> = <Cmat^T k, theta>`, and the product structure of
  `Pi_1 x ... x Pi_r` factorises the sum.
* `Phi(g)` is invariant under column rotations, since `S_j(mu) -> e^(i mu s_j) S_j(mu)`.
* For every `s`, `sum_tau g(Cmat(theta(tau) + s)) <= Phi(g)`, and the term `tau = tau^0` with `s = -theta(tau^0)` equals
  `g(0) >= 1`.
* So Proposition C is correct, and as stated it is close to trivial.

I also checked an extension the write-up does not state, but which is needed for actual counting:
* *Distinct tuples or one shared set.* If the count is restricted to tuples of distinct primes, inclusion-exclusion
  produces diagonal terms like `S(mu_j + mu_l)`. Under a column-wise rotation each term keeps its modulus and gains a
  phase. So the sum of moduli still dominates every rotated count, and the conclusion `Phi >= 1` survives.

**Labels.**
* The "amplification" paragraph is a heuristic, and is labelled so. Its conclusion depends on how the constants
  compare (`(A-1) r/2` against the amplification exponent), and the write-up words it as "for `A` larger than an
  absolute constant".
* The Fejer-positivity inequality is correct: the Fejer kernel is `>= 0`, with value `K+1` at `0`.
* §4.2 (averages) is quoted or heuristic. The write-up says "I know of none". But the summary lists "Averaging
  statements of Kubilius type give exceptional sets, never empty ones" next to the proved items. It should be
  labelled heuristic there.

## 4. Theorem D (lattice cross-gap collapse)

**Lemma 5.2.** Correct: `Im(pi conj alpha)` is a nonzero integer for non-parallel `alpha`. The PCG part is also
correct: a prime parallel to `pi` is `+-pi` (both are primitive).

**Lemma 5.3. Correct, and it holds with `Norm alpha <= n` instead of `n + 1`.**
* Minkowski already gives `Norm alpha <= Q + 1/Q < n + 1`, hence `<= n`, and the contradiction uses only that
  `alpha`.
* The converse is the Lemma 5.2 computation with `beta` in place of `pi`: if `phi = arg beta` with `beta` primitive
  and `Norm beta <= n`, then every non-parallel `alpha` has `|sin| >= 1/(|alpha||beta|) >= (n Norm alpha)^(-1/2)`.
* So, exactly:

```text
LCG_n(phi; n)   <=>   phi = arg beta (mod 2 pi) for some primitive beta in Z[i] with Norm beta <= n.
```

The fake-violated input is therefore *literally* the statement "each angle is the argument of a Gaussian integer of
norm `<= p_j`". This sharpens the write-up's own conclusion ("it reconstitutes the actual problem exactly"):
* LCG is not merely equivalent in effect to the original problem;
* the class "angle models + LCG" is, by definition, the set of angle models that are shadows of Gaussian-integer
  tuples.

Two exact checks:
* `sharp_ls_check.py` Part 1 confirms the `Norm alpha <= n` form on 12,546 composite directions.
* `lcg_fakes.py` finds every fake violator with `Norm alpha <= p_j`.

**Theorem D.** Correct.
* `arg y_x = <2a_x - e, phi> - (pi/2)k_x = theta_0 + delta_x`.
* The `delta_x` are distinct, so the `y_x` are distinct.
* Arc length: `N''^(1/2) Delta <= C N''^(1/2) N^(-1/4) <= C N''^(1/4)`.
* Both directions of the infimum identity hold. Normalising an actual cluster (removing the Gaussian gcd, inert primes
  and `1+i`) lowers `N` and keeps the arc constant, because `L/|g| <= C sqrt(R)/|g| <= C sqrt(R/|g|)`. Actual
  normalised clusters are angle models (Lemma 0.1) and satisfy LCG (Lemma 5.2).

The remark that `outside.md` Lemma 0.1 at `v = e_j` is LCG against `alpha in {1, i, 1 +- i}` is correct.

## 5. Prop D.1 / Conjecture GPA / Lemma D.2, and the sharp large sieve

**Proposition D.1.** Correct as a conditional statement. GPA(`p_j`) and PCG together leave only lattice directions of
norm `<= p_j`. Then follow the proof of Theorem D.

**What GPA asks.**
* `|sin(phi - arg varpi)| < (h Norm varpi)^(-1/2)` is equivalent to `|varpi| |sin(phi - arg varpi)| < h^(-1/2)`, i.e.
  `dist(varpi, R e^(i phi)) < h^(-1/2)`.
* So GPA(`h`) says: *every line through `0` whose direction is not a lattice direction of norm `<= h` passes within
  distance `h^(-1/2)` of a Gaussian prime of norm `!= h`.* The region is a strip of width `2 h^(-1/2)` at all scales at
  once, not a sector at one scale.
* The comparison with GRH in item 1 is right:
  * the half-angle at norm `X` is `(hX)^(-1/2)`;
  * a smoothed sector with `sum |ghat| = O(1)` gives a GRH error `<< X^(1/2) log^2 X` against a main term
    `~ X^(1/2) h^(-1/2)/log X`.
* "GPA is beyond GRH" and "it is not implied by GRH" should read "not known to follow from GRH; the standard explicit
  formula argument falls short by `h^(1/2) (log X)^3` at every scale". Non-implication is not proved.
* Item 2 (lattice directions of norm `> h`) is right. The quadratics `f_m` with `m = +-1` also have no fixed prime
  divisor:
  * an odd `q` dividing every value would divide all coefficients, hence `q^2 | disc = -4`;
  * parity: if all values were even, then `(1+i)` divides both `alpha` and `beta`, so `Im(alpha conj beta)` is even.
* So the Bunyakovsky heuristic applies.

**Independent checks (exact; evidence only, since GPA is a conjecture).**
* `gpa_lemmaD2_check.py` (A): for `h = 29, 53, 197` (values the write-up did not test), every primitive direction with
  `h < Norm beta <= 20h` has a GPA witness. That is 266, 478 and 1784 directions, with largest witness-norm ratio
  77.1.
* (B): the Lemma D.2 hypothesis (primes on both lines `Im(w conj pi) = +-1` with `Re > 0`) holds for **all 8977 split
  `p <= 2*10^5`**, double the write-up's range. The longest walk along a line is 113 steps, the same as the write-up's
  worst case at `p = 91577`. Also checked:
  * every witness has `Norm != p`;
  * for `p < 2000`, an angle just inside the forbidden arc violates PCG.

**Lemma D.2.** Correct.
* `sin(arg w_+ - arg pi) = 1/(|w_+||pi|)` with positive cosine, so PCG excludes the open arc
  `(arg pi, arg pi + 2 rho_+)`, and symmetrically on the other side.
* The claim `Norm w_+- != p` is right: the `Im` values available at norm `p` are `0, +-p, +-(a^2-b^2), +-2ab`, never
  `+-1` for `p >= 5`.

**The sharp large sieve is not PCG (major misclassification).**
* Remark 3.4.2 says that `(R1)` with `Lambda = c(X + K)`, `K >> X^(1/2)`, "is the dual large sieve for points with
  spacing `>> 1/X`", and that "these pointwise statements ... collapse the fakes back to actual clusters". The §6
  table then lists "the sharp `(X + K)` large sieve" with PCG, as "yes under GPA (Prop. D.1)".
* The sharp large sieve only sees spacing. Take `a` on two points at `psi`-distance `t`. The left side is
  `2(2K+1) + 2 D_K(t)`, close to `4(2K+1)` for `K << 1/t`. With the Montgomery-Vaughan constant `2K + 1 + (pi/2) X`,
  this is violated for `X << K << 1/t`.
* So the sharp large sieve at all heights is equivalent, up to constants, to the **max-scale cross-gap**

```text
WCG_c :  |sin(phi_p - arg varpi)| >= c / max(p, Norm varpi)     for all Gaussian primes varpi,
```

  together with the same gap among the grafted angles. This is much weaker than PCG's `(p Norm varpi)^(-1/2)`.
* Proof that WCG does not collapse: let `beta` be primitive with `Norm beta <= p/c^2`. If `Norm varpi <= p`, then
  `|beta||varpi| <= p/c`. If `Norm varpi >= p`, then `|beta| <= |varpi|/c`. Either way
  `|sin(arg beta - arg varpi)| >= 1/(|beta||varpi|) >= c/max(p, Norm varpi)`. So `phi = arg beta` satisfies WCG_c for
  every `beta` of norm in `(p, p/c^2]`. When `beta` is composite, no prime is parallel to it. Yet `phi` violates LCG
  (Lemma 5.3) and, by GPA item 2, PCG.
* Exact check (`sharp_ls_check.py` Part 1, `c = 1/2`, `p = 101, 1009, 10009`), over all 12,546 composite primitive
  `beta` with `p < Norm beta <= 4p`:
  * WCG_{1/2} holds against every Gaussian prime of norm `<= 10^5` (all associates);
  * LCG_p is violated in every case;
  * PCG_p is violated by a prime of norm `<= 10^5` in 9,360 cases.
* Illustration (Part 2: optimal constant `Lambda` as the largest Gram eigenvalue; 232 elements of `G` of norm
  `<= 1500`; graft at `p_0 = 1493`):

| `K` | genuine `Lambda/(2K+1)` | graft at `arg(65+12i)` (norm 4369, composite) | graft at distance `1/(40X)` from another prime |
|---|---|---|---|
| 1500 | 1.230 | 1.230 | 2.002 |
| 6000 | 1.049 | 1.049 | 1.974 |
| 24000 | 1.015 | 1.015 | 1.625 |
| 96000 | 1.002 | 1.002 | 1.018 |

What this shows:
* The sharp large sieve detects near-coincidences, but not a graft onto a lattice direction of norm in `(p, 4p]`, which
  violates LCG.
* So Proposition D.1 does **not** apply to the sharp large sieve.
* The fakes it allows include all shadows of Gaussian-integer tuples with `Norm beta_j <= 4 p_j`. Those cost only
  `O(M)` in `W`, so within lattice-direction fakes it is LCG up to `O(M)`.
* Possibly it also allows non-lattice angles. Admissibility needs, at every dyadic scale `X >= p`, no prime angle within
  `c/X`. The expected number of violations is `sum_X c/log X`, which diverges only like `log log`. So Haar-random fakes
  heuristically violate it, but adversarial Cantor-type angles are not excluded.

Correction: move "the sharp `(X+K)` large sieve" to its own row of the §6 table: "max-scale spacing; status open:
neither capped nor shown fake-killing; reduces to actual clusters only on lattice-direction fakes". There are two
cases, and neither gives a usable mechanism:
* If it kills every non-lattice fake, the surviving class is, up to `O(M)` in `W`, the set of shadows of
  Gaussian-integer tuples (Theorem D with `Norm beta_j <= 4 p_j`), i.e. the original problem.
* If it does not, the surviving non-lattice fakes bar it.

The "equivalent up to constants" claim has two halves:
* Sharp large sieve implies WCG: the two-point computation above.
* WCG implies the sharp large sieve (Montgomery-Vaughan): all points of norm `<= X`, grafted or genuine, are then
  `>= c/X`-separated, so `Lambda <= 2K + 1 + O(X/c)`.

## 6. Fake-violated input claim: tested on round-5 fakes and skew double core

The write-up claims that LCG is true for actual primes and false for every angle model that is not the shadow of an
actual cluster of no larger `W`. By §4 this is correct, and it is even an exact characterisation. I tested it on
explicit fakes.

`lcg_fakes.py` builds explicit realisations. In each, `delta` is a random rational in `[0, Delta]^M` and `phi` is an
exact particular solution of (1.1) plus a random kernel element (Haar-generic), computed with mpmath at
`sum log10(2 p_j^2) + 60` digits. The realisation residual is `<= 10^-167`. The primes are chosen as in the source
theorems:
* `construct.md`: window of log-width `<= 1/(4r)` near `10^7`, `p >= 21 r^2`;
* `lw.md` 6.1: window of log-width `1/(8q)` above `(40M)^4`, `C_0 = 1`, `k = e_1`.

**Evidence that each sample is a fake.** The Lemma 0.1 gaps hold with these minimum ratios of distance to required
distance:

| sample | `M` | `r` | `W/(M log M)` | singles | pairs `e_j +- e_l` | characters `||c||_1 <= 4` | 5000-20000 random monomials |
|---|---|---|---|---|---|---|---|
| `construct.md`, `q = 11`, `b = 8` | 12 | 100 | 54.1 | 26.7 | 31 | 1.2e9 (4554 chars) | 26.7 |
| skew double core `q = 3` | 4 | 6 | 22.0 | 1.9e3 | 2.5e7 | 772 (54) | 772 |
| skew double core `q = 7` | 8 | 14 | 19.4 | 5.5e3 | 2.1e7 | 1.1e3 (868) | 5.5e3 |
| skew double core `q = 11` | 12 | 22 | 18.2 | 3.0e3 | 8.4e7 | 830 (4554) | 3.0e3 |

**LCG at the fake angles.** For each `j` I computed `rho_j = min over 0 < Norm alpha <= p_j` of
`|Im(alpha e^(-i phi_j))| sqrt(p_j)`. LCG fails iff `rho_j < 1`. The search scans the dominant coordinate up to
`sqrt(p_j)` with the nearest integer in the other, and the candidates are re-evaluated in mpmath.

| sample | primes violating LCG | `min rho_j` | `max rho_j` (worst witness) |
|---|---|---|---|
| `construct.md` | **100/100** | 0.0004 | 0.939 (`alpha = 60 - 41i`, norm 5281) |
| skew `q = 3` | **6/6** | 0.048 | 0.707 |
| skew `q = 7` | **14/14** | 0.016 | 0.851 |
| skew `q = 11` | **22/22** | 0.096 | 0.900 |

At the genuine angles `arg(pi_j)` of the same primes, the exact integer `|Im(alpha conj pi_j)|` over the same candidate
set has minimum exactly `1`, so LCG holds there. **So the claimed fake-violated input is confirmed: it holds for the
actual primes and fails at every prime of every tested fake, already with `Norm alpha <= p_j`.**

How to read this:
* LCG's fake-killing power is exactly its being the original problem (§4).
* By Theorem D, a fake of size `O(M log M)` satisfying LCG *is* an actual cluster of size `O(M log M)`. So LCG is not
  "provably false for every fake of size `O(M log M)`": that is equivalent to the open lower-bound question for actual
  circles.
* The write-up states this correctly ("fake-killing, but no reduction").

The skew double core of `exact.md` and `lw.md` 6.1 and the `construct.md` fake also violate the fibre-integrality
input of §1 above. That input is a second, different fake-violated input. It is not in the write-up and it uses no
external primes. It comes in two forms:
* *Per-prime form:* `p_j^(2k) e^(4ik phi_j) in Z[i]`, a Gaussian integer of norm `p_j^(4k)`. Already at `k = 1` this
  pins `phi_j` to finitely many values mod `2 pi`, so it is again the actual problem.
* *Product form for the fibre of `N`:* `prod_j 2 p_j^(2k) cos(4k phi_j) in Z` for all `k`. This is a P3-all-type
  constraint, and whether its residue versions at bounded moduli are capped is the open P3-all question.

## 7. Exact computations (referee_analytic_checks/)

| script | what | result |
|---|---|---|
| `lcg_fakes.py` | explicit `construct.md` fake (`q = 11`, `b = 8`) and skew-double-core fakes (`q = 3, 7, 11`, `lw.md` 6.1 primes): realisation, gap evidence, LCG at fake and genuine angles, fibre integrality | LCG violated at 142/142 fake primes (witness norm `<= p_j`); genuine angles satisfy LCG (exact); fibre sums non-integral at the fakes (§1, §6) |
| `sharp_ls_check.py` | Part 1 (exact): composite `beta`, `p < Norm beta <= 4p`, `p = 101, 1009, 10009`. Part 2 (illustration): optimal large-sieve constants | 12,546 directions: LCG violated in all, WCG_{1/2} holds against all primes of norm `<= 10^5`, PCG violated by such a prime in 9,360. Lattice graft leaves `Lambda` unchanged; near-coincidence graft doubles it for `K <~ 1/t` |
| `gpa_lemmaD2_check.py` | (A) GPA at lattice directions for `h = 29, 53, 197` (exact search, evidence); (B) Lemma D.2 lines for all split `p <= 2*10^5` (exact) | (A) 0 failures in 2528 directions, max witness ratio 77.1. (B) 0 failures in 8977 primes, max 113 steps |
| `composite_graft.py` | illustration: all-ideal large sieve under a multiplicative graft | Lemma-3.2-type bound useless (exceeds `Lambda`); true change `<= 7.5%` at `X = 3000`, `p = 101, 197` |

I did not rerun the write-up's `window_fluct.py` and `gpa_check.py`. Their numbers are labelled evidence, and no
proof uses them. I checked the derived quantities in §7.1 by hand: median `|S|/sqrt(n)` lies between 0.54 and 0.59,
and the Rayleigh median is `sqrt(ln 2) ~ 0.833`. The write-up's `lcg_check.py` and `lines_check.py` are superseded by
the independent checks above. They are consistent with them, including the extremal prime `p = 91577`.

## 8. Corrected statement

**Status: barrier. No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` is still the best
bound.

1. **Proposition A (proved).** A fibre statement (Hecke characters on the norm fibre of `N`: sector counts, `S_k(N)`,
   Riesz products) that is derived from identities and inequalities valid for every real angle vector holds for every
   angle model. Such statements are therefore capped by `outside.md` A.3 and `lw.md` 6.1.
   * *Not covered:* arithmetic fibre inputs, such as the integrality of `N^(2k) S_k(N) = sum_(|z|^2=N) z^(4k)` and
     its residues mod `q`.
   * These are true for actual circles and false for generic fakes. They are explicit on the four fakes of §6.
   * They belong to the P3-all class, which is open, and they need no external primes.
2. **Theorem B (proved, with the corrected ranges).** Fix `c > 0` and `gamma in (0,1]`. Every `(c, gamma)`-robust
   *prime-supported* ensemble input is satisfied up to `1 + o(1)` by the graft of:
   * the `lw.md` 6.1 fakes with `P = M^A`, `A > max(4, 1/gamma)`, `W <= (2A + o(1)) M log M`; or
   * the `sharp.md` P3-char fakes with `A > max(2, 1/gamma)`, `W <= (bA + o(1)) M log M`, `b >= 33`.

   Consequences and limits:
   * Arguments that use these inputs, together with the profile, every monomial gap and Haar-generic properties, with
     `gamma`, `c` and the residue exponent `A` all fixed before the target constant, cannot prove
     `M = o(log R/loglog R)`.
   * Published inputs have absolute `gamma`.
   * Composite-supported inputs (the all-ideal Hecke large sieve, Type I/II sums) need a multiplicative graft. That
     extension is plausible but not proved; Lemma 3.2's bound does not give it.
3. **Proposition C (proved).** The absolute-value Fourier count of prime tuples is invariant under column rotations,
   so it is `>= 1` and can never certify `T = 0`. This also holds with distinct-tuple inclusion-exclusion. The
   amplification and Kubilius statements are heuristics.
4. **Theorem D (proved, sharpened).** `LCG_n(phi; n)` holds iff `phi = arg beta` with `beta` primitive and
   `Norm beta <= n`.
   * Hence the class "angle models + LCG at every prime" has exactly the actual `W_min`.
   * The fake-violated input claimed by the write-up is verified: LCG holds at every actual prime (exact) and fails at
     every prime of explicit realisations of the `construct.md` fake and of the skew double core (`q = 3, 7, 11`;
     142/142 primes).
   * It is fake-killing only because it *is* the actual problem.
5. **Proposition D.1 and Lemma D.2 (proved, conditional).**
   * D.1 holds under GPA. GPA says that every line through `0` that is not a lattice direction of norm `<= h` passes
     within `h^(-1/2)` of a Gaussian prime of norm `!= h`. It is not known to follow from GRH, and at lattice
     directions it is Bunyakovsky-type.
   * D.2's hypothesis holds for every split `p <= 2*10^5`.
6. **Misclassification corrected.** The sharp `(X + K)` large sieve is equivalent, up to constants, to the max-scale
   cross-gap `|sin(phi_p - arg varpi)| >= c/max(p, Norm varpi)`.
   * That cross-gap is satisfied by lattice directions of norm up to `p/c^2`, which violate LCG and, typically, PCG.
   * So it is not covered by Proposition D.1, and its status (capped or fake-killing) is open.
   * Like LCG, it can at best reconstitute the actual problem, up to `O(M)` in `W`.

**Overall.** Narrow-sector distribution of Gaussian primes used through robust prime-supported statements cannot see
the `O(M)` selected primes (Theorem B, Proposition C). Used pointwise at PCG or LCG strength, it reconstitutes the
lattice (Theorem D, Proposition D.1). The dichotomy is not exhaustive, for two reasons:
* arithmetic fibre inputs (P3-all) are open;
* the sharp large sieve is unclassified.

Neither is a new growth mechanism: the first is the known open class, and the second reduces to actual clusters on its
lattice-direction fakes.
