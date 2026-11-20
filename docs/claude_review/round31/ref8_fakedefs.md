# Referee report, block "fakedefs": does the VIII.3 object answer the record's question, and what does it cap?

Scope: Part VIII.3, (VIII.10) and VIII.4 of `astra_free_VIII/answer.md` (lines 1741-1946), read against the
record's definitions. The internal steps of VIII.3 (uncertainty estimate, character margins, residue witnesses,
moment lifts, tail count) are the subject of the blocks `ref8_fakecore.md` and `ref8_faketail.md`. Both report
every step correct, with one correction (fakecore C1: the weight coefficient at fixed `rho` is `4d/rho`). This
report takes those steps as given **only where it says so**. It checks the definitional fit itself.

Scripts (own code, exact integer arithmetic except where floats are marked) are in `round31/ref8_fakedefs/`:
`defs_check.py` -> `defs_check_out.txt`, `window_numbers.py` -> `window_numbers_out.txt`. I read Astra's two
scripts but did not import them. No literature search was needed or done.

## 0. Results

1. **The object is an angle model in the record's exact sense, and it lies inside the record's low-height window,
   for every fixed `0 < delta < 1/12`** (Proved, by checking each clause of overview Def. 5.2, round5
   `construct.md` Def. 1.2 and round21 `fakes.md` §1; the arithmetic is Computed). It has `M = m` points, genuine
   split primes, exponent cap two, normalised norm `N'` with `P_+(N') <= m^(23/12+o(1))`,
   `W = log N' <= (23/3 + o(1)) m log m` along a diagonal sequence, `omega(N') >= (q+1)/2 > M/4`, free dimension
   `m_free = (q-1) - (m-1) = (1 + o(1)) M`, and `Lambda = sum p^(-1/2) >> q^(xi/2) -> infinity`.
2. **It answers question (Q) of round 21 affirmatively for `0 < delta < 1/12`** (Proved, conditional on the
   correctness of VIII.3 as refereed in fakecore/faketail). The set of good angles has positive Haar measure on a
   component of a realisation coset. The arithmetic angles (`arg beta_j`, `Norm beta_j <= p_j`) form a finite set
   there. So non-arithmetic gap-admissible points exist on realisation cosets of barrier-scale window profiles with
   `m_free >= cM` and `Lambda -> infinity`. This is scenario P of `fakes.md` §6.3, for these profiles.
   For `delta >= 1/12` (that is, `P_+ <= M^(23/12)`) (Q) stays open.
3. **No refereed statement of round 21 is contradicted** (Proved, by line-by-line comparison, §3).
   * Theorem A bounds a Haar proportion. Astra's set has positive but tiny Haar measure, far below the bound
     `exp(-2.17 sum_J p_j^(-1/2))`. There is no conflict.
   * Corollaries A.1-A.2 concern a construction that **draws `phi` from Haar measure on the coset**. Astra draws
     from a different, absolutely continuous law: deterministic offsets `t`, a common shift `c_0` of size `~1/H`,
     an integer base-`B` digit series, and a final perturbation of size `B^(-J-2)`. Its density relative to Haar is
     astronomically large. Her union bound is taken under this law. **The hypothesis that fails is "picks `phi` from
     Haar measure on a realisation coset".** The hypotheses `m >= cM` and all primes `<= M^(2-delta)` both hold.
   * J1 concerns `r` primes in one window of multiplicative width `1/(4r)`. Astra's window has width
     `sigma ~ q^(1-d) log q`, which is `q^(2-d) log q` times wider. Its hypothesis fails. The equal weights are
     replaced by a random assignment of primes to columns plus martingale concentration of the pair margins,
     which needs only `sigma = o((q log q)^(-1/2))`.
   * Lemma T does not apply: the angles are non-arithmetic. Theorem C does not apply: there is no complete initial
     segment of split primes, since all primes lie in `[P, 2P]`. Theorem LS holds automatically.
   * Two **heuristics** of round 21 are refuted for `delta < 1/12`: "window fakes cannot come from a
     positive-probability random choice" (README §1.2, listed as heuristic), and fakes.md §6.3 "no
     positive-probability recipe produces them". Both hold only for Haar recipes. Prove's Heuristic 6.3
     ("window fakes plausibly exist") is confirmed in this range.
4. **Weight: 23/3 is a genuine but small improvement on the record's 8, and a large one on the record's 27 for
   point residues** (Proved from the refereed blocks, with the quantifiers of §4).
   * Row 1 (`A_0`): `8 -> 23/3` (factor 0.958). The 8 of `lw.md` Thm 6.1 also certifies Lang-Waldschmidt `LW_J(eps)`
     and Haar-**majority** genericity. The new object certifies neither.
   * Row 2 at `X = 2M` (`A_1 = A_{2M}`): `27 -> 23/3`.
   * At fixed `rho` the coefficient is `4d/rho`. (VIII.8) at `d = 1.99` and fixed `rho in (0.99, 1)` gives
     `7.96/rho in (7.96, 8.04)`, which beats 8 only for `rho > 0.995` (fakecore C1). Only the diagonal (VIII.10)
     reaches `23/3`.
   * The gap to the proved lower bound is `[3, 23/3]` for `A_{2M}` (Theorem 1.1 holds on `A_1`). It is
     `[4 - o(1), 23/3]` for models whose primes all exceed `2M`, by the transfer of `round8/stability.md`
     Thm 5.2 (Prop 7.1). In growth units: no argument from these inputs proves `M <= (6/23 - eta) log R/loglog R`,
     and Theorem 1.1 gives `2/3`.
5. **Capped and uncapped inputs inside the window** (§5). Capped for `0 < delta < 1/12`: row 1 (`A_0` with a.e.
   genericity), row 2 at `X = 2M`, row 3 (via Prop 7.1), row 10. Conditionally capped: row 7 robust inputs with
   `gamma > 1/(2 - xi)`, via `prove.md` Cor 5.4 applied to the spread version, under its Caveat 3. Uncapped in
   every window: point residues at `X = M^A` with `A > 1`; P3-char (row 4); P3-all*; the union large sieve at the
   genuine constant (cross-gap); square-root-quality analytic inputs; Haar-majority and full-torus Haar-a.e.
   inputs; LW rows 11-12 (not checked); the association of each prime norm with its genuine Gaussian argument.
   Uncapped for all inputs: windows with `delta >= 1/12`.
6. **Two record sentences are now false or outdated.** These are the expert-note sentence on the window and the
   overview bullet "every known fake puts its primes in one window `[P, P e^(1/(4r))]`" (§5.5, "Prime-weight
   concentration"). Proposed wording is in §6.

## 1. The record's definitions (verbatim content, with sources)

**Normalised cluster** (`two_thirds/proof.md` Lemma 1.1; `round8/stability.md` §1). Divide `z_1, ..., z_M` by
`G = gcd(z_x)`. Then `N' = N/|G|^2` is odd, every `p | N'` splits, and
`u_x = eps_x prod_p pi_p^(a_xp) conj(pi_p)^(e_p - a_xp)` with `min_x a_xp = 0` and `max_x a_xp = e_p`. The arc
constant does not increase: `C' = C|G|^(-1/2) <= C`.

**Angle model / fake circle** (overview Def. 5.2 = `round4/outside.md` §1; `fakes.md` §1.1). It consists of
* distinct primes `p_j = 1 mod 4`, exponents `e_j`, distinct rows `a_x in prod [0, e_j]`;
* `u_x = 2a_x - e`, `w(v) = sum |v_j| log p_j`, `W = sum e_j log p_j`, `Delta = C e^(-W/4)`;
* a real `phi` that **realises** the profile: `<u_x, phi> = theta_0 + delta_x + (pi/2) k_x`, with
  `delta_x in [0, Delta]` and `k_x in Z`;
* such that `phi` is **gap-admissible**: for every `v in Z^r \ {0}`,
  `dist(<v,phi>, (pi/2)Z) >= arcsin e^(-w(v)/2)` and `dist(<v,phi>, pi/4 + (pi/2)Z) >= arcsin(e^(-w(v)/2)/sqrt 2)`.

An input class *contains a fake* when an angle model of size `W = O(M log M)` satisfies all its inequalities.
"Barrier scale" means `W = O(M log M)`. A fake at `W = c M log M` shows that the class cannot prove
`M <= (2/c - eta) log R/loglog R`.

**Realisation coset, free dimension** (`fakes.md` §1.1). `L = span_Z{u_x - u_1}` (equivalently the row lattice of
`B = 2(a_x - a_1)`), `m_free = r - dim L_R`. `Real(delta)` is empty or a coset of `H = L^perp`, with identity
component `H_0` of dimension `m_free`.

**"Capped"** (overview §5.4). A barrier result shows that no argument using only these inputs proves
`M = o(log R/loglog R)`. Kind F means that the class contains fakes at `W <= cM log M` for a fixed `c`.
"Capped in the window" means the same with the fakes required to have `P_+(N') <= M^(2-delta)`.

**Classes.**
* `A_0` (row 1): exponent profile, pair chords, cut/Plotkin, all monomial gaps, rational characters,
  "Haar-generic angles". `construct.md` Cor 6.2 and `lw.md` Rem. 2 read the last item as "Haar-generic
  properties of the angles". The witness `lw.md` Thm 6.1 has Haar success probability `>= 1/2`, with `delta`
  uniform and `phi` Haar on the realisation torus.
* `A_1(Q_0)` (`construct.md` Def. 1.2, "point residues"): residue data `rho_x^(b) in (Z[i]/b)^*` for every odd
  prime power `b = l^a <= Q_0` with `l` not dividing `N`. They satisfy `rho conj(rho) = N mod b`, compatibility
  across levels, and axiom (R):
  `|sin(<v(c), phi> - psi_u(c))| >= D_u(c) e^(-w(v(c))/2)/sqrt 2`, with `psi_u = arg(u)/2 + (pi/4) c.k`, for
  every integral zero-sum `c` and every unit `u`. Here `D_u(c)` is the actual collision modulus. If `l <= Q_0`
  divides `N`, level-wise cofactor data are added. `A_1 = A_1(2M)`, written `A_{2M}` in the overview's row 2.
  `A_1^>(Q_0)` denotes the models whose primes all exceed `Q_0` (`stability.md` §7).
* Row 3: `A_X` plus near-extremal structure. This means sub-cluster gcds, smoothness, prime counts
  (`omega(N') > M/4`, `e log p < 4W'/M + 8 log^+ C`), Plotkin balance, collision and level stability, rank.
  `stability.md` Prop 7.1 proves that every model in `A_1^>(Q_0)`, `Q_0 >= 5M/4`, satisfies all of them.
* Prime-level residues (row 4, P3-char; P3-all*, overview §6.3): residues that factor through prime residues
  `r_j(l)` with `Norm r_j(l) = p_j` and the forced parity law, imposed on characters (P3-char) or on every
  monomial (P3-all*).

**The low-height window** (`round21/README.md` notation): `Y = P_+(N') <= M^(2-delta)` for a fixed `delta > 0`.
Here `M` is the number of points and `N'` is the normalised norm.

**Question (Q)** (`fakes.md` §6.3; README §1.2, §3). For barrier-scale profiles with all primes `<= M^(2-delta)`
(so `Lambda > M^(delta/2)/4`) and valid character margins: is `G(p) ∩ Real(delta)` non-empty, and if so, is it
non-empty only at arithmetic points? The judge's form (§7) adds `m_free >= cM`, `Lambda -> infinity`.
* "Yes" caps only `A_0` and the row-3 profile inputs (referee correction).
* "Arithmetic" means `phi_j = arg beta_j` with `Norm beta_j <= p_j`.

## 2. Astra's object, clause by clause

Notation as in VIII.3:
* `q = 3 mod 4` prime, `K_ij = chi(j - i)`, `A = {i : K_i0 = +1}` (size `(q-1)/2`);
* `S ⊂ A` random of density `rho` (so `m = |S| = (rho/2 + o(1)) q`);
* columns `I = F_q \ {0}`, `a = 1 + K` on `I`, `e_j = 2`;
* `q - 1` split primes in `[P_0, P_0 e^sigma] ⊂ [P, 2P]`, `P = q^d`, randomly permuted on `I`;
* offsets `eps_i ≈ Delta(1/4 + beta_i/2)` ordered by auxiliary Gaussian prime directions `beta_i`, with
  `Delta = e^(-W/4)/10` in quarter turns.

| clause of the record | Astra's object | label |
|---|---|---|
| distinct primes `= 1 mod 4` | `q-1` distinct split primes in `[P,2P]` | Proved |
| distinct rows, exponents | rows `1 + K_x`, `x in S`. Distinct since pair `l^1` distance on `I` is `q+1` | Proved; Computed (`q = 43, 67, 103, 127, 163, 383`) |
| normalised (`min a = 0`, `max a = e`) | columns lacking an endpoint are truncated or deleted at the end of VIII.3. `W` decreases, a common phase is added, all gaps survive by zero extension, `s_0` is changed | Proved (re-derived). Computed: no column lacks an endpoint in any computed case |
| arc condition, `C` | quarter-turn width `Delta`, physical width `(pi/2)Delta = (pi/20)e^(-W/4)`, so `C = pi/20 <= 1`. After normalisation `W' <= W` and the condition only weakens | Proved |
| realisation, windings | every term of `theta` except `t` has **exactly** constant `K_S`-image (as real vectors), and `2K_S t = eps + const`. So `k_x = 0` and `theta_0` is free. No common phase is set to zero | Proved (re-derived from `KK^T = qI - J`) |
| both sharp gaps for **every** `v != 0` | VIII.3 completion + tail | taken from fakecore/faketail (Proved there) |
| rational characters | every zero-sum rational `lambda` on `S` with `K_S^T lambda` integral is integral: `q lambda = K v_ext` and the Vandermonde step mod `q`, since `|S| <= (q-1)/2` | Proved (re-derived). Computed: rank of `K_S^T` mod `q` equals `m` (`q <= 383`); Smith invariants of the row differences all 1 at `q = 43` |
| barrier scale | `W = (2d + o(1)) q log q = (4d/rho + o(1)) m log m` | Proved; see §4 |
| window `P_+(N') <= M^(2-delta)` | `P_+ <= 2q^d = m^(d+o(1))`, `d = 2 - xi`, any `xi < 1/12`. Inside the window for every fixed `delta < xi` | Proved. Computed: the convergence is slow. At `d = 1.99`, `log P_+/log m = 2.0673` at `q = 10^12` and `1.9945` at `q = 10^200` |
| row-3 prime count `omega(N') > M/4` | `omega(N') >= (q+1)/2 ≈ m/rho` | Proved |
| row-3 prime powers `e log p < 4W'/M + 8 log^+ C` | `2 log(2P)` against `(16d/rho + o(1)) log q` | Proved |
| (R) at all odd prime powers `<= 2m`, all four units | point residues `s_0 omega_i/conj(omega_i)`, auxiliary split primes `2m < l_i <= q^(1.01)`. `D_u <= sqrt2 L_0^(s/2)`, `D_1 <= L_0 abs(beta_i - beta_j)`. All primes of `N` exceed `2m`, so no cofactor data is needed | definitional fit Proved (re-derived: her inequality is (R) with `v(c) = K_S^T c`, `<v(c),phi> = (pi/4) c.eps`, `k = 0`, `psi_u = arg(u)/2`); estimates from fakecore |
| "Haar-generic" | positive Haar measure on a component of a realisation coset of `H_0` (dimension `q - m`) | Proved (faketail 5b). **Positive, not majority, measure** (§4) |
| non-arithmetic | the arithmetic angle vectors are finitely many. A positive-measure set minus a finite set is non-empty | Proved |
| (Q) hypotheses `m_free >= cM`, `Lambda -> infinity` | `m_free = q - m = (2/rho - 1 + o(1)) m`; `Lambda >= ((q+1)/2)(2P)^(-1/2) >> q^(1-d/2) -> infinity` | Proved; Computed (`m_free/m = 1.01` at `rho = 0.995`) |

Small slip (no effect). In "Forced characters" the support bound "at least `q/2` nonzero entries" should be
`(q-1)/2` on `I`. A nonzero polynomial of degree `(q-1)/2` has at least `(q+1)/2` non-roots on `F_q`, and column 0
is always a root. With `(q-1)/2` the bound for `s = 4` reads `h >= 1.2q - 0.4`, and (VIII.9) is unaffected
asymptotically. Computed for `s = 4`: see the table in §7.

## 3. Apparent conflicts with round 21, resolved

**Theorem A and Corollaries A.1-A.2.** Theorem A (`fakes.md` §3.1) states
`nu(G ∩ (g + H_0)) <= prod_J mu_1(G_1(p_j)) <= exp(-kappa sum_J p_j^(-1/2))`. With `|J| = m_free ≈ q/2` and all
primes `<= 2q^d`, the right side is `<= exp(-0.77 q^(1-d/2))`. This is positive. Astra's support has far smaller
Haar measure. All its coordinates lie within `O(1/H)` of `t + c_0`, so on the `(q-m)`-dimensional component its
measure is at most about `(2/H)^(q-m) = exp(-Omega(q log q))`. So both statements hold together (Proved; the
estimate `(2/H)^(q-m)` uses that a coordinate projection of `H_0` onto a basis of `R^r/L_R` is a finite
covering). The corollaries quantify over constructions that "pick `phi` from Haar measure on a realisation coset
(for any rule choosing `delta`, `k` and the coset)" (Cor A.1). Astra's law is the push-forward of
`(c_0, X_1, ..., X_J, zeta)`, uniform on a box, under an affine map into a tiny region of the coset, with
integer digits `D_s` drawn independently. It is absolutely continuous, but it is not Haar. So:
* Cor A.2's hypotheses (`m_free >= cM`, primes `<= M^(2-delta)`) are **satisfied**;
* its conclusion (Haar proportion `<= exp(-kappa c M^(delta/2))`) is **true** for her profile;
* her success probability `1 - o(1)` is a probability under her own law, to which Cor A.2 does not speak.
The round-21 referee had already noted that Cor A.2 "is not a proof that no Haar point succeeds". The design
point is that a union bound under a concentrated law sees densities `O(H)`, `O(B^k)` and `O(B^(J+2))` instead of
`1`, and pays for them with the integer structure of the digits. Neither round-21 write-up considered this.

**J1** (README §4). If `r` split primes lie in `[a, a e^(1/(4r))]` then `a >= 16 r(r-1) e^(-1/(4r))`, which is
about `M^2`. Astra's primes lie in `[P_0, P_0 e^sigma]` with `sigma ≍ q log q/P`. Computed with actual primes
(`window_numbers_out.txt`):

| `q` | `d` | shortest `sigma` for `q-1` split primes in `[P,2P]` | `4r sigma` | `sigma (q log q)^(1/2)` |
|---|---|---|---|---|
| 383 | 1.99 | 3.38e-2 | 51.6 | 1.61 |
| 1019 | 1.99 | 1.50e-2 | 61.2 | 1.26 |
| 1019 | 23/12 | 2.41e-2 | 98.1 | 2.02 |

So the window is 50-100 times wider than J1's at these sizes, and the ratio grows like `q^(2-d) log q`. J1 does
not apply, and no J1 window exists at these `P` (`p_min < 16r(r-1)`). The pair margins survive anyway. At
`q = 383`, with `S = A` and five random permutations, the minimum over all 18,145 pairs of `w(K_i - K_i') - W/2` is
`24.28` to `24.31`, against `2 log p_min = 24.98`. The worst case for an arbitrary assignment would allow
`2 log p_min - (q-1) sigma = 12.08` (Computed, floating point on logarithms). The condition
`sigma (q log q)^(1/2) -> 0` holds only asymptotically (`q^(1.5-d)(log q)^(1.5)`). At the computed sizes the
product is still above 1, but the observed loss is only 0.7. So J1 is correct, and it is not an obstruction:
its hypothesis (equal weights to `O(1/r)`) is a feature of the old designs. `fakes.md` §2 predicted this
heuristically via grouping. Astra achieves it by randomisation instead.

**Lemma T, Corollary T.** These apply only to `phi_j = arg beta_j`, `Norm beta_j <= p_j`. Astra's good set has
positive measure, so it contains non-arithmetic points. If her construction produced arithmetic points it would
give actual clusters and disprove Cilleruelo-Granville. It does not, and she does not claim it. No conflict.

**Theorem C, Proposition B.** They need every split prime below `p` to divide `N` (in particular `5 | N`). All of
Astra's primes are `>= P`. No conflict.

**J2.** "Neither answer to (Q) gives growth." Unchanged. Astra's "yes" caps a class and proves nothing about
actual clusters.

**J3, the pinned case `m_free = 0`.** It remains untouched by measure arguments. For capping this does not
matter: a cap needs one sequence of fakes, and Astra's has `m_free ≈ M`.

## 4. Weight: 8 versus 23/3, and the gap to 3

**Quantifiers.**
* For fixed `rho in (0,1)` and `xi in (0, 1/12)`, and every `q >= q_0(rho, xi, ...)` in the sequence, the
  model has `W/(m log m) <= 4(2 - xi)/rho + o(1)` and `log P_+/log m <= 2 - xi + o(1)`.
* A diagonal sequence `rho_k -> 1`, `xi_k -> 1/12` gives `limsup W/(m log m) <= 23/3` and
  `limsup log P_+/log m <= 23/12` (VIII.10).
* For a fixed window `delta < 1/12`, every member with `xi_k > delta` and `q` large lies in it. So **the window
  `delta` is capped at every coefficient `> 23/3`**, independently of `delta`.
* The thresholds are not uniform in `rho` or `xi`. That is harmless for a cap, which is an existence statement
  for unbounded `M`.
* `C = pi/20` is fixed. For a smaller fixed `C`, shrinking `Delta` by `C` leaves every margin intact, because the
  pair needs only `L_0/P -> 0` and the characters have slack `P^(-0.03q)`. *This is my remark, Heuristic; Astra
  states only constant one.*
* Convergence is slow. Computed at `rho = 0.999`, `d = 23/12`: `W/(m log m) = 8.08` at `q = 10^6`, `7.75` at
  `10^30`, `7.70` at `10^100`. At `d = 1.99`, `rho = 0.995`: `8.43`, `8.08`, `8.02`. Every statement is
  asymptotic.

**Same class?**
* Row 1 (`A_0`). The record's witness (`lw.md` Thm 6.1) has `W <= (8 + o(1)) M log M`, primes `>= (40M)^4`
  (outside every window), Haar success probability `>= 1/2`, and also `LW_J(eps)`. Astra's object satisfies
  every clause of row 1 at `23/3`, provided "Haar-generic" means "outside a Haar-null subset of the realisation
  component". It does **not** have Haar-majority genericity. Its good set has measure `<= exp(-c q^(1-d/2))`
  (Theorem A), and its law is concentrated near one direction (first construction) or near a rational grid
  (spread construction). An input that holds for Haar-most but not Haar-almost-all points is therefore not
  certified. Neither is `LW_J(eps)` (rows 11-12), which Astra did not test. The single-form criterion of row 11
  needs `w(supp v) >= W/4 + x_0` for every character. Pairs give `w(supp v) - W/4 = 2 log P - o(1)` (Proved: the
  pair support on `I` is `(q+3)/2` and `sum(2·1_supp - 1) = 4`). For other characters only `supp >= (q-1)/2`
  is proved, and with unequal weights the margin can then be negative of order `q sigma`. The computed minimum
  support for `s = 4` exceeds `(q-1)/2` (§7), so this is plausible but **Gap** (not proved).
  **Verdict: within row 1 as worded, a genuine improvement `8 -> 23/3`. Within row 1 plus Haar-majority or LW,
  no improvement is shown.**
* Row 2 at `X = 2M`. The refereed witness was `construct.md` Thm 2 at `(27 + o(1))`. `lower.md` gives 48.
  `lw.md` Remark 1's 12 is unrefereed. Astra's object lies in `A_1(2m)` with all primes `> 2m`, so the improvement
  is `27 -> 23/3`, a factor 3.5.

**The gap.**
* Theorem 1.1 uses only `A_1` inputs (`construct.md` §1.2, Cor 6.1), so `W >= 3 M log M - O(M)` on `A_{2M}`.
  For `A_{2M}` the bracket is now `[3, 23/3]` (previously `[3, 27]`). In growth constants: the best constant
  provable from `A_{2M}` lies in `[6/23, 2/3]` (previously `[2/27, 2/3]`).
* For models whose primes all exceed `5M/4` (Astra's), `stability.md` Prop 7.1 transfers Thm 5.2, which then reads
  `W >= (4 - o(1)) M log M`. So within that subclass the bracket is `[4, 23/3]`.
* The architecture has `W ≈ 2·r·log P` with `r ≈ 2M` columns at exponent 2, that is `4 (log P/log M) M log M`.
  It would reach the subclass floor 4 only if `P = M^(1+o(1))`. The method stops at `d > 23/12`, because the
  digit sup norm `E = q^(1.5+a)` must fit below `B < P` (faketail item 6). *Heuristic remark.*
* For `A_0` alone the record proves no `M log M` lower bound that I found. Every lower bound on record uses
  residue pigeonhole.

## 5. What is capped in the window, and what is not

Every statement below assumes VIII.3 and (VIII.10) as refereed in fakecore/faketail. "Window" means
`P_+(N') <= M^(2-delta)` with `0 < delta < 1/12` fixed. "Capped" means kind F at every coefficient `> 23/3`.

**Capped (Proved, conditional only on those blocks):**
* row 1, `A_0`: exponent profile, pair chords, cut/Plotkin and cut-identity consequences, every monomial gap
  (both sharp forms), all rational characters (all integral here), Haar-generic properties in the almost-everywhere
  sense on the realisation component;
* row 2 at `X = 2M` (`A_1 = A_{2M}`): point residues with axiom (R) at every odd prime power `<= 2M`, all four
  units, actual collision moduli, and the round-7 consequences valid on `A_X`;
* row 3: all profile-structure inputs, by `stability.md` Prop 7.1, since the object lies in `A_1^>(2m)` and
  `2m >= 5m/4`;
* row 10 (Baker/Matveev): every character is supported on `>= (q-1)/2 ≈ M` primes;
* row 5 to the extent that it reduces to row 1 (`outside.md` A.1). *I did not re-read A.1.*

**Conditionally capped (Gap: conditional).** Row 7 robust prime-supported inputs with `gamma > 1/(1+eta)`, primes
`>= M^(1+eta)`, `eta < 1 - xi`, via `prove.md` Cor 5.4 applied to the VIII.4 spread models. Astra did not draw this
consequence. Two parts are clear:
* The assigned angles are `L`-spread in the sense of `prove.md` Def. 5.1. The spacing `>= 1/(3q)` gives at most
  `3q eta + 1 <= 8 s eta + 4 (Ls)^(1/2)` points per arc (with `s >= (q+1)/2`), and separation up to conjugation
  `>= 1/(6q) >= M^(-2)` (Proved, given VIII.4).
* The **genuine** angles of the `q-1` primes in the annulus `[P_0, P_0 e^sigma]` must also be spread (Caveat 3).
  The lattice bound gives `O(q log q · eta + 1)` per arc, which is too weak by `log q`. An upper-bound sieve for
  Gaussian primes in polar boxes is needed. **Not checked.**
* Even then, square-root-quality inputs (`gamma = 1/2`) stay uncapped, since `1/(1+eta) > 1/(2-xi) > 1/2`.

**Uncapped (in every window):**
* point residues at moduli `X = M^A`, `A > 1`. Astra proves `X = 2m` only. Her `D_u` bound needs auxiliary primes
  `> X` with `log L_0/log P < 0.51`, so the design cannot pass `X ≈ P^(1/2) ≈ M^(0.96)` without change.
  For `A >= 2 - delta` the "all primes `> X`" design is impossible in the window;
* row 4, P3-char, and P3-all*: no factor-by-factor residue realisation (Astra says so);
* the union large sieve at the genuine constant, i.e. the max-scale cross-gap against external primes (row 7
  below `2D`);
* GRH-quality prime sums and `x^(1/2)` box counts;
* Haar-majority properties on the coset, and Haar-a.e. properties on the full torus `T^r`;
* LW rows 11-12 (not tested; §4);
* additive (archimedean) data, ZE, and the other §5.5 items, which are unchanged;
* the association of each prime norm with its genuine Gaussian argument. VIII.11, a proved planar area count for
  genuine clusters, excludes the **first** construction. It cannot exclude the spread one: there
  `(1 + sqrt P sigma)(1 + sqrt P eta) ≍ q log q >= m` (Proved, arithmetic).

**Uncapped for every input class:** windows with `delta >= 1/12` (`P_+ <= M^(23/12)`). Below them, genuine clusters
are bounded when `P_+ <= (2 - eta) M log M` (VIII.4, refereed in another block). So the part of the window that
is open for `A_0`/`A_{2M}` is `(2 - o(1)) M log M <= P_+(N') <= M^(23/12 + o(1))`.

## 6. Proposed wording

**Expert note §1.** Replace "... and the low-height window $P_+(N')\le M^{2-\delta}$ remain uncapped; here $N'$ is
the normalised norm and $P_+(N')$ its largest prime factor (`round21/README.md`, status **open**), ..." with:

> ... remain uncapped. The low-height window $P_+(N')\le M^{2-\delta}$ ($N'$ the normalised norm, $P_+$ its
> largest prime factor) is **capped for $0<\delta<1/12$, for pair, monomial-gap, rational-character and
> point-residue inputs only**. Legendre (Paley) profiles with exponent cap two and randomly assigned primes in
> $[P,2P]$, $P=M^{2-\xi}$, give angle models with every monomial gap, axiom (R) at every odd prime power $\le2M$,
> a positive-measure set of non-arithmetic angles on a realisation component, $P_+(N')\le M^{23/12+o(1)}$ and
> $W\le(23/3+o(1))M\log M$ (`round31` Part VIII.3, (VIII.10); refereed). This answers round 21's question (Q) for
> $\delta<1/12$ and lowers the fake weight from $8$ (classes $A_0$) and $27$ (point residues at $2M$) to $23/3$.
> It caps nothing at $\delta\ge1/12$. In every window prime-level residues (P3-char, P3-all$^*$), point residues
> at moduli $M^A$ with $A>1$, the genuine-constant union large sieve and square-root-quality analytic inputs stay
> uncapped. Actual clusters have $P_+(N')\ge(2-o(1))M\log M$ (Part VIII.1).

Also add to the fake-construction bullets after the Paley item:

> * Legendre profiles with exponent cap two and primes below $M^{2-\xi}$, $\xi<1/12$:
>   $W\le(23/3+o(1))M\log M$, inside the low-height window, with point residues to $2M$. They have positive, not
>   majority, Haar measure on the realisation coset and are not checked against Lang–Waldschmidt.

**Overview §5.5.** Replace the head of the "low-height window" item with:

> * **The low-height window** $P_+(N')\le M^{2-\delta}$.
>   * **Capped for $0<\delta<1/12$, for $A_0$ (rows 1, 3, 10) and $A_{2M}$ (row 2 at $X=2M$)**, at every
>     coefficient $>23/3$ (`round31` VIII.3 and (VIII.10); refereed in `ref8_fakecore.md`, `ref8_faketail.md`,
>     `ref8_fakedefs.md`). The witness is a Legendre matrix of order $q$, $m\approx\rho q/2$ rows positive on
>     column 0, $a=1+K$, $e=2$, $q-1$ split primes in $[P_0,P_0e^\sigma]\subset[q^d,2q^d]$ randomly permuted
>     ($\sigma\asymp q^{1-d}\log q$, far wider than the equal-weight window of J1). The offsets are ordered by
>     auxiliary Gaussian prime directions, and the angles are a base-$B$ digit series in integral moment lifts.
>   * Round-21 Theorem A and J1 are not contradicted: the angles are drawn from a non-Haar absolutely continuous
>     law of tiny Haar support, and the primes are not in one equal-weight window. This refutes, for
>     $\delta<1/12$, the round-21 heuristic that window fakes admit no positive-probability recipe.
>   * **Open:** windows with $\delta\ge1/12$, and inside every window: P3-char, P3-all$^*$, point residues at
>     $M^A$ with $A>1$, the genuine-constant union large sieve (cross-gap), GRH-quality inputs, Haar-majority and
>     LW properties, and the association of prime norms with their genuine Gaussian arguments. Row 7 robust
>     inputs with $\gamma>1/(1+\eta)$ are capped conditionally (`round21/prove.md` Cor 5.4 applied to the spread
>     version, under its Caveat 3).
>   * Actual clusters: $M\le1+\omega(N')+K_C$, hence $P_+(N')\ge(2-o(1))M\log M$ (VIII.1, VIII.4).

Keep the existing sub-bullets on large sieve, rigidity and "below the window". Replace the sub-bullet "Open (Q)"
by "(Q) is answered yes for $\delta<1/12$ (above) and is open for $\delta\ge1/12$". In the bullet
**Prime-weight concentration**, replace "every known fake puts its primes in one window $[P,Pe^{1/(4r)}]$" with
"the fakes before round 31 put their primes in one window $[P,Pe^{1/(4r)}]$; the round-31 window fakes do not
($\sigma\asymp q^{1-d}\log q$)". In §5.4, update the witness columns of rows 1-3 to cite (VIII.10) at
$(23/3+o(1))M\log M$, keeping `lw.md` Thm 6.1 at 8 as the witness for rows 11-12.

## 7. Computations

`defs_check_out.txt` (exact; `S` random of density 0.99 in `A`):

| `q` | `m` | `K^T=-K`, `KK^T=qI-J` | pair `l^1` on `I` | coefficient sum | pair support | rank `K_S` (mod `10^9+7`, so `<=` rank over Q `<= m`) / mod `q` | `m_free/m` | columns missing an endpoint | `s = 4`: min `h`, `h/q`, min support, `(q-1)/2` |
|---|---|---|---|---|---|---|---|---|---|
| 43 | 21 | yes | 44 | 2 | 23 | 21 / 21 (Smith invariants of the differences: twenty 1's) | 1.048 | 0 | 60, 1.395, 24, 21 |
| 67 | 33 | yes | 68 | 2 | 35 | 33 / 33 | 1.030 | 0 | 94, 1.403, 38, 33 |
| 103 | 51 | yes | 104 | 2 | 53 | 51 / 51 | 1.020 | 0 | 146, 1.417, 59, 51 |
| 127 | 63 | yes | 128 | 2 | 65 | 63 / 63 | 1.016 | 0 | 182, 1.433, 74, 63 |
| 163 | 81 | yes | 164 | 2 | 83 | 81 / 81 | 1.012 | 0 | not run |
| 383 | 187 | yes | 384 | 2 | 193 | 187 / 187 | 1.048 | 0 | not run |

At `q <= 163` density 0.99 selected all of `A`. Every computed `s = 4` minimum lies above `1.2q` and above
`(q-1)/2`. The minimum `s = 4` support exceeds `(q-1)/2` by about `0.09q` to `0.11q`, which is relevant to the
row-11 Gap of §4.

`window_numbers_out.txt` (actual primes from sympy; floating point only for logarithms): the tables of §3 and §4.
