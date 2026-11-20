# Round 31, Part VIII: judge report on low-height primes on integer circles

Source: Part VIII of `astra_free_VIII/answer.md` (heading "# Part VIII. Low-height primes on integer circles",
line 1667, to the end). Also read: Astra's summary `final8.txt`, the coordinator's question `reply8.txt`, and her
scripts `viii_window_check.py` and `viii_erasure_check.py`. Seven referees each took one block. **All seven reports
are present; no referee failed.**

| key | report | scope |
|---|---|---|
| rank | `ref8_rank.md` | VIII.1 (rank theorem), (VIII.4), (VIII.5), field test, real-data test |
| order | `ref8_order.md` | VIII.2, the height-order inequality (VIII.6) and its function-field failure |
| fakecore | `ref8_fakecore.md` | VIII.3 first half: Paley matrix, primes, uncertainty (mod-four) estimate, characters, offsets, residue witnesses |
| faketail | `ref8_faketail.md` | VIII.3 second half: moment lifts, atoms, completion, tail count, Haar measure, (VIII.10) |
| fakedefs | `ref8_fakedefs.md` | definitional fit of the VIII.3 object; what it caps; consistency with round 21 |
| spread | `ref8_spread.md` | VIII.4: area test (VIII.11), spread construction, the COMPUTED paragraphs |
| extend | `ref8_extend.md` | analyst: what VIII.1 buys, sharper constants, whether it can be pushed |

The judge's scripts are in `judge_VIII/`. They use the judge's own exact Gaussian arithmetic and exact rational
ranks, and import nothing from Astra or the referees:
* `viii1_test.py` (and `viii1_test2.py`, which adds a per-deficiency arc table): enumeration of all lattice points
  on every circle `x^2 + y^2 = N`, `N <= NMAX`, all windows of consecutive points with arc constant `<= CMAX`
  starting in the first quadrant, exact normalisation, exact checks of VIII.1. Outputs `viii1_1e6_C9.out`,
  `viii1_3e6_C4.out`, `viii1_1e6_C9_arc.out`.
* `tenpoint.py`: the 10-point record cluster and the constants `K_C`.
* `mod4_check.py`: the mod-four support bound (eigen-identity, exact minimum weights of the reduced eigenvector
  codes, the support transfer). Output `mod4_check.out`.

**Every claim labelled [J] is the judge's own and has not been refereed.**

Definitions used throughout. `N = R^2` is a positive integer and the circle is `|z|^2 = N`. A *cluster* is a set of
`M >= 2` distinct points of `Z[i]` on the circle in an arc of length `<= C sqrt R`. The *normalised* cluster is
obtained by dividing by the Gaussian gcd `g` of the points: `z_i' = z_i/g`, `N' = N/|g|^2 = prod p^{e_p}`,
`W = log N'`, `omega(N')` the number of distinct primes of `N'`, `P_+(N')` the largest of them. Then `N'` is odd and
every `p | N'` is split. Fix `pi_p | p` and put `a_ip = v_{pi_p}(z_i')` (the *exponent row* of mark `i`). The
*chord constant* is `C' = max |z_i' - z_j'| / N'^{1/4} <= C`. The *low-height window* is `P_+(N') <= M^{2-delta}`
for a fixed `delta > 0`. An *angle model* (fake circle) keeps primes and exponent profile but replaces the genuine
prime angles by real angles satisfying the realisation and gap conditions (expert note §1). `A_0` is the class of
pair, cut, monomial-gap and rational-character inputs with Haar-generic angle properties (overview row 1);
`A_1(Q_0)` adds point residues satisfying round 5's axiom (R) at all odd prime powers `<= Q_0`; `A_{2M} = A_1(2M)`.

## 0. Verdict in five lines

1. **Theorem (VIII.1), new and unconditional:** every normalised cluster satisfies `M <= 1 + omega(N') + K_C`, with
   `K_C = C^{4+o(1)}` explicit, and `M <= 1 + omega(N')` when the chord constant is below `sqrt 2` (referee
   sharpening). The judge rechecked the proof line by line and found 0 failures in 1.27 million exact window tests [J].
2. **Consequences for growth are thin:** `P_+(N') >= (2 - o(1)) M log M`; coefficient `1/2` instead of `2/3` when `R`
   is squarefull (every prime of `N` has exponent `>= 4`); **nothing improves for general `R`.** The rank bound is
   function-field valid and holds in `A_1` fakes, so it cannot give the uniform bound.
3. **(VIII.6) is a correct, uniform two-logarithm inequality of low novelty** (record Prop B2 plus Cramer). It bites
   only on finite-order coordinates, gives no growth bound, and has not been tested against any fake.
4. **The window fake is established:** for every fixed `0 < delta < 1/12` the window is capped for `A_0`, `A_{2M}` and
   the row-3 and row-10 inputs, at `W <= (23/3 + o(1)) M log M`. This answers round 21's question (Q) in that range
   and contradicts no round-21 theorem. The judge checked the mod-four support bound [J]. No step is false.
5. **For the integer problem the growth rate `(2/3 + eps) log R/loglog R` and the uniform bound are unchanged.**
   Part VIII's main value is to the map of barriers. Score on the goal: **3/10**.

## 1. Referee-corrected statements per block

Labels: **Proved** means a referee checked a complete proof line by line, or an exact computation was reproduced
with independent code. **Computed** means an exact finite computation, reproduced. **Heuristic**, **Gap**,
**False**, **Withdrawn** as usual. **No referee found an error in any mathematical step of Astra's.** One
quantitative sentence was called False by fakecore; I overrule that verdict in §1.3. Extend's FALSE and GAP verdicts
concern questions the referee posed, not claims of Astra's (§1.7).

### 1.1 rank: VIII.1 and its consequences (referee: CORRECT on every item, one CORRECT WITH CORRECTIONS)

**Proved.**
* **Theorem VIII.1 (corrected statement, quoted in full).** *Let `M >= 2` distinct points of `Z[i]` lie on
  `|z|^2 = N` in an arc of length `<= C sqrt R`, `R = sqrt N`. Normalise as above, and let `C' <= C` be the chord
  constant of the normalised cluster. Let `S` be a finite set of rational primes `l = 3 mod 4` with*
  `(1/2) log 2 + sum_{l in S} log l/(l+1) > log C'`.
  *Append to each exponent row `a_i` the one-hot indicator of the class of `z_i' mod l` in the norm fibre
  `{x in F_{l^2} : x^{l+1} = N' mod l}` (size `l+1`), for each `l in S`. Then the `M` enlarged rows are affinely
  independent. Consequently the exponent rows have affine rank `>= M - 1 - K_C` and `M <= 1 + omega(N') + K_C`, where
  `K_C = sum_{l in S} l`. In particular, if `C' < sqrt 2` then the exponent rows are affinely independent and
  `M <= 1 + omega(N')`.* The `(1/2) log 2` term is the referee correction (Correction 1): `N'` is odd, so
  `(1+i)` divides every `z_i' - z_j'` and no `g_ij`, hence `|b_ij| >= sqrt 2 prod_{colliding l} l`. Astra's
  (VIII.2) omits this term and states "for `C < 1` one may use the empty set"; that is correct but not sharp.
  Further refinement: a block `z' mod 2` (classes `1, i`, affine cost 1) gives `M <= omega(N') + 2` for
  `C' < 2^{3/4}`.
* *Hypothesis to state:* the marks must be distinct (`b_ij != 0`).
* *Constants (Computed, rank; reproduced by the judge, `tenpoint.py`):*

  | `C` | 1.3 | 1.7 | 2 | 3 | 4 | 5.9 | 8 | 8.2394 | 10 |
  |---|---|---|---|---|---|---|---|---|---|
  | `K_C`, Astra's (VIII.2) | 3 | 21 | 21 | 94 | 310 | 1301 | 3885 | 4168 | 8684 |
  | `K_C`, with `(1/2) log 2` | 0 | 3 | 10 | 40 | 94 | 381 | 1150 | 1150 | 2411 |

  Since `sum_{l <= x, l = 3 (4)} log l/(l+1) = (1/2) log x + O(1)`, `K_C = Theta(C^4/log C)`.
* **(VIII.4)** *`M <= pi(P_+(N'); 4, 1) + 1 + K_C`, hence `P_+(N') >= (2 - o_C(1)) M log M`, and for each fixed
  `eta > 0` the cluster size is bounded by a constant depending only on `(C, eta)` on clusters with
  `P_+(N') <= (2 - eta) M log M`.* This improves the record's constants `1/2` and `1` (round8 `stability.md` Thm 3.3(3),
  Prop 3.4). It does not reach the open window `P_+ <= M^{2-delta}`.
* **(VIII.5)** *If every prime of `N'` occurs to exponent `>= k` in the original `N`, then
  `M <= (2/k + o_{C,k}(1)) log R/loglog R`.* It beats `2/3` exactly for `k >= 4`. **Correction 2 (meaning):** if `R`
  is an integer, "every prime of `N` has exponent `>= 4`" (equivalently `>= 3`) means `R` is squarefull, and the
  coefficient is `1/2`; "`N` squarefull" is automatic for integer `R` and gives `1`, no information.
* **Field test.** VIII.1 transfers verbatim to `F_2[t]` and `F_3[t]` with the constant quadratic extension
  (odd-degree irreducibles are inert, fibres of size `p^d + 1`, inert mass diverges). Hence, by the record's filter
  F1, **it cannot by itself prove the uniform bound.** On VI.5 the affine rank is `M - 1`. The VII.1 orthogonality
  remark was not re-derived (not needed).

**Computed (rank).** 0 failures in about 1.42 million windows (exhaustive `N <= 10^7` at `C <= 3`, exhaustive
`N <= 10^6` at `C <= 9`, structured `N <= 10^13`). The 10-point record cluster `N = 1176852625` has affine rank 5,
raised to 6 by the mod-3 block and to 9 = `M - 1` by the mod-7 block. The golden and Pell families keep affine rank 4
and 3 while `omega` grows.

### 1.2 order: VIII.2, the height-order inequality (referee: CORRECT, two CORRECT WITH CORRECTIONS)

**Proved (referee-checked; analytic input quoted, see the caveat).**
* **(VIII.6), corrected statement.** *Let `B` be the `(M-1) x r` matrix with rows `2(a_i - a_1)`, `L = row_Z(B)`,
  `w(v) = sum_p |v_p| log p`, and for `v` in the rational row space let `D_L(v)` be the least `D >= 1` with
  `Dv in L`. There are an absolute effective `c_* > 0` and `M_0(C)` such that every normalised cluster with
  `M >= M_0(C)` and every nonzero finite-order `v` satisfy `(1 + w(v)) log(2 D_L(v)) >= c_* W`.* One may take
  `c_* = 1/(13(2K+3))`, `K` the constant of the two-logarithm Baker–Wüstholz bound.
* **Window corollary (corrected):** *for any fixed `A > 0` and any prime `p <= M^A` whose coordinate vector has
  finite order modulo `L`, `log(2 D_L(e_p)) >= (3c_*/A - o_C(1)) M`.* No window hypothesis is needed on the other
  primes, and Astra's `delta`-dependence of `c_delta` is spurious.
* The minor bound `log H <= (1/6 + o_C(1)) W` (VIII.7), the correction `t*` with `B t* = k` on all rows, `D | d`, the
  margin `1/12`, the height `w(v)/2` and non-torsion of `alpha` are all checked.
* **Function-field failure (Proved):** on VI.5, `w(v) = W/(M-1)` and `D_L(v) = M` (computed for `M = 4..32`), so
  (VIII.6) fails there. The principal-unit explanation is **Heuristic** and not needed.
* *Caveat on the source:* the Baker–Wüstholz statement was checked only through the record's refereed
  `round31/ref_baker.md` and search snippets; the Crelle paper was not read.

**Computed.** On the 10-point cluster: Smith invariants of `B` are `2,2,2,2,24`; coordinate orders `24,8,6,4,24`;
exact phase `t* = (7/24, 3/8, 1/6, 1/4, 11/24)`. **Correction:** this cluster lies outside the lemma's range:
`(s+1) H ||u||_inf = 65.2 > 1`. The exact phase exists because `B` has full column rank there, not by the lemma.

**Novelty: low.** (VIII.6) is `growth/auxiliary.md` Prop B2 (unrefereed; Matveev perfect-power certificates, which
already names "the missing small exponent of `Z^s/Lambda_B`") plus a Cramer bound and the round-31 two-log lemma.

### 1.3 fakecore: VIII.3 first half (referee: CORRECT on every step; one "False as stated", overruled)

**Proved.**
* Paley matrix `K_ij = chi(j - i)`: `K^T = -K`, `KK^T = qI - J`, row distance `q+1`, class `A` of size `(q-1)/2`,
  `a = 1 + K`, `e = 2`, column 0 deleted. Computed at `q = 43, 67, 127, 131, 383`.
* `q-1` split primes in `[P_0, P_0 e^sigma] subset [P, 2P]`, `P = q^d`, randomly permuted; every pair margin
  `w(K_i - K_i') - W/2 >= 2 log P - o(1)`; `W = (2d + o(1)) q log q`. The `o(1)` is about `-1` at computable `q`;
  harmless, since the pair instance of (R) needs only `2 log P - O(1)`.
* **The uncertainty (mod-four support) estimate, corrected statement.** *For primitive integral `v` with
  `sum v = 0` and `h = ||v||_1`, every odd prime `l <= h`, `l != q`, has at least `sqrt q + 1 - h` coordinates of
  `Kv` nonzero mod `l`, and `Kv mod 4` has at least `sqrt q + 1 - h` nonzero coordinates.* (Astra: `sqrt q - h - 1`.
  "Unit coordinate" must be read as "coordinate nonzero modulo the chosen prime".) This is the classical square-root
  bound for extended quadratic-residue codes in Paley-matrix language, valid in every residue characteristic
  other than `q`. Mod 2 alone fails (minimum 2).
* Random `S`: `gcd((Kv)_{S^c}) <= 2`, the norm bound, and the Vandermonde support bound mod `q`.
* Forced characters are integral; `||v||_1 >= 1.2q` for `s >= 4`; `(s+1)||v||_1 >= q||lambda||^2 + sq/2`.
* **Ordered residue witnesses:** Astra's inequality is exactly round 5's axiom (R) (`construct.md` Def 1.2) with zero
  windings, all four units and all odd prime powers `<= 2m`. `D_u <= sqrt2 L_0^{s/2}` and `D_1 <= L_0|beta_i - beta_j|`
  are correct; both sides of the pair instance scale with `|beta_i - beta_j|`, which is why the ordering is needed.
  Computed with exact Gaussian arithmetic at `m = 20, 40, 60`. The auxiliary primes exist only asymptotically.

**Ruling on fakecore C1 (and fakedefs' repetition of it).** Fakecore calls (VIII.8),
`W <= (7.96 + o(1)) m log m`, "False as stated", because at fixed `rho in (0.99, 1)` the coefficient is `7.96/rho`.
**I overrule this.** The displayed statement (VIII.8) does not fix `rho`. Astra's proof ends: "Taking fixed `rho`
tending to one and a diagonal sequence of sufficiently large `q` gives (VIII.8). No uniform starting threshold in
`rho` is being asserted." The diagonal argument is valid: for `rho_k -> 1` choose `q_k` so large that the
construction at `rho_k` exists with `W/(m log m) <= 7.96/rho_k + 1/k`; the `q_k` give the sequence. No uniformity in
`rho`, and no `rho_q -> 1` along one sequence, is needed. **Corrected wording (presentation only):** *for each fixed
`rho in (0.99,1)`, `W <= (7.96/rho + o_rho(1)) m log m`; along a diagonal sequence `rho_k -> 1`,
`W <= (7.96 + o(1)) m log m`.* The same applies to (VIII.10) (`23/(3 rho)` at fixed `rho`, `23/3` on the diagonal),
which Astra already states as a diagonal limit. Label: **Proved**.

### 1.4 faketail: VIII.3 second half and (VIII.10) (referee: CORRECT on every item)

**Proved.**
* Lifts `d^(s) = j^s mod q`, `s <= floor(q/4)`, with `K_S d^(s)` constant and sup norm `O(q^{3/2} sqrt(log q))`.
  Computed at `q = 103, 127`.
* `D_s = d^(s) + q K_{S^c}^T xi_s`: sup norm `<= E` with probability `1 - o(1)`; both atom bounds;
  `Pr(v.D_s = 0) = O(q^{-0.62} h^{-1/2})` (primitive balanced, `h <= q^{0.10}`) and `O(q^{-0.62})` (balanced,
  `h < q/2`). Constants depend on `rho`, roughly `(1-rho)^{-1/2}`.
* Constant `K_S`-images, the Newton step (some `v.D_s != 0` with `s <= h/(2 g_0) <= J`; attained), safety of
  short words on the bounded-digit event.
* Tail: unbalanced words, balanced noncharacters with `h >= q/2`, nonprimitive words.
* **The count of words with vanishing initial moments (Astra's second priority):** at most
  `binom(q+t-1, t) binom(q+t-k, t-k+1)` primitive balanced words of length `h = 2t` have `p_1, ..., p_{k-1}` equal;
  injectivity of `(P, N_1) -> v` by Newton; probability `a^{k-1}` by independence of the `xi_s`; exponents
  `F = 49/1000 - 3gamma/4 <= -11/1000` and `F = 49/1000 - gamma/2 <= -1/1000`, recomputed exactly. Checked exactly
  at `q = 7, 11, 13`.
* Final normalisation; positive Haar measure on the realisation component (tangent space spanned by `1` and
  `K_{S^c}^T`; computed at `q = 43, 59, 103, 127`).
* **(VIII.10), quoted:** *there is a sequence of angle models satisfying every monomial gap and axiom (R) at every
  odd prime power `<= 2m` with `limsup log P_+(N')/log m <= 23/12` and `limsup W/(m log m) <= 23/3`; for every
  `0 < delta < 1/12` its members eventually lie in the window `P_+ <= m^{2-delta}`.* The parameter system closes
  exactly at `xi = 1/12`.
* **Dependency (now settled):** (VIII.10) uses the uncertainty estimate and `gcd((Kv)_{S^c}) <= 2` up to
  `h <= q^tau`, `tau` near `1/4`. The proof of the estimate has no restriction on `h`, and the union bounds close for
  every `tau < 1/2` [J, §1.8].

**Computed.** `F_11` counts `40370, 3520, 440, 0` reproduced (labels in all of `F_11`; over `F_11 \ {0}` the counts are
`22530, 1980, 240, 0`; immaterial). All margins are asymptotic (`q^{0.003}` to `q^{0.005}`); `P_+ <= m^{1.995}` needs
`q > 10^{180}` roughly.

### 1.5 fakedefs: what the object is and what it caps (referee: CORRECT WITH CORRECTIONS; record wording FALSE)

**Proved (conditional only on the blocks of §§1.3–1.4, which passed).**
* The object is an angle model in the record's exact sense (overview Def 5.2, `construct.md` Def 1.2): `M = m` marks,
  exponent cap two, `P_+(N') <= 2q^d = m^{d+o(1)}`, `W <= (4d/rho + o(1)) m log m`, `omega(N') >= (q+1)/2 > M/4`,
  free dimension `m_free = q - m ~ M`, `Lambda = sum p^{-1/2} -> infinity`, zero windings, all rational characters
  integral, arc constant `pi/20`. **Corrections:** it is inside the window only for fixed `delta < 1/12`, not for all
  `delta > 0`; the support slip "`q/2` nonzero entries" should read `(q-1)/2` (no effect); convergence is very slow.
* It answers round 21's (Q) affirmatively for `0 < delta < 1/12` and these profiles. (Q) stays open for
  `delta >= 1/12`.
* No conflict with round 21: Theorem A is an upper bound `exp(-2.17 sum p^{-1/2})` on a Haar proportion; Astra's
  good set has positive Haar measure far below it (about `(2/H)^{q-m}`), so both hold. Cor A.1–A.2 assume `phi` is drawn from Haar
  measure on the coset; Astra's law is absolutely continuous but concentrated, so that hypothesis fails. J1 assumes
  one equal-weight window of width `1/(4r)`; her primes span `sigma ~ q^{1-d} log q`, 50–100 times wider at computed
  sizes. Lemma T and Theorem C have absent hypotheses. **The round-21 heuristic "window fakes admit no
  positive-probability recipe" is refuted for `delta < 1/12`.**
* Weight: within row 1 as worded (`A_0` with almost-everywhere Haar genericity) `8 -> 23/3`; within `A_{2M}`,
  `27 -> 23/3`. It does not certify the Haar-*majority* genericity or the Lang–Waldschmidt condition `LW_J(eps)` of
  the record's 8-witness (rows 11–12 unchecked; the single-form support condition for characters other than pairs
  is a **Gap**). Brackets: `[3, 23/3]` for `A_{2M}`; `[4 - o(1), 23/3]` for models whose primes all exceed `5M/4`.
* **False (record wording):** the expert-note sentence that the low-height window "remains uncapped", and the
  overview §5.5 bullet "every known fake puts its primes in one window `[P, P e^{1/(4r)}]`", are outdated (§4).

### 1.6 spread: VIII.4 (referee: CORRECT; one CORRECT WITH CORRECTIONS)

**Proved.**
* **(VIII.11) area test:** *if all active primes of a normalised genuine cluster lie in `[P, P e^sigma]` and the
  chosen Gaussian factors have arguments in an interval of length `eta`, `0 <= sigma, eta <= 1`, then
  `M <= 1 + K_C + O((1 + sqrt P sigma)(1 + sqrt P eta))` with an absolute `O`.* Novelty minimal (round 21 `prove.md`
  Lemma 3.1 with a radial bin, plus VIII.1).
* The first VIII.3 fake has all angles within `O(q^{-0.243})` of one direction and `sqrt P sigma = q^{0.005} log q`;
  the test allows `O(q^{0.757} log q)` marks against `>= (q+1)/2`. The bare sector count detects it; VIII.1 is not
  needed. With the (VIII.10) parameters the product is `o(q)`.
* **Correction (function-field status):** an ultrametric sector count exists (computed in `F_9[t]`, saturating at
  `sqrt P`); what has no counterpart is a radial bin narrower than one degree, which is exactly what detects the
  first fake.
* **Spread construction:** the models of (VIII.10) can have distinct prime angles separated by `>= 1/(3q)`, and
  `>= 1/(6q)` up to conjugation, at the same weight and height. Every step checked (rational background with
  denominator `12q`, `|T| = O(q^{xi/2} + 1)`, full `F_3` column rank of `K_{S^c,T}`, digit reordering, winding shift
  cancelling in (R), tails).
* **Referee consequence (new, not stated by Astra):** the spread fakes are `L`-spread in the sense of round 21
  `prove.md` Def 5.1. This discharges the fake-side hypothesis of round 21 Cor 5.4. The window would then be capped
  for `A_0` plus robust prime-supported inputs, *once* the genuine angles of a suitable prime set are also spread
  (Caveat 3, **open**).
* Astra's verdict ("the missing input is the association of each prime norm with its genuine Gaussian argument")
  restates the record's diagnosis (round 21 Theorem C, Lemma T; P3-all*; round 31 §1). It is not new.

**Computed (all reproduced exactly):** `q = 383`: 191 positive rows, 17,766 pairs of length 384, 95 lifts (maximum
coefficient 7,558 with her interpolation completion; 9,215 or 8,661 with the proof's own randomised construction,
both below `q^{1.5} sqrt(log q)`), `|T| = 124`, rank 124 mod 3. At `q = 383`, `|T| > sqrt q`, so this rank is an
observation outside the asymptotic range. **Illustration (floating point, one sample):** the max-scale cross-gap
`|sin(phi_p - arg w)| >= c/p` fails for 226 of 382 spread angles at `c = 1`; genuine angles give 0.

### 1.7 extend: analyst (VIII.1 re-derived CORRECT WITH CORRECTIONS; three negative answers to the referee's own questions)

**Proved.**
* VIII.1 with the `(1/2) log 2` sharpening (identical to rank's Correction 1), and the exact barycentre identity
  `E d_ij - W/2 = sum_p e_p log p [2t_p^2 - 2(s_p - 1/2)^2 - 2 m_p]` (notation in `ref8_extend.md` §2).
* **Scope:** VIII.1 uses only pair instances of axiom (R) and the norm-fibre counts. Hence it holds in every
  arc-admissible `A_1(Q_0)` model with `Q_0 >= max S` (up to the normalisation constant of (R), which the judge did
  not recheck). Every `A_1` fake of the record satisfies it; Paley profiles satisfy it with a factor-2 margin
  (`omega = q - 1 ~ 2m`, affine rank `m - 1`).
* **Cleanest corollary:** `M <= 1 + K_C + omega_1(rad N)` (`omega_1` counts primes `= 1 mod 4`), for every `R`. With
  `kappa = log N/log rad(N)`, `M <= (2/kappa + o(1)) log R/loglog R`; this beats `2/3` exactly when `kappa > 3`.
  With Robin's bound (quoted, not rechecked) an explicit bound with no threshold `R_0` follows, but with a worse
  constant than Theorem 1.1.
* Graded Theorem 1.1: `W' >= 3M log M + 2M sum_{p <= M, p not dividing N'} log p/(p-1) - O_C(M)` (minimal novelty).
  The only obstruction between 3 and 4 is small split primes dividing `N'` with balanced levels.
* Robust form: `c^T G' c >= (G_S - log C) ||c||_1^2` on the subspace orthogonal to `1` and the residue indicators.

**FALSE (referee's question (1), not Astra's claim):** "inserting VIII.1 into Theorem 1.1 improves the constant 3".
Prop 3.2: data satisfying the record's inequality (5.1), with `omega >= M` and affinely independent rows, have
`W = 3M log M + O(M)`. So these inputs together give exactly `3`. Astra made no such claim.

**GAP (referee's question (2)):** no strengthening of VIII.1 (weighted rank, sub-clusters, bounded exponents)
gives a constant `> 3`. Negative result, correctly stated.

**FALSE (referee's question (3)):** "an integer-only refinement of the barycentre argument that the Paley or VIII.3
fakes violate". Prop 5.1: every inequality derived from pair instances of (R) at moduli `<= Q_0` holds in
`A_1(Q_0)`. A separating input must use genuine prime arguments or prime-level residues, the archimedean order, or
primes of `N'` above `2m`.

**Ruling.** None of these is a gap or an error in Part VIII. They are negative answers, and they agree with Astra's
own statement that she did not improve the integer growth rate or `2/3`.

**Disagreement resolved [J].** Extend reports the first exponent dependence at arc constant `3.34`; rank reports
`2.83`. My exhaustive enumeration (`viii1_1e6_C9.out`) finds rank-deficient 3-point windows at arc constant
`2.8287` (for example `N = 707281 = 29^4` and `N = 579121`, `N >= 10^4`), confirming rank. Extend's `3.34` is an
artefact of its data set (norms with at least 64 representations). In chord constant the minimum is exactly `2.0`
(`N = 25, 50, 100, 125`), confirming rank's "no deficiency below `C' = 2`".

### 1.8 Judge's checks [J]

**(a) VIII.1, line by line.** I re-derived each step:
1. Normalisation is a similarity; chords scale by `1/|g|` and `sqrt R` by `1/sqrt|g|`, so `C' <= C/sqrt|g| <= C`.
2. If an inert `l` divided `N'` it would divide every `z_i'` (it is prime in `Z[i]`), contradicting gcd 1; likewise
   `1+i`. So `N'` is odd and split-only.
3. `v_{pi}(g_ij) = min(a_ip, a_jp)` and `v_{pibar}(g_ij) = e_p - max(a_ip, a_jp)`, so
   `log|g_ij| = (W - d_ij)/2`. With `|g_ij||b_ij| <= C' e^{W/4}` this is (VIII.3).
4. An affine dependence of the enlarged rows splits into probability vectors `P, Q` with disjoint supports and
   equal barycentres in every coordinate, including the indicators, so the residue laws mod each `l` agree.
5. `|x - y| = x + y - 2xy - 2 min(x,y)(1 - max(x,y))` on `[0,1]^2`; with `i ~ P`, `j ~ Q` independent,
   `E|X - Y| <= 2mu(1-mu) <= 1/2`, so `E d_ij <= W/2`.
6. The reductions mod `l` lie in the norm fibre (conjugation is Frobenius on `F_{l^2}`; `N' mod l != 0`), of size
   `l + 1`. The collision probability is `sum_t P_l(t)^2 >= 1/(l+1)`.
7. A collision gives `l | g_ij b_ij`, and `l` does not divide `g_ij` because `N(g_ij) | N'`; so `l | b_ij`. Distinct
   inert primes and `1+i` are pairwise coprime, so `|b_ij| >= sqrt2 prod_{colliding} l`.
8. Hence `E log|b_ij| >= (1/2) log 2 + sum_S log l/(l+1) > log C'`, and (VIII.3) gives `E d_ij > W/2`, a strict
   contradiction.
9. An indicator block lies in the hyperplane "coordinates sum to 1", so it has affine rank `<= l`; affine rank is
   subadditive under concatenation of coordinates. So `M - 1 <= omega(N') + K_C`.

I find no gap. The proof is complete for distinct marks.

**(b) VIII.1 on integer circles** (`viii1_test.py`, exact). Every circle `N <= 10^6` at arc constant `<= 9`
(663,125 windows with `M >= 3`, plus the record norms `1176852625, 32988780417125, 235370525, 10414625`), and every
circle `N <= 3*10^6` at arc constant `<= 4` (603,990 windows). In every window the script checked exactly: the
identity `N(gcd(z_i', z_j')) = prod_p p^{e_p - |a_ip - a_jp|}` for every pair; `(1+i) | b_ij` for every pair;
affine independence of the enlarged rows with Astra's `S` and with the sharpened `S`; `M <= 1 + omega(N') + K_C`; and
affine independence of the bare exponent rows when `C' < sqrt 2`. **Failures: 0 in every test, over 1,267,115 window tests (the two runs overlap on `N <= 10^6`, arc `<= 4`).** The largest deficiency `M - 1 - omega(N')` is 0 for chord constant `< 2` and at most 6 for chord constant
`<= 9`; on the 10-point cluster the affine rank is 5, then 6 with the mod-3 block and 9 with the mod-7 block
(`tenpoint.py`), as Astra states. First arc constants of each affine-rank deficiency `d` on circles `N >= 10^4`
(`viii1_1e6_C9_arc.out`): `d = 1, 2, 3, 4, 5` first occur at arc constant `2.8287, 4.0046, 5.1825, 6.6426, 8.7484` (witnesses `N = 855625, 21025, 21125, 60125, 27625`), reproducing rank's table exactly. Any valid `K_C` must be at least these deficiencies; the theorem's `K_C` is 3 to 4 orders of magnitude larger.

**(c) The step Astra asked to have prioritised that the referees leave least settled: the mod-four support bound.**
Fakecore proved it, spread re-derived its mod-3 form, and faketail left one open dependency: its use up to
`h <= q^tau` with `tau` near `1/4` in (VIII.10). I checked it line by line:
1. `C^2 = -qI`; for balanced `v` (zero at infinity) `Cv = (Kv, 0)` and `Cv = v mod 2`.
2. With `omega = (1+g)/2`, `g = sqrt(-q)`, the vector `y(v) = (Cv - v)/2 + v omega` is integral and satisfies
   `C y = g y`; `y_-` is the same with `1 - omega`, and `y_+ - y_- = gv`. Since `v` is primitive and `g` is a unit at
   every prime above `l != q`, one of `y_+, y_-` is nonzero modulo any prime `L` above `l`.
3. For `l = 2`: if `v_i = 0` and `(Kv)_i = 0 mod 4`, then `y_i = (Kv)_i/2 = 0 mod 2`. So
   `supp(y mod L) subset supp(v) cup supp(Kv mod 4)` (the infinity coordinate is `-(sum v)/2 = 0`). For odd `l` the
   same holds with `Kv mod l`. Mod 2 alone this inclusion fails, which is why mod 4 is needed.
4. `C_ab = chi(det(s_a, s_b))` with `s_i = (1,i)`, `s_inf = (0,1)`; for `gamma in SL_2(F_q)` with
   `gamma s_a = c_a s_{sigma(a)}`, the signed permutation `(Py)_{sigma a} = chi(c_a) y_a` commutes with `C`. It is
   defined over `Z` and transitive on `P^1(F_q)`, so it moves a coordinate nonzero mod `L` to infinity while
   keeping an integral eigenvector.
5. For that eigenvector, row infinity gives `sum f = -lambda y_inf`, nonzero mod `L`. Over `Q(g, zeta_q)` the
   Fourier transform gives `(chi(xi) G - lambda) fhat(xi) = 0` for `xi != 0` with `G^2 = -q`, so `fhat` is supported
   on `{0}` and one quadratic class. Hence `f * (f o n) = (sum f)^2/q` identically for a nonsquare `n`. This
   identity has coefficients integral away from `q` and reduces mod `L` to a nonzero constant on all `q` places. A
   convolution of two vectors with support `b` has support `<= b^2`, so `b >= sqrt q`, and the full support is
   `>= sqrt q + 1`.
6. Hence `|supp(Kv mod 4)| >= sqrt q + 1 - h`, and the same for odd `l <= h`. **The proof uses no bound on `h`.**
   The union bound `(2q+1)^{q^tau + 1} (q^tau + 1) rho^{sqrt q - q^tau - 1} -> 0` and the norm test
   `(2q+1)^{h} exp(-c_rho q/h) -> 0` both hold for every fixed `tau < 1/2`. So the dependency faketail flagged is
   discharged.

Computed (`mod4_check.out`, exact): the identity `C y(b) = g y(b)` for random integral `b` at `q = 7..47`; the
reduction mod a prime above 2 of the lattice `{y(b)}` has dimension `(q+1)/2` and exact minimum weight 4, 8, 8, 12 at
the split primes `q = 7, 23, 31, 47` (both primes above 2), and, over `F_4`, 6 and 8 at the inert `q = 11, 19`; all
are `>= sqrt q + 1`. (These are the extended quadratic-residue codes, including the `[48, 24, 12]` code at
`q = 47`.) On 300 random primitive balanced words per length `h = 2, ..., 10` at `q = 103, 127, 167, 199`, the
support inclusion of step 3 and the bound `>= sqrt q + 1` for the nonzero reduction held every time; the minimum
support of `Kv mod 4` was 39 to 101, far above `sqrt q + 1 - h`. The numerics do not stress the bound at these
sizes; the proof carries it.

### Withdrawn or corrected wording

Nothing of Astra's mathematics is withdrawn. Corrected wording: (VIII.2) and "for `C < 1` the empty set" (sharpened
to `C' < sqrt 2`); "four-full integer norms" (means `R` squarefull); "`c_delta`" (no `delta`-dependence); the
COMPUTED 10-point illustration of (VIII.6) (outside the lemma's range); "no verbatim counterpart at ultrametric
infinity" (only for radial bins finer than one degree); the uncertainty bound (`sqrt q + 1 - h`); the fixed-`rho`
coefficient (`7.96/rho`, `23/(3 rho)`); the support slip `q/2` to `(q-1)/2`; the broken file links in `final8.txt`
(the scripts are in `round31/astra_free_VIII/`).

## 2. New relative to the record and the literature, and correct

1. **VIII.1, `M <= 1 + omega(N') + K_C` unconditionally, and `M <= 1 + omega(N')` for `C' < sqrt 2`.** The record had
   `omega(N') > M/4` and `omega(N') >= (M/2)(1 - M^{-1/3})^2` for `C <= 1` (round8 `stability.md` Thm 3.3, Prop 3.4),
   a cut-space Gram bound `M <= Omega(N') + 1` (Prop 4.1(b)), and `M <= omega(N') + O_C(1)` only as the *ceiling* of
   conditional or method-limited inputs (round 20 `wild.md` Prop 6.8, (UI-l); lit18 Props A, B). The new ideas are
   the barycentre inequality `E|X - Y| <= 2mu(1-mu)` in exponent space and the inert residue blocks. Literature
   novelty is **unverified**: two referees ran one web search each and read snippets only (Cilleruelo–Granville,
   Cilleruelo–Córdoba, Bourgain–Rudnick); none showed a bound of this shape.
2. **(VIII.4)–(VIII.5):** `P_+(N') >= (2 - o(1)) M log M`; coefficient `2/kappa`, hence `1/2` for squarefull `R`.
3. **(VIII.6)** in uniform form with absolute `c_*` (new packaging of record Prop B2; low novelty).
4. **The first fakes inside the low-height window** (VIII.3, (VIII.10)), with point residues at `2M`, positive Haar
   measure and non-arithmetic angles, at `W <= (23/3 + o(1)) M log M`; the first point-residue fakes below 27. New
   devices: random assignment of primes in a wide interval with martingale control of the Paley margins (bypassing
   J1); the mod-four square-root support bound used as an erasure estimate; integral moment lifts with base-`B`
   digits and the count of words with vanishing initial moments; residue witnesses ordered by auxiliary Gaussian
   prime directions.
5. **Spread window fakes** (VIII.4), and the referee observation that they are `L`-spread (round 21 Def 5.1).
6. The area test (VIII.11): correct, minimal novelty.
7. Referee additions: the `(1/2) log 2` sharpening and the 2-adic block; `K_C = Theta(C^4/log C)`; the `kappa`
   corollary; the exact barycentre identity; no-additivity with Theorem 1.1; the window bracket `[3, 23/3]` for
   `A_{2M}`.

## 3. Rulings

**(a) Is `M <= 1 + omega(N') + K_C` a THEOREM about integer circles? Yes.** Exact statement: the corrected
Theorem VIII.1 quoted in §1.1, with `S` any finite set of primes `l = 3 mod 4` such that
`(1/2) log 2 + sum_{l in S} log l/(l+1) > log C'`, and `K_C = sum_{l in S} l`. The smallest admissible `K_C` (prefix
`3, 7, 11, 19, ...`) is 0 for `C' < sqrt 2`, 3 for `C' < 1.73`, 10 at `C' = 2`, 40 at 3, 94 at 4, 1150 at
`C' = 8.24`, and in general `Theta(C^4/log C)`. It is **new to the record** (it turns a conditional ceiling into a
theorem) and, as far as the referees' snippet-level searches show, not in the literature. Three qualifications:
(i) `K_C` is very lossy: the observed deficiency at `C' = 8.24` is 4, against `K_C = 1150`; (ii) it is valid verbatim
over `F_2[t]` and `F_3[t]`, and holds in every arc-admissible `A_1(Q_0)` model with `Q_0 >= max S`, so every fake of
the record satisfies it; (iii) consequently it cannot prove the uniform bound by itself.

**(b) Consequences for the growth rate.**
* *Subwindow:* clusters with `P_+(N') <= (2 - eta) M log M` have `M <= M_0(C, eta)` (Proved). This is a vanishing
  sliver of the window (`log P_+ ~ log M`); it does not touch `P_+ <= M^{2-delta}`. It is also an `A_1` statement.
* *k-full case:* if every prime of `N'` has exponent `>= k` in `N`, then `M <= (2/k + o(1)) log R/loglog R` (Proved);
  improvement on `2/3` only for `k >= 4`, i.e. for integer `R`, exactly when `R` is squarefull, with coefficient
  `1/2`. More generally `M <= (2/kappa + o(1)) log R/loglog R`, `kappa = log N/log rad N`.
* *General `R`:* **nothing improves.** VIII.1 gives only `W >= (1 + o(1)) M log M` in general, and combined with
  Theorem 1.1's inputs the constant 3 stays exactly 3 (extend Prop 3.2).

**(c) The height-order inequality (VIII.6): a theorem, of little use.** It is proved (given the Baker–Wüstholz
two-logarithm bound, quoted from the record's refereed `ref_baker.md`). It applies only to coordinates of finite
order modulo the row lattice, which are generic only in pinned or partly pinned profiles. It gives no growth bound,
alone or with VIII.1: the trivial upper bound on orders is `exp(W/6)`, far above `exp(cM)`. Its only new content is
that, in a pinned cluster, every coordinate at a prime `<= M^A` has order `>= exp(cM)`, where Liouville gives about
`M/4`. Its function-field failure (VI.5) shows it uses characteristic 0, but it has not been shown to fail on any
integer fake, so it has no demonstrated separating power.

**(d) The low-height-window fake: established.** For every fixed `0 < delta < 1/12` there are angle models in the
window satisfying every monomial gap, every rational character, axiom (R) at all odd prime powers `<= 2M` and all
units, with a positive-Haar-measure set of non-arithmetic angles on a realisation component, at
`W <= (23/3 + o(1)) M log M` (diagonal sequence). Every internal step passed its referee; the two steps Astra asked to
have checked first that a referee alone covered (the moment count, faketail; the residue witnesses, fakecore) are
proved, and the third (the mod-four bound) is proved by fakecore and rechecked by the judge, including its range in
`h`. **Capped in the window for `delta < 1/12`:** row 1 (`A_0` with almost-everywhere Haar genericity), row 2 at
`X = 2M` (`A_{2M}`), row 3 (via `stability.md` Prop 7.1), row 10 (Baker–Matveev, characters supported on `>= M`
primes), each at every coefficient `> 23/3`. **No conflict with round 21:** Theorem A (Haar proportion; positive),
Cor A.1–A.2 (Haar-drawn constructions; not this law), J1 (equal-weight windows; not used), Lemma T and Theorem C
(hypotheses absent) all stand. What falls is round 21's *heuristic* that window fakes have no positive-probability
recipe, and the record sentence that the window is uncapped.

**(e) What remains uncapped.**
* Every input in the deeper window `delta >= 1/12`, i.e. `(2 - o(1)) M log M <= P_+(N') <= M^{23/12}`.
* In every window: prime-level residues (P3-char, P3-all*); point residues at moduli `M^A`, `A > 1`; the union
  large sieve at the genuine constant (max-scale cross-gap), which genuine angles satisfy with constant 1 and the
  spread sample fails for 226 of 382 primes; square-root-quality (GRH) prime sums and box counts; Haar-majority and
  full-torus Haar-almost-everywhere properties; Lang–Waldschmidt rows 11–12; (VIII.6) (untested on fakes); the
  association of each prime norm with its genuine Gaussian argument. Row 7 robust prime-supported inputs are capped
  only conditionally (round 21 Cor 5.4 plus Caveat 3, open).
* Outside the window, the record's caps are unchanged (expert note §5).

**(f) Consequences for the integer problem.** **The best growth bound for general `R` stays
`M <= (2/3 + eps) log R/loglog R`. The uniform bound `M(C) < infinity` is untouched.** Part VIII adds one clean
unconditional structural theorem (VIII.1), which improves the growth coefficient to `1/2` only when `R` is
squarefull, and one barrier result (the window fake) that removes the low-height window, for `delta < 1/12`, from
the list of places where local inputs might still prove growth. Bluntly: the low-height restriction is not a cheaper
sub-problem for `A_0`/`A_{2M}` inputs; any proof of growth in the window must use one of the uncapped items in (e).
VIII.1 also shows that the record's "at best `M <= omega(N')`" ceilings for S-unit, Subspace and subgroup methods
(expert note §5 items 5 and 8) are already reached unconditionally, so those methods can add something only by
beating `omega(N')`.

## 4. Recommended edits to the expert note (recommendations only; nothing was written to the repository)

1. **§1, after Theorem 1.1, add:**
   "> **Theorem 1.2** (**proved**, round 31 Part VIII.1; `round31/README_VIII.md`). For a normalised cluster with chord
   constant $C'\le C$, let $S$ be a finite set of primes $\ell\equiv3\ (4)$ with
   $\frac12\log2+\sum_{\ell\in S}\frac{\log\ell}{\ell+1}>\log C'$. Then the exponent rows, enlarged by the residue
   classes of the marks mod each $\ell\in S$, are affinely independent, so $M\le1+\omega(N')+K_C$ with
   $K_C=\sum_{\ell\in S}\ell=\Theta(C^4/\log C)$; for $C'<\sqrt2$, $M\le1+\omega(N')$. Consequently
   $P_+(N')\ge(2-o(1))M\log M$, and $M\le(\frac2\kappa+o(1))\log R/\log\log R$ with $\kappa=\log N/\log\mathrm{rad}N$,
   which is $\frac12$ when $R$ is squarefull. The theorem holds verbatim over $\mathbb F_2[t]$, $\mathbb F_3[t]$ and in
   every $A_1$ angle model, so it does not by itself give the uniform bound, and it does not improve $\frac23$ for
   general $R$."
2. **§1, replace** "and the low-height window $P_+(N')\le M^{2-\delta}$ remain uncapped; here $N'$ is the normalised
   norm and $P_+(N')$ its largest prime factor (`round21/README.md`, status **open**)," **by:**
   "remain uncapped. The low-height window $P_+(N')\le M^{2-\delta}$ ($N'$ the normalised norm, $P_+$ its largest
   prime factor) is **capped for $0<\delta<1/12$ for pair, monomial-gap, rational-character and point-residue
   inputs at moduli $\le2M$**: Legendre (Paley) profiles with exponent cap two and randomly assigned primes in
   $[P,2P]$, $P=M^{2-\xi}$, give angle models with $W\le(23/3+o(1))M\log M$ and a positive-measure set of
   non-arithmetic angles (round 31 Part VIII.3, (VIII.10); refereed, `round31/README_VIII.md`). It is open for
   $\delta\ge1/12$, and in every window for prime-level residues, point residues at moduli $M^A$ with $A>1$, the
   genuine-constant union large sieve and square-root-quality inputs. Actual clusters have
   $P_+(N')\ge(2-o(1))M\log M$ (Theorem 1.2),"
3. **§1, fake bullets, add after the Paley bullet:** "* Legendre profiles with exponent cap two and primes below
   $M^{2-\xi}$, $\xi<1/12$, give fakes at $W\le(23/3+o(1))M\log M$ inside the low-height window, with point residues
   at moduli $\le2M$ (round 31 Part VIII.3; **proved**). They have positive, not majority, Haar measure on the
   realisation component and are not checked against Lang–Waldschmidt."
4. **§5 item 2, replace** "Exceptions inside the low-height window are listed in `round21/README.md`." **by:**
   "Inside the low-height window the record's earlier exceptions (`round21/README.md`) stand only for
   $\delta\ge1/12$ or for the inputs listed as uncapped in `round31/README_VIII.md` §3(e)."
5. **§5 item 5, append:** "Since $M\le1+\omega(N')+K_C$ holds unconditionally (Theorem 1.2), these methods add
   nothing unless they beat $\omega(N')$." **§5 item 8, after "capped at $M\le\omega(N')+O_C(1)$", insert:** "(a
   bound now proved unconditionally, Theorem 1.2)".
6. **§5 item 4, append:** "A uniform two-logarithm inequality $(1+w(v))\log(2D_L(v))\ge c_*W$ for finite-order
   coordinates holds on actual clusters (round 31 Part VIII.2; **proved**, low novelty, cf. `growth/auxiliary.md`
   Prop B2); it gives no growth bound and has not been tested against fakes."
7. **Pointer table, replace** "round 21: low-height window (open)" **by** "round 21 / round 31 Part VIII: low-height
   window (capped for $\delta<1/12$; open otherwise) | `round21/README.md`, `round31/README_VIII.md`".
8. **Sentences to avoid:** "the low-height window is uncapped"; "every known fake puts its primes in one window
   $[P,Pe^{1/(4r)}]$"; "window fakes have no positive-probability recipe"; "VIII.1 improves the constant $2/3$";
   "$M\le\omega(N')+1$ for $C<1$ only"; "fakes violate the rank bound"; "(VIII.6) separates fakes from circles";
   "coefficient $1/2$ for integer radii" (only squarefull radii).

## 5. Open sub-questions, ranked by honest probability of progress on the integer problem

All probabilities are the judge's estimates of affecting the integer problem (growth rate or uniform bound).

1. **Do pair inputs of `A_1` force more than `3M log M`?** P ≈ 4%. The `A_{2M}` bracket is now `[3, 23/3]`, and
   Theorem 1.1 uses only `A_1`-type inputs, so a proof that the per-pair inequalities (not their average) force
   `W >= (3 + eta) M log M` would improve `2/3` to `2/(3 + eta)` for every `R`. Extend shows the constant moves only
   through small split primes dividing `N'` with balanced levels. *First step:* solve the exact cut-cone linear
   program (column generation) for `M <= 40`, with adversarial levels at split `p <= M` dividing `N'` and adversarial
   residue classes; in parallel, build pair-only `A_1` models with balanced small primes and geometric residues. If
   `W = (3 + o(1)) M log M` is reached, `2/3` is the limit of pair inputs.
2. **Interpolating Theorem 1.1 and VIII.1 for `2 < kappa < 3`.** P ≈ 2%. Exponents `>= 2` give slack at middle
   levels in the barycentre identity. *First step:* test whether the chain lemma at large primes with `e_p >= 2`
   plus VIII.1 gives `W >= c(kappa) M log M` with `c(kappa) > max(3, kappa)`.
3. **The max-scale cross-gap against the spread fakes.** P ≈ 1% (high map value). This is the one uncapped input in
   the window that genuine angles satisfy with constant 1 and that involves all Gaussian integers, not only the
   cluster's. *First step:* either build a positive-measure modification of the spread fake satisfying
   `|sin(phi_p - arg w)| >= c/max(p, N w)` for a fixed `c`, or prove a Theorem-A-type scarcity bound under the
   construction's law.
4. **Caveat 3 (genuine angles of a short-interval prime set are `L`-spread).** P ≈ 0.5%. It would turn round 21
   Cor 5.4 into an unconditional cap of robust inputs in the window. *First step:* an upper-bound sieve for Gaussian
   primes in polar boxes of angle `>= P^{-1/2}` in norm intervals with `sigma sqrt(q log q) = o(1)`.
5. **Smaller `K_C`, and independence without residues for `C' < 2`.** P ≈ 0.3% (no growth effect; `K_C` polynomial
   in `C` cannot give the uniform bound by itself). Data: no deficiency below chord constant 2, deficiency `<= 6` for
   `C' <= 9`. *First step:* show that a minimal circuit with `E log|b| < log C` has most cross chords
   `b = unit * (1+i)`, and count such configurations.
6. **Window fakes for `delta >= 1/12`.** P ≈ 0.2%. The method stops because the digit sup norm `E = q^{1.5+a}` must
   fit below `B < P`. *First step:* replace the Hoeffding lifts by lifts of sup norm `q^{1+o(1)}`.
7. **(VIII.6) against pinned window fakes.** P ≈ 0.2%. *First step:* compute `min (1 + w(v)) log(2 D_L(v))/W` over
   finite-order `v` in the VIII.3 profiles and in any pinned fake.

**Recommendation on what to ask the colleague next.** Ask item 1: the per-pair `A_1` problem, with the explicit goal
of either an unconditional constant `> 3` (hence a growth coefficient below `2/3` for every `R`) or an `A_1` fake at
`(3 + o(1)) M log M` closing the bracket `[3, 23/3]`. This is the only item on the list whose success would change
the integer growth rate. Ask her to record VIII.1 with the `(1/2) log 2` sharpening and the `kappa`-corollary, and
to state (VIII.8) with explicit quantifiers in `rho`. Do not ask for further window constructions beyond item 3.

## 6. Scores (0–10, on the integer goal)

| block | correctness | score on goal |
|---|---|---|
| VIII.1 rank theorem and (VIII.4), (VIII.5) | correct; constant sharpened by referees; new | 3.5 |
| VIII.2 height-order inequality | correct; low novelty; no growth use | 1 |
| VIII.3 window fake, (VIII.8), (VIII.10) | correct in every step; one quantifier wording clarified | 3 |
| VIII.4 area test and spread fakes | correct; area test minimal novelty; spread version feeds round 21 Cor 5.4 | 1.5 |
| **Part VIII overall** | **no error found; careful, with honest scope statements** | **3** |

The quality as mathematics is high (about 8/10): the window fake is an intricate construction that survived four
independent referees, and VIII.1 is short and clean. The score reflects that neither the growth rate for general `R`
nor the uniform bound moved.
