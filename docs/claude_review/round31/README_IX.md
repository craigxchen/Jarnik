# Round 31, Part IX: judge report on projecting genuine marks and factoring residue witnesses at all heights

Source: Part IX of `astra_free_IX/answer.md` (heading "# Part IX. Projecting genuine marks and factoring residue
witnesses at all heights", line 1955, to the end, line 2315). Also read: Astra's summary `final9.txt`, the
coordinator's question `reply9.txt`, her work notes (`ix_atom_work.md`, `ix_crossgap_work.md`, `ix_midpoint_work.md`,
`ix_rigidity_work.md`, `ix_sharp_work.md`) and her scripts `ix_projection_check.py`, `ix_factor_check.py` (read,
not imported). Record: `expert_note/expert_note.md` §§1 and 5, `round31/README_VIII.md`,
`round5/round6_sharp/sharp.md` §1.1, `round5/lower.md` Def. 1.1, `round5/construct.md`.

Eight referees each took one block. **All eight reports are present; no referee failed.**

| key | report | scope |
|---|---|---|
| proj | `ref9_proj.md` | IX.2 (projection theorem (IX.1), certificate (IX.2), corollaries, field check, Paley scope), first paragraph of IX.6 |
| restate | `ref9_restate.md` | IX.1 (restatements), the rational-direction fake of IX.6, the Gram paragraph of IX.6 |
| factor | `ref9_factor.md` | IX.3 first half: Paley set-up, integral right inverse, witnesses (IX.6), same-level deletion, pair estimate, offsets |
| factor2 | `ref9_factor2.md` | IX.3 second half: nonpair slack (IX.7), union bound, real correction, spread variant, e-all failure, q = 11 computation |
| height | `ref9_height.md` | IX.4: support property, atom bound, parameter system closing at 1/6, quantifiers |
| crossgap | `ref9_crossgap.md` | IX.5: (IX.8), Schottky measure, conditional fibre measures, replicated profile, finite test, twelve-copy variant |
| rigidity | `ref9_rigidity.md` | IX.7: rigidity for rational prime witnesses, negative control |
| meaning | `ref9_meaning.md` | analyst: what Part IX changes in the map of inputs, what remains, what to ask next |

The judge's scripts are in `judge_IX/`. They use the judge's own exact Gaussian arithmetic and import nothing from
Astra or the referees:
* `proj_test.py`: enumeration of all lattice points on every circle `x^2 + y^2 = N`, `N <= NMAX`; all maximal
  windows of consecutive points with arc length `<= C sqrt R` starting in `[0, pi/2)`, `M >= 3`; exact normalisation
  by the Gaussian gcd; exact exponent rows; the larger-sieve chain of (IX.1) for every primitive `b` with
  `|b|^2 <= 50` and every certificate direction; the certificate identity, angle claim and derivative bound of
  (IX.2) for all pairs, all `2e_i` and (for `M <= 7`) every `c in {-1,0,1,2}^M` with `sum c = 2`, `||c||_1 = 4`;
  the pair identity `T(e_i + e_j) Norm(gcd(z_i', conj z_j')) = N'` and the bound `T <= C'^2 sqrt N'`. Every fibre
  count `nu_p` used is checked by brute force inside the run. Outputs `proj_2e5_C6.out`, `proj_1e6_C2.5.out`.
* `tenpoint_IX.py`: the record 10-point cluster under (IX.1). Output `tenpoint_IX.out`.
* `samelevel_test.py`: the IX.3 witnesses and the same-level deletion, built from scratch (own unimodular column
  reduction for the right inverse). Output `samelevel_test.out`.
* `rigidity_check.py`: controls for the IX.7 mechanism. Output `rigidity_check.out`.

**Every claim labelled [J] is the judge's own and has not been refereed.**

Definitions used throughout. `N = R^2` is a positive integer and the circle is `|z|^2 = N`. A *cluster* is a set of
`M >= 2` distinct points of `Z[i]` on the circle in an arc of length `<= C sqrt R`. The *normalised* cluster is
`z_i' = z_i/g`, `g` the Gaussian gcd of the points, `N' = N/|g|^2 = prod p^{e_p}`, `R' = sqrt N'`, `W = log N'`,
`s = omega(N')`, `P_+ = P_+(N')`. `N'` is odd and split-only; `a_ip = v_{pi_p}(z_i')` is the *exponent row*. The
*chord constant* is `C' = max |z_i' - z_j'| / N'^{1/4} <= C`. The *low-height window* is `P_+ <= M^{2-delta}`. An
*angle model* (fake) keeps primes and profile and replaces the genuine prime angles by real angles satisfying the
realisation equation and the monomial gaps (G). `A_{2M}` adds round 5's point-residue axiom (R) at odd prime powers
`<= 2M`. *P3-char* and *P3-all* are the classes `F^char_X`, `F^all_X` of `sharp.md` Def. 1.1 (factored residues:
per-prime residues `r_j(q)` of norm `p_j`, point residues the monomials in them; (e-char) tests zero-sum mark
characters at the four units, plus the same-level pair clause of `lower.md` Def. 1.1(e) at parent primes `p_l <= X`;
(e-all) tests every monomial `v` against collisions of `prod (r_j/conj r_j)^{v_j}` with `i^s`). *P3-all\** is the
record's open question whether `W_min^all = O(M log M)` at cutoff `X = 3M` (or `3M^A`). A *diagonal sequence*
lets the density `rho -> 1` (and here `xi -> 1/6`) while the Paley order `q` is taken large enough for each step.

## 0. Verdict in five lines

1. **One new unconditional theorem about integer circles, restricted to orientation:** for every `0 != b in Z[i]`,
   the integers `Re(conj(b) z_i)` of `M` lattice points on a circle have diameter `X_b` with `2 + X_b >> M^{3/2}`;
   hence an arc within angle `kappa R^{-1/2}` of `arg b` holds `<< (|b|(C kappa + C^2))^{2/3} + 1` points. Correct
   (two referees, and the judge on 24,568 windows and 1.68 million exact sieve tests, 0 failures) [J].
2. **It does not touch the growth rate for general `R`.** It needs a rational direction of height `o(s^{3/2})`
   essentially on the arc, which generic arcs do not have (the best direction has height about `R'^{1/4}`), and
   which the record's extremal data never have. It improves the classical near-axis/diagonal observation from
   exponent 1 to 2/3; the sieve step is very probably folklore.
3. **The factored-residue fakes are established:** one growth-scale angle model satisfies P3-char at **every**
   cutoff at once, with global witnesses in `Q(i)^*` of exact norm `p_j`, including the same-level clause at every
   parent prime (judge rechecked the deletion: 0 failures in 22,992 same-level identities) [J]; with IX.4 the window
   is capped for `0 < delta < 1/6` at coefficient `22/3`; with IX.5 the single-label max-scale cross-gap at constant
   one is capped outside the window.
4. **The rigidity theorem IX.7 is correct and is a characterisation of genuineness, not a bound.** Its hypothesis
   (one rational lift per prime) is empty at every finite cutoff, so it says nothing about P3-all\*, which remains
   the open middle of the residue column together with twisted collisions at `M^A`.
5. **For the integer problem the growth rate `(2/3 + eps) log R/loglog R` and the uniform bound are unchanged.**
   Part IX is careful, correct and mostly a further shrinking of the map. Score on the goal: **2.5/10**.

## 1. Referee-corrected statements per block

Labels: **Proved** means a referee checked a complete proof line by line, or an exact computation was reproduced
with independent code. **Computed**, **Heuristic**, **Gap**, **False**, **Withdrawn** as usual. **No referee found
an error in any mathematical step of Astra's.** The FALSE verdicts in the structured summaries (meaning, three
items) concern questions the referee posed or framings of the coordinator's question, not claims of Astra's (§1.8).

### 1.1 restate: IX.1 (referee: CORRECT on every item; the IX.6 fake CORRECT WITH CORRECTIONS, §1.6)

**Proved.**
* **VIII.1 restated** (Theorem 1.2 of the note), exactly as in `README_VIII.md` §1.1, with the `(1/2) log 2` term:
  *for `S` a finite set of primes `l = 3 mod 4` with `(1/2) log 2 + sum_{l in S} log l/(l+1) > log C'`, the exponent
  rows enlarged by the residue classes mod each `l in S` are affinely independent; `M <= 1 + s + K`,
  `K = sum_{l in S} l`; `M <= 1 + s` for `C' < sqrt 2`; `M <= 2 + s` for `C' < 2^{3/4}` (block mod 2).* Distinct
  marks are required. The mod-2 block was already in README_VIII §1.1: not new.
* **Referee addition (Proved, Computed):** the 2-adic tower. The unit classes mod `(1+i)^k` meeting a fixed odd norm
  number `f_k = 1, 2, 4, 4, 8, 8, 16, 16` (`k = 1..8`). *Appending the class of `z_i'` mod `(1+i)^k` gives affine
  independence whenever `(1/2) log 2 * sum_{j<=k} 1/f_j + sum_{l in S} log l/(l+1) > log C'`, at extra affine cost
  `f_k - 1`.* In particular `M <= 4 + s` for `C' < 2`. The 2-adic contribution is capped at `(5/4) log 2`. 0 failures
  on 108,320 windows (`N <= 10^6`, arc `<= 4`).
* **Radical corollary:** *`M <= 1 + K_C + omega_1(rad N)`, hence `M <= (2/kappa + o_C(1)) log R/loglog R` with
  `kappa = log N/log rad N`, the `o_C(1)` uniform in `kappa`.* Not new (README_VIII §1.7).
* **(VIII.8) with quantifiers:** *for fixed `0.99 < rho < 1` and all primes `q = 3 mod 4`, `q >= q_0(rho)`: an
  `A_{2m}` window model with `m = (rho/2 + o_rho(1)) q`, `P_+ <= 2q^{1.99} <= m^{1.995}` (needs `q >~ 10^184`),
  `W <= (7.96/rho + o_rho(1)) m log m`; `7.96 + o(1)` only along a diagonal; no threshold uniform in `rho`.* This is
  the judge's wording of README_VIII §1.3.

### 1.2 proj: IX.2 and the orientation criterion (referee: CORRECT on every item; two CORRECT WITH CORRECTIONS)

**Proved (referee line by line; judge line by line and computed [J], §1.9(a)).**
* **(IX.1), the integer projection theorem, quoted in full (corrected wording):** *Let `M >= 3` distinct points of
  `Z[i]` lie on `|z|^2 = N` (no arc hypothesis). For `0 != b in Z[i]` let `X_b` be the diameter of the integers
  `x_i = Re(conj(b) z_i)`, and `D_b = N|b|^2`. For odd `p` let `nu_p = (p - chi_4(p) + 1 + (D_b/p))/2` if `p` does not
  divide `D_b`, `nu_p = 1` if `p = 3 mod 4` divides `D_b`, `nu_p = p` if `p = 1 mod 4` divides `D_b`; `nu_2 = 2`. Then*
  `log(2 + X_b) >= sum_{p <= M} log p/nu_p - 2 theta(M)/M`,
  *and consequently, with absolute constants,*
  `2 + X_b >= c_0 M^2 exp(-sum_{p <= M, p = 1 (4), p | D_b} log p/p) >= c_1 M^{3/2}`.
  The same holds for the normalised cluster with `N'` in place of `N`. The exponent `3/2` is the method's limit in
  the worst case (all split `p <= M` dividing `N`); if no split `p <= M` divides `D_b` it gives `M^2`.
* **(IX.2), the exponent-row certificate:** *for every `c in Z^M` with `sum c_i = 2`, put `v_p = sum_i c_i a_ip - e_p`
  and `T(c) = prod_{p | N'} p^{|v_p|}`. Then `prod z_i'^{c_i}/N' = u alpha/conj(alpha)` with `|alpha|^2 = T(c)`, `u`
  a unit; `b = eta alpha` (`eta in {1, i, 1+i, 1-i}`, `eta/conj eta = u`) has `arg b` within `||c||_1 Delta/4` of the
  arc midpoint modulo `pi`, and `X_b <= C^2 |b| (||c||_1/4 + 1/2)`. Hence*
  `(1 + ||c||_1) sqrt(T(c)) >= c_2 C^{-2} M^{3/2}` *with `c_2` absolute.*
* **Corollaries:** *(a) along any sequence with `s -> infinity`, `min_{sum c = 2} (1 + ||c||_1)^2 T(c) = o(s^3)`
  implies `M = o(s)`, hence `M = o(log R/loglog R)`; (b) if two distinct rows satisfy
  `sum_p |a_ip + a_jp - e_p| <= 2`, then `M <<_C 1 + P_+^{2/3}`; (c) for fixed `eta > 0`, such clusters with
  `P_+ <= M^{3/2-eta}` have `M <= M_0(C, eta)`; (d) for a squarefree `N'` the hypothesis of (b) says that the two
  rows agree in at most two columns.*
* **Orientation and single-mark criteria (IX.6, first paragraph):** *if some `0 != b` has
  `|b|(1 + sqrt(R') |sin(theta - arg b)|) = o(s^{3/2})`, `theta` the arc midpoint, then `M = o(s)`; for `c = 2e_i`,
  `T = N'/gcd(Re z_i', Im z_i')^2`, so a mark of primitive length `o(s^{3/2})` suffices.*
* **Paley scope:** *the VIII.3-4 fakes have `log T(c) >> q log q` for every `c` with `sum c = 2`* (at least `(q-1)/2`
  nonzero `v_j`, sharp; Computed `q = 7..43`).
* **Scope sentence to add (referee, Proved):** *a pair certificate means `z_j' = (b/conj b) conj(z_i')`, i.e.
  `z_j'` is the reflection of `z_i'` in the line through `b`; always `T(e_i + e_j) = N'/Norm(gcd(z_i', conj z_j'))
  <= C'^2 sqrt(N')`. Generic large circles have pair `T ~ N'^{0.2-0.3}`, and the extremal data of Theorem 1.1
  (`ref8_extend.md` Prop 3.2) admit no `c` with `T = o(M^3)`.*

**Corrected wording.**
* *Function-field check:* proved for constant quadratic extensions `F_{p^2}(t)/F_p(t)` (`p^{deg T} >>_p M^3`).
  "A direct test on every field example, including VI.5 and VII.1" must read "on every field example over a
  constant quadratic extension (VI.4, VI.5, the first model of VII.1)"; VII.1's ramified model `eta^2 = Delta(t)` is
  not covered as written (**Gap**, minor, unused).

**Computed (proj, reproduced).** 6,079 fibre counts `p <= 251`; 10-point cluster: five `T = 1` pairs, projections
48503, 48505, 48509, 48511, 48515 (in Astra's frame), sieve chain `13.44 <= 32.80 <= 49.70`; D = 12 ternary field
witness.

### 1.3 factor: IX.3 first half (referee: CORRECT on every item, one CORRECT WITH CORRECTIONS)

**Proved (referee; judge for the deletion and witnesses, §1.9(b)).**
* *Paley profile.* `K_ij = chi(j - i)` on `F_q`, `q = 3 mod 4`, column 0 deleted, `a_ij = 1 + K_ij`, `e_j = 2`, marks
  `S` sampled at density `rho` from `{i : K_i0 = 1}`. On the full class column `j` has `(q-3)/4` entries 2 and
  `(q-1+2chi(j))/4` entries 0, so every column survives normalisation with probability `1 - o(1)`; mark `i` has an
  odd exponent exactly at label `j = i`.
* *Right inverse.* `D` (rows `a_i - a_{i0}`) satisfies `D = H mod 2` (`H` rows `e_i - e_{i0}`) and has an integral
  right inverse `U` (all Smith invariants 1). `V = U(H - D)/2` is integral and `F = I + 2V` gives `DF = H` (IX.4).
* *Witnesses (IX.6).* With windings having even differences and `u = U((k_i - k_{i0})/2)`,
  `r_j = i^{u_j} pi_j prod_l gamma_l^{V_jl}`, `gamma_l = pi_l/conj pi_l`, has `r_j conj(r_j) = p_j` and every mark
  ratio at every odd exterior prime power is `gamma_i/gamma_{i'}`; every zero-sum character gives
  `prod gamma_i^{c_i}`, independently of the choice of `U`.
* *Same-level deletion at parent primes.* At `p_l` omit `gamma_l`; for a same-level pair the missing coefficient
  `K_il - K_{i'l}` is 0, the owner `l` is the unique mark at level 1 and is never in a same-level pair, so the ratio
  is again `gamma_i/gamma_{i'}` at every power. The record (lower.md Def. 1.1(e)) requires nothing at `p_l` for
  pairs at different levels or for non-pair characters.
* *Global collision cost.* For every character and unit `u_0`, the product of all colliding powers, at every cutoff,
  divides `A_c - u_0 conj(A_c) != 0`, so it is `<= sqrt2 L_0^{n/2}`; for pairs at unit one the odd part divides
  `Im(pi_i conj pi_{i'})`, `1 <= D_1 <= L_0 |beta_i - beta_{i'}|`.
* *Narrow sector, (IX.5), offsets.* Primes in one sector of width `min(1,C)/1000` and one radial bin
  (`sigma = O_C(q log q/P)`); `w(a_i - a_j) - W/2 = 2 log P_0 + o(1)` for every pair when `d > 3/2`; the unit-one
  pair inequality holds with a factor of about 35 to spare. All bounds are uniform in the cutoff.

**Corrected wording.** (i) `m/omega(N') -> rho/2` at fixed `rho`, `-> 1/2` only along a diagonal. (ii) The
record's same-level clause is unit-one; "every Gaussian unit" is stronger than required and is also proved.
(iii) "parent primes below `2m`" is vacuous: every `p_j >= q^d > 2m`. (iv) The witness heights depend on `U` and are
enormous for a generic `U`; irrelevant for P3-char, decisive for e-all.

### 1.4 factor2: IX.3 second half (referee: CORRECT, two CORRECT WITH CORRECTIONS)

**Proved.** The nonpair bounds `h >= 1.2q`, `(n+1)h >= q||c||_2^2 + nq/2`, `h >= qn/m` (for `n >= 4`); (IX.7)
`h - b_0 n >= (1 + a_rho) q` and `h - b_0 n >= 0.49h` with `a_rho = (1 - rho)/2 - o(1)` (it tends to 0 as `rho -> 1`,
which is why `m/omega -> 1/2` needs a diagonal); the union bound with radii `2 L_0^{n/2} P_0^{-h/2}` for all units
and both gaps; the real correction `t_j = (2q)^{-1} sum K_ij eps_i` with an unrestricted common phase; the spread
variant's windings `2Kb` have even differences.

**IX.3, corrected statement, quoted in full (factor2):**
> *Fix `C > 0`, `rho in (0.99, 1)` and a moment-construction parameter set as in VIII.3, VIII.10 or IX.4 with prime
> exponent `d` (any `d` admissible there). For every prime `q = 3 mod 4` with `q >= q_0(C, rho, d)` there is an
> angle model with `m = (rho/2 + o(1)) q` marks, `omega(N') = q - 1`, `N' = prod p_j^2`, `p_j in [q^d, 2q^d]`, and
> global witnesses `r_j in Q(i)^*`, `r_j conj(r_j) = p_j`, such that the model and witnesses lie in `F^char_X(m, C)`
> of `sharp.md` Def. 1.1 for **every** `X >= 3`, with the same tuple `(p, e, a, phi, theta_0, delta, k)` and with
> `r_j(q)` the reductions of `r_j` (singular factor deleted at parent primes). `W/(m log m) = 4d/rho + o(1)`.
> Diagonalising `rho -> 1` gives `m/omega(N') -> 1/2` and, with IX.4, `limsup log P_+/log m <= 11/6`,
> `limsup W/(m log m) <= 22/3` (with VIII.10 alone: `23/12` and `23/3`). For `X >= X_0(model)` the same model and
> witnesses are not in `F^all_X`.*

**Corrected wording.** "These witnesses fail e-all" is proved for these witnesses and, more generally, for every
witness family in the S-unit group of the parent other than the genuine primes up to associates; it is not "no
witnesses for this model satisfy e-all". The bad cutoff is `|Im alpha^n|` (astronomical). Any even `n` works. "Max
entry 9" is a property of FLINT's particular `U` (another `U` has 18). "d > 3/2" must read "d admissible for the
moment construction used".

**Computed (factor2, reproduced).** `q = 11`: 330 exterior and 90 same-level identities, failure of the singleton
test at `(23, 173, i^3)`, and smaller failures, e.g. `(13, 101, 0)`.

### 1.5 height: IX.4 (referee: CORRECT on every item, quantifier paragraph CORRECT WITH CORRECTIONS)

**Proved (given IX.3, now established).** The support property for `h <= q^tau`, `tau < 1/2`; the Halász-type atom
bound `sup_b P(sum w_i X_i = b) <= q^{o(1)} (A sqrt(qh))^{-1}` for `A = q^a`, any `a > 0`; the parameter system
closes exactly at `xi = 1/6` (two constraints bind there simultaneously: `tau < 1/2` from the square-root
uncertainty support, and `1.5 + a + kappa < d` from lifts of sup norm `q^{3/2}`).

**IX.4, corrected statement, quoted in full (height):**
> *Fix `C > 0`, `0 < xi < 1/6` and `rho in (0.99, 1)`. There is `q_0 = q_0(C, xi, rho)` such that for every prime
> `q = 3 mod 4`, `q >= q_0`, the construction of IX.3 + VIII.3 has a realisation with positive probability, and
> with positive Haar measure on its realisation component, which is a P3-char fake at every cutoff, with
> `m = (rho/2 + o(1)) q`, `P_+(N') <= 2q^{2-xi}` and `W/(m log m) <= 4(2 - xi)/rho + o_{C,xi,rho}(1)`. Consequently
> for each fixed `0 < delta < 1/6` and each `C` these fakes eventually lie in the window `P_+ <= m^{2-delta}`. Along a
> diagonal `xi_k -> 1/6`, `rho_k -> 1`: `limsup log P_+/log m <= 11/6`, `limsup W/(m log m) <= 22/3`. The endpoint
> `delta = 1/6` is not reached.*

### 1.6 restate (continued): the IX.6 rational-direction fake and the Gram paragraph

**Proved (corrected statement, quoted in full):**
> *Fix `b = 2 + i`, `beta = arg b`. Let `a > 0`, `K >= 1` satisfy `dist(k beta, (pi/4)Z) >= a(2 + |k|)^{-K}` for all
> integers `k != 0` (two-logarithm Baker–Wüstholz; the constants depend on `b`). Fix `A > 8(K+1)`. For all large
> primes `q = 3 mod 4`, `m = q + 1`, there is an angle model in `A_1(2m)` (all monomial gaps, axiom (R) at all four
> units at odd prime powers `<= 2m`, zero windings) with `s = 9m - 8` primes in `[P, P e^{1/(4s)}]`,
> `m^A <= P <= 2m^A`, `W <= (9A + o(1)) m log m`, arc constant 1, and all `m` angles in
> `[beta + m^{5/4} Delta, beta + (m^{5/4} + 1) Delta]`, `Delta = e^{-W/4}`. The real numbers `Re(conj(b) z_x)` lie in
> an interval of length `(sqrt5 + o(1)) m^{5/4}`.*

Meaning (restate, Proved): (IX.1) and the IX.6 criterion are not consequences of the `A_1(2m)` axioms; any proof of
them must use integrality of the projected coordinates (which Astra's proof does). The fake does **not** separate
(IX.1) from P3-char or the IX.3 axioms, says nothing about genuine clusters, and lies far outside the window.

**Gram paragraph (Proved, minimal novelty).** In IX.3, `d_ij - W/2 = 2 log P_0 + o(1)` two-sided, so the real
cofactor magnitudes `B_ij = l_ij exp((d_ij - W/2)/2)` satisfy `max log B_ij = O_C(log m)`. The mean statement is a
triviality: the all-ones vector gives `mean log B_ij <= W/(4(M-1)) + log C'` in every model. What stays outside:
integrality of `|b_ij|^2`, the complex cocycle `g_ij b_ij + g_jk b_jk + g_ki b_ki = 0`, and the coupling of
`arg b_ij` to the actual Gaussian primes.

### 1.7 crossgap: IX.5 (referee: CORRECT, with corrections; one citation not read in full)

**Proved.**
* *Schottky measure.* There are absolute `p_0`, `c_* > 0` and a probability `mu` on quarter-turns with
  `|mu^(k)| <= C_0 (1 + |k|)^{-epsilon}` such that every `y` in its support satisfies
  `|sin((pi/2) y - arg varpi)| >= 1/max(p, Norm varpi)` for every Gaussian prime `varpi` and every `p >= p_0`, and
  `|sin((pi/2) y - arg w)| >= c_*/Norm w` for every `0 != w in Z[i]`. The support is a set of signed continued
  fractions with partial quotients in a sparse set `= 0 mod Q`; good approximants have norms divisible by a chosen
  `l`. The power Fourier decay is **cited** from Bourgain–Dyatlov (GAFA 27 (2017), Theorem 2, nonlinear phase);
  the referee checked the hypotheses but could not read the full text (arxiv blocked). Label: Proved modulo that
  citation.
* *Conditional measures on prescribed sums* (`g > g_0 + 2/epsilon` suffices).
* *Application, corrected:* *for some fixed integer `g >= g_*(mu)`, fixed `d > max(2, 2/beta)` (unspecified, since
  `epsilon` is not explicit), `rho in (0.99, 1)`, and Paley primes `q >= q_0(C, rho, g, d)`: angle models with
  `m = (rho/2 + o(1)) q`, `s = g(q - 1)`, `W/(m log m) -> 4gd/rho`, satisfying P3-char at every cutoff (as in IX.3)
  and lines 2 and 3 of (IX.8).*
* *Line 1 on the same model* along an unbounded sequence of parents; and the **elementary twelve-copy variant**
  (`g = 12`, any `d > 2`, line 1 only, no Fourier-decay theorem). Fully proved.
* *Referee addition (Proved, new):* these fakes satisfy round 21's sharp union large sieve with the genuine constant
  `D(K, X) = 2K + pi X/2` (Theorem LS item 2 had `2D`).

**Corrected wording.** As written, (IX.8) includes the owner's own direction `w in Z[i]^* Z pi_j`, where genuine
angles give sine 0; genuine angles satisfy all three lines with constant 1 for every other `w`
(`|Im(pi_j conj w)| >= 1`). The fakes satisfy the stronger version, so they satisfy every record form a fortiori.
Lines 2 and 3 are Lebesgue-null conditions (no Haar-positive fake can satisfy them).

### 1.8 rigidity and meaning (rigidity: CORRECT, two CORRECT WITH CORRECTIONS; meaning: see below)

**IX.7, corrected statement, quoted in full (rigidity; judge rechecked, §1.9(c)):**
> *Assume the monomial gaps (G). Fix a label `j`. Suppose the residues of `p_j` at all odd prime powers outside a
> finite set `E_j` are reductions of one `r_j in Q(i)^*` (hypothesis (L)), and that they satisfy the norm condition
> of Def. 1.1(d) there (then `r_j conj(r_j) = p_j` follows). Suppose the e-all inequality (1.2) at unit `-1` holds for
> `v = n e_j`, `n = 1, 2, ...`, at every finite cutoff. Then `r_j = u pi_j` with `u` a unit and `pi_j` a Gaussian
> prime of norm `p_j`, and `phi_j = eps_j arg r_j (mod pi)`. If the same tests hold for `v = n(e_i + e_j)` on the
> edges of a spanning tree, all `eps_j` are equal. In every case the model, with `phi_j` replaced by `arg rho_j`
> (`rho_j = r_j` or `conj r_j`), is an actual normalised cluster with common phase `theta_0 - pi sum_j e_j m_j` and
> the same `delta`, `k`. The cutoff needed depends on the instance (`E_j`, `J`, the Diophantine type of `phi_j`) and
> the angular step is ineffective.*

(G) cannot be dropped (`phi_j = 0`, `r_j = pi_j` passes every ray test). **Finite-cutoff lifts (rigidity Prop. 3,
meaning §3.2; Proved, and the judge agrees [J]):** *at every finite cutoff `X`, every factored datum (d) at odd prime
powers `<= X` is the reduction of `r_j = pi_j gamma/conj(gamma)` for some `gamma in Z[i]` (local Hilbert 90 plus
CRT), of exact norm `p_j` and height `e^{O(X)}`.* So (L) is a hypothesis on the projective limit only.

**meaning.** Not a re-referee of proofs. Its FALSE verdicts are on (a) the referee's own question whether Dirichlet
directions make the projection usable for generic arcs (it is vacuous: `X_b >> R'^{1/4}`), (b) the framing "every
residue-type input is either capped or the problem itself" (false: P3-all\* and twisted collisions sit in the
middle), (c) whether Astra's closing question is a route to the rate (it is a restatement). None is a claim of
Astra's; she says herself that she cannot produce certificates for general clusters and that IX.7 says nothing for
finite cutoffs. I uphold all three as answers to those questions [J]. Its **Theorem O** (orientation form of (IX.1)
in original coordinates) is correct; §3(a) adopts it.

**Computed (meaning).** IX.3-type witnesses with one choice of `V`: `log J_j ~ 1.2 q log P`; 6/10 labels at `q = 11`
and 16/18 at `q = 19` fail the singleton e-all test already at `X = 3m`. Evidence only (one `V`, small `q`).

### 1.9 Judge's checks [J]

**(a) (IX.1) and (IX.2), line by line, and on genuine circles.**
1. `x_i^2 + y_i^2 = D_b` with `y_i = Im(conj(b) z_i)`, so each value of `x` comes from at most two marks
   (`z = w/conj b`), `n >= M/2`.
2. Mod odd `p` not dividing `D`, the conic has `p - chi_4(p)` points; `x` has two `y` when `D - x^2` is a nonzero
   square and one when `D = x^2`; counting gives `nu_p = (p - chi_4(p) + 1 + (D/p))/2`. Inert `p | D`: `x^2 + y^2 = 0`
   forces `x = y = 0`. Split `p | D`: all `x`.
3. Cauchy–Schwarz: ordered pairs in the same class `>= n^2/nu_p - n`; each such difference is divisible by `p`;
   distinct primes dividing a nonzero integer contribute at most its log. Summing gives
   `sum_{p <= M} log p (n^2/nu_p - n) <= sum_{i != j} log|x_i - x_j| <= n(n-1) log X_b`.
4. Mertens with the progression `1 mod 4` for bad split primes gives the two bounds; `theta(M)/n = O(1)`.
5. Certificate: `z_i' = eps_i prod pi_p^{a_ip} conj(pi_p)^{e_p - a_ip}`, `z^{-1} = conj(z)/N'`, so
   `prod z^c/N' = u prod (pi_p/conj pi_p)^{v_p}`; `b/conj b = e^{2i arg b}` gives `2 arg b = sum c_i arg z_i (mod 2pi)`
   for any lifts, and with `sum c_i = 2` the deviation from the midpoint is `<= ||c||_1 Delta/4`.
6. Mean value theorem: `X_b <= |b| R' Delta (|sin(theta - arg b)| + Delta/2)` and `R' Delta^2 <= C^2`.
I find no gap.

Computation (`proj_test.py`, exact): all circles `N <= 2*10^5` at arc constant `<= 6` (19,108 windows) and all
`N <= 10^6` at arc constant `<= 2.5` (5,460 windows). **1,682,031 sieve-chain tests and 502,767 certificate tests
(identity with a unique unit, angle claim, derivative bound), 0 failures**; `n >= M/2` always; every fibre count used
agreed with brute force; the pair identity `T(e_i + e_j) Norm(gcd(z_i', conj z_j')) = N'` and the bound
`T(e_i+e_j) <= C'^2 sqrt N'` held in every case, and `T(2e_i) gcd(Re z_i', Im z_i')^2 = N'` always. Minimum of
`(2 + X_b)/M^{3/2}` is `0.265` (`N = 5525`, `M = 8`, `b = 1`, `X_b = 4`; a window symmetric about the axis), minimum of
`(1 + ||c||_1) sqrt(T) C'^2/M^{3/2}` is `1.21`. On the record 10-point cluster (`tenpoint_IX.py`): the cluster is
exactly symmetric under `x <-> y`, arc constant 8.2393, projections `x + y = 48503, 48505, 48509, 48511, 48515`,
chain `13.44 <= 32.80 <= 49.70`, reproducing proj and meaning.

**(b) The same-level deletion in IX.3 against the record's definition.** The record (`sharp.md` Def. 1.1(d), (e-char)
with "the same-level pair clause at the primes `p_l <= X`, exactly as `lower.md` Def. 1.1(e)") asks, at a parent
prime `p_l`: level residues `rho_x(p_l) = i^{-k_x} prod_{j != l} r_j(p_l)^{a_xj} conj(r_j(p_l))^{e_j - a_xj}` with
`Norm r_j(p_l) = p_j`, and, in the pair inequality for `c = e_x - e_y` with `a_xl = a_yl` only, the extra factor
`p_l^{a}` with `a = max{a : rho_x(p_l) = rho_y(p_l) mod p_l^a}`. Line by line:
1. `r_j^{(l)} = i^{u_j} pi_j prod_{l' != l} gamma_{l'}^{V_{jl'}}` is a unit mod `p_l` (all `p_j`, `j != l`, differ
   from `p_l`), of norm `p_j`, and `r^{(l)}_j/conj = (-1)^{u_j} prod_{l' != l} gamma_{l'}^{F_{jl'}}`.
2. With `e = 2`, `r^a conj(r)^{2-a} = p_j (r/conj r)^{K_ij}`, so
   `rho_x/rho_y = i^{-(k_x-k_y)} (-1)^{sum_{j != l}(K_xj - K_yj) u_j} prod_{l' != l} gamma_{l'}^{sum_{j != l}(K_xj - K_yj) F_{jl'}}`.
3. Full sums are `(Du)_x - (Du)_y = (k_x - k_y)/2` and `((D_x - D_y)F)_{l'} = (e_x - e_y)_{l'}`; the omitted `j = l`
   terms carry the factor `K_xl - K_yl = 0` for a same-level pair. So the sign cancels the winding and the exponents
   are `(e_x - e_y)_{l'}`, `l' != l`.
4. `a_{il} = 1` iff `i = l`, so `l` is not in a same-level pair and the `l`-coordinate of `e_x - e_y` is 0: the ratio
   is exactly `gamma_x/gamma_y`, a unit mod `p_l`.
5. Hence `p_l^a | pi_x conj(pi_y) - conj(pi_x) pi_y = 2i Im(pi_x conj pi_y)`, a nonzero rational integer coprime to
   the exterior collision moduli; the total collision modulus of the pair at all cutoffs divides its odd part, and
   the pair inequality of §1.3 covers it. Different-level pairs and non-pair characters impose nothing at `p_l`.
The deletion is correct. The record's clause is unit-one; Astra's version at all four units also holds via
`D_{u_0}`.

Computation (`samelevel_test.py`, own exact code, own unimodular column reduction; `q = 11, 19, 23`, random `S` at
`rho = 1, 0.8`, 17 instances, primes `> 100` randomly assigned, nonconstant windings with even differences):
`DU = I`, `V` integral, `DF = H`; norms of all witnesses and level residues as required; **45,888 exterior pair
identities and 16,320 random zero-sum characters at moduli `l^a`, `l < 60`, `a <= 3`; 22,992 same-level identities
at `p_l, p_l^2, p_l^3`; 3,832 comparisons of the record's same-level collision exponent with
`min(3, v_{p_l}(Im(pi_x conj pi_y)))`; 0 failures.** Owner unique at level 1 in every case. Negative control (one odd
winding difference): 1,532 failures. My right inverse has entries up to `2.8 * 10^8` at `q = 23`: the witness
heights are huge, as factor says.

**(c) The rigidity proof of IX.7, line by line.**
1. *Structure.* The norm of `r_j = eps (1+i)^t prod q^{c_q} prod pi_p^{a_p} conj(pi_p)^{b_p}` is `p_j`, so `t = c_q = 0`
   and `a_p + b_p = [p = p_j]`; `d_p = a_p - b_p = a_p + b_p mod 2`. So `r_j/conj r_j = u alpha/conj alpha` with
   `u = eps^2 = +-1`, `J = |alpha|^2 >= p_j` odd, and `arg alpha/pi` irrational (disjoint supports, `J > 1`).
2. *Clean real parts.* For `l in E` not dividing `J`, `alpha` is a unit mod `l` (including `l = 2`); with `L` a
   multiple of 4 and of these orders, `alpha^n = 1 mod l`, so `Re alpha^n = 1 mod l`; for `l | J`, exactly one prime
   above `l` divides `alpha`, so `2 Re alpha^n = conj(alpha)^n` is nonzero mod it. So `Re alpha^n` is odd and prime to
   every parent prime.
3. *Charging.* If `q^k || Re alpha^n`, `q` odd, outside `E`: `alpha^n = -conj(alpha)^n mod q^k`, `conj alpha` a unit,
   `u^n = 1`, and (L) makes `r_j mod q^k` the residue; so the collision with `-1` has exponent `>= k`. Taking
   `X >= |Re alpha^n|`: `|cos(n phi_j)| >= p_j^{-n/2} |Re alpha^n| = (J/p_j)^{n/2} |cos(n arg alpha)|` (IX.9).
4. `|cos(kL arg alpha)| >= 1/2` infinitely often gives `J <= p_j`, so `J = p_j` and `r_j = eps pi_j` or `eps conj pi_j`.
5. *Angle.* With `z = r_j`, step 3 gives `|cos(kL phi)| >= |cos(kL beta)|`. `a_k = cos(2kL phi) - cos(2kL beta) >= 0`
   has Cesàro mean 0 (`2L phi`, `2L beta` not in `2 pi Z`, by irrationality from (G) and from disjoint supports);
   `a_k^2 <= 2a_k` gives mean `a_k^2 = 0`; expanding,
   `a_k^2 = cos^2 A + cos^2 B - cos(A+B) - cos(A-B)` has mean 1 unless `L(phi +- beta) in pi Z`. So
   `phi = eps beta + pi t/L (mod pi)`.
6. *Removing `t`.* On `n = 1 mod L'` the `E`-part of `Re z^n` is frozen at `B`, so `|cos(n phi)| >= |cos(n beta)|/B`;
   `n phi mod pi` is dense on the progression; as `n phi -> pi/2`, `|cos(n beta)| -> |sin(pi t/L)|`, forcing `L | t`.
7. *Edges.* If the single signs and the edge sign are not all equal, one of `2 beta_i`, `2 beta_j`,
   `2(beta_i + beta_j)` is `0 mod pi`, impossible. *Recovery:* `2a_xj m_j pi = 0 mod 2pi`, so the `pi`-multiples shift
   every mark by the same `pi sum e_j m_j`, absorbed in `theta_0`; the drifting common phase is kept.
No gap. Controls (`rigidity_check.py`, label 13, `pi = 3 + 2i`): fake lift `pi (2+i)^2/(2-i)^2` with the genuine
angle fails at `n = 2` (cutoff 19); fake lift with `gamma = 4 + i` fails at `n = 1` (cutoff 29); genuine lift with
angle `-arg pi + pi/3` fails at `n = 4` (cutoff 17); genuine lift and either genuine sign pass for `n <= 80`. In
practice the failures come at small `n`; the theorem's mechanism needs no such luck but gives no rate.

**(d) Referee disagreements.** There are none of substance. Factor and factor2 agree on the diagonal reading of
`m/omega -> 1/2`; height and factor2 agree that `d` must be admissible for the moment construction; rigidity and
meaning give the same finite-cutoff lift proposition by two routes (CRT with local Hilbert 90; weak approximation on
the norm-one torus) [J: both are correct]. Restate's open dependency (two-sided pair surplus for the `g`-column
version) is discharged by crossgap §4 (all primes in one bin).

### Withdrawn or corrected wording

Nothing of Astra's mathematics is withdrawn. Corrected wording: "every field example, including VII.1" (constant
extensions only); `m/omega(N') -> 1/2` (diagonal only); "parent primes below `2m`" (vacuous); "these witnesses fail
e-all" (these and S-unit witnesses only; bad cutoff astronomical); "d > 3/2" (admissible `d`); "max entry nine"
(depends on `U`); IX.4 threshold `q_0(C, xi, rho)` over Paley primes; IX.5 "for fixed `g`" (`g >= g_*(mu)`; `d`
unspecified; coefficient `4gd/rho`); (IX.8) including the owner's direction (genuine angles fail it there);
IX.6 "absolute `a, K`" (depend on `b`); IX.7 norm hypothesis redundant; "M <= 2 + s for C' < 2^{3/4}" (already in
the record).

## 2. New relative to the record and the literature, and correct

* **(IX.1)/Theorem O (new to the record; literature novelty unverified, sieve step probably folklore).** Round 4
  `hardness.md` Lemma 4.2(b) gave `M <= 2|w|(kappa + C)^2 + 2` near a rational direction `w`, and the trivial
  projection count gives `M <= 2 + 2|b|(C kappa + C^2/2)`. Gallagher's larger sieve on the projected norm fibre
  replaces the linear dependence by the power `2/3` (`1/2` when no small split prime divides `N|b|^2`). A brief
  web search by the judge and snippet searches by two referees found no source for this application
  (Cilleruelo–Córdoba and Cilleruelo–Granville treat short arcs by other means, e.g.
  [Close lattice points on circles](https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/close-lattice-points-on-circles/E14DBCBF07B4B98FECF3282E2BA672B5),
  [lattice points on short arcs](https://arxiv.org/pdf/2012.10784)) [J].
* **The exponent-row certificate (new, elementary):** a sum-two row combination with small `T` is exactly a
  rational direction of height `sqrt T` on the arc; for pairs, a near reflection symmetry.
* **All-height P3-char fakes (new, barrier):** first cap of a residue class uniform in the cutoff, with parent
  primes below the cutoff and global witnesses; supersedes `sharp.md` Thm 5.1's `33 max(A,2)` for large `A`.
* **The window to `delta < 1/6` at `22/3` (quantitative improvement of a fake, 1/12 -> 1/6).**
* **Cross-gap fakes (new barrier):** constant one against all Gaussian primes, absolute constant against all lattice
  directions, finite sharp test; plus the referee's union-sieve corollary at the genuine constant. Outside the window.
* **IX.6 fake (new, map entry):** integrality of projections is not implied by `A_1(2m)`.
* **IX.7 (new as an elementary residue statement):** the residue-language analogue of `composite.md` Prop 3.1 and the
  round-31 sieve rigidity; makes `sharp.md` §6.6 rigorous at all heights. A characterisation, not a bound.
* **Referee additions:** 2-adic tower (`M <= 4 + s` for `C' < 2`); finite-cutoff rational lifts; union large sieve
  at constant `D(K, X)` for IX.5 fakes; pair bound `T(e_i + e_j) <= C'^2 sqrt N'`.

## 3. Rulings

**(a) The projection theorem is a THEOREM about integer circles. Exact statement (Theorem O):** *Let `M` distinct
lattice points lie on `|z| = R` in an arc of length `<= C sqrt R` with midpoint direction `theta`. For every
`0 != b in Z[i]`,* `M^{3/2} <= c_3 (2 + |b|(C sqrt R |sin(theta - arg b)| + C^2/2))` *with `c_3` absolute; in
particular, if `|sin(theta - arg b)| <= kappa R^{-1/2}` then `M << (|b|(C kappa + C^2))^{2/3} + 1`, and
`M << (|b|(C kappa + C^2))^{1/2} + 1` when no split `p <= M` divides `N|b|^2`.* The normalised version, with the
row certificate (IX.2), is the form quoted in §1.2.
* **New or repackaging?** Both, honestly. The *hypothesis* (arc near a rational direction of small height) is the
  classical one (round 4 Lemma 4.2(b), CG07's axis cases, Chan's bound near axis points). The *conclusion* is new to
  the record: exponent `2/3` instead of `1`, which changes the admissible range from `|b| kappa = o(log R/loglog R)`
  (trivial projection count) to `|b| kappa = o((log R/loglog R)^{3/2})` [J]. The restricted rate result
  "`M = o(omega(N'))` under a certificate with `(1 + ||c||)^2 T = o(s^3)`" is that statement in row language.
  "Bounded `M` for nearly complementary rows with `P_+ <= M^{3/2-eta}`" is genuinely beyond the trivial count
  (which would only give `P_+ >> M`, weaker than Theorem 1.2's `(2 - o(1)) M log M`), so in that class it is a real
  improvement; but the class is near-reflection-symmetric clusters.
* **For which clusters does it give anything?** Only clusters whose arc lies within angle `O(T^{1/2} kappa R'^{-1/2})`
  of a Gaussian direction of height `o(s^{3/2})`: a proportion about `(6/pi^2) T^2/sqrt R'` of arc positions
  (meaning, computed), i.e. `(log R)^{O(1)}` positions out of `sqrt R`. Exactly symmetric clusters (`T = 1`, all
  record clusters with `M >= 8`) are covered, where the classical bound already gives `O_C(1)`. Generic clusters have
  certificates with `T ~ N'^{0.2-0.3}` and best directions of height `~ R'^{1/4}` (Dirichlet; vacuous, meaning
  §2.5); the extremal data of Theorem 1.1 admit no small certificate. It gives no parent-only statement, so nothing
  for `P_+ <= M(log M)^A` as a class of radii.

**(b) The all-height factored-residue fakes are established.** IX.3 passed both of its referees and the judge's
check of the step Astra asked for first (the same-level deletion). Exactly capped now, at `W = O(M log M)`:
* **P3-char (`sharp.md` Def. 1.1 (a)-(d), (e-char), same-level clause) at every cutoff simultaneously**, on one
  model with global witnesses `r_j in Q(i)^*` of norm `p_j`, with parent primes below the cutoff: **yes**. Constant
  `23/3` with VIII.10, `22/3` with IX.4 (diagonal). The record's scope limit "all primes of `N` above `X`" is removed
  for this class.
* **Point residues at moduli `M^A` with `A` not fixed in advance** (`A_X`, axiom (R), four units, odd prime powers),
  via `F^char_X subset A_X`: **yes for zero-sum characters**, uniformly in `A`, also in regime (1.4). **Not**: the
  2-adic clause, and twisted point collisions (`sum c != 0`), which Def. 1.1 does not contain and IX.3 does not treat.
* **The window for `0 < delta < 1/6` at coefficient `22/3`:** **yes**, for `A_{2M}` inputs and all-height P3-char,
  as a limsup along a diagonal; open for `delta >= 1/6`.
* **The constant-one max-scale cross-gap:** **yes for single labels**: constant one against every Gaussian prime at
  every norm, absolute `c_*` against every Gaussian integer, constant one against lattice directions of norm
  `<= p_j` (the latter also elementarily, twelve copies), jointly with all-height P3-char, but only at prime heights
  `m^d`, `d > 2` large and unspecified: **outside the window**, and modulo the Bourgain–Dyatlov citation for lines 2
  and 3. Not capped: the cross-gap inside the window, the composite and product-scale cross-gaps of round 20, the
  union large sieve in full operator form.

**(c) The rigidity theorem is established**, with the corrections in §1.8. Meaning: within witnesses that come from
one rational number per prime, the all-height e-all test at the single unit `-1` on rays already forces an actual
cluster. Because a rational lift exists at every finite cutoff, the hypothesis is pure uniformity in the cutoff. So
**it is not true that every residue-type input is now either capped or the problem itself.** The two ends are
settled: character-tested residues at any height are capped (IX.3); all-monomial residues at all heights with
rational lifts are the problem (IX.7). The middle is open.

**(d) What remains uncapped and not equivalent to the problem.**
* **P3-all\***: e-all (every monomial, four units) at polynomial cutoffs `X = 3M^A`, `A` chosen after the weight
  constant. IX.3's witnesses fail e-all only at an astronomical, unquantified cutoff (proved); numerically they fail
  at `X = 3m` for one choice of `V` (computed); no fake and no lower bound.
* **Twisted point collisions** (`sum c != 0`) at moduli `M^A`, open since round 5.
* e-all at all cutoffs with adelic witnesses having no rational lift (no easier than P3-all\* for lower bounds).
* The window for `delta >= 1/6`: `(2 - o(1)) M log M <= P_+ <= M^{11/6}`.
* The cross-gap inside the window; composite and product-scale cross-gaps; the operator-form union large sieve.
* Square-root-quality prime sums and box counts; Lang–Waldschmidt rows 11-12; Haar-majority genericity;
  (VIII.6) (untested on fakes); the Gaussian-integer cocycle without the norm condition (open both ways since
  Part IV); the 2-adic residue clause.
* On the integer side: integrality of projections is a genuine separating input (IX.6), but it carries a rate only
  in orientation windows.

**(e) Consequences for the integer problem.** **The best growth bound for general `R` stays
`M <= (2/3 + eps) log R/loglog R`. The uniform bound `M(C) < infinity` is untouched.** Part IX adds one restricted,
position-dependent theorem (Theorem O) that gives `M = o(log R/loglog R)` only for arcs near Gaussian directions of
height `o((log R/loglog R)^{3/2})`, and three barrier results that remove further input classes from the list of
places where a proof of growth might hide. Bluntly: the thread has now capped every local and character-tested input
at every height, and the residue route to the rate reduces to one question the record has had open since round 6
(P3-all\*), or to the problem itself.

## 4. Recommended edits to the expert note (recommendations only; nothing was written to the repository)

1. **§1, after Theorem 1.2 and its bullets, add:**
   "> **Theorem 1.3** (**proved**; round 31, Part IX.2; `round31/README_IX.md` §§1.2, 3(a)). Let $M$ distinct lattice
   > points lie on $|z|=R$. For $0\ne b\in\mathbb Z[i]$ the integers $\mathrm{Re}(\bar b z_i)$ have diameter $X_b$ with
   > $2+X_b\gg M^{3/2}$ (absolute constant; $\gg M^2$ if no split $p\le M$ divides $N|b|^2$). Hence an arc of length
   > $\le C\sqrt R$ whose midpoint direction $\theta$ satisfies $|\sin(\theta-\arg b)|\le\kappa R^{-1/2}$ holds
   > $\ll(|b|(C\kappa+C^2))^{2/3}+1$ points.

   This is position-dependent: it gives $M=o(\log R/\log\log R)$ only when $|b|(1+\kappa)=o((\log R/\log\log R)^{3/2})$
   (the trivial projection count needs $o(\log R/\log\log R)$), and nothing for generic arcs, whose best Gaussian
   direction has height about $R^{1/4}$. In row language: a combination of exponent rows with coefficient sum two and
   $T(c)=\prod p^{|\sum_ic_ia_{ip}-e_p|}$ gives $(1+\|c\|_1)\sqrt{T(c)}\gg_CM^{3/2}$; two nearly complementary rows force
   $P_+(N')\gg_CM^{3/2}$. The record's ten-point cluster is symmetric about the diagonal. Angle models satisfying
   all of $A_{2M}$ violate the projection bound in a rational-direction window (Part IX.6), so it uses integrality."
2. **§1, replace** the paragraph from "Actual clusters satisfy $W\ge3M\log M-O(M)$. Not everything is capped." to
   "... square-root-quality prime sums and box counts (`round31/README_VIII.md` §3(e))." **by:**
   "Actual clusters satisfy $W\ge3M\log M-O(M)$. Not everything is capped. Residues that factor through prime residues
   of the genuine norms and are tested through characters (P3-char, `round5/round6_sharp/sharp.md` Def. 1.1 with the
   same-level clause) are capped **at every modulus at once**: round 31 Part IX.3 gives one angle model, with global
   witnesses in $\mathbb Q(i)^*$ of norm $p_j$ and parent primes below the cutoff, at $W\le(23/3+o(1))M\log M$, and
   $(22/3+o(1))M\log M$ with $P_+(N')\le M^{11/6+o(1)}$ (Part IX.4; diagonal sequences; refereed,
   `round31/README_IX.md`). Hence point residues at moduli $M^A$ with $A$ chosen after the cluster are capped for
   zero-sum characters at odd prime powers. The single-label max-scale cross-gap at constant one is capped outside
   the low-height window (Part IX.5). Still uncapped: P3-all$^*$ (prime-level residues on every monomial at moduli
   $3M^A$), twisted point collisions at moduli $M^A$, the sharp-constant large sieves in operator form and the
   composite cross-gaps, every cross-gap inside the window, square-root-quality prime sums and box counts,
   Lang–Waldschmidt and the Gaussian cocycle. With one rational lift per prime and all moduli at once, the
   all-monomial residue test holds only on actual clusters (Part IX.7); since such lifts exist at every finite
   modulus, this is the problem itself and says nothing about P3-all$^*$. The low-height window
   $P_+(N')\le M^{2-\delta}$ ($N'$ the normalised norm, $P_+(N')$ its largest prime factor) is **capped for
   $0<\delta<1/6$** for pair, cut, monomial-gap, rational-character, point-residue and character-tested prime-residue
   inputs (Parts VIII.3, IX.4), and stays **open** for $\delta\ge1/6$ and for the uncapped inputs above."
3. **§1, fake bullets, append to the Legendre bullet:** "With factored witnesses at every modulus and primes below
   $M^{11/6+o(1)}$, at $W\le(22/3+o(1))M\log M$ (round 31 Part IX.3–4; **proved**)."
   **In the residue-design bullet, replace** "Referee scope limits: all primes of $N$ must exceed the modulus bound
   $X$, and" **by** "Referee scope limits (removed for zero-sum characters by round 31 Part IX.3):".
4. **§5 item 1, replace the caveat by:** "Caveat: caps uniform in $A$ are proved for zero-sum characters and
   factored prime residues tested through characters (round 31 Part IX.3), not for twisted collisions or the
   2-adic clause; P3-all$^*$ is uncapped."
5. **§5 item 2, replace the last sentence by:** "Inside the low-height window the earlier exceptions
   (`round21/README.md`) stand only for $\delta\ge1/6$ or for the inputs listed as uncapped in
   `round31/README_IX.md` §3(d). The single-label max-scale cross-gap at constant one is capped outside the window
   (round 31 Part IX.5)."
6. **§5 item 12 (round 31 history), add a bullet:** "*Part IX* (`round31/README_IX.md`): Theorem 1.3 (orientation,
   restricted); all-height P3-char fakes; window capped for $\delta<1/6$; cross-gap fakes outside the window;
   rigidity for rational witnesses (a characterisation of genuineness, ineffective in the modulus). No change to the
   growth rate for general $R$."
7. **Pointer table:** replace the window row by "round 21, round 31 Parts VIII–IX: low-height window (capped for
   $\delta<1/6$; open otherwise) | `round21/README.md`, `round31/README_VIII.md`, `round31/README_IX.md`"; add
   "Theorem 1.3, orientation (projection) theorem | `round31/README_IX.md` §§1.2, 1.9(a), 3(a)"; add `README_IX.md`
   to the round-31 row.
8. **Sentences to avoid:** "P3-all is capped" or "P3-all$^*$ is answered"; "every residue-type input is capped or is
   the problem itself"; "point residues at moduli $M^A$ are capped uniformly in $A$" (zero-sum only); "Part IX
   proves growth in the window" or "for a class of radii"; "the max-scale cross-gap is capped" without "single-label,
   outside the window"; "the window is capped for $\delta\le1/6$"; "rational lifts are a restriction at a finite
   modulus"; "IX.3's witnesses fail P3-all$^*$" (computed for one $V$ and small $q$ only); "Theorem 1.3 is a new
   method" (the sieve step is probably folklore).

## 5. Open sub-questions, ranked by honest probability of progress on the integer growth RATE

All probabilities are the judge's estimates of a growth-rate statement for general `R` (in brackets: a useful
restricted or barrier result).

1. **P3-all\* in regime (1.4).** P ≈ 2% (15%). It is the only residue class neither capped nor the problem. A lower
   bound there would be a theorem about actual clusters (`sharp.md` §6.7). *First step:* a finite-cutoff,
   height-quantified IX.7: for witnesses `r_j = pi_j gamma/conj(gamma)` with `Norm gamma <= H`, find the least cutoff
   `X(H, E_j)` at which the ray tests on `n e_j`, `n <= n_0`, at the four units force `gamma/conj gamma = +-1`; in
   parallel, minimise the witness heights over right inverses `U` in IX.3 (a closest-vector problem) and test e-all at
   `X = 3m^A`. Either a fake passing e-all at `W = O(M log M)` or a forced bound `J_j <= p_j^{1+o(1)}`.
2. **The Gaussian cocycle without the norm condition.** P ≈ 1% (5%). *First step:* integral models on the Legendre
   profile with all divisibilities `g_ij | z_i - z_j` exact and norms only approximately `N`.
3. **Twisted point collisions at `M^A`.** P ≈ 0.5% (10%). *First step:* extend IX.3's object to words
   `v = K_S^T c`, `sum c != 0`, at `X = 2M`, where the free common phase enters.
4. **The deep window `delta >= 1/6`.** P ≈ 0.5% (5%). Needs both lifts of sup norm `o(q^{3/2})` and a support bound
   past `sqrt q`; or a restricted theorem at `P_+ <= M(log M)^A`.
5. **The cross-gap inside the window.** P ≈ 0.3%. *First step:* replace the union bound `(24q)^h P^{-h/2}` of the
   twelve-copy variant by a word count adapted to the coarse structure, at `P = q^{2-xi}`.
6. **Producing small-height directions for non-symmetric clusters (to use Theorem 1.3).** P ≈ 0.2%. Dirichlet gives
   only `R'^{1/4}`; this would be a near-reflection-symmetry theorem, which no known mechanism supplies.

**Recommendation.** Ask the colleague for one more pass, on item 1 only, posed as an explicit alternative: either a
factored-residue fake passing e-all at moduli `3M^A` (`A` after the weight constant) at `W = O(M log M)`, or a
height-quantified finite-cutoff version of IX.7 followed by a lower bound `W >= f(M) M log M`, `f -> infinity`, for
near-genuine witnesses. Twisted collisions on IX.3's object as a side item. Ask her to record Theorem O, the scope
sentence, and the corrected quantifiers. **If that pass returns another fake, stop the residue route for internal
work**: the thread will then have capped every local input class the record knows, and the remaining items are the
problem itself or have been open since round 6 with no mechanism proposed.

## 6. Scores (0–10, on the integer goal)

| block | correctness | score on goal |
|---|---|---|
| IX.1 restatements | correct; nothing new (mod-2 block already in the record) | 0.5 |
| IX.2 projection theorem, certificate, corollaries | correct; new to the record, restricted, probably folklore method | 3 |
| IX.3 all-height factored fakes | correct in every step; same-level deletion rechecked | 3 |
| IX.4 window to `1/6` | correct; quantifiers corrected | 2 |
| IX.5 cross-gap fakes | correct modulo an unread citation; outside the window | 2 |
| IX.6 rational-direction fake, Gram remark | correct; map entry; Gram remark near-tautological | 1 |
| IX.7 rigidity | correct; characterisation, not a bound | 1.5 |
| **Part IX overall** | **no error found; honest scope statements throughout** | **2.5** |

As mathematics the quality is high (about 8/10): IX.3–IX.5 are intricate constructions that survived independent
referees and the judge's own checks, and IX.2 is short and clean. The score reflects that neither the growth rate
for general `R` nor the uniform bound moved, and that the one integer theorem bites only on exceptional arc
positions.

Sources consulted by the judge's literature search: [Close lattice points on circles (Cilleruelo–Granville)](https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/close-lattice-points-on-circles/E14DBCBF07B4B98FECF3282E2BA672B5),
[Discrete fractional integrals, lattice points on short arcs, and diophantine approximation](https://arxiv.org/pdf/2012.10784).
