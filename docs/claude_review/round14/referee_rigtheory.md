# Referee report: round14/rigtheory.md (Rig(K))

Status: complete (referee: adversarial, independent code in `referee_rigtheory_checks/`)

## 0. Verdict

**The result survives in its headline form only after corrections: "Rig(K) undecided; no growth, no cap".** The
lattice theory and the one-modulus data are correct and were reproduced. There are three substantive problems.

1. **(M1) "Strict", RigS and Theorem 4.1 ignore a symmetry of the class.** The per-column associate twist
   `pi_j -> i pi_j` maps P3-all\* members to members with the same `W` (Prop. T below). On labels it acts by
   `ell -> ell + (m/2) t`.
   * As stated, RigS forces every member *and* its twist to be strict. That is impossible at any `q = 3, 5 mod 8` with
     `(p_j/q) = 1` (Lemma T'), so RigS implies the class is empty in the regime, i.e. it contains Theorem 4.1's
     conclusion.
   * At the level of Rig(K), the twisted genuine system is good at every `K` and non-strict at 14/17, 37/40 and 61/65
     forced moduli (`K = 3, 4, 5`). So Rig(K) can never yield RigS as stated.
   * The fix is easy: strictness up to `(eps, t)`. Theorem 4.1 then holds as a class-level transfer, with H^gen for
     the twisted primes.
2. **(M2) The claim "`q_core >= M` is reached at bounded `K`" needs `r >= M^gamma`**, which regime (1.4) does not supply.
3. **(M3-M5) The joint-core evidence is weaker than presented.**
   * Joint rigidity claims rest on HiGHS MILP runs (floating point, no certificate) yet are labelled "exact". §4a
     gives what I could certify exactly.
   * The `K = 4, 5` joint cores and the `K = 4` coherence are conditional on fixing genuine labels at 3 and 7, and at
     `K = 5` also on the whole flexible top above 163. Joint growth from `K = 3` to `K = 4` is therefore not
     established unconditionally: `37 <= q_core^(J,eps)(4) <= 61`.
   * The "anticode-threshold explanation" of the round-10 data is a numerical coincidence, not a consequence of
     Theorem 3.4. The rigid range extends well past the covered range (e.g. 137-157 at `K = 5`), and the onset at 163
     is F1-only.

**Verified independently with exact computation:**
* the class (parity, per-unit budgets, units, moduli) matches P3-all\* restricted to norm `<= K`;
* the genuine collision identity;
* forced sets for `K = 2..5`;
* one-modulus cores for `K = 2..5`, including every round-10 count compared;
* all 22 joint witnesses (genuine, all levels, 0 violations).

**Proofs checked:** Theorems 2.5 and 3.4, Lemmas 3.1 and 3.3, Cor. 3.5, Props 3.6(1-2) and 3.7.

## 1. Class definition check (P3-all\*, parity, budget, monomials, moduli, units)

Checked against `sharp.md` Def. 1.1 (e-all), Prop. 2.1, Cor. 2.2, Prop. 6.1 and `round13/README.md` §7.2.

* **Parity.** `ell_j(q) = nu_j(q) mod 2`, `nu_j(q) = [(p_j/q) = -1]`, nothing else constrained at prime level: this is
  exactly Cor. 2.2(1). Re-checked independently (`referee_rigtheory_checks/t_genuine.py`, own residue-field code, not
  the author's random-generator tables): genuine labels have this parity at every odd prime `q < 400`, `q` not dividing
  `N`, `ps = {5,13,17,29}`.
* **Per-unit budget.** `B(v,s) = w(v)/2 + [s odd] (log 2)/2` is (1.2) of `sharp.md` with `nu_s = 2^(-1/2)` for odd `s`.
  Collision of `(v, s)` at `q^a` means `prod_j (r_j/conj r_j)^(v_j) = i^s mod q^a`, i.e. `<v, ell(q^a)> = s log(i)`.
  Genuine identity "collision iff `q^a | D_s(Z_v)`, `D = (Y, X-Y, X, X+Y)`": re-checked at level 1 for all `v` with
  `||v||_1 <= 5`, all units, all `q < 400` (198,560 checks, 0 mismatches).
* **Monomials and moduli.** K-good = budget for every `v` with `1 <= ||v||_1 <= K` and every unit, demand summed over
  all `q in S_K = {q <= Q_K}` odd, not dividing `N`, at all levels. That is Rig(K) of round 13 §7.2 verbatim; it is
  **weaker** than P3-all\* (which charges every `v`, every `q <= X`), and the write-up says so (P3-all\* = Rig(inf) on
  `[3, X]`). Level data at `q = p_l` are not charged by (1.2), and the write-up does not charge them: correct.
* **Units / normalisation (minor slip, no effect).** §1 fixes "a generator with `log i = m/4`" and asks labels to be
  "compatible under reduction". For inert `q` these two conventions are incompatible at even levels: if `g_a` reduces to
  `g_1` and `log_(g_1) i = m_1/4`, then `log_(g_a) i` must reduce to `m_1/4` mod `m_1`; since `m_a/4 = q^(a-1) m_1/4`, this
  forces `log_(g_a) i = 3 m_a/4` whenever `q^(a-1) = 3 mod 4`. The unit of a collision must be read off at level 1 (it is the same at every level). The
  author's code does exactly this (`sigc2.c` reads `s` at level 1 and only tests "some unit" at higher levels), and the
  resulting `s <-> -s` relabelling at 3 and 7 (top-level normalisation) is harmless because the set of all labels at a
  free modulus is closed under `ell -> -ell`; at fixed (G) moduli the code uses `D_s` directly. Wording fix only.
* **The symmetry group of the class is larger than the write-up's (major, §3.6 below).** Besides per-column
  conjugation `eps`, the class is invariant under **per-column associate twists** `pi_j -> i pi_j`. The write-up's
  notion of *strict* (and hence RigS and Theorem 4.1) is not invariant under them.

**Verdict on the class:** the label class, budgets, monomials, moduli and units of Rig(K) are exactly the P3-all\*
budget data restricted to norm `<= K` (neither stronger nor weaker), up to the wording slip on `log i`. The coherence
notions "genuine class" and "genuine up to conjugation" are invariant under the full symmetry group (they live mod `n`);
"strict" is not (see 3.6).

## 2. Regime check (1.4)

* Rig(K), Def. 2.1-2.4 and Theorem 2.5 do not involve `M`, `c_0`, `A` except through `r`, the primes and
  `Q_K <= X = 3M^A`, which holds once `A >= K + o(1)`; since `K = K(c_0)` and `A` is chosen after `c_0`, the order of
  quantifiers in RigS ("there are `K(c_0)`, `A_1(c_0)` such that for `A >= A_1`") is right.
* **All computations are at `r = 4` (some `r = 5, 6`) with the five smallest split primes.** These are formally
  compatible with `p_j <= M^(h_0)`, `r <= rho M` for large `M`, but regime (1.4) is an asymptotic statement with
  `r` growing; the data are (correctly) labelled evidence. Round 13 asked for `r = 4..8` and primes `<= 400`; `r >= 7`
  and non-minimal primes were not run.
* **`r >= M^gamma` is an extra hypothesis, not part of regime (1.4)** (which only gives `r <= rho M`; near-extremality
  `3 M log M - O(M) <= W <= c_0 M log M` with `p_j <= M^(h_0)` forces `sum e_j >~ 3M/h_0`, not many columns). §4.2
  item 3 states it; the round summary ("Only `q_core >= M` is needed, which the anticode range reaches at bounded `K`")
  drops it. By the write-up's own Theorem 2.5(3), `q_core(K) <= 4A_K + 5 = O_K(r^(floor(K/2)))`; if `r = M^(o(1))`
  no bounded `K` reaches `Q_eta = (M^(1-eta), M]`, so RigS cannot come from Rig(K) at bounded `K` in that part of the
  regime. Correct statement: "reached at bounded `K >= 2/gamma + 1` **when `r >= M^gamma`**; for `r = M^(o(1))` it needs
  `K -> infinity` with `M`, where the bounded-`K` theory (Thm 2.5, Thm 3.4 constants, Prop. 3.6 asymptotics `r -> inf`,
  `K` fixed) does not apply as stated."

## 3. Proof checks

### 3.1 Theorem 2.5 (reductions): correct
(1) Restriction lowers demands. For the extension, take a collision-free parity-respecting label at a free `q`; every
lift has no collision, because a collision at `q^a` reduces to one at `q` (the unit changes consistently). The bound
`max F_K^ex <= 4A_K + 5` follows from `L4K.md` Thm 1(b). (2) Monotonicity of `Core_J`, of the exact and relaxed `Core_1`,
and of `F_K^ex`: re-derived, correct. The same argument proves monotonicity of `Core^eps_J`. **Caveat:** only the
*cores* are monotone. `q_core(K)` (the largest `Q` with `F_K ∩ [3,Q]` inside the core) is not monotone a priori,
because `F_(K+1)` gains new forced moduli below `Q`. The write-up only claims monotone cores, which is fine. (3) Correct:
`Core_J ⊆ Core_1^exact ⊇ Core_1^relaxed`. The computed `Core_1` is the relaxed one, which is a lower bound for the exact
one; `Core_J ⊆ Core_1^relaxed` is not proved, but nothing uses it.

### 3.2 Lemma 3.1, Cor. 3.2: correct
`|D_s| = nu_s^(-1) |Z_v| |sin(arg Z_v - s pi/4)|` is the last line of the proof of `sharp.md` Prop. 1.2. "Not a
genuine collision at `q_0`" in Cor. 3.2 should read "`(v, s)` is not a genuine collision with unit `s` at `q_0`". The
two readings agree here, because `q` cannot divide two of the `D_s` (otherwise `q | Z_v` for an inert `q`, or both primes
above a split `q` divide `Z_v`).

### 3.3 Lemma 3.3 and Theorem 3.4 (anticode rigidity): correct
The path induction in Lemma 3.3 is valid (`f(z + e_j) - f(z)` is a generator; `rho(0) = 0` gives `f(0) = 0`). In Thm
3.4, the map `g` is onto `Z/n`, so `|c'(D)| = n`, so `g` is a bijection, so `Lambda_(c') ∩ (D-D) = Lambda_c ∩ (D-D)`.
In (2), `Lambda_c ⊆ Lambda_(c')`, so `c' = u c` with `u` a unit. In (1), the Lemma 3.3 generators of `Lambda_(c')` lie in
`D + S_1 - D`, so equal index forces equality. Translation invariance holds.

Random test (`t34.py`; 3,000 random non-class `c'` at each of `q = 41, 79, 157, 263`): no counterexample to the
conclusion. The test is weak, since the hypothesis fails at all four moduli (see below).

### 3.4 Cor. 3.5 and Prop. 3.6: correct as stated; Prop. 3.6(3) and the "explanation" overreach
* Cor. 3.5(1)-(3) re-derived. The transport inequality is just the budget rearranged.
* Prop. 3.6(1)-(2) are correct and elementary. (2) is pigeonhole.
* **The covering hypothesis fails exactly where it would be needed.** `cover.py` uses balls and double balls (the
  anticodes of the write-up), with `ps = {5,13,17,29}`:
  * `K = 3`: covered forced moduli are `{3, 7, 11, 19, 23, 37}`. 31 and 41 are not covered, yet the one-modulus core
    reaches 41 (`n = 10 > |S_1| = 9`).
  * `K = 5`: covered moduli run up to 131, with gaps at 73, 83, 97, 103. **137-157 are not covered**, yet they are
    rigid.
  * At the first flexible `K = 5` modulus 163, my own scan finds 80 non-class survivors, **all F1** (F2 = 0).
  
  Theorem 3.4 only excludes F2 on covered moduli. It therefore says nothing about the onset of F1 flexibility at
  `n = 41`, and it does not reach 137-157. So "the norm-5 behaviour ... is exactly the anticode threshold
  `a(4,4) = |S_2| = 41`" (§0, §5.3 item 2) is a numerical coincidence, recorded as a successful prediction at `K = 6`.
  It is not a consequence of the theorem. It fails at `K = 7`, and at `r = 5` the core exceeds the anticode size, as
  the author notes. The word "explained" must be downgraded to "consistent with (evidence)".
* Prop. 3.6(3) needs `r >= M^gamma`, which is not part of regime (1.4) (§2).

### 3.5 Prop. 3.7: correct
(1): one collision costs `log q`. (2): `|D_s(Z_u)| <= sqrt2 |Z_v||Z_(v')| < q`.

### 3.6 Theorem 4.1 and RigS: the implication is valid, but RigS as stated is not a rigidity statement (MAJOR)

**Proposition T (referee; proved).** Let `X = (p, e, (a_x), phi, theta_0, delta, k, r) in F^all_X(M, C)` and
`t in {0,1}^r`. Put `r'_j(q) = i^(t_j) r_j(q)` at every level, `phi' = phi + (pi/2) t`,
`theta_0' = theta_0 - (pi/2)<e, t>`, `k'_x = k_x + 2<a_x, t>`, and leave everything else unchanged. Then
`X_t := (p, e, (a_x), phi', theta_0', delta, k', r') in F^all_X(M, C)`, with the same `p, e, a, W, r`, and its labels are
`ell'(q^a) = ell(q^a) + (m_(q^a)/2) t`.

*Proof.* (a) is unchanged.
* (b): `<2a_x - e, phi'> = <2a_x - e, phi> + (pi/2)(2<a_x,t> - <e,t>) = theta_0' + delta_x + (pi/2) k'_x`.
* (c): `<v, phi'> = <v, phi> mod pi/2`, and Lemma 0.1 only sees `<v, phi>` mod `pi/2`.
* (d): `Norm(i r_j) = Norm r_j`. Also `(i r)/conj(i r) = -r/conj r`, which shifts the label by `m/2`.
* (e-char): the point residues become
  `rho'_x = i^(-k'_x) prod_j (i r_j)^(a_xj) conj(i r_j)^(e_j - a_xj) = i^(-<e,t>) rho_x`. This is a common unit, and it
  cancels in `prod rho_x^(c_x)` since `sum c = 0`.
* Level residues at `p_l`: the factor is `i^(-2 a_xl t_l - ...)`, common to every pair with `a_xl = a_yl`. Only such
  pairs enter the same-level clause (`lower.md` Def. 1.1(e)).
* (e-all): `prod (r'_j/conj r'_j)^(v_j) = (-1)^(<v,t>) prod (r_j/conj r_j)^(v_j)`. So `m'_(v,s) = m_(v, s - 2<v,t>)`,
  `nu_s` keeps its parity, and `|sin(<v,phi'> - s pi/4)| = |sin(<v,phi> - (s - 2<v,t>) pi/4)|`. QED

For an actual cluster, `X_t` is the same cluster written with the associates `i^(t_j) pi_j` (Prop. 1.2 with that choice
of Gaussian primes).

**Lemma T' (referee; proved).** Suppose `n_q` is odd (`q = 3, 5 mod 8`), `(p_j/q) = 1`, `t_j = 1`, and `ell(q)` is
strict for some `eps` (`ell = u (eps∘L)` mod `m_q`, `u = 1 mod 4`). Then `ell + (m_q/2) t` is strict for **no** `eps'`.

*Proof.* Reduce mod 4. Since `n` is odd, `m/2 = 2n = 2 mod 4`. `L_j` is even, so `eps L_j = L_j mod 4` for both signs.
Hence `ell_j = L_j` and `ell'_j = L_j + 2 mod 4`, while strictness of `ell'` would force `ell'_j = L_j mod 4`. QED

**Consequences.**
1. Let `X` be any member and `t = e_1`. RigS applied to `X` and to `X_t` forces both to be strict at every `q in Q_eta`.
   By Lemma T', this forces `(p_1/q) = -1` at every `q in Q_eta` with `q = 3, 5 mod 8`. Whenever the primes allow such a
   `q` with `(p_1/q) = 1`, RigS as stated therefore implies **that the class has no member in the regime**, which is
   already the conclusion of Theorem 4.1. Such a `q` exists:
   * unconditionally if `p_1 <= (log M)^C` (Siegel-Walfisz);
   * under GRH for every `p_1 <= M^(h_0)` and large `M`. (Any column `j` can play the role of column 1.)
   
   So Theorem 4.1 is a valid implication, but in its present form its hypothesis is not a rigidity statement. It
   contains the conclusion.
2. At the level of Rig(K) the effect is unconditional and computable. The twisted genuine system `L + (m/2) e_1` is
   K-good for every `K` and every level, because it is the genuine system of `(i pi_1, pi_2, ...)`. It is in the
   genuine class at every `q`, since only mod `n` matters there. Yet it is **non-strict** (no `eps`, no
   `u = 1 mod 4`) at 14 of the 17, 37 of the 40 and 61 of the 65 forced moduli for `K = 3, 4, 5` (`twist.py`, exact).
   * Lemma T' predicts the subset with `n` odd and `(5/q) = 1`, and that subset is contained in the computed list.
   * Hence the "strict" joint core of §4.1 is essentially empty at every `K`, and **Rig(K) can never yield RigS as
     stated**, whatever the growth of the class-level core.
3. **Fix.** Replace "strict for `eps`" by "strict for `(eps, t)`":
   `ell(q^a) = u (eps∘L(q^a) + (m_(q^a)/2) t)`, with `u = 1 mod 4` and one pair `(eps, t)` for all `q in Q_eta`
   (coherent). The proof of Theorem 4.1 then goes through verbatim with `pi^(eps,t)_j = i^(t_j) pi_j^(eps_j)`. Pair
   collisions become `q | D_(k_x - k_y)(Z^(eps,t)_(a_x - a_y))`, and H^gen must be read for these primes. The coherence
   notion (mod `n`) and all `Core^eps_J` data are unaffected.
4. **Two further corrections.**
   * The symmetry group of good systems acting on labels is `{+-1}^r x (Z/2)^r` (`4^r` systems, up to `u`), not
     `{+-1}^r`.
   * The one-modulus "non-strict class survivors" at 41, 79 and 251 should be re-examined modulo local twists. Even so,
     they do not bear on the joint strict core, which item 2 settles.

**Smaller points in §4.2.**
* RigS names a `K = K(c_0)` that plays no role in its statement. It should either say "every `K`-good system (on
  `S_K`)" or drop `K`.
* Item 4 asserts `|sin((delta_x - delta_y - arg eta)/2)| >= 1/sqrt 2` for `eta != 1`. For `eta = +-i` this is
  `sin(pi/4 - |delta_x - delta_y|/2)`, i.e. only `>= 1/sqrt 2 - O(Delta)`. This is harmless.

## 4. Independent exact computations (referee code only; `referee_rigtheory_checks/`)

Everything below was written from scratch: residue-field discrete logs instead of the author's random-generator
tables, my own C enumerators, and my own budget arithmetic in Python integers. No author code was imported. The
author's witness pickles were read as data, after checking with `pickletools` that they contain only data opcodes.

| check | script / output | result |
|---|---|---|
| genuine collision identity, parity | `t_genuine.py` | 198,560 checks (`q < 400`, `||v||_1 <= 5`, all units), 0 mismatches |
| forced sets `F_K^ex`, `K = 2..5` | `fdfs.c`, `forced.py`, `forced_K4.txt`, `forced_K5.txt` | `#F = 7, 17, 40, 65`, max `41, 103, 199, 349`: **agree**. Exhaustive DFS. Every free modulus has an explicit verified free label. Scanned to 3000 / 8000 / 12000, which exceeds `4A_K + 5` (`A_K` from `L4K.md` Thm 1(b) with the Sidon set `{1,2,4,8}`, `B = 8`: 6549 for `K = 4`, 11745 for `K = 5`; 1625 for `K = 3`) |
| one-modulus relaxed core, `K = 2..5` | `onemod.c`, `onemod.py`, `onemod_K4.txt`, `onemod_K5.txt` | `q_core^1 = 7, 41, 79, 157`, first flexible `11, 43, 83, 163`. Isolated rigid moduli 89-103 (`K = 4`) and 181, 191 (`K = 5`). Round-10 numbers 47/83/97/101/103/107 (`K = 4`) and 167/199/211 (`K = 5`: 120/96, 120/40, 884/832) reproduced. First F2 modulus 19, 47, 127, 197 for `K = 2, 3, 4, 5`. All 96/40 survivors at 167/199 are F1: **agree** |
| all-norm one-modulus reference, `K = 4` | `onemod.py ... full` | rigid for every forced `q_0 <= 107`; first flexible 109 (72 survivors, all F1): **agrees** with §5.7 |
| 22 joint witnesses (`K = 3`, `G = {3,7}`; `K = 4`, `G = {3,7,11}`) | `vwit.py` | every witness: label set = my `F_K^ex` exactly; parity, lifts and absence of level-2 collisions OK; 0 budget violations in integers (512 resp. 1280 pairs `(v,s)`, both signs). The label at `q_0` lies outside every `eps`-class at 41-103 (`K = 3`) and at 67, 71, 79, 83, 89, 97, 101, 103, 107, 127, 151, 199 (`K = 4`): **agree**. These are genuine explicit non-genuine survivors (all levels, unconditional) |
| Theorem 3.4 hypothesis (covering by balls / double balls) | `cover.py` | covered forced moduli: `K = 3`: {3,7,11,19,23,37}; `K = 4`: up to 53; `K = 5`: up to 131 with gaps 73, 83, 97, 103; `K = 6`: up to 193 with gaps. The rigid one-modulus range extends past the covered range at every `K` (§3.4) |
| associate twist | `twist.py` | twisted genuine system non-strict at 14/17, 37/40, 61/65 forced moduli (`K = 3, 4, 5`); Lemma T' subset contained (§3.6) |

## 5. Implications for growth, checked against the barrier logic

* **No growth and no cap is claimed. Correct.** The `2/3` bound and the P3-all\* bracket are untouched.
* **Fake test (§1), correct.** Every fake with genuine primes, given the genuine labels of its primes (or any of its
  `4^r` twisted or conjugated versions), is K-good for every `K` and meets every rigidity conclusion. So rigidity of
  Rig(K) type can at most route P3-all\* to statements about genuine labels (P3-G-type pair statistics, H^gen, L1). That
  is where a fake would have to fail. P3-G fakes are known only for `log p >= 2.08 X`, so in regime (1.4) the routed
  question is not known to be capped. The write-up's §4.2 item 2 says exactly this.
* **Theorem 4.1 in its present form is not a usable reduction** (§3.6): RigS as stated is not invariant under the
  associate twists, which are symmetries of the class. It implies emptiness of the class in the regime (given a `q` as
  in Lemma T'). So it cannot be the output of any rigidity argument. After the `(eps, t)` fix the transfer is valid, with
  H^gen read for the primes `i^(t_j) pi_j^(eps_j)`.
* **The size claim needs `r >= M^gamma`** (§2). Without it, bounded `K` gives `q_core <= O(r^(K/2))`, which for
  `r = M^(o(1))` never reaches `Q_eta`.
* **Stall does not by itself cap P3-all\* (§4.3). Correct.** A cap needs all norms, pair demands and slab completion.
  The flexible tops computed here use F1 slack from uncharged moduli above `F_K^ex` (large cofactors). That slack is
  absent in P3-all\*, as the author's all-norm one-modulus test (reproduced at `K = 4`) indicates.
* **Status of "growth is supported".**
  * One-modulus core: reproduced exactly at `K = 2..5`. It is a stable power of the forced range at `r <= 6` with the
    smallest primes. This is evidence, correctly labelled.
  * Joint core: the only unconditional facts are
    * `Core^eps_J(3) = F_3 ∩ [3, 37]` (rigidity at 37 from the author's HiGHS run; my certification is in §4a);
    * `37 <= q_core^(J,eps)(4) <= 61`, from monotonicity and the verified witness at 67.
    
    The `K = 4` lower bound 61 and coherence are conditional on genuine labels at all levels of 3 and 7. The `K = 5`
    statement is conditional on genuine labels at 3, 7 **and on the entire flexible top (every forced modulus above
    163)**. That condition removes the moduli from which a joint deformation would draw slack (Cor. 3.5(3)), so it is
    close to a one-modulus test. "The joint core grows with `K` (`n = 9, 15`)" compares an unconditional with a
    conditional number: **joint growth is not established even at `K = 4`.**

## 6. Issues list

**Major**
* **M1 (strictness / RigS / Thm 4.1).** The class P3-all\* is invariant under per-column associate twists
  `pi_j -> i pi_j` (Prop. T). These act on labels by `ell -> ell + (m/2) t`.
  * The write-up's "strict" (`u = 1 mod 4`, `eps` only) is not twist-invariant.
  * RigS as stated implies that the class is empty in the regime (Lemma T', whenever a suitable `q in Q_eta`
    exists). So Theorem 4.1's hypothesis contains its conclusion.
  * At the level of Rig(K), the twisted genuine system is K-good for every `K` and non-strict at 14/17, 37/40 and
    61/65 forced moduli. So **no Rig(K) result can ever give RigS as stated.**
  * Fix: strictness up to `(eps, t)` (`4^r` genuine systems), with H^gen for the primes `i^(t_j) pi_j^(eps_j)`.
* **M2 (size claim).** "The size route 3 needs (`q_core >= M`) is reached at bounded `K`" requires `r >= M^gamma`,
  which regime (1.4) does not provide. For `r = M^(o(1))`, Theorem 2.5(3) gives `q_core(K) = M^(o(1))` at every bounded
  `K`.
* **M3 (certification standard).** Every rigidity (infeasibility) claim for joint cores rests on HiGHS MILP runs: `K = 3`
  at `q_0 = 37`; `K = 4` on `[11, 61] ∪ {73}`; `K = 5` on `[11, 157]`; `K = 4` coherence.
  * These runs use floating-point logs, a `1e-9` tolerance, and no checkable infeasibility certificate.
  * They are labelled "exact". The round standard asks for exactly certified infeasibility. They are
    solver-verified, not certified. §4a gives the referee's exact certification where it was feasible.
  * By contrast, the one-modulus rigidity claims are exhaustive enumerations and are exact; I reproduced them.
* **M4 (overclaimed explanation).** Theorem 3.4 excludes only F2 deformations, and only on covered moduli. The covered
  range is strictly smaller than the rigid range at every `K` computed: at `K = 5`, 137-157 are not covered; at `K = 3`,
  31 and 41 are not covered. The onset of flexibility at 163 is F1-only. So "round-10 data explained ... exactly the
  anticode threshold" is a numerical coincidence. It fails at `K = 7` and at `r = 5`. It is not a consequence of the
  theorem.
* **M5 (joint growth not established).** The `K = 4` and `K = 5` joint cores and the `K = 4` coherence are conditional
  on fixing genuine labels at 3 and 7. At `K = 5` they are also conditional on the whole flexible top above 163.
  * Unconditionally: `Core^eps_J(3) = F_3 ∩ [3, 37]` (37 per M3), and `37 <= q_core^(J,eps)(4) <= 61`.
  * "The joint core grows with `K` (`n = 9, 15`)" and "tracks the one-modulus core" are not established.

**Minor**
* m1. The `log i = m/4` convention is incompatible with reduction-compatible generators at even levels of inert `q`.
  Units must be read at level 1. The code does this correctly, and the swap at 3 and 7 is harmless by `ell -> -ell`
  symmetry.
* m2. RigS names a `K = K(c_0)` that does not occur in its statement.
* m3. `|sin((delta_x - delta_y - arg eta)/2)| >= 1/sqrt 2` holds only up to `O(Delta)` for `eta = +-i`.
* m4. Cor. 3.2 should say "not a genuine collision **with unit `s`**". This is harmless, since `q` divides at most one
  `D_s`.
* m5. §3.4 "local certificates provably fail" is exact only at `r = 4`, `K = 5`. In the regime it is heuristic.
* m6. Scope versus round 13 §7.2 (`r = 4..8`, primes `<= 400`): only `r <= 6` with the smallest primes was run.
  Untested: `r >= 7`, larger primes, and the joint core beyond `K = 4`.
* m7. One-modulus flexibility is reported at the level-1 relaxation. Its all-level validity is evidence (`lift2` was not
  run for every survivor). The flexibility claims used in the joint section are backed by all-level witnesses, which I
  verified.
* m8. The non-strict one-modulus class survivors at 41, 79 and 251 should be re-examined modulo local twists. This is
  only relevant to the old notion of strictness.

## 7. Corrected statement

**Class.** Rig(K) concerns parity-respecting label systems on the odd primes `q <= Q_K` not dividing `N` (all
levels). They meet the per-unit budget `m_(v,s) <= nu_s^(-1) e^(w(v)/2)` for every `v` with `1 <= ||v||_1 <= K` and
every unit `s`. This is exactly the norm-`<= K` part of P3-all\*. The symmetry group acting on good systems contains
`{+-1}^r x (Z/2)^r` (column conjugations and associate twists), together with `u = 1 mod 4` scalings.

**Proved (all quantifiers checked):**
* Theorem 2.5: only `F_K^ex` matters, the cores are monotone, and `q_core <= 4A_K + 5 = O_K(r^(floor(K/2)))`.
* Lemma 3.1 (slack identity) and Cor. 3.2.
* Lemma 3.3, Theorem 3.4 and Cor. 3.5: there is no F2 deformation at any modulus covered by an anticode of diameter
  `<= K-1`, and the joint transport inequality holds.
* Prop. 3.6(1)-(2) and Prop. 3.7.
* Theorem 4.1, **after correction**: "strict and coherent" must mean `ell(q^a) = u_q (eps∘L(q^a) + (m_(q^a)/2) t)`, with
  `u_q = 1 mod 4` and one `(eps, t)` on `Q_eta`; H^gen is read for the primes `i^(t_j) pi_j^(eps_j)`.
  * As written (`eps` only), RigS implies that the class is empty in the regime whenever some `q in Q_eta` has
    `q = 3, 5 mod 8` and `(p_j/q) = 1` (Prop. T, Lemma T'). Rig(K) can never yield it: the twisted genuine system is
    good at every `K` and non-strict at most moduli.
  * The size statement "`q_core >= M` at bounded `K`" holds only under the extra hypothesis `r >= M^gamma`
    (`K >= 2/gamma + 1`).
  * For actual circles alone the corrected theorem is vacuous. Its content is the class-level transfer of H from
    designed to genuine labels.

**Exact, reproduced independently (`ps = {5,13,17,29}`):**
* Forced sets for `K = 2..5`.
* One-modulus relaxed cores `q_core^1 = 7, 41, 79, 157` (`n = 2, 10, 20, 39`) for `K = 2..5`, together with every
  round-10 count compared.
* All 22 joint witnesses (all levels, unconditional). They give `q_core^(J,eps)(3) <= 37` and
  `q_core^(J,eps)(4) <= 61`.
* `Core^eps_J(3) ⊇ F_3 ∩ [3, 31]` by exhaustion.
* The `K = 6, 7` one-modulus figures were not re-checked here, except the spot check in §4a.

**Solver-verified only (HiGHS, float logs, no certificate), unless §4a says otherwise:**
* `37 in Core^eps_J(3)`;
* given genuine labels at all levels of 3 and 7: `Core^eps_J(4) ⊇ F ∩ ([11, 61] ∪ {73})`, and coherence on
  `F ∩ [11, 61]`;
* given genuine labels at 3, 7 and every forced modulus above 163: `Core^eps_J(5) ⊇ F ∩ [11, 157]`.

**Evidence only:**
* the anticode-threshold "law". It matches at `K = 5, 6` for `r = 4`, fails at `K = 7` and at `r = 5`, and is not
  implied by Theorem 3.4 because the rigid range is not covered;
* growth of the one-modulus core in `K` and `r` (`r <= 6`, smallest primes);
* the all-norm reference being more rigid (one-modulus only).

**Answer to Rig(K): undecided.** The one-modulus core grows on the tiny instances computed. Growth of the joint core
is **not** established even from `K = 3` to `K = 4` unconditionally (`37 <= q_core^(J,eps)(4) <= 61`). The strict joint
core as the write-up defines it does not grow at all: it is essentially empty at every `K`, because of the twists. The
corrected strict notion was not tested. No growth and no cap. `M <= (2/3 + eps) log R/loglog R` and the P3-all\*
bracket are unchanged.
