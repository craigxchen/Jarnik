# Round 29: Judge's README ("abstract, creative, crazy ideas")

## 0. Verdict (short)

* **Growth rate: unchanged.** `M <= (2/3+eps) log R/loglog R` (two referees) stands. No lane proved, and nothing here
  implies, any `M = o(log R/loglog R)` statement.
* **Uniform bound: unchanged (open).**
* **No idea passes both filters for the arithmetic problem.** Every arithmetic input proposed this round is either
  angle-only or residue-only, so fakes satisfy it (fails F1), or it is a single exactly critical inequality (fails F2).
  The only structure that does pass F1 is the moment surface `Sigma_L`, and that structure lives only in the
  function-field model, where it is open and has no arithmetic transfer.
* **One killed headline.** The clothes lane claimed that smooth moduli `q` defeat the 2/3 engine (Prop 4.4). That claim
  is false (kill_clothes, Lemma K). The residue-class twin obeys the full 2/3 bound for every odd `q` coprime to `N`.
* **Genuinely new and correct:**
  * the exact 8-point "cube" cluster over `C(t)`/`Q(i)(t)`/`F_{q^2}(t)` at `c = 0`, with its 8-point `Z[i][u]` template;
  * the self-similar decomposition of the neighbouring-circle pair count `E^W`;
  * the character-vanishing kill of the Li–Liu device;
  * the twin form of Theorem 1.1;
  * several exact small-case values in the function-field model.
* **Cross-lane correction (judge, unrefereed, but checked exactly in §6).** kill_moon's revised question ("is
  `sup_D M_R(0)` finite, and does it equal 6?") and its "about 6" heuristic are **refuted by ff's Theorem C**. The two
  lanes use the same model, and the cube family gives 8 points at `D = 8`, `c = 0`.
* **Dedicated follow-up: at most one,** the function-field decision problem (R_4) (§4). It is a calibration of the
  model, not a route to the growth rate.
* **Round score: 2/10 on the goal** (§5). The round was honest and well checked, and the creativity was real. Nothing
  moved.

Labels: **proved** (proof checked by the lane's killer), **exact** (finite exact computation), **evidence**,
**heuristic**, **open**, **unrefereed** (the judge's own items, marked in place and collected in §6).

## 1. Inputs and what was checked

| lane | write-up | killer report | killer verdict | novelty (killer) |
|---|---|---|---|---|
| ff (function fields) | `ff.md` | `kill_ff.md` | survives (obstruction), with corrections | repackaged barrier + one new exact object |
| spectral (Li–Liu template) | `spectral.md` | `kill_spectral.md` | survives (obstruction), with corrections | repackaged barrier |
| clothes (other guises) | `clothes.md` | `kill_clothes.md` | **does not survive** (headline Prop 4.4 false) | repackaged barrier |
| moon (18 moonshots + FF) | `moon.md` | `kill_moon.md` | survives (obstruction), with corrections | repackaged barrier + small calibration |

**What I read.**
* Every kill report: in full.
* The lane write-ups: verdicts, the key theorem statements and the filter sections.
* The expert note §1. The round-28 README §0.

**My own checks.**
* An exact re-check of Theorem C (`judge_checks/cube_judge.py`, §6).
* A comparison of the ff and moon model definitions (§6).

**The uploaded paper.** The Li–Liu conditional Goldbach paper was read in full by spectral and checked against the PDF
by kill_spectral. Two lanes found by WebSearch that it is not available online; the same authors' arXiv 2606.05224 is
a different paper.

## 2. Killer-corrected statements per lane

### 2.1 ff — the problem over `C(t)` (score 4/10)

Setting: monic `z` of degree `D` with `z zbar = N`, and `deg(z_i - z_j) <= K = floor(D/2+c)`. Put `L = D-1-K`.

1. **(proved; already in the record as `families.md` Lemma P)** The arc condition is `L` exact real-linear equations
   `sum_k (a_k - b_k) Im rho_k^m = 0`, for `m <= L`. Distinct points are at `l1`-distance `>= L+1`.
2. **(proved)** `M(D,c) <= 2^{2c+2}(D+2)`.
   * Function-field fakes (arbitrary columns in `R^L` in place of the moment vectors) reach `>= D+2c+2` points for
     `D = 2^k` and `0 <= c <= D/2-k-2`.
   * The boundary case `c = D/2-k-1`, which includes `D = 8`, `c = 0`, is **unproved**: sign normalisation makes
     columns coincide there.
3. **(proved, exact; new)** The cube family. Take `y in F^4` with `sum y = 0` and `prod y = 1`, and put
   `r_k = eps_k/y_k + i y_k`.
   * It gives exact 8-point clusters at `c = 0` for every `D ≡ 0 (mod 8)`, over `Q(i)` and `F_{q^2}`.
   * Hence `sup_D M(D,c) >= 8` for every `c >= 0`.
   * It also gives the first 8-point polynomial template in `Z[i][u]`. Its constant `C ≈ 3.18e9` makes it a new
     mechanism, not a record.
4. **(proved; the record's free-curve argument, `gauss.md` Thm 4.3(ii))** A dominating family of rank-`r` clusters
   satisfies `r(L+1)+1 <= 2D`. The proof needs complexification and the corrected full-rank case.
5. **(elementary; the twin of expert-note §5 item 8)** `sup M < infinity` iff the code rank is bounded.
   * **(R_4)**, rank `<= 4` for large `D`, would give `M <= 8` via an exact LP. (R_4) is **open**.
   * `limsup M = 8` is a **heuristic**.
6. **(exact)** Squarefree values:
   * `M(D,0) = 4, 3, 6, 3, 5` for `D = 2..6`;
   * `M(7,0) >= 4` (from the record);
   * `M(7,1) >= 8` and `M(8,1) >= 14` (killer, vertical-line roots).
7. Arithmetic transfer: **negative**. The `D/2` independent equations come from the infinite residue field. A `p`-adic
   arc gives one congruence, and Fermat quotients (Buium delta-characters) are residues that fakes already satisfy.

### 2.2 spectral — the Li–Liu template (score 3/10)

1. **(checked against the PDF)** Li–Liu's Hypotheses S and R are, by the paper's own (2.28), (2.30) and Thm 6.1,
   operator-norm bounds on `T_N`. They are strictly stronger than the Goldbach asymptotic. GRH supplies only
   holomorphy, and the GRH norm bound `(log N)^2` is weaker than the trivial `log N`.
2. Transplant **H** (Hecke angle characters) is **fake-blind** (proved). Here
   `lambda_k(N) = 4 prod_p U_{e_p}(cos 4k phi_p)`, with a factor `(-1)^{k e_2}` at 2.
3. Transplant **Q** (auxiliary modulus on a strip of width `W`) is Parseval for the pair count `E^W >= M(M-1)` on
   neighbouring circles.
   * For `W >= 5R`, the supremum over non-trivial positions is `>= 2 floor(min(C/(2 sqrt5),1) sqrt R)` (proved).
   * The position average is of order `(W/R) C^2 log R` (exact counts; heuristic constant `2/pi^2`).
   * **(proved; new, killer)** Self-similarity: `E^W` is exactly a sum over `P'` of short-arc counts of the same problem
     on circles of radius `R/|P'|`, with constant `C/sqrt|P'|`. So any hypothesis bounding `E^W` uniformly is provably
     a restatement.
4. Transplant **A** (Li–Liu verbatim on `a^2+b^2 = N`) is exactly diagonal on primitive points for `R > 16C^2`
   (proved).
5. The scale obstruction (Prop 2.0) holds **only for statistics of bounded degree in `M`**. An `eps`-loss that is
   uniform in a growing degree would not be vacuous, but no such input is known.
6. The Li–Liu residual-kernel numerics are thin evidence: only four points lie inside the Hardy strip.

### 2.3 clothes — the problem in other clothing (score 2/10; killed)

1. **(proved; dictionary)** Take `q` odd with `gcd(q,N) = 1`.
   * A residue class `{z : z ≡ z0 mod q}` on the circle is a coset of `ker(a -> prod u_p^{a_p})` in the exponent box.
   * `T(N,q) <= max_f r_f(N) <= 2^{omega(q)} T(N,q)`, where `f` runs over the primitive forms of discriminant
     `-4q^2`. The killer proved the converse direction.
2. **(proved; corrected)** For **every** such `q >= lambda N^{1/4}`, a twin class of size `M` satisfies
   `3M log M <= log N + (14 + 4 log^+(2/lambda)) M`.
   * Prop 4.4 ("smooth `q` leaves the pigeonhole step with nothing") is **false**, by the killer's Lemma K: confinement
     to one class mod `l^t` restricts the residues mod `l^{t+b}` to `l^b` classes.
   * The lane's Lead (the smooth-conductor regime) therefore has a false premise.
3. Twin boundedness is **open**. It is blocked by F1 (group fakes) and by F2 (`alpha = 1/4` is the Coppersmith /
   Lenstra critical point). No implication between the twin and the arc problem is known. The twin as "Lenstra at
   `alpha = 1/4`" is already in `round4/hardness.md`.
4. Peeling gives only `M < 2W/(log P_+ - 4 log^+ C) + O(1)`, which is the record's Thm 3.10. Every fake is
   `M^{O(1)}`-smooth, so the smooth/rough dichotomy route is closed (proved).
5. Small errors: the CC exponent at `k = 4` is 2/5, not 3/8, and the twin evidence is thinner than stated.

### 2.4 moon — 18 moonshots and a function-field development (score 3/10)

1. Each of the 18 screens reaches its stated conclusion, with these corrections:
   * the Erdős–Turán "vacuous for every `N`" proof is invalid; the conclusion is shown only on the five record clusters;
   * `s_m` carries a factor 4;
   * idea 17 is misreported: the 7- and 5-point record clusters sit **one collision above** the forced minimum at
     `l = 3`.
2. **(proved)** Over `F_q(t)`, `q ≡ 3 (mod 4)`, under hypothesis (H) (no small-degree prime factors of `N`):
   `M <= (1/2+o(1)) log_q R / log_q log_q R`. This is a direct transfer of Theorem 1.1's mechanism.
3. **(proved)** Over `R(t)` at `c = 0`:
   * `M <= 2D` for even `D` and `M <= D+1` for odd `D` (Plotkin, squarefree). Hadamard profiles make this sharp for
     local data.
   * `M = 6` exactly at `D = 4`, via de Bruijn–Tengbergen–Kruyswijk.
   * Sperner blocks give 6 points at `D = 8` (proved) and at `D = 12` (evidence).
4. **(proved; new, killer)** The Goldbach-device kill is sharper than moon stated. Every lattice point `z` divides `N`
   in `Z[i]`, so every character mod `N` vanishes identically on the circle.
5. **Judge's cross-lane correction (unrefereed; §6).** moon's "generic count predicts 5", and kill_moon's revised "about
   6, does it equal 6?", are both **false at `D = 8`**. ff's Theorem C gives an 8-point cluster in the identical model.
   So: `sup_D M_R(0) >= 8`, the record's integer-box value at `D = 8` (4 points) was far from optimal, and the open
   question is the finiteness of the supremum (ff's (R_4)), not its being 6.

## 3. The filters, and whether anything passes

| idea | F1 (fakes) | F2 (criticality) | arithmetic relevance |
|---|---|---|---|
| FF moment-surface constraint (ff, moon §4) | passes **in the FF model** (fakes = arbitrary columns; `Sigma_L` separates them) | single configuration: unproved; families: Thm D (non-critical only because it concerns families) | none: the transfer is negative (residue field) |
| Li–Liu / spectral (H) | fails (proved fake-blind) | Cauchy–Schwarz, exactly critical | none |
| Li–Liu (Q), neighbouring circles | not fake-blind, but uses only `E^W >= M(M-1)` | self-similar restatement (proved) | restatement |
| Li–Liu (A), characters mod `N` | degenerate (diagonal; characters vanish) | — | none |
| residue twin (clothes) | fails (group fakes) | critical (`alpha = 1/4`) | parallel problem, no implication |
| peeling / smooth–rough | fails (fakes are smooth) | critical | Thm 3.10 again |
| moon's 18 screens | 15 fail F1 or F2, idea 10 (Rédei) has no testable input, idea 17 is weak data | — | none |

**Answer: no idea passes both filters for the arithmetic problem.** Structures the record had never used were all
visited:
* the circle outside the cluster: the whole circle is a product set, i.e. one-point Liouville;
* neighbouring circles: self-similar;
* derivations: Thm D, which concerns families;
* non-archimedean twins: the residue twin, which carries the same 2/3 bound;
* group structure: kernel cosets;
* L-functions of `N`: fake-blind, or vanishing.

Each unused structure either reduces to pair or residue data or reproduces the problem itself. This confirms
expert-note §5 from new directions, but it adds no new mechanism.

## 4. Follow-up (at most two; I name one)

**Follow-up 1 — decide (R_4) in the function-field model.** Over `C(t)` at `c = 0`, is the rank of the 0/1 kernel code
of a cluster bounded by 4 for large `D`, and hence `M <= 8`? Or do rank-5 configurations exist?

* *Why this one.* It is the only question this round where the filters bite in our favour. Fakes are arbitrary columns,
  so any answer must use the moment surface `Sigma_L`, which is real non-fake structure. It is also exactly the
  function-field shadow of expert-note §5 item 8 (bounded rank ⇔ uniform). A rank-5 example would suggest that the
  conjecture is not "geometric". A proof would be the first exceptional-set theorem in any model of the problem.
* *Concrete first step.*
  * Run an exhaustive or homotopy search at `D = 9, 10, 12`, `c = 0`, for 5 independent 0/1 vectors of weight
    `>= L+1` in the moment kernel, over the cube locus and over Sperner-block / mixed-sign generalisations of it.
  * In parallel, try to prove `M <= 8` on the 2-parameter cube locus alone, by excess intersection on that surface.
  * Pre-register a stop rule: no rank-5 example at `D <= 12` and no proof on the cube locus ⇒ stop.
* *Honest probabilities.* (I would give this follow-up 2/10 if judged on the arithmetic goal.)
  * About 15% that a dedicated round settles (R_4) in either direction for the model.
  * Below 2% that it changes anything about the arithmetic growth rate. The transfer is blocked by the
    residue-field obstruction (§2.1 item 7), and nothing this round weakens that block.

**Not recommended:**
* the residue twin: its premise is killed, and it is parallel to the arc problem;
* bounds on the spectral `E^W`: provably a restatement;
* the `F_q(t)` (H)-removal: a routine transfer that cannot affect the integers;
* further Li–Liu-type conditional templates: S and R are norm restatements, and characters mod `N` vanish on the
  circle.

## 5. What it means, and scores

* **Growth rate:** unchanged, `(2/3+eps) log R/loglog R`.
* **Uniform bound:** open, and unchanged. The conditional routes remain Vojta-type (Questions A and B, expert note §§2–3).
* **The uploaded Goldbach paper:** a correct but formal reformulation. Its hypotheses are operator-norm bounds that are
  stronger than the conclusion. Its device of characters to the modulus `N` has no non-trivial analogue on
  `x^2 + y^2 = N`. No transplant of it improves anything.

| lane | score | reason |
|---|---|---|
| ff | 4/10 | One genuinely new exact object (the cube family) and a clean open question (R_4). Thm D and the linearisation were already in the record, and Thm B has a boundary gap. |
| spectral | 3/10 | A careful, correct reading of Li–Liu, and three decisive negatives. Prop 2.0 was overstated. The best new item (self-similarity) came from the killer. |
| moon | 3/10 | Broad and mostly correct screening, plus useful function-field calibration. One false data claim, one invalid proof, a wrong lead heuristic. |
| clothes | 2/10 | A good dictionary, but the headline new negative is false and the rest is Thm 3.10 and `round4/hardness.md` again. |
| **round** | **2/10** | Creative, honest and well killed. No growth improvement, no idea passes both filters, and one follow-up of model interest only. |

## 6. Judge's own items (unrefereed)

1. **Exact re-check of Theorem C** (`judge_checks/cube_judge.py`, sympy over `Q(i)`).
   * Input: `y = (-9/4, -2/3, 1/4, 8/3)`, with `sum y = 0` and `prod y = 1`.
   * Over all 256 transversals of the 8 root pairs, the ball of radius `deg <= 4` around `Q` contains **exactly 8**
     points.
   * The maximum pairwise difference degree is 4.
   * `N = Q Qbar` is real and squarefree, of degree 16.
   * This agrees with ff and kill_ff.
2. **Model identity.** moon §4.1 and ff §1 both use:
   * monic `z` of degree `D` with `z zbar = N`;
   * pairwise `deg(z_i - z_j) <= floor(D/2+c)`;
   * conjugation on coefficients.

   Over `R(t)`, moon's `B = C[t]`, which is ff's `C(t)` setting. Hence ff's `M(8,0) >= 8` applies to moon's
   `M_R(0)`. That refutes kill_moon's "about 6" heuristic and the "= 6?" form of the open question. It is consistent
   with Plotkin (`M(8,0) <= 16`) and with ff's LP cap of 8 under (R_4).
3. **Consistency of small values.** moon's `M = 6` at `D = 4` equals ff's `M(4,0) = 6`. moon's Sperner-block 6 at
   `D = 8` lies below ff's 8.

Files:
* this README: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round29/README.md`;
* the judge's script: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round29/judge_checks/cube_judge.py`.
