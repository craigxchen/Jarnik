# Round 31, Part X: judge report on finite cutoff and the coordinate-ray tests on one parent

Source: Part X of `astra_free_X/answer.md` (heading "# Part X. Finite cutoff and the coordinate-ray tests on one
parent", line 2317, to the end, line 2498). Also read: Astra's summary `final10.txt`, the coordinator's question
`reply10.txt` (question (B), P3-all\*, was the one she took), her scripts `x_height_check.py`, `x_digit_check.py`,
`x_local_check.py` (read, not imported). Record: `round31/README_IX.md` (format and the previous rulings),
`round5/round6_sharp/sharp.md` §1.1 (Def. 1.1, classes `F^char_X`, `F^all_X`), `expert_note/expert_note.md` §1 and
§5, `round14/c1construct.md` Lemma P0 and Proposition P.

Six referees each took one block. **All six reports are present; no referee failed.**

| key | report | scope |
|---|---|---|
| finite | `ref10_finite.md` | X.1 (finite rigidity, direction pinning, dense families, scope of the height parameter); X.5 least-cutoff table |
| rays1 | `ref10_rays1.md` | X.2 first half: statement, coarse profile, class encoding, local construction, digit covering (X.1), index bound (X.2), angular completion, word count |
| rays2 | `ref10_rays2.md` | X.2 second half: which words are tested, the deliberate failure, auxiliary heights (X.3), limit of the order argument; X.5 negative control |
| inverse | `ref10_inverse.md` | X.3 (finite-lift dichotomy, mod-four bound on owner twists) |
| weight | `ref10_weight.md` | X.4 (weight bound, restricted rate, lcm lemma) and the calibrations paragraph of X.5 |
| meaning | `ref10_meaning.md` | analyst: what Part X settles; genuine collision structure; exact computations on the explicit parents `q_P` |

The judge's scripts are in `judge_X/` (own exact code; nothing of Astra's or the referees' imported):
* `x1_judge.py` → `x1_judge.out`: order lemma brute force; the X.1 mechanism on 5,529 random nonintegral lifts on two
  parents; the X.5 row `P = 100` rebuilt from scratch (cases, kill test by exact rational rectangles, survivor
  verified at `phi = arg(2+i)` for every exponent).
* `x3_judge.py` → `x3_judge.out`, `x3_judge_q31.out`: the X.3 parity lower bound on the IX.3 Legendre design, and
  explicit integral `V` with `D V = (H - D)/2` attaining it, by an own Hermite-form integer solver.
* `x4_meaning_judge.py` → `x4_meaning_judge.out`: residue tuples of genuine factors at inert primes (`D(x)`,
  `F_2`-rank, divisibility of `Im(pi_j conj pi_k)`); a bounded-rank genuine parent; full enumeration of the lattice
  points of `q_P` for `P = 100, 200`.

**Every claim labelled [J] is the judge's own and has not been refereed.**

Definitions used throughout. As in `README_IX.md`. In particular: a *label* is a physical prime `p_j = 1 mod 4` of
the parent `N = prod p_j^{e_j}`, `W = log N`, `s` = number of labels, `m` (or `M`) = number of marks (points). A
*witness* (factored residue datum) at label `j` is a family of units `r_j(q) mod q^{A_q}` of norm `p_j` at the odd
primes `q <= X`, `q != p_j` (`sharp.md` Def. 1.1(d)). For `v in Z^s \ {0}` and `s' in Z/4` the **e-all test** is
(1.2): `|sin(<v, phi> - s' pi/4)| >= nu_{s'} m_(v,s') e^{-w(v)/2}`, where `nu = 1` (`s'` even) or `2^{-1/2}` (`s'`
odd), `w(v) = sum |v_j| log p_j`, and `m_(v,s')` is the product over odd `q <= X`, `q` not dividing `N`, of the
largest `q^a <= X` with `prod_j (r_j/conj r_j)^{v_j} = i^{s'} mod q^a`. The *ray* of `v` is the set of tests on all
`n v`, `n >= 1`. **P3-char**/**P3-all** are `F^char_X`/`F^all_X`. **P3-all\*** is the record's open question whether
`W_min^all = O(M log M)` at cutoff `X = 3M`, or `3M^A` with `A` chosen after the weight constant. A *rational lift*
of a witness at label `p` is `r = pi gamma/conj(gamma)` (`Norm pi = p`, `gamma in Z[i]`) whose reductions agree with
the witness at the primes considered; its *twist norm* is `Norm gamma`. A *signed word of mass `h`* is `v` with
`||v||_1 = h`. For a partition of the labels into classes, the *class-sum vector* of `v` is
`sigma(v)_kappa = sum_{class(j) = kappa} v_j`.

## 0. Verdict in five lines

1. **No mathematical error in Part X.** Every proof step passed its referee; the judge rechecked X.1 line by line,
   rebuilt the X.5 row `P = 100` (least cutoff 439, survivor `985 - 1158i`) and the X.3 bound (exact minima 20, 72,
   110 at `q = 11, 19, 23`) [J]. Corrections are of wording and scope: one sentence of X.2 is false as worded (§1.2).
2. **P3-all\* is not settled in either direction.** X.2 is a partial fake: one parent sequence at
   `W ~ (1536H^2/rho) m log m` passes P3-char and e-all at every cutoff on every power of every word of mass `<= H`
   whose class-sum vector (for `b(m) -> infinity` constructor-chosen classes) is nonzero. Open: the balanced words
   inside a class, in particular the label differences `e_j - e_k` that a pigeonhole argument picks from the data.
3. **The height route to a lower bound is dead.** X.1 (proved, new, effective): on a parent with `T = max v_2(l-1)`,
   two ray tests at `n = 2^{T-1}` force any lift of twist norm `<= H` to be genuine once `X >= sqrt2 (pH^2)^{n/2}`
   (polynomial on 5 mod 8 parents). But X.2 (5 mod 8 version) shows passing data need not have such lifts, and X.3
   shows the natural exact witnesses have log twist `>= (1/2 - o(1)) sum log p`. Lift height is never charged to `W`.
4. **X.4 is a reformulation without content for genuine circles** (its bounded-rank hypothesis fails for every
   genuine parent with polynomially bounded labels; computed rank `= pi_3(x)` [J]). The residue route is closed for
   internal work, on the merits and in substance by the previous judge's rule.
5. **For the integer problem the growth rate `(2/3 + eps) log R/loglog R` and the uniform bound are unchanged.**
   New data (meaning, spot-checked [J]): on `q_P` (all split primes in `[P, 2P]`) no two lattice points lie within
   `8 q_P^{1/4}` for `200 <= P <= 450`. Score on the goal: **2/10**.

## 1. Referee-corrected statements per block

Labels: **Proved** = a referee checked a complete proof line by line, or an exact computation was reproduced with
independent code. **Computed**, **Heuristic**, **Gap**, **False**, **Withdrawn** as usual. **No referee found an
error in any proof step of Astra's.** The verdicts FALSE in the structured summaries concern one sentence of X.2
(rays2), framings of X.4 (weight) and the question "does Part X answer P3-all\*" (meaning), which Astra does not
claim. The one GAP (finite, item 8) is the absence of any consequence for `W` from X.1, which Astra also states.

### 1.1 finite: X.1 and the X.5 table (referee: CORRECT on X.1; two CORRECT WITH CORRECTIONS)

**Proved (referee line by line; judge line by line and computed [J], §1.7(a)).**

* **X.1, finite rigidity, quoted in full (referee wording):** *Fix a parent `S` of split primes, a label `p in S`, a
  height bound `H`. Put `T = max_{l in S} v_2(l - 1)` and `n = 2^{T-1}`. Let `r = pi gamma/conj(gamma)`,
  `Norm pi = p`, `Norm gamma <= H`, and suppose the witness at `p` equals the reduction of `r` at every odd exterior
  prime power `<= X` prime to the reduced denominator `b` of `r` (arbitrary residues at the primes of `b`). If
  `X >= sqrt2 (pH^2)^{n/2}` and the two e-all tests on `n e_p` at the units `i` and `-i` hold, then `b = 1` and
  `r in {u pi, u conj(pi) : u in {1, i, -1, -i}}`, for every real `phi_p`.* The same holds with `n = 2^k` for any
  `k` such that no `l in S` prime to `J` has `ord_l(alpha/conj alpha) = 4n`; some `k <= |S|` always works (lift-
  dependent), so `X >= sqrt2 (pH^2)^{2^{|S|-1}}` is a uniform cutoff.
* *Mechanism.* `r = alpha/b` in lowest terms, `J = Norm alpha = p b^2`, `b <= H`, `b` and `J` odd, `alpha`,
  `conj alpha` of disjoint support. For `alpha^n = x_n + i y_n`, the integers `x_n - y_n` (charged at unit `i`) and
  `x_n + y_n` (unit `-i`) are odd, nonzero, prime to `J`, of size `<= sqrt2 J^{n/2}`, and prime to every parent
  prime `l`: a divisor `l` would give `(alpha/conj alpha)^n = +-i` in the cyclic group of order `l - 1`, which for
  `n` a power of two means order exactly `4n = 2^{T+1}`, impossible. So both are fully charged, and squaring and
  adding the two sine inequalities gives `1 >= J^n/p^n = b^{2n}`.
* **Direction pinning (corrected precision):** *if `min S > sqrt(2 max S)` and `X >= sqrt(2p)`, a genuine lift
  `r = x + iy` passing the singleton tests `e_p` at units `+-i` forces `phi = theta (mod pi)` or
  `phi = pi/2 - theta (mod pi)`, `theta = arg r`.* Astra's "`phi = +-arg r (mod pi/2)`" is true but weaker. This is
  round 14 Lemma P0 / Proposition P in its unpolluted case. For a nonintegral candidate with `sqrt(2p) b < min S`
  the same `n = 1` argument applies at `X >= sqrt(2p) b`.
* **Dense families:** *on the product of all primes `5 mod 8` in `[P, 2P]` (`T = 2`, `n = 2`), for every label and
  every `H <= P^h`, the cutoff `X >= 2 sqrt2 P^{1 + 2h}` forces every lift of twist norm `<= H` genuine and pins
  every angle (`P > 4`). With `v_2(l - 1) <= T_0` the relative density among split primes is `1 - 2^{1-T_0}`
  (asymptotic in `P`) and the cutoff exponent is `2^{T_0-2}(1 + 2h)`.*
* **Scope of the height parameter:** the coordinator's literal conclusion `gamma/conj(gamma) = +-1` is false even on
  genuine circles (`gamma = conj pi` gives `conj pi`; `gamma = 1 + i` gives `i pi`); the correct conclusion is
  `r in {u pi, u conj pi}`. This corrects the question, not Astra.
* **Conditional consequence (referee, Proved):** *on a comparable 5 mod 8 parent, a member of P3-all at cutoff
  `X >= 2 sqrt2 P^{1+2h}` in which every label's witness is the reduction of a lift of twist norm `<= P^h` is an
  actual cluster* (Proposition P of round 14, whose hypothesis (i) holds because `|x +- y| <= sqrt(2p) < min S`),
  so `W >= 3M log M - O(M)`.

**Contrapositive and meaning (finite, Proved; status of the GAP).** A nongenuine datum passing the two tests at `X`
has every rational lift with `b > ((X/sqrt2)^{2/n}/p)^{1/2}`. Nothing follows for `W`: lift height is not part of
the weight, and every finite datum has a lift of height `e^{O(X)}` (README_IX §1.8). X.1 excludes only fakes whose
witnesses have polynomial-height lifts, on parents with bounded `T`. The GAP is "X.1 gives a lower bound for
P3-all\*", which Astra does not claim.

**Computed (finite, reproduced; row `P = 100` reproduced by the judge [J]).**

| parent in `[P, 2P]` | `P` | labels | nonintegral cases | least cutoff | survivor at `X - 1` | proved bound (*) |
|---|---:|---:|---:|---:|---|---|
| all split | 100 | 10 | 190 | 439 | label 181, `b = 113`, `985 - 1158i`, `phi = arg(2+i)` | about `10^111` |
| all split | 500 | 36 | 2,556 | 3,881 | label 997, `b = 557` | about `10^576` |
| all split | 1,000 | 67 | 8,911 | 8,681 | label 1601, `b = 1381` | about `10^317` |
| 5 mod 8 | 1,000 | 36 | 2,556 | 5,783 | label 1621, `b = 1493` | about `1.1 * 10^10` |

**Corrected wording (table).** (a) *"Least cutoff" means: the least `X` at which no nonintegral lift with
`Norm gamma <= 2P` whose reductions are defined at every exterior prime `<= X` (so `b` is a parent prime) passes
(1.2) on all `n e_p` at all four units for any real `phi`.* That class is narrower than X.1's hypothesis, which
allows exterior primes in `b`. With zero charge at the primes of `b` (most favourable to survival, not shown
realisable), two lifts on the first parent (`b = 185`, `b = 145`) survive at 439 and die at 769 and 443 (Computed,
finite). (b) Almost every case dies at `n = 1`; the table is consistent with X.1 but does not exercise its dyadic
mechanism, and (*) is not sharp (Astra: "sufficient, not optimal").

### 1.2 rays1 + rays2: X.2 (rays1: CORRECT WITH CORRECTIONS; rays2: one sentence FALSE as worded)

**X.2, corrected statement, quoted in full (rays1 wording, with the rays2 correction of the last paragraph):**
> *Fix `C > 0`, an integer `H >= 2` and `rho in (0.99, 1)`. Put `B = 64H^2`, `g = 2B = 128H^2`, `d = 3`. Fix any
> `b(q) -> infinity` with `U(q) := 100 g B^{T(q)} = o(log log q)`, where `T(q) = b(q)` (one-hot rows) or
> `ceil(2H log_2(2b(q))) + 2` (random detecting rows); this forces `b <= (log log q)^{1/(2H log_2(64H^2)) + o(1)}`,
> e.g. `(log log q)^{1/32}` at `H = 2`. Then there is `q_0 = q_0(C, H, rho, b(.))` such that for every prime
> `q = 3 mod 4`, `q >= q_0`, there is a tuple `(p, e, a, phi, theta_0, delta, k, r)` with:*
> * *profile: the `g`-fold replicated Paley profile of IX.3/IX.5, `s = g(q-1)` physical primes in `[P, 2P]`,
>   `P = q^3`, in one sector and one radial bin, `m = (rho/2 + o(1)) q` marks, `W/(m log m) -> 1536 H^2/rho`,
>   `m/s -> rho/(256 H^2)`;*
> * *the tuple is in `F^char_X(m, C)` for **every** `X >= 3` (same-level clause included), and satisfies (G);*
> * *for a partition of the labels into `b(q)` classes (balanced, or with boundedly many reserved singletons), every
>   `v` with `1 <= ||v||_1 <= H` and `sigma(v) != 0`, every `n != 0`, every unit and every `X`: (1.2) holds for `nv`.
>   This includes every `n e_j`, every `e_j + e_k`, every `e_j +- e_k` across classes, every word of odd mass;*
> * *the angles are nongenuine, and by IX.7 no label's witness is cofinitely the reduction of one element of
>   `Q(i)^*`.*
>
> *The words of mass `h <= H` not covered are exactly those with `sigma(v) = 0`, the short vectors of the lattice
> spanned by same-class differences: a share `O_H(1/b + 1/s)` (about `b^{-h/2}` for even `h`, none for odd `h`);
> among differences, `sum_c C(n_c, 2)/C(s, 2) <= 1/b + o(1)`. The tuple is not claimed to lie in `F^all_X`.
> **Capped** are arguments that use only P3-char, (G) and e-all on a set `E_m` of words of mass `<= H` **fixed
> before the data**, such that some partition into `b(m)` classes gives every word of `E_m` a nonzero class sum.
> For difference tests this means a `b(m)`-colourable test graph, which by Turán has at most a `1 - 1/b(m)` share
> of the pairs.*

**False as worded (rays2; rays1 concurs in substance).** Astra's sentence "This is a counterexample to a growth
argument requiring even this fraction of the tests" (line 2356) is false: a test graph using a share
`1 - o(1/b(m))` of the differences, or containing a clique on `b(m) + 1` labels, is not capped, and a pigeonhole
argument that picks equal-residue pairs from the data (X.4, `sharp.md` §6.7) is not capped at all, because the
construction places every uncontrolled coincidence inside a class (rays1 computed 67,706 untested words with
`D | sum v_j L_j`). The corrected statement is the last paragraph above. [J: I uphold this. Density is the wrong
measure; the relevant property is that the word set is data-independent and `b(m)`-separable.]

**Proved (rays1, line by line, and Computed).** The digit-sum covering (X.1): attainment of every block sum
(digit sums of `g >= B` slots fill an interval of `g + 1` integers; rounding loss `< 4g` against slack `>= 100g`),
digit windows `B^{t-1}L/D mod 1 in [a_t/B, (a_t + 2)/B]` (needs `z_u <= 16H + 1 <= B - 1`), and `D` not dividing
`sum v_j L_j` for tested `v` (error `2H/B < 1/(4H)`; `B > 8H^2` would do). The index bound (X.2):
`gcd(sum v_j L_j, d_a) | G = lcm(1..ceil U)`, so a collision of `n v` with any unit at `l^a` forces
`d_a | 4nG`, `l^a <= 6nG`, colliding product `e^{O(nG)}` uniformly in `X`. The angular completion
(Haar on block-sum fibres; union bound `O_H(q^{-1/2+o(1)})`). Single labels pass already with one class. Computed:
2,631,680 digit-interval and 6,985,602 tested-word checks; 8,068 tested words x 12 powers x 4 units with actual
Gaussian residues mod `l^a` (`l = 3, 5, 7, 13`); 0 failures. Astra's script reproduced her counts.

**Quantifiers (rays1).** `C, H, rho` first; then `g, B, d` and the weight constant `12g/rho`; then `b(.)`; then per
`q`: marks, primes, coarse data, classes, local data at each exterior `l` independently (compatible levels by
König's lemma), fine angles. `q_0` is not effective (`log q` must exceed `e^U`, `U >= 3.4 * 10^9` at `H = 2`). The
constant grows like `H^2`: words of unbounded length are not capped by one sequence.

**The deliberate failure (rays2: CORRECT WITH CORRECTIONS), corrected statement:** *the X.2 instance **can be
chosen** so that one same-class difference `e_j - e_k` (labels from distinct blocks with equal Legendre parity at
an exterior prime `l_* in (4P, 8P)`) has equal logs throughout the tower at `l_*`; then the unit-one test on
`e_j - e_k` fails for every angle at every cutoff `X >= l_*`, while all positive statements of X.2 still hold.*
Summary line 7 "can be forced to fail" must read "can be chosen to fail". The failing cutoff is `>= 4q^3`, of order
`m^3`: regime `3M^A` with `A >= 3`, not the P3-all\* cutoff `3M`. *Heuristic (rays2, not checked by the judge):*
ties at four exterior primes in `(q, 1.48q)` should move the failure below `3M`. It is a choice, not a forced
failure, so it says nothing about whether the omitted edges can be repaired.

**Auxiliary heights (X.3 of the text; rays2 CORRECT, sharper form), quoted in full:** *fix the 5 mod 8 instance
(all physical primes `5 mod 8`). For every label `p_j`, every `X >= sqrt(2 P_max)` and every `r` with
`r conj r = p_j` agreeing with the witness at all exterior odd prime powers `<= X` prime to its reduced denominator
`b_j`: `b_j > (X/(sqrt2 p_j))^{1/2}`. Hence `Norm gamma > (X/(sqrt2 p_j))^{1/2}` for every representation
`r = pi gamma/conj gamma`.* Since such lifts exist at every `X`, the least lift height diverges with `X` on one fixed
parent with fixed weight. Computed (rays2): 19,325 nonintegral lifts on the 36 primes `5 mod 8` in `[1000, 2000]`,
all X.1 identities exact.

**Limit of the order argument (rays2, Proved):** a common index bound `G` for all differences forces
`log G >= psi_odd(s - 1) ~ s`. It blocks only the common-index method, not pair-by-pair constructions.

**Computed (rays2, reproduced): negative control** 101, 109 at 1511: 1511 prime, `3 mod 4`, first prime `>= 1500`
with `(101/l) = (109/l)`; equal ratio classes; the test needs `|sin| >= 14.40`. Astra's five-mark bookkeeping
(90 / 480 / 56,440 checks) was not re-run (her script imports a missing module `ix_factor_check`): **not
reproduced**.

### 1.3 inverse: X.3 (referee: CORRECT WITH CORRECTIONS, scope only)

**X.3, corrected statement, quoted in full (inverse wording):**
> *Fix `q = 3 mod 4`, the Legendre profile (`K_ik = chi(k - i)` on `I = F_q^*`, `a_ik = 1 + K_ik`), `S` in the
> positive class of column 0, `m = |S|`, split primes `p_k in [P, 2P]` with factors `pi_k`, units `u_{ii'}`, and
> then a cutoff `X`. Let `r_k in Q(i)^*` satisfy (Loc_X): at every odd prime `l <= X` not a parent prime, all `r_k`
> are `l`-units, `r_k conj r_k = p_k`, and `prod_k (r_k/conj r_k)^{a_ik - a_i'k} = u_{ii'} gamma_i/gamma_{i'} mod l`
> for all `i, i' in S`. Write `r_k/pi_k = u_k G_k/conj G_k` (`G_k` conjugate-primitive), `H_k = |G_k|^2`,
> `T = sum log p_k`, `B = sum (log p_k)^2`, `L_X = sum log l` over those `l`. Then either (1)
> `sum_k |d_k| log H_k >= L_X - (q+3)/2 log 2P - log 2` for some pair `d = a_i - a_{i'}`, or (2) the pair identities
> hold in `Q(i)^*`, the parent valuations `V_jk = v_{pi_k}(r_j/pi_j)` satisfy `D(I + 2V) = H`, and the **owner**
> twists satisfy `sum_{i in S} log H_i >= (1/2)(mT - sqrt(mqB) - m log P_max)`, with all but
> `qB/(2 eta T - log P_max)^2 = O(eta^{-2})` owners having `log H_i >= (1/2 - eta) T`. The same bound holds for
> `|V| + |E|` if the targets are altered by even parent valuations `2E` (Remark R). It concerns owner rows only and
> says nothing about e-all passing or failing, nor about witnesses that are local at the cutoff.*

**Proved (inverse, line by line).** `||d||_1 = q + 1`, `d = h mod 2`; the case-1 bound is "a nonzero integral
element of bounded norm divisible by every exterior `l <= X`" with no gcd discarded; case 2 uses only parent
valuations (exterior twist factors allowed); the mod-four relation `V_ik - V_i'k = (h_k - d_k)/2 mod 2` gives
`(-1)^{V_ik} = sigma_k K_ik` off the diagonal; Paley `K^T K = qI - 11^T` and Cauchy–Schwarz give the bound.

**Computed (inverse; judge reproduced [J], §1.7(b)).** Unit weights, full positive class: exact minimum of the total
owner weight over all integral `V` = parity bound `m(q-3)/2` at `q = 11, 19, 23` (20, 72, 110); strictly above it at
`q = 31` (240 > 210) and `q = 43` (504 > 420). Astra's Cauchy–Schwarz form (10.8, 48.8, 78.2, 158.4, 333.1) is
weaker with the same leading term.

**Consequence (inverse, Proved).** README_IX §5 item 1, second half ("minimise the witness heights over right
inverses `U`"), is answered negatively for heights: every exact global witness family of the IX.3 design has
`log J_i = log(p_i H_i^2) >= (1 - 2eta) T` for all but `O(eta^{-2})` owners. The previous judge's single-`V` value
`log J ~ 1.2 q log P` is within a factor about 1.2 of optimal. **Not checked (inverse):** e-all on the
height-optimal witnesses at `X = 3m^A`.

### 1.4 weight: X.4 (referee: CORRECT as logic; FALSE as a description)

**Proved.** *(X.4), stronger form: for a member of `F^all_X` and `x <= X`, at most `D(x)` labels have
`log p_j < theta_3(x)`; hence `W >= (s - D(x)) theta_3(x)`. Here `D(x)` = number of distinct tuples
`(r_j(l)/conj r_j(l))_{l <= x, l = 3 (4)}` and `theta_3(x) = sum_{l <= x, l = 3 (4)} log l`.* One line from (1.2) at
unit one on `e_j - e_k`, `nu_0 = 1`, `|sin| <= 1`.

**Restricted rate, corrected statement, quoted in full (weight):** *Fix `C`. For a sequence of members of
`F^all_X(m, C)` with `m -> infinity` and `X >= K log m` for every fixed `K` (e.g. `X = 3m`), for each fixed `K > 0`:
`liminf W/(m log m) >= (K/2) liminf (s - D(K log m))/m`. In particular `D(K log m) = o(m)` for every fixed `K`
implies `W/(m log m) -> infinity`; Astra's bounded-rank hypothesis is one sufficient condition. Contrapositive:
`W <= A m log m` forces `D(K log m) >= s - (2A/K + o(1)) m`.*

**Corrected wording.** In the proof of `log L(x) = o(x)` (`L(x) = lcm_{l <= x, l = 3(4)} (l + 1)`) the bound
"`psi(x/B)`" for small primes is wrong (powers up to `x + 1` occur); the true contribution is
`<= pi(x/B) log(x+1) = (1 + o(1)) x/B`, so the order and the lemma stand. Computed: `log L(x)/x` = 0.145, 0.109,
0.097, 0.083, 0.073 at `x = 10^2..10^6`.

**False (as a description): "a restricted rate theorem with content for genuine circles."** For genuine factors,
equal tuples at all inert `l <= x` force `prod l | Im(pi_j conj pi_k) != 0`, so `e^{theta_3(x)} <= sqrt(p_j p_k)`
for **any** two split primes, with no cluster. So on genuine data (X.4) is a statement about the prime set; on
natural parents (all split primes in `[P, 2P]`; first `s` split primes) `D(x) = s` once
`e^{theta_3(x)} > 2P`, and (X.4) reads `W >= 0`; their tuple rank is `pi_3(x)`, unbounded (PNT in fixed moduli).
Genuine parents meeting the hypothesis exist (primes `a^2 + (Q_y b)^2`, rank 0), but there all but `m^{o(1)}` primes
exceed `m^{K/2 - o(1)}` for every `K` and the conclusion follows from the prime sizes. Meaning proves the converse
cleanly: the hypothesis fails on every genuine parent with `P <= m^{O(1)}`. [J: confirmed, §1.7(c).] For fakes X.4
only asks for nearly injective tuples, which generic witnesses have: it excludes no known fake.

**Calibrations of X.5 (weight).** Function-field exponent arithmetic: Proved. The field families' growth and the
"T5" sentence: Not checked (the weight referee could not locate T5; it is item (T5) of `reply10.txt`, the twin
relaxation of the other thread [J]).

### 1.5 meaning (analyst; not a re-referee)

* *Genuine witnesses pass every residue test on every parent (Proved, meaning §2.1):* for genuine `pi_j` and any
  `v`, the collision modulus of `v` at unit `i^u` is the truncated odd exterior part of a fixed nonzero integer
  `D_u(prod pi^{v+} conj pi^{v-})`, with equality in (1.2) before truncation. Computed: 360 pair-unit cases at
  `P = 100`, 0 mismatches. So residue tests can contribute to a growth proof only by excluding nongenuine models,
  i.e. through rigidity.
* *Cluster words are the X-smooth shadow of (T1) (Proved step):* on a genuine cluster the mark-difference character
  `a_i - a_j` is tested in P3-char, and its collision modulus is the X-smooth part of the (T1) cofactor `t_ij`.
* *Reduction of the residual question to (A) (Heuristic).* Any P3-all\* lower bound would apply to genuine circles
  and so must in particular give `M = o(s)` on `q_P`. [J: this direction is trivially true of any growth theorem;
  the heuristic content is that the residue machinery offers no route to it other than rigidity plus (T1).]
* **Computed (meaning; P = 100 and 200 reproduced by the judge [J]):** on `q_P`, minimal arc distance between two
  lattice points divided by `q_P^{1/4}`: 0.76 (`P = 50`), 3.13 (100), 6.66 (150), 468 (200), 4.6e3 (250), 3.3e4
  (300), 6.6e5 (350), 4.8e6 (400), 8.1e7 (450). `M_max(q_P, C)` for `C = 1, 4, 8`: (2,2,3), (1,2,2), (1,1,2), then 1
  for `200 <= P <= 450` (certified by meet-in-the-middle with fixed-point error bounds and 80-digit recomputation
  of the ten best candidates); `P = 500` heuristic only. Primes `5 mod 8`: 25.6, 739, 1.3e4, 8.6e5 at
  `P = 200, 400, 500, 600`. The gaps follow the random model `(pi/4) 3^{-s}` within a factor 7.
* *Feasibility of (A1) on `q_P` (Heuristic):* a lower bound `|sum c_p theta_p - k pi/4| >> q^{-1/4}` for all
  `c in {0, +-1}^s`; Liouville is exactly critical (`|D| <= s/2` only), Matveev loses a tower. Estimates: 3% for
  `M = o(s)` on `q_P`, 5% for `M <= (1/3 - c)s`.

### 1.6 Withdrawn or corrected wording (collected)

Nothing of Astra's mathematics is withdrawn. Corrected: "a counterexample to a growth argument requiring even this
fraction of the tests" (false; data-independent `b(m)`-separable word sets only); "an omitted difference can be
forced to fail" ("can be chosen to fail"; cutoff of order `m^3`); "`phi = +-arg r (mod pi/2)`" (precisely `theta`
or `pi/2 - theta` mod `pi`); "least cutoffs ... every local-unit rational lift" (only lifts with parent-prime
denominators; X.1's wider class not certified, at most 769 on the first parent under the zero-charge model); "almost
every owner" (all but `O(eta^{-2})`, owner rows only); X.4 "restricted rate" (reformulation; `D`-form); `psi(x/B)`
(misworded); summary line 10 (genuine circles: no content). The weight constant of X.2 is `1536H^2/rho`, not
uniform in `H`.

### 1.7 Judge's checks [J]

**(a) X.1 line by line, and one table row.**
1. *Lowest terms.* `r = pi gamma^2/Norm gamma`, so `b | Norm gamma <= H`. A rational prime dividing `alpha` either
   divides `b` (not minimal) or has square dividing `p` (impossible). `b` even would make `4 | J = x^2 + y^2` with
   `gcd(x, y) = 1`: impossible, so `b` and `J = p b^2` are odd. An inert prime cannot divide `J` (it would divide
   `x` and `y`); a split `l | J` has exactly one factor dividing `alpha`. Disjoint supports.
2. *Integers.* `x_n^2 + y_n^2 = J^n` odd, so `x_n, y_n` of opposite parity and `x_n +- y_n` odd; a common prime of
   `x_n, y_n` would split and put both factors in `alpha^n`. `(x_n - y_n)(1 - i) = alpha^n - i conj(alpha)^n` and
   `(x_n + y_n)(1 + i) = alpha^n + i conj(alpha)^n`; nonzero because `alpha/conj alpha` is not a root of unity; a
   prime of `J` dividing either would put the factor of `alpha` into `conj alpha^n`.
3. *Parent primes.* `l | x_n -+ y_n`, `l` not dividing `J`, gives `u^n = +-i` for `u = alpha/conj alpha` in the
   cyclic norm-one group of order `l - 1`; `ord(u^n) = 4` with `n = 2^k` forces `ord u = 4n` (an odd part of
   `ord u` survives in `ord u^n`). With `4n = 2^{T+1}` this cannot divide `l - 1`.
4. *Charging.* The units `i`, `-i` charge exactly `x_n - y_n`, `x_n + y_n` (unit `i^s` at `q^a` iff
   `q^a | alpha^n - i^s conj alpha^n`); every prime power factor is odd, exterior, prime to `b`, `<= X`. With
   `nu = 2^{-1/2}`: `sin^2(n phi - pi/4) + sin^2(n phi + pi/4) = 1 >= ((x_n - y_n)^2 + (x_n + y_n)^2)/(2p^n) =
   J^n/p^n = b^{2n}`.
No gap. Computed (`x1_judge.py`): order lemma on all `u in Z/N`, `N < 300`, `n = 2^k <= 32` (840 cases, 0 failures);
2,635 nonintegral lifts on the split primes in `[100, 200]` (`T = 6`, `n = 32`, `H = 400`) and 2,894 on the primes
`5 mod 8` in `[1000, 2000]` (`n = 2`, `H = 4000`): oddness, nonvanishing, coprimality to `J` and to every parent
prime, size bound and the sum-of-squares identity, 0 failures. **Row `P = 100` rebuilt from scratch:** 10 labels,
190 cases (`alpha = pi_p pi_l^2` or `pi_p conj(pi_l)^2`, `b = l`); kill test by the exact rational rectangle for
`p^n(cos 2n phi, sin 2n phi)` for each `n <= 40`; at `X = 437, 438` exactly one case survives every single-exponent
test, label 181, `b = 113`, `alpha = 985 - 1158i`; at `X = 439` none; the survivor (and its conjugate) passes every
exponent at `phi = arg(2+i)` (checked to `n = 211`, after which the inequalities hold automatically). **Least
cutoff 439 confirmed.**

**(b) X.3 bound.** The lower bound is the parity count: each odd `V_ik` costs at least `log p_k`. `x3_judge.py`
computes, column by column, the forced parity pattern of the owner entries (one free bit per column), the minimal
number of odd entries, and then searches for an integral `V` with owner entries in `{0, +-1}` on exactly those
positions, non-owner entries free, solving `D V_{.k} = ((H - D)/2)_{.k}` by an own column-Hermite solver (every
solution verified). Results: `q = 11, 19, 23`: parity bound 20, 72, 110 = `m(q-3)/2`, attained, so these are the
exact minima. `q = 31`: parity bound 210, not attained in 15 of the 30 columns (a vector of weight equal to the
number of forced odd entries must be `+-1` there and 0 elsewhere, so the exact minimum exceeds 210), consistent with
the referee's 240. The bound is correct, intrinsic to the design (all `V`, all `U`), and nearly sharp at small `q`.

**(c) X.4 for genuine circles.** `x4_meaning_judge.py`: split primes in `[100, 200]` and `[1000, 2000]`, genuine
factors, tuples at the first inert primes `3, 7, 11, 19, 23, 31, 43, 47`. `D(x)` = 4, 10 (`P = 100`; full from
`x = 7`) and 4, 29, 63, 67 (`P = 1000`; full from `x = 19`), matching weight's table; the `F_2`-rank of the Legendre
vectors equals `pi_3(x)` = 1, ..., 8 at every step, so the tuple rank is maximal; all 612 equal-tuple pairs satisfy
`prod l | Im(pi_j conj pi_k) != 0` and `(prod l)^2 <= p_j p_k`. Bounded-rank genuine parent: 12 primes
`a^2 + (231 b)^2` (smallest 213,613) have one tuple at `l = 3, 7, 11` (`omega = 1`), rank 0. So: the joint-rank
hypothesis **can** hold on genuine parents, but only on specially built ones whose small inert residues are
trivial, and there (X.4) carries nothing beyond the sizes of the primes; on the natural dense parents it fails and
(X.4) is vacuous for `x` past `(2 + o(1)) log P` scale. I uphold weight and meaning.

**(d) `M_max` on `q_P` [J].** Full enumeration of all `4 * 2^s` lattice points (`P = 100`: 4,096; `P = 200`:
262,144; double precision, rounding `<= 10^-14` rad against gaps `>= 6 * 10^-8` rad): minimal arc distance
`/ q^{1/4}` = **3.125** (`P = 100`) and **468.4** (`P = 200`); `M_max` for `C = 1, 4, 8` = **(1, 2, 2)** and
**(1, 1, 1)**. Both rows of meaning's table reproduced.

**(e) Referee disagreements.** (1) rays1 calls the "fraction of tests" sentence a matter of interpretation, rays2
calls it false as worded: same substance; I rule **false as worded**, corrected as in §1.2. (2) finite says X.1
"closes the low-height coherent-witness route to fakes on such parents", meaning says X.1 is "not progress toward a
lower bound": both true and compatible (X.1 restricts how fakes can be built; X.2/X.3 show the restriction is not
binding on fakes that use local witnesses). (3) meaning estimates a full P3-all\* fake at 70%; nobody disputes the
direction. No disagreement on any proof.

## 2. What is new relative to the record, and correct

* **X.1, effective finite rigidity (new, Proved).** The record had IX.7 (ineffective, all cutoffs, rational lifts)
  and round 14's Proposition P (genuine residues at the pinning moduli). X.1 is the height-quantified bridge: lifts
  of twist norm `<= H` are forced genuine at `X >= sqrt2 (pH^2)^{n/2}`, polynomial in `pH^2` on 5 mod 8 parents (and
  for bounded `v_2(l-1)`). The use of a dyadic power to make the two diagonal coordinates prime to every parent prime
  is the new idea; it is elementary.
* **X.2, partial P3-all\* fake (new barrier, Proved).** First fake passing e-all, at every unit and every cutoff, on
  all powers of all singleton rays and all short words with a nonzero class sum, jointly with P3-char at every
  cutoff. Mechanism: block-sum-preserving local resampling plus a positional digit code that keeps tested logs away
  from all high divisors of the cyclic local orders (a failure of the product formula that no global data can have).
* **Lift heights diverge on a passing parent (new, Proved).** On the 5 mod 8 version, every finite-cutoff rational
  lift has twist norm `> (X/(sqrt2 p_j))^{1/2}`. This is what kills the height route (§3(a)).
* **X.3 (new, Proved, negative).** No short exact witnesses for the IX.3 owner design, over all right inverses and
  all exterior twists: owner log twist `>= (1/2 - eta) T` for all but `O(eta^{-2})` owners; referee's Remark R and
  exact minima.
* **Meaning's data on `q_P` (new to the record, Computed).** First exact values of `M_max` on the explicit dense
  parents that are extremal for Theorem 1.1; truth `M = 1` from `P = 200` against a known bound of order `s/3`.
* **Minimal novelty:** X.4 (one line from (1.2); a reformulation), the lcm lemma (standard sieve), the limit of the
  order argument, the function-field transfer.

## 3. Rulings

**(a) Is P3-all\* settled? No, in neither direction.** Exactly:
* **Capped now** (at `W = O_{C,H}(m log m)`, one parent sequence, one witness family valid at every cutoff):
  P3-char at every cutoff (IX.3) **together with** e-all at every unit, every cutoff and all powers on every word of
  mass `<= H` whose class-sum vector is nonzero for a constructor-chosen partition into `b(m) -> infinity` classes
  (X.2). This contains every singleton ray, every pair sum, every cross-class difference, every odd-mass word. Any
  argument using only these tests on a data-independent, `b(m)`-separable word set cannot prove
  `W/(m log m) -> infinity`. `H` is fixed before the weight constant (`~1536H^2/rho`).
* **Shown to be "the problem itself" at finite cutoff:** on comparable parents with bounded `v_2(l - 1)` (e.g. all
  primes `5 mod 8` in `[P, 2P]`), a member of P3-all at a polynomial cutoff `X >= 2 sqrt2 P^{1+2h}` whose witnesses
  are reductions of lifts of twist norm `<= P^h` is an actual cluster (X.1 + round 14 Prop. P). This moves the
  "problem itself" end of the residue column (IX.7, limit `X -> infinity`) down to polynomial cutoffs, **for
  polynomial-height witnesses only**.
* **Not shown either way:** e-all on the balanced words inside classes (the lattice spanned by same-class
  differences), in particular data-chosen label differences `e_j - e_k`, at `3M^A`; words of length `H(m) -> infinity`.
* **The height route is closed (Proved):** X.1's hypothesis is not implied by passing the tests (X.2, 5 mod 8: least
  lift heights diverge with `X` on a fixed passing parent), the natural exact witnesses have exponentially large
  height (X.3), and every finite datum has a lift of height `e^{O(X)}`. A lower bound for P3-all\* cannot be
  obtained by bounding lift height and applying X.1.

**(b) Is the residue route closed for internal work? Yes.** By the previous judge's rule ("if that pass returns
another fake"): in substance, not literally. Part X returned a fake for every test except the data-chosen balanced
words, and showed the other arm of the posed alternative (height-quantified IX.7 followed by a lower bound) cannot
be completed. On the merits: genuine witnesses pass every residue test on every parent (meaning §2.1), so residue
tests can help only by forcing rigidity; rigidity at finite cutoff requires the lift-height hypothesis, which passing
data do not satisfy; what is left (same-class pair collisions across primes) has no proposed mechanism, and any
lower bound from it would include a rate on the explicit parents `q_P`, where the problem is a linear-forms-in-
logarithms question at the Liouville threshold. **Do not commission a further P3-all\* pass.**

**(c) X.4: not a restricted rate theorem with content for genuine circles.** It is a correct one-line consequence of
(1.2); for genuine factors it reduces to the divisibility `prod l | Im(pi_j conj pi_k)`, a property of any two split
primes. Its hypothesis fails on every natural parent and on every genuine parent with polynomially bounded labels;
where it holds the conclusion is trivial from the prime sizes. For fakes it is a mild necessary condition (nearly
injective inert tuples). Record it in the `D`-form of §1.4, labelled a reformulation.

**(d) What remains uncapped and not equivalent to the problem.**
* The **same-class balanced words** of P3-all\* at `3M^A`: whether non-liftable local witnesses can make every
  label difference `e_j - e_k` pass, i.e. a design of pair coincidences across the exterior primes `<= 3M` with
  colliding product `<= sqrt(p_j p_k)` for every pair, compatible with X.2's singleton data. Words of unbounded
  length.
* **Twisted point collisions** (`sum c != 0`) at moduli `M^A` (open since round 5); the 2-adic clause.
* The **Gaussian cocycle (T1) on a fixed parent** (question (A) of `reply10.txt`, not taken by the colleague):
  either `M = o(s)` on one explicit dense family such as `q_P`, or an integral-cofactor fake.
* The window for `delta >= 1/6`; the cross-gap inside the window; the operator-form union large sieve; square-root-
  quality prime sums and box counts; Lang–Waldschmidt.
* Minor: the X.5 least cutoffs for X.1's full class (exterior denominator primes).

**(e) Consequences for the integer problem.** **None for the rate: `M <= (2/3 + eps) log R/loglog R` and the
uniform bound `M(C) < infinity` are unchanged.** Part X proves a correct effective rigidity statement that applies
only to witnesses no fake is obliged to use, and a barrier that removes the last residue class with a plausible
mechanism (singleton rays plus height) from the list of places a growth proof might hide. The one new fact about
actual circles is computational: on the explicit dense parents `q_P` the truth is `M_max = 1` for `200 <= P <= 450`
with margin `468` to `8 * 10^7` (in units of `q_P^{1/4}`), against a known bound of order `s/3`. Bluntly: the residue
route has been mapped to its end and carries no rate.

## 4. Recommended edits to the expert note (recommendations only; nothing was written to the repository)

1. **§1, "Not everything is capped" paragraph: replace** "Not capped: P3-all$^*$ (prime-level residues tested on
   every monomial at moduli $3M^A$), twisted point collisions" **by:**
   "Within P3-all$^*$ (prime-level residues tested on every monomial at moduli $3M^A$), one parent sequence at
   $W=O_H(M\log M)$ passes, at every modulus at once and together with P3-char, every power of every single-label
   ray, every pair sum and every signed word of length $\le H$ whose coefficient sum is nonzero on one of
   $b\to\infty$ label classes (round 31, Part X.2; **proved**, `round31/README_X.md` §1.2). Not capped: the balanced
   words inside a class of P3-all$^*$, in particular label differences $e_j-e_k$ chosen from the data, twisted
   point collisions"
   (rest of the sentence unchanged).
2. **Same paragraph, after** "since such lifts exist at every finite modulus, this is the problem itself and says
   nothing about P3-all$^*$." **insert:**
   "Effectively (Part X.1; **proved**): if all primes of $N$ are $\equiv5\pmod8$ and lie in $[P,2P]$, a witness that is
   the reduction of $\pi_j\gamma/\bar\gamma$ with $\mathrm{Norm}\,\gamma\le P^h$ and passes the two diagonal tests on
   $2e_j$ at a modulus $\ge2\sqrt2P^{1+2h}$ is a genuine prime factor, and with genuine factors the model is an
   actual cluster (round 14, Prop. P). Passing data need not have such lifts: on one fixed fake every lift at label
   $p_j$ has $\mathrm{Norm}\,\gamma>(X/(\sqrt2p_j))^{1/2}$ (Part X.2), so bounded-height rigidity yields no lower
   bound for P3-all$^*$. Genuine witnesses pass every residue test on every parent, so residue tests can serve a
   growth proof only through rigidity."
3. **§1, Legendre bullet, append:** "A $128H^2$-fold replicated version also passes all single-label rays and all
   words of length $\le H$ with nonzero class sum, at $W\le(1536H^2/\rho+o(1))M\log M$ (Part X.2; **proved**). For
   the IX.3 witness design no right inverse gives short witnesses: almost every owner twist has log norm
   $\ge(\frac12-o(1))\log\mathrm{rad}\,N$ (Part X.3; **proved**)."
4. **§5 item 1, replace the caveat's last clause** "P3-all$^*$ is uncapped." **by:** "P3-all$^*$ is uncapped only on
   the balanced words inside label classes, in particular data-chosen label differences (round 31, Part X.2); the
   residue route is closed for internal work (`round31/README_X.md` §3(b))."
5. **§5 item 12 (round 31 history), add after the Part IX bullet:**
   "* ***Part X: finite cutoff and coordinate rays*** (`round31/README_X.md`; six referee blocks; no error found;
   judge's score on the goal 2/10). Effective finite rigidity for bounded-height rational lifts (polynomial modulus
   on the $5\bmod8$ dyadic parent); a partial P3-all$^*$ fake (all single-label rays and all short words with nonzero
   class sum, every power, unit and modulus, with P3-char); on it the least lift heights diverge with the modulus;
   no short exact witnesses on the IX.3 design (mod-four bound, exact minima 20, 72, 110 at $q=11,19,23$); X.4 is a
   reformulation whose hypothesis fails on genuine parents with polynomially bounded primes. Computed: on
   $q_P=\prod_{p\equiv1\,(4),\,P\le p\le2P}p$ no two lattice points lie within $8q_P^{1/4}$ for $200\le P\le450$
   (closest pair $468$ to $8\cdot10^7$ times $q_P^{1/4}$). *Judge's reading*: the residue route is closed for
   internal work; the rate is unchanged. Do not write "P3-all$^*$ is capped", "Part X gives a finite-cutoff lower
   bound", "X.4 proves a growth rate", "the omitted difference must fail", "the twist must be $\pm1$",
   "$M_{\max}(q_P,C)=1$ is proved"."
6. **Pointer table:** add "round 31 Part X: finite rigidity, partial P3-all$^*$ fake, data on $q_P$ |
   `round31/README_X.md`, `round31/ref10_*.md`"; add `README_X.md` to the round-31 row.
7. **Sentences to avoid** (also in §5): "P3-all$^*$ is capped" or "answered"; "X.2 defeats any argument using a
   fraction $1-o(1)$ of the tests" (only data-independent, class-separable word sets); "X.1 gives a lower bound for
   P3-all$^*$"; "X.4 is a restricted growth rate for circles"; "an omitted difference is forced to fail" (chosen, at
   modulus of order $M^3$); "the twist must be $\pm1$" (false on genuine circles); "the least cutoff for bounded lifts
   is 439" without "for lifts with parent-prime denominators"; "$M_{\max}(q_P,C)=1$ is proved" (computed for
   $200\le P\le450$).

## 5. Should this thread stop?

**The residue route stops now. The thread gets one final question, labelled final, and then stops whatever the
answer.** The reason not to stop outright: the one lead both the coordinator and the previous judges called live,
question (A) (the Gaussian cocycle (T1) on one fixed parent), was offered and not taken, and the meaning referee's
data make it concrete: on `q_P` the truth is `M = 1` with margin `((1.4P)^{1/4}/3)^s`, so a fixed-parent question
is now well calibrated.

**The one question (final).** *On the explicit parents `q_P` (the product of all primes `p = 1 mod 4` in `[P, 2P]`,
or of all primes `5 mod 8` there), with `C` fixed before `P`, do one of the following. (A2, preferred) Build an
integral-cofactor fake on `q_P`: Gaussian integers `z_1, ..., z_M` with `M` larger than any fixed constant (ideally
`M >> s`), whose exponent rows are rows of `q_P`, with `z_j - z_i = 2i eta H_ij t_ij` as in (T1), `t_ij` nonzero
rational integers with `|t_ij| <= (C/2) sqrt(Q_ij) q_P^{-1/4}`, the triangle identities exact in `Z[i]`, all
divisibilities exact and the signed residues kept, dropping only `|z_i|^2 = q_P` (state the norm defect); or prove
that none exists and say what the proof uses. (A1, if she sees a route) Prove `M_max(q_P, C) <= (1/3 - c) s`, or
`M = o(s)`, or state exactly which Diophantine inequality on the `theta_p` would suffice and show (by a fake on the
same prime set) that (T1), the cocycle and the signed residues do not imply it.* Full proofs or exact reproducible
computations; test against the meaning referee's data (`ref10_meaning/`), the twin relaxation (T5) and the function-
field families.

**Honest probabilities.** Moving the worst-case growth rate for general `R`: **<= 1%** (a theorem on `q_P` is about
one explicit family and moves the worst case only if its method extends to all dense profiles). A restricted theorem
on `q_P` (`M = o(s)` or a constant below `1/3`): about 3%. A useful barrier (an integral-cofactor fake capping the
cocycle input on fixed parents): about 25%. The uniform bound: < 0.5%. After this question the thread should stop.

## 6. Scores (0–10, on the integer goal)

| block | correctness | score on goal |
|---|---|---|
| X.1 finite rigidity, pinning, dense families | correct; new and effective; applies only to bounded-height witnesses | 1.5 |
| X.2 partial P3-all\* fake, lift-height divergence, deliberate failure | correct; one sentence false as worded; strongest new barrier of the part | 2.5 |
| X.3 no short exact witnesses | correct; closes one sub-question negatively | 1 |
| X.4 weight bound under joint rank | correct logic; reformulation, no content for genuine circles | 0.5 |
| X.5 computations and calibrations | reproduced (one row by the judge); scope narrower than stated | 0.5 |
| meaning's data on `q_P` (referee addition) | computed, spot-checked | 1 |
| **Part X overall** | **no error found; careful and honest scope statements, one overstatement** | **2** |

As mathematics the quality is high (about 7.5/10): X.2's digit-code construction is intricate and survived two
referees with millions of exact checks, X.1 is short and clean, X.3 is tight at small `q`. The score reflects that
neither the growth rate nor the uniform bound moved, that the rigidity statement binds no fake, and that the one
claim beyond its proof (the "fraction of the tests" sentence) concerned exactly the tests that remain open.
