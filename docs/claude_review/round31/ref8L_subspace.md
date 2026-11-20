# Referee report on Astra Part VIII.4: the absolute-threshold Subspace hypothesis

Block: section VIII.4 of `astra_lemma_VIII/answer.md` (lines 1837-1868), with the scope remark at line 1890 and the summary items in the opening paragraph (line 1746) and the twelve-line summary item 8 (line 1925).
Scripts and outputs: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref8L_subspace/` (`check_viii4.py`, `run_eq.py`, `part1_out.txt`, `part2_out.txt`). I wrote the code myself. I read Astra's `viii_check.py` only to confirm that it does not test VIII.4. It does not.

## 1. Verdicts

| Claim | Verdict | Label |
|---|---|---|
| VIII.4 main theorem: in every dimension `k >= 2`, for every `0 < delta < k-1` and every fixed `H_0`, covering the solutions of the round-8 system takes `>= exp(c s^theta log s)` proper subspaces for every `theta < 1 - delta/(k-1)` | correct; it refutes round 8's `QS(m, delta, f, H_0)` for every `m`. Two attribution corrections apply (§3) | Proved |
| VIII.4 five-mark equal-norm theorem (`k = 5`, `0 < delta < 2`) | correct | Proved |
| Referee extension: the equal-norm obstruction holds in every dimension `2 <= k <= 6` over the whole arc range `0 < delta < (k-1)/2` | new here | Proved (§5) |
| Scope remark (line 1890): both constructions lie below the relative cutoff `(prod_{p in P} p)^gamma` | correct. It agrees with the record's open `QS^prod` | Proved |
| Small instances, exact arithmetic | all pass | Computed (§6) |

**Most important point.** The record already contains a refutation of the fixed-threshold hypothesis: lit18 Prop A gives `f(s) >= exp(s^(1/2-o(1)))`. That proposition is still listed as **unrefereed** (`overview/overview.md` row 13 and §6.7). Astra's proof is shorter and fully elementary, and its exponent is never worse than Prop A's: it is strictly better when `delta < (m-1)/2` and matches Prop A in the limit `delta -> (m-1)/2`. Once this report is accepted, the fixed-threshold refutation of `QS` can be recorded as **proved and refereed**. Lit18 Prop B, which handles thresholds `H_0(s)`, is not covered by this report and stays unrefereed.

The equal-norm version is new to the record. It shows that "equality of norms" alone does not let an `S`-integer input escape at a fixed threshold. The record's wording at `overview/overview.md:1730`, `lit18/README.md:135` and `lit18/check_tools.md:40, 321` must change (§7).

## 2. The round-8 hypothesis, quoted exactly

**The system.** `round8/weakest.md` §5.1, lines 663-672:

> Work over `K = Q(i)` with normalised absolute values: `|x|_inf = |x|^2` at the complex place and `|x|_w = N(w)^(-v_w(x))` at a Gaussian prime `w`. [...] For a finite set `P` of split rational primes, let `S_P` consist of `inf` and the places above `P`, with `s = |P|`.
> **The system `Sigma_m(P)`.** At `inf` use the forms `X_1` and `X_i - X_1` (`i = 2..m`). At each finite `w in S_P` use `X_1, ..., X_m`. [...]
> `prod_(v in S_P) prod_(i=1..m) |L_(i,v)(x)|_v / ||x||_v  <=  H(x)^(-m-delta),      x in K^m \ {0}.   (5.1)`

**The hypothesis.** `round8/weakest.md` lines 690-691:

> **Hypothesis `QS(m, delta, f, H_0)`.** For every finite set `P` of split primes, with `s = |P|`, the solutions of (5.1) with `H(x) >= H_0` lie in the union of at most `f(s)` proper linear subspaces of `K^m`.

The facts about this hypothesis that matter for VIII.4 are these.
* **Cutoff.** The lower height cutoff `H_0` is a fixed constant that does not depend on `P`. It is an absolute threshold.
* **Label.** The record calls it a "Hypothesis". `round8/README.md:331` calls it "non-standard" and, for `m >= 4` and `1 < delta < (m-1)/2`, "neither refuted nor supported". `round8/README.md:392` lists it among the "New simultaneous hypotheses [...] non-standard and unsupported".
* **What depends on it.** Only a conditional result: `weakest.md` Theorem 5.2 (lines 693-698), "Assume `QS(m, delta, f, H_0)` [...] carries `M <= f(omega(N')) + m - 1`". This yields growth when `f(s) = o(s)` and the uniform bound when `f` is bounded.
* **Round 8 did not ask for polynomial `f`.** The record's target was sublinear `f`; see Theorem 5.2 and Remark 5.4 at line 763: "sublinear `S`-dependence of the number of exceptional subspaces of `Sigma_m` for some `m >= 4` and `delta in (1, (m-1)/2)`".
* **Where "polynomial in `|S|`" comes from.** The phrase appears in the coordinator's question (`reply8.txt`: "a quantitative Subspace conjecture with polynomial dependence on |S|") and in `lit18/tools.md` §2.2 Consequence 1: "No polynomial-in-`|S|` quantitative Subspace [...] theorem exists in the standard (fixed-threshold) form".
* **Status in the record.**
  * `round8/weakest.md` Prop 5.3 (refereed) refutes `QS(m, delta, o(s))` for `m <= 3`.
  * lit18 Prop A (`lit18/tools.md` §2.2; recounted in `lit18/check_tools.md` §2.1) refutes it for every `m` with `f(s) >= exp(s^(1/2-o(1)))`. Prop A is unrefereed per `overview/overview.md:937, 1059, 1798`.
  * lit18 Prop B (`check_tools.md` §2.2) covers every threshold `H_0(s)` but forces only `s - m + 1` subspaces. It is also unrefereed.

Astra's counterexample uses exactly this system, these normalisations and an absolute threshold. It therefore refutes **that** statement, not a variant. It also refutes every weakening with larger `f`: any `f` below `exp(c s^theta log s)`, so in particular every polynomial `f`.

## 3. The main theorem, corrected statement

**Theorem (VIII.4; Proved).**
* **Data.** Fix integers and reals with `k >= 2`, `0 < delta < k-1`, `H_0 >= 1`, `0 < theta < 1 - delta/(k-1)` and `0 < c < 1 - theta - delta/(k-1)`.
* **Primes.** For `X >= 2`, let `P_X` be the set of primes `p ≡ 1 (mod 4)` with `X <= p <= 2X`, and put `s = |P_X|`.
* **Covering bound.** There is `X_0 = X_0(k, delta, H_0, theta, c)` such that for every `X >= X_0` the following holds. Any family of proper linear subspaces of `Q(i)^k` whose union contains every solution of (5.1) (with `m = k`, `P = P_X`) of height `H(x) >= H_0` has at least `exp(c s^theta log s)` members.
* **Every `s` occurs.** As `X` increases by 1, `s(X)` changes by at most 1, because `2X+2` is even. Since `s(X) -> infinity`, every sufficiently large integer occurs as `s(X)` with `X >= X_0`.
* **Consequence.** `QS(k, delta, f, H_0)` forces `f(s) >= exp(c s^theta log s)` for **all** sufficiently large `s`. Hence `QS` is false for every `m >= 2`, every `0 < delta < m-1`, every `H_0`, and every `f` with `f(s) = exp(o(s^(1-delta/(m-1)) ))`. This includes every polynomial `f` and every `f = o(s)`.

Two corrections apply to Astra's framing.
1. **What round 8 proposed.** Round 8 proposed sublinear `f`, not polynomial `f`. Refuting polynomial `f` is stronger, so the conclusion stands. "The specific [...] system proposed in round 8, a polynomial bound" should read "the system of round 8 with any sub-exponential bound `exp(o(s^(1-delta/(m-1))))`".
2. **Novelty.** The refutation itself is not new. lit18 Prop A already refutes the fixed-threshold form with `exp(s^(1/2-o(1)))`, uniformly in `delta < (m-1)/2`. The new content is the following.
   * The exponent improves to `s^(1-delta/(m-1)-o(1))`. For `delta < (m-1)/2` this is strictly larger than `1/2`, and as `delta -> 0` it tends to `1`. The record's upper bound is `exp(O(s))` (Evertse 1996; the record quotes it from memory, and I did not read it). The gap therefore narrows to `s^(1-delta/(m-1)-o(1)) <= log f(s) <= O(s)`.
   * The proof uses no counting of smooth `S`-integers in squares.

**Remark (Proved).** The same tuples have `H >= T^(2b) = (2X)^(2bj)`. The theorem therefore also holds for every threshold `H_0(s) <= exp(2b (1-o(1)) s^theta log s)`, with `b = delta/(k-1)`. It does not reach the sparse thresholds `exp(A s log s)` of lit18. Those are handled only by the unrefereed Prop B, with a linear count.

### 3.1 Step-by-step check of the proof

All steps are Proved. I checked each one.

1. **`s ≍ X/log X`.** This is the prime number theorem for the fixed progression `1 mod 4`. Then `log(2X) = log s + log log s + O(1) = (1+o(1)) log s`.
2. **Products.** Put `j = floor(s^theta)` and `T = (2X)^j`. The `C(s,j)` products of `j` distinct primes of `P_X` are distinct rational integers in `[X^j, T]`, by unique factorisation.
3. **Bins.** Put `b = delta/(k-1) < 1`. Cover `[X^j, T]` by intervals of multiplicative ratio `1 + T^(-b)`. This needs `ceil(j log 2 / log(1+T^(-b))) <= 2 j T^b` intervals, since `log(1+y) >= y log 2` for `0 < y <= 1`.
4. **Pigeonhole.** Some bin holds a set `E` with `log|E| >= log C(s,j) - log(2jT^b)`. Here `log C(s,j) >= j log(s/j) >= (1-theta) j log s` and `log(2jT^b) = b j log(2X) + O(log j) = (b + o(1)) j log s`. So `log|E| >= (1-theta-b-o(1)) j log s`, and `j log s >= (s^theta - 1) log s`. The claimed exponent follows, with every `c < 1-theta-b`.
5. **Gcd and heights.** Take distinct `u_1, ..., u_k` in `E` and `g = gcd(u_i)` over `Z`.
   * The Gaussian gcd equals `g`, because the ideal generated in `Z[i]` by rational integers is `gZ[i]`. So `x = (u_i/g)` is primitive in `Z[i]^k`.
   * Then `H(x) = prod_v max_i |x_i|_v = (max u_i/g)^2`, since `max_i |x_i|_w = 1` at every finite place. So `H <= T^2`.
   * `g` divides the nonzero integer `u_i - u_1`, and `|u_i - u_1| <= (min_l u_l) T^(-b)`. Hence `max u/g >= T^b` and `H >= T^(2b)`.
   * Consequently `H >= H_0` once `X` is large.
6. **Exact left side.**
   * Every `x_i` divides `u_i`, so it is a unit at every finite place outside `S`.
   * At a finite place `w in S`: the forms are the coordinates and `max_i |x_i|_w = 1`. The contribution is `prod_i |x_i|_w`.
   * Over all finite `w in S`, the product formula (with `|.|_inf = |.|^2`) gives `prod_i |x_i|^(-2)`.
   * At infinity: `|x_1|^2 prod_{i>=2} |x_i - x_1|^2 / H^k`.
   * The product is exactly `H^(-k) prod_{i>=2} ((u_i - u_1)/u_i)^2`. Each factor is `<= T^(-2b)`, so the left side is `<= H^(-k) T^(-2b(k-1)) = H^(-k) T^(-2delta) <= H^(-k-delta)`, because `H <= T^2`.
   * This matches Astra's displayed equality. The normalisations agree with `weakest.md` (5.1): there the norm `||x||_v` is the max norm, raised to the `k`-th power through the product over `i`.
7. **Count.** Let `V` be a proper subspace and `sum c_i X_i = 0` a hyperplane containing it, with `c_{i0} != 0`.
   * On a raw grid point `(u_1, ..., u_k)` in `E^k`, the other coordinates determine `u_{i0}`. So `V` contains at most `|E|^(k-1)` raw points.
   * Membership is invariant under the scaling by `1/g`, so covering the solutions covers the raw ordered distinct tuples.
   * Hence the cover has at least `(|E|)_k / |E|^(k-1) >= |E|/2` members once `|E| >= 2k^2`.
8. **Range.** The construction needs `b < 1 - theta`, that is, `theta < 1 - delta/(k-1)`. This contains the arc range `0 < delta < (k-1)/2`, where `theta` can be taken up to `1/2`.

**Consistency check.** For fixed `delta`, `exp(c s^theta log s)` with `theta < 1` is below the recalled Evertse upper bound `exp(O_{k,delta}(s))`, which holds above the threshold `max(k^(4k/delta), H(L))`. There is no contradiction.

## 4. The five-mark equal-norm theorem

**Theorem (VIII.4, equal norms; Proved).**
* **Data.** Fix `k = 5`, `0 < delta < 2`, `H_0`, `0 < theta < 1 - delta/2`, `delta/2 < b < 1 - theta`, and `0 < c < 1 - theta - b`.
* **Statement.** For `X >= X_0`, with `P_X` as above, the following holds. Any family of proper subspaces of `Q(i)^5` whose union contains every solution of (5.1) with `H >= H_0` **and five coordinates of equal complex modulus** has `>= exp(c s^theta log s)` members.

### 4.1 Step-by-step check of the proof

All steps are Proved.

1. **Directions.** Fix a Gaussian prime `pi_p` above each `p in P_X`. For each `j`-subset `A` put `G_A = prod_{p in A} pi_p` and `u_A = G_A / conj(G_A)`.
   * The `u_A` are distinct. Suppose `u_A = u_B` and `p in A \ B`. Then `pi_p` divides `G_A conj(G_B)`. But `pi_p` does not divide `G_B conj(G_A)`, because `p` is split and `pi_p` is not associate to `conj(pi_p)`.
   * `N(G_A) <= T`.
   * By pigeonhole over `<= 2 pi T^b + 1` arcs of width `T^(-b)`, some arc holds a set `E` with `log|E| >= (1-theta-b-o(1)) j log s`.
2. **Coordinates.** Take parameters `(u_1, ..., u_4)` in `E^4` and five of the six pair monomials `u_i u_j`. Multiplying by `prod_nu conj(G_nu)` gives `y_ij = G_i G_j prod_{nu != i,j} conj(G_nu)`.
   * These are Gaussian integers with common norm `prod_nu N(G_nu) <= T^4`, supported above `P_X`.
   * Remove the Gaussian gcd. The coordinates keep equal moduli, the vector is primitive, and `H = |y|^2 <= T^4`.
   * Their arguments lie in an interval of length `2T^(-b)`.
3. **Lower height.** If the coordinates are distinct, `1 <= |y_i - y_1| <= sqrt(H) |e^(i phi) - 1| <= sqrt(H) 2T^(-b)`. So `H >= T^(2b)/4`.
4. **Left side.** As in §3.1 step 6, the left side equals `H^(-5) prod_{i>=2} |(y_i - y_1)/y_i|^2 <= H^(-5) (2T^(-b))^8`. Then `(2T^(-b))^8 <= T^(-4delta) <= H^(-delta)` as soon as `T^(8b-4delta) >= 2^8`. This holds for large `X`, since `b > delta/2`.
5. **Count.**
   * A coordinate equality `u_i u_j = u_k u_l` is a nonzero quadratic in the parameters (or, when the pairs share an index, the linear condition `u_j = u_l`). There are at most 10 such equalities. By the Schwartz-Zippel bound over `C` on the grid `E^4`, each removes `<= 2|E|^3` parameter tuples.
   * A hyperplane `sum c_ij X_ij = 0` pulls back to the quadratic `sum c_ij u_i u_j`. This is nonzero because distinct monomials are linearly independent. So the hyperplane holds `<= 2|E|^3` parameter tuples.
   * Hence the cover needs `>= (|E|^4 - 20|E|^3)/(2|E|^3) >> |E|` members.
   * The parent norm `prod_nu N(G_nu)` varies with the tuple, but every parent is supported in `P_X`. This is the meaning of Astra's "the parent may vary within the ambient `P`".

**Comparison with the record.** Round 8's Prop 5.3 already uses equal-norm tuples, but only for `m <= 3` and with polynomial counts `s^(2-delta)`. lit18 Prop A uses unequal norms. The equal-norm exponential obstruction in dimension 5 is therefore new.

## 5. Referee extension: equal norms in dimensions 2 to 6 (Proved)

**Statement.** The construction of §4 with `r = 4` parameters covers the whole arc range in every dimension `k` from 2 to 6.
* **`k <= 4`.** Use `k` parameters and the linear coordinates `y_i = G_i prod_{nu != i} conj(G_nu)`. Then `H <= T^k`, the archimedean product is `<= T^(-2b(k-1))`, and the hypothesis needs `delta < 2b(k-1)/k`. Since `2(k-1)/k >= (k-1)/2` exactly when `k <= 4`, this covers `delta < (k-1)/2`.
* **`k = 5, 6`.** Use five, respectively all six, of the pair monomials of four parameters. Then `H <= T^4`, the archimedean product is `<= (2T^(-b))^(2(k-1))`, and the hypothesis needs `delta < b(k-1)/2`. This covers `delta < (k-1)/2` as `b -> 1`.
* **Count.** The count is as in §4 (linear forms for `k <= 4`, quadratics for `k = 5, 6`).

**Limit of the method.** For `k >= 7` this family needs `r >= 5` parameters (squarefree monomials of a fixed degree in `r` variables number at most `C(r, floor(r/2))`). It then covers only `delta < 2b(k-1)/r <= 2(k-1)/5`, which is not the full arc range. Whether the equal-norm obstruction holds for `k >= 7` near `delta = (k-1)/2` is open here.

## 6. Exact computations (Computed)

`check_viii4.py` computes the left side of (5.1) and `H` **place by place**: the archimedean factor, then both Gaussian primes above every prime dividing a coordinate norm. It uses exact `Fraction`s and does not use the product formula. It checks `LHS <= H^(-k-delta)`, with `delta = a/q`, exactly as `LHS^q H^(kq+a) <= 1`. Throughout, `X = 1000` and `s = |P_X| = 67`.

### 6.1 Rational squarefree products (`part1_out.txt`)

| `j` | `k` | `delta` | `b` | `C(s,j)` | fixed-partition bins | `|E|` (best window) | tuples checked | all solve (5.1) | min `delta_eff` | lower bound `(|E|)_k/|E|^(k-1)` |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 3 | 1/2 | 1/4 | 2211 | 63 | 78 | all 456,456 | yes | 0.524 | 75.0 |
| 2 | 5 | 1 | 1/4 | 2211 | 63 | 78 | 20,000 sampled | yes | 1.070 | 68.4 |
| 3 | 3 | 1/2 | 1/4 | 47,905 | 623 | 197 | 20,000 sampled | yes | 0.521 | 194.0 |
| 3 | 5 | 3/2 | 3/8 | 47,905 | 10,756 | 25 | 20,000 sampled | yes | 1.586 | 16.3 |

* In every tuple checked, the place-by-place value equals Astra's closed form `H^(-k) prod ((u_i - u_1)/u_i)^2` exactly, and `H = (max u/g)^2` exactly.
* `T^(2b) <= H <= T^2` held in every tuple checked.
* For every row, the exact all-tuple sufficient condition `((max E - min E)/min E)^(2(k-1)) (max E)^(2delta) <= 1` holds. This is the argument of §3.1 step 6, so **every** ordered distinct tuple of `E` solves (5.1) in all four rows.
* The last row has `delta = 3/2`, inside the arc range `delta < 2` for `k = 5`.
* At these sizes the forced count (75 and 194) already exceeds the linear bound `s - k + 1` of lit18 Prop B.

### 6.2 Equal norms (`part2_out.txt`)

* **Setup.** `j = 3`, `T = 8·10^9`, `b = 3/10`, angular window `1.07e-3`. The best window holds 45 directions. I used 16 of them, to keep the run time manageable.
* **Parameter tuples.** There are 65,536, of which 21,856 have a coordinate coincidence. These are exactly the tuples with a repeated input, `16^4 - (16)_4`. The remaining 43,680 are "good".

| version | `k` | `delta` | all 43,680 good tuples solve (5.1) | equal moduli after gcd removal | min `delta_eff` | `T^(2b)/4 <= H <= T^4` | lower bound `good/(2|E|^3)` |
|---|---|---|---|---|---|---|---|
| five pair monomials (Astra) | 5 | 1/2 | yes | yes | 0.748 | yes | 5.33 |
| six pair monomials (§5) | 6 | 1/2 | yes | yes | 0.919 | yes | 5.33 |
| linear (§5) | 4 | 1/3 | yes | yes | 0.559 | yes | 5.33 |

* **Near the endpoint.** A run at `delta = 3/2`, `b = 9/10` found only 2 directions in a window of width `1.2e-9`, so it tested nothing. Values of `delta` near the endpoint 2 need much larger `X`. This is a finite-size effect; the asymptotic proof is unaffected.

The computations confirm the identities and inequalities on small instances. The asymptotic counts rest on the proofs in §§3-5.

## 7. What in the record depended on the refuted hypothesis, and what must be reworded

**Nothing unconditional depended on `QS`.** The only result using it is `round8/weakest.md` Theorem 5.2. That theorem remains a correct implication, but with a fixed `H_0` its hypothesis is false. With thresholds depending only on `s`, it is capped by the unrefereed Prop B. The relative-threshold form `QS^prod` (`lit18/check_tools.md` §2.5; `overview/overview.md` §6.7) is untouched: by Astra's line 1890, both constructions have `log H <= 4j log(2X) = o(s log X)`. Astra's own relative-height hypotheses (VIII.5-VIII.6) are refereed in other blocks.

The following wording must change.
1. **`round8/README.md:331` (row `QS`) and `:392` (item 2).** Replace "For `m >= 4`, `delta in (1, (m-1)/2)` it is neither refuted nor supported" and the listing as an open "new simultaneous hypothesis" with: "False for every `m >= 2`, `0 < delta < m-1` and fixed `H_0`. Any admissible `f` satisfies `f(s) >= exp(s^(1-delta/(m-1)-o(1)))` (Astra VIII.4, refereed in `round31/ref8L_subspace.md`; also lit18 Prop A)." The same applies to `round8/weakest.md` lines 37, 725-726, 759-765 (Remark 5.4) and line 948. Those lines say that the bound is false only for `m <= 3` and open for `m >= 4`.
2. **`overview/overview.md` row 13 (line 937), lines 1059 and 1798, and §6.7 (lines 1335-1338); `expert_note/expert_note.md:676`.** The fixed-threshold part of the cap becomes **proved, refereed** (this report), with the sharper exponent `exp(s^(1-delta/(m-1)-o(1)))`. The part for thresholds `H_0(s)` (Prop B, `>= s - m + 1`, hence `M <= omega(N')`) stays unrefereed.
3. **`overview/overview.md:1730`, `lit18/README.md:135`, `lit18/check_tools.md:40` and `:321`.** Each says that an escaping `S`-integer input "must use equality of norms, or a threshold tied to `prod p`". Equality of norms **alone** does not escape. At a fixed threshold, equal-norm tuples supported on `O(s^theta)` of the primes of `P` still force `>= exp(c s^theta log s)` subspaces, for every `theta < 1 - delta/2` when `k = 5` (§4) and for a positive `theta` throughout the arc range when `2 <= k <= 6` (§5). What escapes the known constructions is **full support** (every `p in P` dividing the common norm) or, equivalently for the cluster application, a threshold tied to `prod_{p in P} p`. The parenthetical "(full support)" in `check_tools.md:40` carries the content and should become the main clause.
4. **Astra's own text.**
   * At lines 1746 and 1839, "a polynomial bound [...] proposed in round 8" should read "round 8's `QS` (which asked for `f = o(s)`), and every sub-exponential weakening `exp(o(s^(1-delta/(m-1))))` of it".
   * The section should cite lit18 Prop A as the earlier refutation with exponent `1/2`.

## 8. Sources

* **Read in full or in the relevant parts:**
  * `round8/weakest.md` §§5.1-5.2;
  * `round8/README.md` lines 190-200, 325-345 and 388-395;
  * `lit18/tools.md` §2.2;
  * `lit18/check_tools.md` §§1-2.5 and §4;
  * `lit18/README.md` §2.2;
  * `overview/overview.md` row 13, §6.7 and §7 item 5;
  * `expert_note/expert_note.md` lines 670-680;
  * Astra `answer.md` Part VIII and `reply8.txt`.
* **Not consulted.** No literature was needed for this block, and I made no web requests. Evertse's `exp(O(s))` upper bound is cited only as the record quotes it from memory (`round8/weakest.md` §5.2). I did not verify it.
