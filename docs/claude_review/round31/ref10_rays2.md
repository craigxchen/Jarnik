# Referee report, block rays2: Part X.2, second half (classes of labels, word count, deliberate failure, lift heights)

Refereed text: `astra_free_X/answer.md`, lines 2348-2404 (section X.2, read in full), the second COMPUTED
paragraph of X.5 (line 2479), summary lines 5-8 (lines 2491-2494), and `x_finite_work.md` (which contains no X.2
material; it holds the X.1 and X.4 drafts). The construction of the parent, of the block-sum-preserving local towers and
of the angles is taken as described, as instructed; its internal consistency (attainment of block sums, compactness,
P3-char preservation, ordinary gaps) is assigned to another referee and is marked "assumed" below.

Referee scripts (own code, exact integers): `round31/ref10_rays2/negcontrol.py`, `round31/ref10_rays2/digits.py`,
`round31/ref10_rays2/lift.py`. Astra's `x_local_check.py` was read; it cannot run here (it imports a missing module
`ix_factor_check`), so the negative control was reproduced independently.

## 0. Verdict

1. **The mechanism is correct.** I re-derived it. A short word whose class sums are not all zero has its exterior log
   kept away from every high divisor of the group order by a positional digit code. So every collision of every power
   `n v` happens at a prime power at most `6nG`, with `log G = o(log log q)`. Proved, given the construction (assumed).
2. **The count is correct, but it measures the wrong thing.** The untested short words are exactly the short vectors
   of the sublattice `Lambda` spanned by the same-class differences. Every word of odd mass is tested. Among the
   differences, the untested share is `≈ 1/b`. The phrase "a counterexample to a growth argument requiring even this
   fraction of the tests" is **false as worded**. What the construction defeats is any test set on which some
   partition into `b(m)` classes leaves every word with a nonzero class sum. For differences this means a `b(m)`-partite
   test graph. By Turán such a graph covers at most a `1 - 1/b` fraction of the pairs, and here `b(m) <= (log log q)^{1/32}`
   at `H = 2`.
3. **The deliberate failure is correct as an existence statement.** The instance *can be chosen* so that one same-class
   difference fails at every cutoff `X >= l_* in (4P, 8P)`. That is not "must fail" and not "cannot be repaired", as
   Astra herself says. The failing cutoff is `≍ m^3`, so the failure lies in the regime `3M^A` with `A > 3`. It is not
   at the P3-all* cutoff `3M`.
4. **The lift-height bound (X.3) is correct.** A sharper form holds: the reduced denominator `b_j` of any finite-cutoff
   rational lift satisfies `b_j > (X/(sqrt2 p_j))^{1/2}`, and `Norm gamma >= b_j` for every Hilbert-90 representation.
   Proved, given that the fake passes the singleton rays (assumed construction plus item 1).
5. **The negative control 101, 109, 1511 is reproduced** (Computed).
6. **What this half caps.** It caps "P3-char at every cutoff + ordinary gaps + e-all on all powers of all words of
   mass `<= H` with a nonzero class sum, at every unit and every cutoff", for a constructor-chosen partition into
   `b(m) -> infinity` classes. It is **not** a P3-all* fake. The input left uncapped is the same-class pair
   collisions: the pair-pigeonhole input of X.4 and of `sharp.md` §6.7.

## 1. Definitions and quantifiers (as reconstructed)

* Labels: the `s = g(q-1)` physical primes `p_j in [P, 2P]`, `P = q^3`, `q = 3 mod 4` prime, `g = 2B`, `B = 64H^2`.
  Blocks are the IX.5 replication groups of size `g`. Inside each block, the product of the exterior local ratios is
  prescribed by the IX.3/IX.5 data.
* A **signed word** is `v in Z^s \ {0}` with mass `||v||_1 <= H`. Repeated coordinates are allowed (`2e_j` has mass 2).
  The **ray** of `v` is `{n v : n >= 1}`. Its test is (1.2) of `sharp.md` Def. 1.1 at all four units `i^t` and at
  every cutoff `X`.
* **Classes**: a partition of the label set into `b` classes, chosen by the constructor. Any partition works for the
  positive statements; balanced sizes are used only for the count. The class-sum vector of `v` is
  `sigma(v) in Z^b`, with `sigma_c = sum_{class(j)=c} v_j`. The tested set is
  `V_H = {v : 1 <= ||v||_1 <= H, sigma(v) != 0}`.
* Number of classes: `T = ceil(2H log2(2b)) + 2` rows and `U = 100 g B^T`. The proof requires
  `U = o(log log q)`. Since `B^T >= (2b)^{2H log2 B}`, this forces
  `2b <= (log log q)^{1/(2H log2(64H^2)) + o(1)}`, for example `(log log q)^{1/32}` at `H = 2`. So `b -> infinity`,
  but extremely slowly.
* Order of quantifiers: first `C` and `H`; then `g, B` (functions of `H`) and the weight constant `O_{C,H}(1)`, which is
  `≍ g`, i.e. `≍ H^2`. Then `q -> infinity`, with `b = b(q)` as above. Then one parent, one compatible family of local
  factors valid at all cutoffs, and the angles (a set of positive fibre measure). All cutoffs `X` and all `n` are
  quantified last. The constant is not uniform in `H`. So words of length `H(m) -> infinity` are not capped.

## 2. Item (1): why the words of `V_H` pass, with all powers. **Proved** (given the assumed construction)

*Detection rows.* Rows `c_1, ..., c_T in {0,1}^b` exist such that every nonzero `sigma in Z^b` with
`||sigma||_1 <= H` has `c_t . sigma != 0` for some `t`. For a uniformly random row and fixed `sigma != 0`, pick `k` with
`sigma_k != 0`. Conditional on the other bits, `c . sigma` takes two values that differ by `sigma_k`, so it vanishes
with probability `<= 1/2`. The number of such `sigma` is below `(2b+1)^H`, and `2^{-T}(2b+1)^H < 1`. Correct. Computed:
`digits.py` found detecting rows for `(H,b) = (2,3), (2,4), (3,3)`, which detect all 24, 40 and 62 class-sum vectors.

*Digit code.* Fix an exterior odd prime `l`, a level `a` and a prime power `r^e || d_a`, where
`d_a = (l - chi_4(l)) l^{a-1}` is the order of the cyclic norm-one group. Let `k_r` be minimal with `r^{k_r} >= U`.
If `e >= k_r`, put `D = r^{k_r}` and require the log `L_j` mod `D` to lie in a cell
`[Dz/B^T, D(z+1)/B^T]` whose base-`B` digits are `z_t in {a_t, a_t + 1}`, with `a_t = 16H (c_t)_{class(j)}`. Then
`B^{t-1} L_j / D mod 1` lies in `[a_t/B, (a_t+2)/B]`. For `v in V_H`, choose a detecting `t` and put
`u = c_t . sigma(v)`, so `1 <= |u| <= H`. Then `B^{t-1} (sum v_j L_j)/D` is congruent mod 1 to
`u/(4H) + err`, with `|err| <= 2H/B = 1/(32H) < 1/(4H) <= dist(u/(4H), Z)`. So `D` does not divide `sum v_j L_j`,
and the `r`-part of `gcd(sum v_j L_j, d_a)` is at most `r^{k_r - 1} < U`. If `e < k_r`, then `r^e < U` trivially.
Hence `gcd(sum v_j L_j, d_a) | G = lcm(1..ceil U)`, which is (X.2). Computed: `digits.py` verified the cell window for
every label and digit, and on about 8,400 random tested words it found `D` never dividing the sum. It also checked the
quantitative distance bound `>= 1/(4H) - 2H/B`, with `D` a power of 3 that is not a multiple of `B^T`.

*Powers.* A collision of `n v` with `i^t` at `l^a` means `n sum v_j L_j = t d_a/4 mod d_a`, which requires
`4 | d_a` (true for `l` odd). So `d_a | 4n sum v_j L_j`, hence `d_a | 4n gcd(sum v_j L_j, d_a) | 4nG`, and
`l^a <= (3/2) d_a <= 6nG`. The colliding prime powers are distinct-prime powers, each `<= 6nG`, so their product
divides `lcm(1..floor(6nG))` and has `log = O(nG)` uniformly in the cutoff. This holds for each unit separately.

*Angles.* `||v||_1 <= H < g`, so on each block `v` is either zero or nonconstant. A nonconstant character of a block
torus pushes Haar measure forward to Haar measure, also after multiplying by `n`. The forbidden measure for `n v` is
`O(e^{O(nG)} P^{-n h/2})`. Summing over the `(2gq)^h` words and over `n` gives `o(1)`, because
`e^{O(G)} = q^{o(1)}` and `P^{1/2} = q^{3/2}`. Correct.

The mechanism therefore controls exactly the words with `sigma(v) != 0`. For `v in Lambda := ker sigma`, the
same-class differences, nothing is controlled. Two labels of one class have the same cells, and their log difference
can be divisible by any `D`.

## 3. Item (2): the count, and which words matter

**Proved (exact form).** A word is untested iff `sigma(v) = 0`, i.e. iff `v in Lambda`, the lattice spanned by the
`e_j - e_k` with `class(j) = class(k)`. Since `sum_j v_j = ||v||_1 mod 2`, **every word of odd mass is tested**. For
mass 2 the untested words are exactly `±(e_j - e_k)` within a class. Their share of all `2s + 4 C(s,2)` signed words of
mass 2 is `≈ 1/(2b)`. Their share of the `C(s,2)` difference rays is `sum_c C(n_c,2)/C(s,2) <= 1/b + o(1)` for
balanced classes. For even `h` the untested share is `≍_h b^{-h/2}`. All of this is consistent with Astra's
`O(1/b + 1/s)`.

**Computed** (`digits.py`, exact generating-function count over all signed words):

| `s`, `b` | h=2 | h=4 | h=6 | odd h | same-class share of differences |
|---|---:|---:|---:|---:|---:|
| 240, 4 | 0.12292 | 0.03966 | 0.01857 | 0 | 0.24686 |
| 240, 8 | 0.06042 | 0.01027 | 0.00272 | 0 | 0.12134 |
| 240, 16 | 0.02917 | 0.00247 | 0.00034 | 0 | 0.05858 |
| 960, 32 | 0.01510 | 0.00067 | 0.00005 | 0 | 0.03024 |

**False as worded; corrected.** The sentence "This is a counterexample to a growth argument requiring even this
fraction of the tests" (line 2356) does not hold for arbitrary test sets of that density. Corrected statement:
*for fixed `H`, let a growth argument use the e-all tests on a set `E_m` of words of mass `<= H`. If there is a
partition of the labels into `b(m)` classes for which every word of `E_m` has a nonzero class sum, the argument is
capped by X.2 (with `b(m)` as in §1).* For differences, `E_m` must then be a `b(m)`-partite graph. By Turán it has at
most `(1 - 1/b(m)) s^2/2` edges. So an argument using a fraction `1 - o(1/b(m))` of the differences, or any test graph
with chromatic number `> b(m)`, is not capped. Since `b(m) <= (log log q)^{1/32}` at `H = 2`, the uncapped fraction
`1/b` decays more slowly than any power of `log log m`.

**Which words matter.** The mark-character words (row differences) are block-constant, of mass `≍ q > H`, and covered
by P3-char, not by `V_H`. A pair-pigeonhole argument picks the pairs whose residue tuples coincide, as X.4 does. Those
pairs depend on the data, and the fake places all uncontrolled coincidences inside classes. So the density of tested
words is not the relevant measure. The relevant question is whether the data-dependent pigeonhole pairs can be forced
across classes, and nothing here bears on it.

## 4. Item (3): the deliberate failure. **Proved as an existence statement**

The statement checks out. The prime `l_*` lies in `(4P, 8P)` (Bertrand), so it is exterior, because the parent primes
lie in `[P, 2P]`. A largest class has size `>= s/b > 2g`, so it meets at least three blocks, and two of its labels
`j, k` from distinct blocks have equal required parity at `l_*`. The parity of the ratio log is the Legendre symbol
`(p_j/l_*)`, for split and inert `l_*` alike. Computed: equal ratios mod `l` are attainable iff the Legendre symbols are
equal, with 0 mismatches over 11 small primes and 28 label pairs (`negcontrol.py`). Equal cells allow a common value.
The remaining `g - 1 >= B` slots in each affected block still attain the block sums (assumed, as in the construction).
At unit one, the collision of `e_j - e_k` contains `l_* > 2P >= sqrt(p_j p_k)`. Test (1.2) would need
`|sin(phi_j - phi_k)| >= l_*/sqrt(p_j p_k) > 1`, so it fails for every angle at every cutoff `X >= l_*`. The positive
tests are unaffected, because `e_j - e_k` lies in `Lambda` and every word of `V_H` uses only cell membership.

**What is proved and what is not.** It is proved that *some* instance of the construction fails e-all on one omitted
difference while passing everything listed in X.2. It is not proved that every instance fails, nor that the
construction cannot be modified so that all same-class differences pass. Whether they can is the P3-all* question
itself. Summary line 7 ("can be forced to fail") should read **"can be chosen to fail"**.

**Correction on regime.** The failure is seen only at `X >= l_* > 4q^3 ≍ m^3`. So the failing instance fails e-all at
cutoff `3M^A` only when `A >= 3` (plus a constant). It is *not* shown to fail at P3-all*'s cutoff `3M ≈ 1.5 rho q`.
**Heuristic (referee):** the same tie argument at four exterior primes in `(q, 1.48 q)` should move the failure below
`3M`. Their product exceeds `q^4 > 2P`. A largest class meets more than `16` blocks, so pigeonhole on the 16 parity
patterns finds a usable pair. I have not checked this in detail.

## 5. Item (4): finite-cutoff lift heights on the 5 mod 8 version. **Proved, with a sharper form**

Logic. Suppose `r = pi gamma/conj gamma = alpha/b` (lowest terms, `Norm pi = p`) agrees with the data at all exterior
odd prime powers `<= X` that do not divide `b`. Since all parent primes are `5 mod 8`, X.1 applies with `n = 2`: no
parent prime `l` can have `ord_l(alpha/conj alpha) = 8`, because `8` does not divide `l - 1`. The integers
`x_2 ∓ y_2` (from `alpha^2`) are odd, nonzero, prime to `J = p b^2` and to every parent prime, and
`|x_2 ∓ y_2| <= sqrt2 p b^2`. If `X >= sqrt2 p b^2`, the two e-all tests on `2 e_p` at the units `i` and `-i` charge
them fully. Squaring and adding gives `1 >= b^4`, so `b = 1`. The fake passes these two tests at every cutoff, because
singletons lie in `V_H`. Then `r` is an associate of `pi` or of `conj pi`. The singleton tests at `n = 1` (cutoff
`>= sqrt(2p)`, with `|x ± y| < sqrt(2p) < P`) give `sin 2phi = sin 2arg r`, so `phi` is a genuine direction up to units
and conjugation. The angles were chosen off that finite set, which is a contradiction.

**Corrected statement.** Fix the 5 mod 8 instance (parent, local data and angles). For every label `p_j`, every cutoff
`X >= sqrt(2 P_max)` and every rational `r` with `r conj r = p_j` that agrees with the data at all exterior odd prime
powers `<= X` prime to its reduced denominator `b_j`, we have `b_j > (X/(sqrt2 p_j))^{1/2}`. Since `b_j <= Norm gamma`
for every `gamma` with `r = pi gamma/conj gamma` (because `r = pi gamma^2/Norm gamma`), (X.3) follows:
`Norm gamma > (X/(sqrt2 p_j))^{1/2}`. Such lifts exist at each `X` (CRT and local Hilbert 90), so the least lift height
tends to infinity with `X` on the fixed instance.

Computed (`lift.py`): on all 36 primes `5 mod 8` in `[1000, 2000]` I tested 20,000 random lifts (19,325 of them
nonintegral). In every case `x_2 ∓ y_2` was odd, nonzero, prime to `J` and to all parent primes, and bounded by
`sqrt2 p b^2`, and the two forced sine bounds square-summed to exactly `b^4`.

## 6. Item (5): what this half establishes about P3-all*

* **Not a P3-all* fake.** The construction does not establish e-all on `Lambda` at any cutoff. The deliberately failing
  instance fails it at `X >= l_*`. Astra says this herself (lines 2356, 2395).
* **Capped class (Proved, given the construction).** Fix `C` and `H`. Then there are fakes with `W = O_{C,H}(m log m)`
  satisfying all of the following:
  * P3-char at every cutoff, including the same-level clause;
  * all ordinary gaps;
  * e-all at every unit and every cutoff, for all powers `n v` of all `v` with `||v||_1 <= H` and
    `sigma(v) != 0`, for a partition into `b(m) -> infinity` classes chosen by the constructor.

  This includes all singleton rays, all pair sums, all cross-class differences, all differences to reserved singleton
  classes, and all odd-mass words.
* **Is the class natural?** Partly. The singleton rays at all cutoffs are exactly the input of the finite-height
  rigidity route (X.1 and IX.7; README_IX §5, item 1, first step). Together with (X.3), this half shows that the route
  needs the lift height as an extra cost: the singleton rays alone force nothing about `W`. That is the genuinely
  useful content. The cross-class short words do not enter any known growth argument.
* **What is left.** The same-class differences and the lattice `Lambda` they span. These are the pair collisions of
  labels across many primes, the input of X.4 and the "pair pigeonhole" of `sharp.md` §6.7. Long words that are not
  block-constant are also untested. The "limit of the order argument" paragraph (line 2404) is correct (Proved). The
  odd prime power `r^k < s` divides the order of the norm-one group mod `r^{k+1} < s^2`, which is exterior since
  `P > s`. Pigeonhole on `s` logs gives a pair that agrees mod `r^k`, so any common `G` has
  `log G >= psi_odd(s - 1) ~ s`. This blocks only the common-index method, not a pair-by-pair construction.

## 7. Item (6): the X.5 negative control. **Computed, reproduced**

The script `negcontrol.py` checked the following.

* `1511` is prime and `3 mod 4`, so it is inert and is not one of the 80 fine labels (the primes `1 mod 4` from 101 to
  1117).
* It is the first prime `>= 1500` with `(101/l) = (109/l)` (both equal `+1`).
* An exhaustive enumeration mod 1511 shows that the norm-101 and norm-109 units have the same 756 ratio classes, i.e.
  the same coset of index 2 in the norm-one group of order 1512.
* Explicitly, `r_101 = 143 + 143i` and `r_109 = 351 + 351i` mod 1511 satisfy
  `(r_101/conj r_101)/(r_109/conj r_109) = 1`.
* Hence at every cutoff `X >= 1511` the unit-one test on `e_101 - e_109` requires
  `|sin| >= 1511/sqrt(11009) = 14.40`, which is impossible for any angles.

Astra's assertion `ell*ell > p_j p_k` is equivalent to `ell > sqrt(p_j p_k)`. The five-mark design bookkeeping in the
same paragraph (90 pair identities, 480 level identities, 56,440 order checks) was not re-run, because the script's
dependency is missing.

## 8. Labels of the refereed statements

| Statement | Label |
|---|---|
| Words in `V_H` pass with all powers at every unit and every cutoff (mechanism (X.1)/(X.2)) | Proved (given the construction) |
| Detection rows exist with `T = ceil(2H log2 2b) + 2` | Proved; Computed for small `(H,b)` |
| `b -> infinity` with `U = o(log log q)` | Proved; quantitatively `b <= (log log q)^{1/(2H log2 64H^2)}` |
| Fraction `1 - o(1)` of the words at each mass; `>= (1 - 1/b - o(1)) C(s,2)` difference rays | Proved; Computed exactly; odd mass fully tested |
| "Counterexample to a growth argument requiring even this fraction" | False as worded; corrected in §3 (partition / `b`-partite form) |
| An omitted difference can fail below `8P` | Proved as existence ("can be chosen"); at cutoff `≍ m^3`, not `3M` |
| Lift height bound (X.3) on the 5 mod 8 version | Proved; sharper denominator form in §5; Computed transfer |
| Limit of the order argument (line 2404) | Proved |
| X.5 negative control 101/109/1511 | Computed, reproduced |
| X.2 is not a P3-all* fake; P3-all* open | Correct (Heuristic status statement) |
