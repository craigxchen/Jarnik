# Round 21: fake circles in the low-height window P_+(N) < M^(2-delta)

Round 21, target "close the window from the other side". Scripts and outputs: `fakes_checks/`. Labels as in the
overview (proved / exact / evidence / heuristic). Nothing here has been refereed.

## 0. Summary and verdict

**Verdict: partial. No fake circle with all primes `<= M^(2-delta)` was constructed, and nonexistence was not proved.
What is proved is why every known construction must fail there, which construction-independent inputs cannot
separate window fakes from actual clusters, and what a window fake would have to be.** No growth statement for
actual circles is claimed.

1. **The operative obstruction is the Haar measure of the gap-admissible set, not the window
   `[P, P e^(1/(4r))]`.** The capped prime-count input (`round8/stability.md` Thm 3.3: `omega(N') > M/4`) forces
   `Lambda = sum_{p|N} p^(-1/2) > M^(delta/2)/4` in the window (§1.2). **Theorem A (proved):** on every realisation
   coset of every profile, the Haar proportion of gap-admissible angle vectors is at most
   `exp(-2.17 sum_{j in J} p_j^(-1/2))`, `J` any set of coordinates independent modulo the character span (size up
   to the free dimension `m`). So Haar completion, the second step of every barrier-scale fake on record, needs
   `P_+ > 9.8 m^2` to succeed with probability `1/2` (`9.8 (M-1)^2` for the double core, against the `84 (M-1)^2` used
   by `outside.md` A.1(ii)), and in the window it succeeds with probability `<= exp(-c M^(delta/2))` (Cor. A.2).
   The equal-weight window is a design feature that grouping removes (§2; heuristic remark).
2. **Arithmetic points do not help (proved).** Angles `arg beta_j` of primitive Gaussian integers of odd,
   pairwise coprime norms `<= p_j` are always gap-admissible, but such a "transplant" model realises a profile iff
   its points form an **actual** cluster on a circle of norm `N' <= N` (Lemma T). A barrier-scale transplant fake
   for unbounded `M` would disprove Cilleruelo–Granville (Corollary T).
3. **The large sieve is implied by the gap axioms (proved).** Every angle model, whatever its prime sizes, satisfies
   the dual large sieve for its own prime angles with the genuine constant `2K + pi X/2`, and for all Gaussian
   prime angles (cluster angles replaced) with constant `2(2K + pi X/2)` (Theorem LS). So the large sieve at
   `K ~ P^(1/2)`, although non-trivial for the cluster in the window, cannot exclude a window fake that the
   monomial gaps do not already exclude; the sharp (`c < 2`) form is the unclassified cross-gap item of overview §5.5.
4. **Complete prime segments are rigid (proved); sparse ones are not (evidence).** `G_1(5) = {arctan(1/2), arctan 2}`
   (Proposition B, from `arctan(1/2) + arctan(1/3) = pi/4`), and more generally **Theorem C**: if every split prime
   `< p` divides `N` and their angles are genuine, the gap inequalities force `phi_p = ± arg pi_p (mod pi/2)` (a
   Minkowski argument plus two tight lattice lines; multi-prime form of round-10 Theorem D). So an angle model is
   genuine on every complete initial segment of split primes dividing `N` (Corollary C; checked independently by a
   50-digit certificate for `13, ..., 53` in all 63 conjugation branches). For **sparse** prime sets, as in the window,
   the gap-admissible set has **positive** measure in every case tested, obeying a Poisson law
   `mu(G(p)) ≈ exp(-c(Pi - 1))`, `c ∈ [1.1, 1.8]` (doubly exponential in `Lambda`); its points are biased towards, but
   not concentrated on, arithmetic points; and a sequential (coordinate-by-coordinate) construction finds points,
   verified on all 360,050 monomials of weight `<= 34`, at `Pi = 17.4` (Haar sampling would need `~e^30` trials), and
   at `Pi = 26.4` (all 1,043,472 monomials of weight `<= 33`; `~e^45` trials); it dies at `Pi = 39` (§5). On the full torus, positivity is **proved** up to `Lambda ~ (1/4) log P` (Proposition A'), a `log^2` gain over
   the `P >~ r^2` of the union bound, but far from the window.

**What decides the window (§6.3).** Everything reduces to one question (Q): on realisation cosets of barrier-scale
profiles with `Lambda -> infinity`, is the gap-admissible set non-empty away from arithmetic points? If yes
(*scenario P*, what the Poisson law extrapolates to), window fakes exist and the window is capped, though no
positive-probability recipe produces them. If no (*scenario R*, what provably happens on complete prime segments, Theorem C), window angle
models are transplant-like, hence actual clusters (Lemma T), and "no window fakes" is a restricted growth theorem
for actual clusters with `P_+(N) < M^(2-delta)`. The evidence does not decide between them. In both scenarios the
large sieve (Theorem LS) and measure bounds (Theorem A) are useless for closing the window: a proof must decide
(Q) from the gap structure of the cluster's Gaussian primes. Proposed wording for overview §5.5: *the low-height
window is uncapped; the record's fake constructions provably fail there (Theorem A, the `r^2` law is intrinsic to
Haar completion); large-sieve inputs are implied by monomial gaps (Theorem LS); arithmetic fakes are actual
clusters (Lemma T); whether non-arithmetic gap-admissible points survive on realisation cosets as
`Lambda -> infinity` is the open question (Q).*

## 1. Setting: angle models, capped input classes, the window

**1.1 Angle models (overview §5.1).** A profile is `P = (p, e, (a_x)_{x<=M})`: distinct primes
`p_1, ..., p_r = 1 mod 4`, exponents `e_j >= 1`, distinct rows `a_x in prod_j [0, e_j]`. Put
`u_x = 2a_x - e in Z^r`, `w(v) = sum_j |v_j| log p_j`, `W = sum_j e_j log p_j`, `Delta = C e^(-W/4)`.
An angle model is a profile with a real vector `phi` that

* *realises* it: `<u_x, phi> = theta_0 + delta_x + (pi/2) k_x` with `delta_x in [0, Delta]`, `k_x in Z` (overview (5.1));
* is *gap-admissible*: for every `v in Z^r \ {0}`,
  `dist(<v,phi>, (pi/2)Z) >= arcsin(e^(-w(v)/2))` and `dist(<v,phi>, pi/4 + (pi/2)Z) >= arcsin(e^(-w(v)/2)/sqrt 2)`
  (overview Lemma 5.1).

Throughout, `T := R/(pi/2)Z` with normalised Haar measure `mu_1`, and `G(p) ⊂ T^r` is the set of gap-admissible
vectors. Both conditions descend to `T^r`: replacing `phi_j` by `phi_j + pi/2` changes `<v, phi>` by
`(pi/2) v_j`, which changes neither distance, and changes `<u_x - u_1, phi>` by an element of `(pi/2)Z`, which is
absorbed into `k_x - k_1`. Conversely, if `<u_x - u_1, phi> = delta_x - delta_1 (mod pi/2)` for all `x`, put
`theta_0 = <u_1, phi> - delta_1`; then (5.1) holds with suitable `k_x`. So, for fixed `delta`, the realisation set is

```text
Real(delta) = { phi in T^r : <u_x - u_1, phi> = delta_x - delta_1  (mod pi/2),  2 <= x <= M }.
```

Let `L = span_Z{u_x - u_1}`, `L_R = L ⊗ R`, `l = dim L_R`, and `m = r - l` (the **free dimension**). `Real(delta)` is
empty or a coset of the closed subgroup `H = {phi : <lam, phi> = 0 mod pi/2 for all lam in L}`, whose identity
component `H_0` is an `m`-dimensional subtorus; `H/H_0` is finite.

**1.2 What "every capped input class" forces.** The fakes asked for must satisfy every row of overview §5.4 that
is capped. Row 3 includes the prime-count statement of `round8/stability.md` Thm 3.3(3): `omega(N') > M/4` for
`C <= 1` (and `> W'/(4W'/M + 8 log^+ C)` in general), and Thm 3.3(2): every prime power has `e log p < 4W'/M +
8 log^+ C`. So a fake in the class has `r > M/4` primes. If all of them are `<= P := M^(2-delta)`, then

```text
Lambda := sum_j p_j^(-1/2) > (M/4) M^(-1 + delta/2) = M^(delta/2)/4,
Pi     := prod_j (1 + 2/(sqrt(p_j) - 1)) >= exp(Lambda) > exp(M^(delta/2)/4).                 (1.1)
```

(`log(1 + 2/(sqrt p - 1)) >= 1/sqrt p` for `p >= 5`, since `2/(s-1) >= 2/s` and `log(1+y) >= y/(1+y)`
with `y = 2/(s-1)`, `y/(1+y) = 2/(s+1) >= 1/s` for `s >= 1`.) Without the prime-count input a loophole "few primes
with large exponents" would have to be examined: the cut-slack identity together with a Gram-rank count only gives
`r >~ sqrt(M)` in my estimate (*heuristic remark*, not used), and with `r <= M^(1-delta/2)` primes the union bound of
§2 would work in the window. The prime-count input closes this loophole; it is a capped input (row 3), so the
fakes asked for must satisfy it.

`Pi - 1` is exactly the quantity whose smallness the known constructions need (§2), and `Lambda` is the quantity
that the obstruction of §3 sees. **The window is precisely the regime in which `Lambda` is a positive power of
`M`.** In the known fakes `Lambda = O(1)` (primes `>= 21 r^2`), or even `O(1/M)` (`lw.md` 6.1, `P >= (40M)^4`).

**1.3 The two meanings of "no fake is known in the window".** Round 10 §6 item 1 says no fake with all primes
`<= M^(2-delta)` is known, and explains this by the equal-weight windows `[P, P e^(1/(4r))]` of the known fakes,
which need `P >~ r^2 log r`. §2 shows that this explanation is incomplete: the equal-weight window is a feature
of the particular combinatorial designs and can be removed; the obstruction that cannot be removed by changing the
design is the measure of `G(p)`, which `Lambda` controls (§3).

## 2. Why the known fakes need P >~ r^2 log r

Every barrier-scale fake on record (`round8/lw.md` Thm 6.1, `round5/construct.md`, `round5/lower.md`,
`round5/round6_sharp/sharp.md`, the grafts of `round10/analytic.md` Thm B) is built in two steps.

1. **Combinatorics.** A design (skew double core, fully flipped Paley core, parity-class design) with all cut
   weights equal up to `O(1/r)`, so that pair and character margins computed with equal weights survive. This
   is where the window `[P, P e^(1/(4r))]` comes from: `r` split primes in an interval of multiplicative width
   `e^(1/(4r))` need, by the prime number theorem in progressions, `P/(8 r log P) >~ r`, i.e. `P >~ r^2 log r`.
2. **Haar completion.** Fix `delta`, draw `phi` from Haar measure on the realisation coset (`round4/outside.md`
   Thm A.1(i)), and bound the probability of a gap failure. Character monomials (`v in L_R`) are handled
   deterministically by the margins. Non-character monomials are handled by the union bound

   ```text
   P(some non-character v fails) <= 1.33 (Pi - 1),     Pi = prod_j (1 + 2/(sqrt(p_j) - 1))        (2.1)
   ```

   (`outside.md` Thm A.1(ii); `lw.md` Thm 6.1(b)), which is `< 1` only if `sum_j p_j^(-1/2) = O(1)`, i.e.
   `P >~ r^2`. `outside.md` uses `p_j >= 21 r^2`; `lw.md` 6.1 uses `P >= (40M/C_0)^4`.

**The two requirements are of different nature.**

* Step 1 is a property of the designs, not of the problem. Equal weights can be traded for **grouped**
  weights: give each cut a block `G_j` of `k` primes with identical columns, and choose the blocks so that the
  block weights `sum_{p in G_j} log p` agree to within `1/(8r)`. The number of `k`-subsets of the split primes
  `<= Y` is `~ (Y/(2 log Y))^k/k!`, while the block weights lie in an interval of length `k log Y`; so for `k >= 2`
  and `Y = M^(1+eps)` the pigeonhole principle offers far more than `r` blocks in a common weight window.
  (*Heuristic remark*: I did not prove that `r` *disjoint* blocks can be extracted, and nothing below uses it.)
  Characters of the grouped profile are images of characters of the ungrouped design, so the margin
  computations of step 1 carry over verbatim. **The window structure is therefore not the obstruction.**
* Step 2 is not a property of a design. It fails in the window for every profile whose free dimension `m` is
  large compared with `sqrt(P)` (the double core has `m = M - 1`), and Theorem A below shows that the failure is not an artefact of the union bound: the Haar probability of success
  is itself at most `exp(-2.17 sum_{j in J} p_j^(-1/2))`. This is the real content of "`P >~ r^2`".

So of the three explanations offered in the task (counting of Gaussian prime angles, the large sieve, the window
structure), none is the operative one for the known method: the operative one is **the Haar measure of the
gap-admissible set**, which collapses once `Lambda = sum p_j^(-1/2)` is large (§3, §5). The large sieve is
implied by the gap axioms and adds nothing (§4). Counting of Gaussian prime angles enters only through the
arithmetic ("transplant") points of `G(p)`, which are actual clusters in disguise (§3.3).

## 3. The obstruction: the gap-admissible set in the window

### 3.1 Theorem A: Haar scarcity on realisation tori (proved)

For a prime `p = 1 mod 4` let

```text
G_1(p) = { t in T : dist(n t, (pi/2)Z) >= arcsin(p^(-n/2)),  dist(n t, pi/4 + (pi/2)Z) >= arcsin(p^(-n/2)/sqrt 2)  for all n >= 1 },
```

the set of admissible values of one coordinate when only the monomials `n e_j` are imposed.

> **Theorem A.** Let `P` be any profile, `L_R`, `H`, `H_0`, `m` as in §1.1, and let `J ⊂ {1, ..., r}` be a set of
> indices such that `{e_j : j in J}` is linearly independent modulo `L_R`. Let `nu` be the Haar probability
> measure on a coset `g + H_0` of `H_0` in `T^r` (every non-empty `Real(delta)` is a finite union of such
> cosets). Then
>
> ```text
> nu( G(p) ∩ (g + H_0) )  <=  prod_{j in J} mu_1(G_1(p_j))  <=  exp( - kappa sum_{j in J} p_j^(-1/2) ),
> kappa = (4/pi)(1 + 1/sqrt 2) = 2.1735...
> ```
>
> One can take `|J| = m`, so the right side is at most `exp(-kappa m P_+(N)^(-1/2))`.

*Proof.* (a) *The projection `pi_J : H_0 -> T^J`, `phi -> (phi_j)_{j in J}`, is onto.* Its image is a compact
connected subgroup of `T^J`. If it were proper, some non-trivial character `psi -> exp(4i sum_{j in J} n_j psi_j)`
(`n in Z^J \ {0}`) of `T^J` would be trivial on it, so the character `chi_v(phi) = exp(4i <v, phi>)` of `T^r` with
`v = sum_{j in J} n_j e_j` would be trivial on `H_0`. The characters of `T^r` are the `chi_v`, `v in Z^r`.
By Pontryagin duality the annihilator of `H = L^perp` is `L`. `H/H_0` is finite (a closed subgroup of a torus
has finitely many components), so the annihilator of `H_0` contains `L` with finite index; hence every `v` with
`chi_v` trivial on `H_0` has a non-zero multiple in `L`, i.e. `v in L_R`. This contradicts the independence of
`{e_j}_{j in J}` modulo `L_R`.

(b) *Under `nu`, the coordinates `phi_j`, `j in J`, are independent and `mu_1`-uniform.* The image of Haar
measure on `H_0` under the continuous surjective homomorphism `pi_J` is a translation-invariant probability
measure on `T^J`, hence Haar. On the coset, `phi_J = g_J + pi_J(h)` with `h` Haar on `H_0`: a translate of Haar.

(c) *Reduction.* The gap conditions for the monomials `v = n e_j` (`n >= 1`) read `phi_j in G_1(p_j)`, because
`w(n e_j) = n log p_j`. So `G(p) ∩ (g + H_0) ⊂ {phi : phi_j in G_1(p_j), j in J}`, whose `nu`-measure is
`prod_{j in J} mu_1(G_1(p_j))` by (b).

(d) *One coordinate.* `G_1(p)` misses the two open windows of `n = 1`: `{dist(t, (pi/2)Z) < arcsin p^(-1/2)}` and
`{dist(t, pi/4 + (pi/2)Z) < arcsin (2p)^(-1/2)}`. They are disjoint in `T` because
`arcsin p^(-1/2) + arcsin (2p)^(-1/2) <= arcsin 5^(-1/2) + arcsin 10^(-1/2) = arctan(1/2) + arctan(1/3) = pi/4` for
`p >= 5`. Their total `mu_1`-measure is `(4/pi)(arcsin p^(-1/2) + arcsin (2p)^(-1/2)) >= kappa p^(-1/2)`
(`arcsin y >= y`). Hence `mu_1(G_1(p)) <= 1 - kappa p^(-1/2) <= exp(-kappa p^(-1/2))`.

(e) *Size of `J`.* The images of `e_1, ..., e_r` span `R^r/L_R`, of dimension `m`, so some `m` of them form a
basis. Each term is `>= P_+^(-1/2)`. QED

**Remarks.** (1) The bound uses only the `r` one-prime monomials `n e_j`; the evidence of §5.1 shows that the
true measure is much smaller (it behaves like `exp(-c(Pi - 1))`, doubly exponential in `Lambda`). (2) A second
moment argument (sketch, not used) over all non-character `v` (pairwise independent under `nu` unless parallel modulo `L_R`) gives
`nu(G ∩ coset) <= 1/E`, `E` the expected number of failures, which is of the same order `exp(-O(Lambda))`; I record
only the product form because it is exact and short. (3) `delta` and the windings play no role: the bound holds
for every coset, hence also after averaging over any distribution of `(delta, k)`.

> **Corollary A.1 (the `r^2` law is intrinsic to Haar completion).** Suppose a construction picks `phi` from Haar
> measure on a realisation coset (for any rule choosing `delta`, `k` and the coset) and produces a gap-admissible
> vector with probability `>= eta`. Then for every `J` as in Theorem A,
> `sum_{j in J} p_j^(-1/2) <= kappa^(-1) log(1/eta)`; in particular `m <= kappa^(-1) log(1/eta) P_+^(1/2)`.
> For the skew double core (`r = 2(M-1)`, `l = M - 1`, so `m = M - 1`) and `eta >= 1/2` this gives
> `P_+ >= (kappa (M-1)/log 2)^2 > 9.8 (M-1)^2`, against the `84 (M-1)^2 = 21 r^2` that `outside.md` A.1(ii) uses.

So the method of the record needs primes of size `>= c M^2` for a structural reason, with the constant
determined to within a factor `~ 9`, and the equal-weight window is irrelevant to it.

> **Corollary A.2 (the window).** Let all primes be `<= M^(2-delta)` and let the free dimension satisfy
> `m >= c M`. Then every realisation coset carries gap-admissible vectors of Haar proportion at most
> `exp(-kappa c M^(delta/2))`. In particular Haar completion with any success probability
> `>= exp(-o(M^(delta/2)))` is impossible in the window.

*Proof.* Theorem A with `|J| = m`, `p_j^(-1/2) >= M^(-1 + delta/2)`. QED

**The matching lower bound, sharpened.** The union bound (2.1) proves `mu(G(p)) > 0` only when `Pi - 1 < 0.75`, i.e.
`Lambda <~ 1/2`. Integrating one coordinate at a time does better on the full torus:

> **Proposition A' (sequential completion; proved).** Order the primes `p_1 < ... < p_r` (all `>= 13`) and put
> `Pi_{<k} = prod_{i<k} (1 + 2/(sqrt(p_i) - 1))`, `kappa_1 = 2.203`. If `kappa_1 Pi_{<k} < sqrt(p_k) - 1` for every `k`,
> then `mu(G(p)) >= prod_k (1 - kappa_1 Pi_{<k}/(sqrt(p_k) - 1)) > 0`. In particular, if all primes lie in
> `[P/2, P]` with `P >= 200`, then `mu(G(p)) > 0` as soon as
> `2 Lambda (1 + 2 sqrt(2/P)) < log((sqrt(P/2) - 1)/kappa_1)`, i.e. for `Lambda <= (1/4 - o(1)) log P`.

*Proof.* Classify `v != 0` (up to sign) by its last non-zero coordinate `k`, normalised so that `v_k >= 1`. The
conditions for such `v` involve only `phi_1, ..., phi_k`; let `A_k(phi_{<k})` be the set of `phi_k` satisfying all of
them. By Fubini, `mu(G) >= prod_k inf_{phi_{<k}} mu_1(A_k(phi_{<k}))`. For fixed `phi_{<k}` and `v = (v', n, 0, ...)`, the
set of bad `phi_k` is `{t : n t + <v', phi_{<k}> in window}`, whose `mu_1`-measure equals that of the windows,
`(4/pi)(arcsin e^(-w/2) + arcsin(e^(-w/2)/sqrt 2)) <= kappa_1 e^(-w/2)` (for `e^(-w/2) <= 13^(-1/2)`,
`arcsin y <= 1.0133 y`, and `(4/pi)(1 + 1/sqrt 2) 1.0133 < 2.2025 < 2.203`), whatever `<v', phi_{<k}>` is. Summing over
`n >= 1` and `v' in Z^(k-1)`: `sum_n p_k^(-n/2) = 1/(sqrt(p_k) - 1)` and `sum_{v'} e^(-w(v')/2) = Pi_{<k}`. For the last
claim, `Pi_{<k} <= Pi <= exp(sum_j 2/(sqrt(p_j) - 1))` and, for `p_j >= P/2 >= 100`,
`2/(sqrt(p_j) - 1) <= (2/sqrt(p_j))(1 + 2 p_j^(-1/2)) <= (2/sqrt(p_j))(1 + 2 sqrt(2/P))`, so
`kappa_1 Pi_{<k} <= kappa_1 exp(2 Lambda (1 + 2 sqrt(2/P))) < sqrt(P/2) - 1 <= sqrt(p_k) - 1`. QED

So the threshold for *provable* existence of Haar-thick gap-admissible sets moves from `Lambda ~ 1` to
`Lambda ~ (1/4) log P`, i.e. from `P >~ r^2` to `P >~ (r/log r)^2`. This is still far from the window, where
`Lambda > M^(delta/2)/4`. (Proposition A' is stated on the full torus. On a realisation coset the same argument
runs over coordinates of `H_0`, and the sums are over monomials grouped by their image in `Z^m`; I did not work
out the resulting condition.)

**What Theorem A does not say.** It bounds a measure, not a cardinality. `G(p)` contains the genuine angle
vector (measure zero in general), and Theorem A cannot exclude that a realisation coset meets `G(p)` in a
non-empty null set, or in a set of tiny positive measure (which §5.4–5.5 suggest is the typical situation). Nor does it say anything when `m` is small: if `rank L = r` (`m = 0`), `Real(delta)` is finite
and the angles are determined by `(delta, k)` up to a finite group; such "rigid" realisations make the angles
rational multiples of `pi/2` (denominator `[Z^r : L]`) up to `O(||lambda|| Delta)`, and their existence is a
discrete problem about characters of `Z^r/L` that no measure argument decides (§6.2).

### 3.2 Proposition B and Theorem C: complete prime segments are rigid (proved)

> **Proposition B.** `G_1(5) = { arctan(1/2), arctan 2 }` in `T = R/(pi/2)Z`. Consequently, in **every** angle
> model with `5 | N`, the angle of the prime 5 is genuine: `phi_5 = ± arg(2 + i)  (mod pi/2)`, i.e. `phi_5` is
> the argument of a Gaussian prime above 5 up to a unit.

*Proof.* `arcsin 5^(-1/2) = arctan(1/2)` and `arcsin 10^(-1/2) = arctan(1/3)`, and `arctan(1/2) + arctan(1/3) = pi/4`.
The open `n = 1` windows in `[0, pi/2)` are `[0, arctan(1/2))`, `(pi/4 - arctan(1/3), pi/4 + arctan(1/3)) =
(arctan(1/2), arctan 2)` and `(arctan 2, pi/2)`. Their union is `[0, pi/2)` minus `{arctan(1/2), arctan 2}`. Both points
lie in `G_1(5)`: `arctan(1/2) = arg(2 + i)`, `arctan 2 = arg(1 + 2i)`, and for every `n` the Gaussian integer
`(2+i)^n = X + iY` is primitive of odd norm `5^n` (5 is split and `(2+i), (2-i)` are coprime), so `|X|, |Y|, |X ± Y| >= 1`,
which is the gap inequality (overview Lemma 5.1). `arctan 2 = pi/2 - arctan(1/2) = -arctan(1/2) (mod pi/2)`. QED

For any other single prime, `mu_1(G_1(p)) > 0` (`0.3019` at `p = 13`, `0.3927` at `p = 17`; `gap1d.py`; evidence),
so Proposition B is special to `p = 5` **as a one-prime statement**. The mechanism behind it is not special: the
windows of a *complete* set of smaller primes tile the circle, exactly as the two `n = 1` windows do at 5.

> **Theorem C (initial-segment rigidity).** Let `p` be a split prime and `S` a set of split primes containing
> **every** split prime `< p`. Let `phi in R^(S ∪ {p})` satisfy `phi_s = arg pi'_s (mod pi/2)` for each `s in S`, where
> `pi'_s` is one of the two Gaussian primes above `s`, and the gap inequalities of Lemma 5.1 for every
> `v = (v', ±1)` with `v' in Z^S` and `w(v') < log p`. Then `phi_p = ± arg pi_p (mod pi/2)`.
>
> **Corollary C.** If every split prime `<= y` divides `N`, then in every angle model the angles of all these primes
> are genuine: `phi_s = ± arg pi_s (mod pi/2)` for each split `s <= y`. In particular, if an angle model has a non-genuine angle at a
> prime `p | N`, then some split prime `< p` does not divide `N`.

*Proof of Theorem C.* Write `theta = phi_p`. For `v' in Z^S` let `A'_(v') = prod_s pi'_s^(v'_s+) conj(pi'_s)^(v'_s-)`, so
`<v', phi_S> = arg A'_(v') (mod pi/2)` and `Norm A'_(v') = e^(w(v'))`. As `v'` runs over `Z^S`, the `A'_(v')` run, up to
units, over all primitive Gaussian integers of odd norm whose prime factors lie above primes of `S`; in particular
over **all** primitive Gaussian integers of odd norm `< p`.

*Step 1 (gaps give separation from every short lattice direction).* For every primitive `z in Z[i]` with
`Norm z < p`: `|Im(z e^(-i theta))| >= p^(-1/2)`, i.e. `|sin(theta - arg z)| >= (p Norm z)^(-1/2)`. Indeed a primitive `z`
is `u (1+i)^eps z_0` with `u` a unit, `eps in {0, 1}` and `z_0` primitive of odd norm, and `conj(z_0) = u' A'_(v')` for
some `v'` with `w(v') = log Norm z_0 < log p`. For `v = (v', 1)`, `<v, phi> = theta - arg z_0 (mod pi/2)` and
`e^(-w(v)/2) = (p Norm z_0)^(-1/2)`. If `eps = 0` the axis gap gives `dist(theta - arg z_0, (pi/2)Z) >= arcsin((p Norm z)^(-1/2))`;
if `eps = 1` then `arg z = arg z_0 + pi/4 (mod pi/2)` and the diagonal gap gives
`dist(theta - arg z_0, pi/4 + (pi/2)Z) >= arcsin((2p Norm z_0)^(-1/2)) = arcsin((p Norm z)^(-1/2))`. In both cases
`|sin(theta - arg z)| = sin dist(theta - arg z, pi Z) >= sin dist(theta - arg z, (pi/2)Z) >= (p Norm z)^(-1/2)`.
(`z = 1` and `z = 1 + i` are the cases `v' = 0`.)

*Step 2 (Minkowski: `theta` is close to `arg pi_p` or `arg conj(pi_p)`).* Suppose
`dist(theta - arg beta, (pi/2)Z) > arcsin(1/p)` for both `beta = pi_p, conj(pi_p)`. For `0 < eps < 1/(3p)` the closed
convex symmetric set `K = {z : |Re(z e^(-i theta))| <= (1+eps) sqrt p, |Im(z e^(-i theta))| <= (1 - eps/2)/sqrt p}` has
area `4(1+eps)(1-eps/2) > 4`, so it contains a non-zero `z in Z[i]`; dividing by the gcd of its coordinates keeps it
in `K`, so `z` may be taken primitive. Then `Norm z <= (1+eps)^2 p + 1/p < p + 1`, so `Norm z <= p`, and
`|Im(z e^(-i theta))| < p^(-1/2)`. If `Norm z < p` this contradicts Step 1. If `Norm z = p`, `z` is an associate of `pi_p` or
`conj(pi_p)` and `|sin(theta - arg z)| < 1/p`, so `dist(theta - arg beta, (pi/2)Z) < arcsin(1/p)` for that `beta`:
contradiction. Hence `theta = arg pi~ + eta` with `pi~` an associate of `pi_p` or `conj(pi_p)` and `|eta| <= arcsin(1/p)`.

*Step 3 (two tight directions pin `eta = 0`).* Let `pi~ = x + iy`; `gcd(x, y) = 1`, so the lines
`L_± = {z in Z[i] : Im(z conj(pi~)) = ±1}` are non-empty and invariant under `z -> z + pi~`, which shifts
`R(z) = Re(z conj(pi~))` by `p`. Choose `z_± in L_±` with `R_± = R(z_±) in (0, p]`; `R_± = p` is impossible
(it would give `Norm z = (p^2 + 1)/p`), so `0 < R_± < p` and `Norm z_± = (R_±^2 + 1)/p < p`. Each `z_±` is primitive
(an integer combination of its coordinates equals `±1`). Since `Im(z e^(-i theta)) = (Im(z conj pi~) cos eta - R(z) sin eta)/sqrt p`,
Step 1 for `z_±` reads `|± cos eta - R_± sin eta| >= 1`. For `z_+`: `cos eta - R_+ sin eta` decreases from `1` at `eta = 0` to
`-1` at `eta = 2 arctan(1/R_+)`, so the inequality fails for `0 < eta < 2 arctan(1/R_+)`; for `z_-` it fails for
`-2 arctan(1/R_-) < eta < 0`. As `2 arctan(1/R_±) > 2 arctan(1/p) > arcsin(1/p)` (for `p >= 5`), Step 2 leaves only
`eta = 0`. QED

*Proof of Corollary C.* Induction over the split primes `5 = s_1 < s_2 < ... <= y`, using Theorem C with
`S = {s_1, ..., s_(k-1)}` at step `k`. QED

**Independent check (`cascade_cert.py`; exact up to the correctness of `mpmath` at 50 digits).** For the backgrounds
`{5}, {5,13}, ..., {5,...,41}` and the next primes `13, 17, 29, 37, 41, 53`, in **all** `1 + 2 + 4 + 8 + 16 + 32` conjugation
branches (with `phi_5 = arg(2+i)` fixed by the global symmetry `phi -> -phi`), the windows of all monomials with
`w(v) <= 20`, each shrunk by `1e-40`, cover `[0, pi/2)` except two intervals of length `2e-40` at `± arg pi_p`, and
exact Gaussian-integer arithmetic exhibits tight monomials on both sides of each (`cascade_cert_out.jsonl`).
For example at `phi_13 = arg(2+3i)`, `phi_5 = arg(2+i)`: `2+3i` is diagonal-tight (`|X - Y| = 1`) on one side and
`(2+3i)(2+i) = 1+8i` axis-tight (`|X| = 1`) on the other. Floating-point runs of `next_angle.py` agree:

| background (genuine) | new prime | `T` | admissible measure | admissible set |
|---|---|---|---|---|
| `{5}` | 13 | 34 | 0 | only `± arg(3+2i)` (lost to boundary rounding) |
| `{5}` | 17 | 34 | 0.0196 | 36 intervals; contains `arg(2+3i)` (13 is missing from the background) |
| `{5, 13}` | 17 | 30 | 0 | only the genuine directions |
| `{5, 13, 17}` | 29 | 26 | `3.5e-17` | one point-like interval at `arg(5+2i)` |
| `{17}` | 13 | 34 | 0.1216 | 151 intervals (5 is missing) |

The rows with an incomplete background (`{5} -> 17`, `{17} -> 13`) show that completeness matters: remove one
smaller prime and the admissible set has positive measure again.

**Relevance to the window.** Corollary C says that angle models are exactly as rigid as actual clusters on any
*complete* initial segment of split primes; the rigidity is a Minkowski (geometry of numbers) phenomenon, the
multi-prime form of `round10/analytic.md` Theorem D. In the window the cluster's primes are a **sparse** subset of the
split primes up to `P`: near-extremal clusters have `r = O(M)` primes (`round8/stability.md` Cor 3.6(3)), a fraction
`O(M log P/P) = O(M^(-1+delta) log M) -> 0` of the split primes `<= P`, so Step 1 is available only for the few `z` composed of cluster primes and the
Minkowski point of Step 2 is typically not among them. That is exactly where §5 finds positive measure again.

### 3.3 Lemma T: the arithmetic points of `G(p)` are actual clusters (proved)

> **Lemma T (transplants).** Let `beta_1, ..., beta_r` be primitive Gaussian integers (not divisible by any rational
> integer `> 1`) of odd norms `n_j = Norm beta_j` with `5 <= n_j <= p_j` and `gcd(n_i, n_j) = 1` for `i != j`. Then
> `phi_j := arg beta_j` is gap-admissible for the primes `p_1, ..., p_r`. It realises a profile `(p, e, (a_x))` with
> constant `C` iff the `M` Gaussian integers
> `Z_x = prod_j beta_j^(a_xj) conj(beta_j)^(e_j - a_xj)` lie on an arc of angular width `C e^(-W/4)` (up to units);
> in that case they are `M` distinct lattice points of the circle of norm `N' = prod_j n_j^(e_j) <= N`, on an arc of
> length `<= C (N'/N)^(1/4) sqrt(R') <= C sqrt(R')`, `R' = sqrt(N')`.

*Proof.* For `v in Z^r` put `A_v = prod_j beta_j^(v_j^+) conj(beta_j)^(v_j^-)`, so `arg A_v = <v, phi>` and
`Norm A_v = prod_j n_j^(|v_j|) <= e^(w(v))`. *`A_v` is primitive of odd norm.* A rational prime `l` dividing
`A_v` divides `Norm A_v`, hence exactly one `n_j`. Since `beta_j` is primitive, `l = 1 mod 4` splits as
`l = pi conj(pi)` and `beta_j` is divisible by only one of `pi, conj(pi)`, say `pi`. The `l`-part of `A_v` is the
`l`-part of `beta_j^(v_j)` (if `v_j > 0`) or of `conj(beta_j)^(-v_j)` (if `v_j < 0`), a power of a single one of
`pi, conj(pi)`; so `l` does not divide `A_v`. `Norm A_v` is odd since all `n_j` are odd. *Gap.* Write
`A_v = X + iY`. If `Y = 0` then `A_v = X` is a primitive integer, so `X = ±1` and `Norm A_v = 1`, forcing `v = 0`
(all `n_j >= 5`); the same for `X = 0`. If `X = ± Y` then `(1 ∓ i) | A_v`, so `Norm A_v` is even: impossible. Hence
`|X|, |Y|, |X ± Y| >= 1` for `v != 0`, and
`sin dist(<v,phi>, (pi/2)Z) = min(|X|,|Y|)/|A_v| >= Norm(A_v)^(-1/2) >= e^(-w(v)/2)`; the diagonal inequality is
the same with `|X ± Y|/(sqrt 2 |A_v|)`. *Realisation.* `<u_x, phi> = arg Z_x` (with the lift fixed by `phi`), so (5.1)
says exactly that the `Z_x` lie in a sector of width `Delta = C e^(-W/4)` modulo units. *Distinctness.*
`Z_x conj(Z_y) = A_(a_x - a_y) * (positive integer)`; if `Z_x = Z_y` this is real, so `a_x = a_y`. *Arc length.* The
`Z_x` have modulus `sqrt(N') = R'`, so their chord lengths are `<= R' Delta = C R' N^(-1/4) = C (N'/N)^(1/4) sqrt(R')`.
QED

> **Corollary T.** If, for some fixed `c` and `C`, transplant angle models with `W <= c M log M` exist for
> arbitrarily large `M`, then actual circles carry `M` lattice points on arcs of length `C sqrt(R')` with
> `log R' <= (c/2) M log M`; hence `M_infinity(C) = infinity` and the Cilleruelo–Granville conjecture fails, at the
> rate `log R/loglog R`. Conversely the conjecture implies that no transplant model is a barrier-scale fake.

So the arithmetic part of `G(p)` cannot supply fakes that cap anything: a transplant model at barrier scale is
an actual near-extremal cluster. (The genuine angles are the case `beta_j = pi_j`.) **A fake in the window must
therefore be a non-arithmetic point (Corollary T) that no positive-probability Haar recipe finds (Theorem A when
`m` is large): a point of a doubly-exponentially thin positive-measure part of `G(p) ∩ coset`, if that part is
non-empty (§5.4, §5.5, §6.3).**

## 4. Large sieve in the window: implied by the gap axioms

No fake was constructed, so the requested check ("do the fakes satisfy the large sieve statements available in
the window?") is answered in the only form that makes sense here: **every angle model, whatever its prime sizes,
satisfies the large sieve statements in question, because they follow from the pair gaps.** So the large sieve
cannot be the input that excludes window fakes; an argument using it is an argument using class `A_0`.

**Setting.** The Hecke characters of `Q(i)` are `lambda^k(alpha) = (alpha/|alpha|)^(4k)`; the angle of a split prime
ideal `(pi)` is `4 arg pi in R/2piZ` (well defined: units change `arg pi` by multiples of `pi/2`), and its conjugate has
angle `-4 arg pi`. An angle model replaces `± 4 arg pi_j` by `± 4 phi_j` for the cluster primes `p_j`. The large
sieve used below is the Montgomery–Vaughan/Selberg form: if `x_1, ..., x_R in R/Z` are `delta`-separated, then for any
`N` consecutive integers `n` and any coefficients, `sum_n |sum_r c_r e(n x_r)|^2 <= (N + delta^(-1) - 1) sum_r |c_r|^2`.

> **Lemma LS1 (spacing).** Let `phi` be gap-admissible for primes `p_1, ..., p_r <= X`. Then the `2r` points
> `± 4 phi_j` are `(4/X)`-separated in `R/2piZ`. The same holds for the genuine angles `± 4 arg varpi` of all split
> Gaussian prime ideals of norm `<= X`.

*Proof.* `dist(4(phi_i ∓ phi_j), 2piZ) = 4 dist(phi_i ∓ phi_j, (pi/2)Z) >= 4 arcsin((p_i p_j)^(-1/2)) >= 4/X`, by the gap
for `v = e_i ∓ e_j` (`w(v) = log p_i p_j <= 2 log X`). For the pair `± 4 phi_j` of one prime:
`dist(8 phi_j, 2piZ) = 4 dist(2 phi_j, (pi/2)Z) >= 4 arcsin(1/p_j) >= 4/X` (`v = 2e_j`). The genuine statement is the same
computation with `phi_j = arg pi_j`, valid for all split primes (it is Lemma 5.1 for the genuine angles). QED

> **Theorem LS.** Let `phi` be gap-admissible for primes `p_1, ..., p_r <= X`, and let `Theta^ext` be the genuine
> angles of all split prime ideals of norm `<= X` other than the `2r` cluster ideals. Put `D(K, X) = 2K + pi X/2`.
> 1. *(cluster angles)* For all `K >= 0` and coefficients,
>    `sum_{|k|<=K} |sum_j (a_j e^(4ik phi_j) + b_j e^(-4ik phi_j))|^2 <= D(K,X) sum_j (|a_j|^2 + |b_j|^2)`:
>    the same inequality, with the same constant, as for genuine Gaussian primes.
> 2. *(all prime angles, cluster angles replaced)* For the multiset `Theta^ext ∪ {± 4 phi_j}` the large sieve holds
>    with constant `2 D(K, X)`.
> 3. Hence every input of the form "the large sieve over Gaussian primes of norm `<= X` holds with constant
>    `c D(K, X)`" with `c >= 2` (union) or `c >= 1` (cluster primes only) is satisfied by every angle model with
>    all primes `<= X`. No condition on `M`, `P` or the regime is needed.

*Proof.* 1: Lemma LS1 and the large sieve with `x = theta/2pi`, `delta = 2/(pi X)`, `N = 2K + 1`. 2: both families are
`(4/X)`-separated (Lemma LS1); write the exponential sum as `S_ext + S_cl` and use
`|S_ext + S_cl|^2 <= 2|S_ext|^2 + 2|S_cl|^2` and 1 for each family. 3: immediate. QED

**What this settles and what it leaves.**

* In the window (`r > M/4 > P^(1/(2-delta))/4`), the cluster angles are a positive-power fraction of all prime
  angles up to `P`, and the large sieve at `K ~ P^(1/2)` is non-trivial for them (`r > P^(1/2)` makes
  `D(K,P) < (2K+1) r`). Theorem LS shows that this non-triviality is irrelevant: the inequality is a consequence of
  the pair gaps `v = e_i ± e_j` and `v = 2e_j`, so it holds for genuine and fake angles alike. A large-sieve
  argument in the window is an argument from class `A_0` (overview §5.4 row 1) and inherits that class's status
  there. To force anything, the large sieve would have to force angular concentration of the cluster primes,
  and the realisation condition (5.1) constrains only the sums `<u_x, phi>`, not the individual `phi_j`
  (Theorem A(b): on a realisation coset the free coordinates are exactly uniform).
* Round 10 Thm B needed grafted primes `>= M^A`, `A > 2`, to get the factor `1 + o(1)`. Theorem LS removes every
  size condition at the price of the factor 2. The range `1 <= c < 2` for the union is the "sharp `(X + K)` large
  sieve" of overview §5.5, equivalent up to constants to the max-scale cross-gap
  `|sin(phi_j - arg varpi)| >= c'/max(p_j, Norm varpi)` against **external** primes; it is not implied by the gap
  axioms and stays unclassified. It is satisfied by the transplant points of Lemma T
  (`|sin(arg beta_j - arg varpi)| = |Im(beta_j conj(varpi))|/(n_j Norm varpi)^(1/2) >= 1/X` unless `beta_j` is
  associate to `varpi` or its conjugate), i.e. by the only gap-admissible points known to exist in the window.

## 5. Computations (scripts in `fakes_checks/`; all evidence unless stated)

All scripts run under `prlimit --as=3500000000`; outputs are streamed to `*.jsonl` next to them.

### 5.1 The Haar measure of `G(p)`: a Poisson law (`mc_gapset.py`, `mc1_out.jsonl`, `mc2_out.jsonl`)

Uniform `phi in T^r`, all monomials with `w(v) <= T` tested (both gap conditions), `n` samples. Truncation makes
the surviving fraction an **upper** estimate of `mu(G(p))`. The genuine vector was checked to survive in every run.

| primes | `Lambda` | `Pi` | `T` | `n` | surviving fraction | `-ln(frac)/(Pi - 1)` | `prod_j mu_1(G_1(p_j))` |
|---|---|---|---|---|---|---|---|
| 13 | 0.277 | 1.77 | 30 | 1e5 | 0.302 | 1.55 | 0.302 |
| 13, 17 | 0.520 | 2.90 | 26 | 1e5 | 0.0361 | 1.75 | 0.119 |
| 13, 17, 29 | 0.706 | 4.22 | 24 | 1e5 | 0.00315 | 1.79 | ~0.06 |
| 13, 17, 29, 37 | 0.870 | 5.88 | 22 | 1e5 | 0.00024 | 1.71 | ~0.035 |
| 13, ..., 41 | 1.026 | 8.06 | 21 | 1e5 | 0 | > 1.6 | ~0.02 |
| 6 primes 97–149 | 0.558 | 3.06 | 22 | 1e5 | 0.0805 | 1.22 | |
| 9 primes 97–181 | 0.788 | 4.86 | 20 | 1e5 | 0.0102 | 1.19 | |
| 8 primes 193–277 | 0.523 | 2.85 | 21 | 1e5 | 0.123 | 1.13 | |
| 12 primes 193–317 | 0.753 | 4.52 | 19 | 6e4 | 0.0213 | 1.09 | |
| 5 | 0.447 | 2.62 | 40 | 2e5 | 0 | — | 0 (Prop. B) |

(The product column uses `mu_1(G_1(13)) = 0.3019`, `mu_1(G_1(17)) = 0.3927` from `gap1d.py` and rough values for
the larger primes.) **Reading.** `mu(G(p))` behaves like `exp(-c (Pi - 1))` with `c` between 1.1 and 1.8, i.e. like
the Poisson prediction `exp(-E)` with `E` the expected number of failing monomials, and far below the single-prime
product that Theorem A proves. If this law persists, (1.1) gives `mu(G(p)) <= exp(-c exp(M^(delta/2)/4))` in the
window, which is smaller than `e^(-K W)` for every fixed `K` once `M^(delta/2) >> log(M log M)`. In words: **in the
window the gap-admissible set is thinner than any resolution that a barrier-scale construction can address**,
so a fake there cannot come from a positive-probability random choice; it has to be found by structure (a
sequential search, §5.5) or be arithmetic. Thinness is not emptiness: the law predicts positive measure.

### 5.2 One more angle on a fixed background (`next_angle.py`, `next_angle_out.jsonl`)

Admissible set of `phi_p` (`p = 101`) when the background angles are fixed, all monomials with `w(v) <= T`.

| background | type | `Pi_bg` | `T` | admissible measure | Gaussian directions of norm `<= 101` inside |
|---|---|---|---|---|---|
| none | — | 1 | 24 | 0.773 | all 48 |
| 13 | genuine | 1.77 | 24 | 0.629 | 44 (those not parallel to `arg pi_13`-monomials) |
| 13, 17 | genuine | 2.90 | 24 | 0.442 | 40 |
| 13, 17, 29 | genuine | 4.22 | 24 | 0.305 | 36 |
| 13, ..., 41 | genuine | 8.06 | 22 | 0.075 | 26 |
| 13, ..., 73 (8 primes) | genuine | 17.4 | 20 | 0.0072 | 20 |
| 13, ..., 73 (8 primes) | Haar | 17.4 | 20 | 0.0255 | 2 |

`ln(measure/0.773) ≈ -0.29 (Pi_bg - 1)` throughout (slope 0.27–0.33): the windows cover like independent intervals,
with no overlap saving, the union-bound sum being `≈ 0.22 Pi_bg`. With a genuine background the Gaussian
directions of odd norm `<= p` coprime to the background survive, as Lemma T predicts; with a Haar background they
have no special status. **Consequence for hybrid constructions** (genuine angles on a supercritical set of primes,
Haar angles on the rest): a free prime `p` survives with probability `~ exp(-c Pi_bg/sqrt p)`; in the window
`Pi_bg >= exp(M^(delta/2)/4)` while `sqrt p <= M`, so hybrids fail as well (evidence; the single-coordinate version
is not covered by Theorem A, whose coordinates must be jointly Haar).

### 5.3 Exact one-prime sets and the cascade (`gap1d.py`, `probe5.py`, `next_angle.py`)

* `G_1(5)`: interval subtraction to `n = 60` at 60 digits leaves one interval of length `4e-62` at `arctan(1/2)`
  (the other point is lost to boundary rounding); `probe5.py` shows that `t = arctan(1/2) ± eta` violates `n = 1`
  for every tested `eta in {1e-3, ..., 1e-40}`. This is the computational shadow of Proposition B.
* `G_1(13)`: measure 0.3019, 684 components (`n <= 60`); `G_1(17)`: 0.3927, 916 components.
* Cascade: table in §3.2.

### 5.4 Are the Haar survivors arithmetic? (`mc_survivors.py`, `refine_surv.py`, `*_out.jsonl`)

Haar survivors of the truncated test were saved, re-tested at larger truncations, and compared with transplant
vectors: for each coordinate, the distance to the nearest direction `arg beta` (`beta` primitive, odd norm
`<= p_j`, norm `> 1`) in units of the mean spacing of those directions (uniform points give median `~0.28`).

| primes | survivors at `T` | still surviving at larger `T` | median normalised distance | all coordinates within 0.1 |
|---|---|---|---|---|
| 13, 17, 29 | 644 of 2e5 (`T = 24`) | 644 at `T = 28, 32, 36` | 0.120 (uniform 0.278) | 61/644 = 9.5% (uniform 1.2%) |
| 13, 17, 29, 37 | 325 of 2e6 (`T = 22`) | 321 at `T = 26`, 319 at `T = 30` | 0.124 (uniform 0.276) | 6/319 = 1.9% (uniform 0.2%) |
| 13, …, 41 | 21 of 6e6 (`T = 21`) | not refined | 0.115 (uniform 0.276) | 0/21 |
| 9 primes 97–181 | 2934 of 3e5 (`T = 20`) | not refined | 0.211 (uniform 0.275) | 0/2934 (32–56 directions per prime) |

**Reading.** (1) The survivors are stable under refinement, so at `Lambda = 0.71` and `0.87` the set `G(p)` has a
genuine positive-measure part (`≈ 3e-3`, `≈ 1.6e-4`). (2) That part is enriched near transplant vectors (median
distance halved, joint closeness 8–9 times the uniform rate) but is **not** concentrated on them: most survivors
are far from every transplant vector in at least one coordinate. With larger primes (many more
transplant directions) the bias is weak (median `0.21` against `0.275`). So a strong rigidity statement ("every point of
`G(p)` is near an arithmetic point") is **false** at moderate `Lambda`; what the data show is a positive-measure,
Poisson-thin, arithmetically biased set.

### 5.5 Sequential construction beyond the union-bound regime (`seq_build.py`, `seq_lookahead.py`, `verify_deep.py`)

Coordinates are chosen one at a time, each uniformly from its exact admissible set given the previous ones (all
monomials with last non-zero coordinate `k` and `w <= T`; `next_angle.allowed_set`); the product of the step measures
is an unbiased estimate of `mu(G_T)`. Completed vectors are re-verified against **all** monomials up to a larger weight.

| primes | `Lambda` | `Pi` | `T` | complete / tries | verified | log of path measure |
|---|---|---|---|---|---|---|
| 13, …, 41 (5) | 1.03 | 8.06 | 21 | 4/4 | 3 at `w <= 27` | −12.6 … −14.2 |
| 13, …, 73 (8) | 1.41 | 17.4 | 20 | 4/4 | **4 at `w <= 34` (360,050 monomials)** | −27.9 … −31.4 |
| 13, …, 97 (10) | 1.62 | 26.4 | 18 / 23 / 27 | 3/3 each | 0 at 22 / 0 at 27 / **1 at `w <= 33` (1,043,472 monomials)** | −38 … −47 |
| 13, …, 109 (12) | 1.81 | 39.0 | 17 / 22 / 26 | 3/3, 3/3, **0/3** | 1 at 20 / 0 at 25 / — | −52 … −65 |
| same, one-step lookahead (6 candidates) | 1.81 | 39.0 | 26 | 1/2 | 0 at 29 | — |

The Poisson-expected number of failures of a "random" point of `G_T` beyond the verified weight is small: the
shell sums `2.2 sum e^(-w/2)` decay geometrically (8 primes: `0.885, 0.19, 0.084, 0.035` on `(20,25], (25,28], (28,31],
(31,34]`, so `≈ 0.03` beyond 34; 10 primes: `1.11, 0.46, 0.11` on `(23,27], (27,31], (31,33]`, so `≈ 0.2` beyond 33;
`verify_deep.py`). So these are, with good confidence (not proof: infinitely many monomials remain), points of `G(p)`: gap-admissible vectors with `Pi = 17.4` and
`26.4`, where `outside.md` A.1(ii) needs `Pi < 1.75` and Haar sampling would need `~e^28` resp. `~e^45` trials. They
are not arithmetic: each has a coordinate at distance `>= 0.25` mean spacings from every transplant direction (maximal
normalised distance `0.25–0.63`; medians `0.08–0.17`, against `0.28` for uniform points). At `Pi = 39` the construction dies: with `T = 26` the admissible set of the 9th–12th coordinate is
empty for every background tried. This is what the Poisson picture predicts once the expected number of surviving
gaps per step drops below 1 (step mass `kappa_1 Pi_{<k}/(sqrt(p_k) - 1) ≈ 7.5` at `p = 109`), so it is not by itself
evidence of tiling; but it marks the end of the experimentally accessible range. **Points of `G(p)` that are
neither Haar-findable nor arithmetic exist at least up to `Pi ≈ 26`.**

## 6. What is proved, what is evidence, what is open

### 6.1 Ledger

| statement | status |
|---|---|
| (1.1): capped prime counts + window give `Lambda > M^(delta/2)/4`, `Pi > exp(M^(delta/2)/4)` | proved (from `round8/stability.md` Thm 3.3) |
| The equal-weight window `[P, P e^(1/(4r))]` is a design feature, removable by grouping | heuristic remark (disjoint blocks not proved) |
| Theorem A: Haar proportion of gap-admissible points on any realisation coset `<= exp(-2.17 sum_J p_j^(-1/2))` | proved |
| Corollary A.1: Haar completion with success `>= 1/2` needs `P_+ > 9.8 m^2` (double core: `9.8 (M-1)^2`) | proved |
| Corollary A.2: in the window with `m >= cM`, Haar completion succeeds with probability `<= exp(-2.17 c M^(delta/2))` | proved |
| Proposition B: `G_1(5) = {arctan(1/2), arctan 2}`; the angle of 5 is genuine in every angle model | proved |
| Theorem C / Corollary C: complete initial segments of split primes are rigid (angles genuine) | proved |
| Cascade certificates `13, ..., 53`, all conjugation branches (independent check of Theorem C) | exact (mpmath, 50 digits) |
| Lemma T: transplant vectors are gap-admissible; transplant models are actual clusters on `N' <= N` | proved |
| Corollary T: barrier-scale transplant fakes for unbounded `M` would disprove CG | proved |
| Theorem LS: large sieve over cluster primes (constant `2K + pi X/2`) and over all primes (constant twice that) holds for every angle model | proved |
| Poisson law `mu(G(p)) ≈ exp(-c(Pi - 1))`, `c ∈ [1.1, 1.8]` | evidence (§5.1) |
| Next-angle measure `≈ exp(-0.29(Pi_bg - 1))` at `p = 101`; hybrids fail | evidence (§5.2) |
| Proposition A': `mu(G(p)) > 0` whenever `kappa_1 Pi_{<k} < sqrt(p_k) - 1` (e.g. `Lambda <= (1/4 - o(1)) log P`) | proved (full torus only) |
| Haar survivors form a positive-measure set, biased towards but not concentrated on transplant vectors | evidence (§5.4) |
| Sequential construction finds non-arithmetic gap-admissible vectors at `Pi = 17.4` (verified to `w <= 34`) and `Pi = 26.4` (to `w <= 33`), Haar measure `~e^(-30)`, `~e^(-45)`; dies at `Pi = 39` | evidence (§5.5) |

### 6.2 What was not achieved, and why

**No fake was constructed in the window, and nonexistence was not proved.** The two known sources of
gap-admissible points are:

* *Haar-generic points.* They exist in profusion when `Pi - 1 < 0.75` (outside.md A.1(ii)), with positive measure on the full torus up to
  `Lambda ~ (1/4) log P` (Proposition A'), and Theorem A shows
  that their proportion on a realisation coset is `<= exp(-2.17 sum_J p_j^(-1/2))`. In the window, for profiles with
  free dimension `m >= cM`, this is `exp(-c' M^(delta/2))`: the method cannot work, and the evidence of §5.1 says
  the true proportion is doubly exponentially small.
* *Arithmetic points* (Lemma T, including the genuine angles). They exist for every prime set, but a transplant
  model at barrier scale **is** an actual cluster (Corollary T), so it cannot cap an input without disproving the
  conjecture.

A window fake would therefore have to come from one of two places, neither settled:

1. *The Poisson-thin positive-measure part of `G(p) ∩ coset`.* Haar completion cannot find it (Theorem A), but
   it need not be empty: at moderate `Lambda` it is non-empty, stable under refinement, and mostly far from
   arithmetic points (§5.4), and a sequential construction finds its points where Haar sampling would need
   `~e^30` trials (§5.5). If it stays non-empty on realisation cosets as `Lambda -> infinity`, window fakes exist
   (non-constructively from the point of view of union bounds). Proposition A' proves non-emptiness on the full
   torus only up to `Lambda ~ (1/4) log P`, far below the window.
2. *Rigid realisations* (`m = 0` or small): when `rank L = r`, the angles are determined by `(delta, k)` modulo a
   finite group and are, up to `O(||lambda|| Delta)`, rational multiples of `pi/2` with denominator `D = [Z^r : L]`.
   The gap conditions then become conditions on a character `chi` of the finite group `Z^r/L`: its kernel must have
   minimum weight `>= W/2 - O(log M)`, and `||chi(v)|| >~ e^(-w(v)/2)` off the kernel. A pigeonhole count does not
   obviously exclude this (a contradiction needs `#{v : w(v) <= W/4 - O(log M)} > D`, and by a rough Hadamard count `D`
   can be as large as `e^((r/2) log(4M))`),
   and a random character fails for the same union-bound reason as a Haar vector. Whether a structured character
   exists is a finite problem for each `M` that I did not attack. Note that the prime-count input gives only
   `r > M/4`, so profiles with `r <= M - 1` (where rigid realisations live) are not excluded by the capped classes.

### 6.3 Two scenarios for the window

Everything above reduces the window question to one statement about the gap-admissible set:

> **(Q)** For barrier-scale profiles with all primes `<= M^(2-delta)` (so `Lambda > M^(delta/2)/4`), and valid
> character margins, is `G(p) ∩ Real(delta)` non-empty for some `delta`, and if so, is it non-empty only at
> arithmetic (transplant-like) points?

* **Scenario P (Poisson positivity).** `nu(G ∩ coset) ≈ exp(-c(Pi_nc - 1)) > 0`, `Pi_nc` the non-character part of
  `Pi`, as the Poisson law of §5.1 extrapolates. Then window fakes exist, at barrier scale if the weighted character
  margins can be arranged (§2 suggests they can), and the window is **capped** for every input class that
  Haar-almost-every point satisfies; the `r^2` law would be an artefact of union bounds. Nothing would be
  constructible in the sense of the record (no positive-probability recipe), but capping needs only existence.
  *For:* positivity and stability at moderate `Lambda` (§5.4), sequential constructions far beyond the union-bound
  regime (§5.5), the near-independence of the windows (§5.2, the Poisson slope).
* **Scenario R (supercritical rigidity).** For `Lambda -> infinity` the windows tile the realisation cosets except
  at arithmetic points, as they provably do on complete initial segments of primes (Proposition B,
  Theorem C). Then window angle models are transplant-like, hence actual clusters (Corollary T): "window fakes" would disprove the conjecture, and "no window fakes" would be a restricted
  growth theorem for actual clusters with `P_+(N) < M^(2-delta)`. *For:* exact tiling is a theorem whenever the
  background is a complete initial segment (Theorem C); the enrichment of survivors near transplant vectors (§5.4).

The data do not decide between P and R: exact tiling is seen (and proved) only for complete initial segments,
which the window's sparse prime sets are not, and the Poisson behaviour is measured only up to `Pi ≈ 40`,
`Lambda ≈ 1.8`, where the sequential construction already starts to die (§5.5). Both scenarios agree on three things, which are the
firm conclusions of this round: (i) no construction method of the record works in the window (Theorem A);
(ii) the large sieve does not distinguish fakes from actual clusters there (Theorem LS); (iii) arithmetic fakes are
actual clusters (Lemma T). **So an argument proving growth in the window must decide (Q) itself, i.e. use the
supercritical gap structure of the cluster's Gaussian primes (a sparse-prime version of Theorem C, which is the
multi-prime Theorem D of `round10/analytic.md` for complete segments), not measure and not the large sieve.**

### 6.4 Open points, smallest first

1. A sparse form of Theorem C: for which prime sets `S` (density in the split primes `< p`, or a covering condition on
   the Minkowski points of Step 2) does a genuine background on `S` still pin `phi_p`? Theorem C needs density 1; §5
   shows that density `<= 1/2` already fails for small primes.
2. Extend Proposition A' to realisation cosets (sequential integration over coordinates of `H_0`) and find the
   largest `Lambda` for which it proves `nu(G ∩ coset) > 0`. Reaching `Lambda -> infinity` would prove scenario P,
   i.e. cap the window.
3. Prove a doubly exponential upper bound `nu(G ∩ coset) <= exp(-c(Pi - 1))` (the upper half of the Poisson law).
   Higher moments of the failure count, with linearly dependent tuples counted, look like the natural route. This
   would not decide (Q) but would quantify how thin window fakes must be.
4. The simplest test of scenario R without the prime 5: does some set of split primes `>= 13` have `mu(G(p)) = 0`?
   (All sets tested have positive measure; a single-coordinate version: for a genuine background `S` with
   `Pi_S >= C p^(1/2) log p`, is the admissible set of `phi_p` null?)
5. Decide the rigid-realisation problem (`m = 0`) for small `M` and primes `<= M^2`.
