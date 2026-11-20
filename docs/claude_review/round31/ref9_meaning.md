# Round 31, Part IX: analyst referee ("meaning")

Scope. This block does not re-referee the proofs of Part IX (other referees do). It reads Part IX in full
(`astra_free_IX/answer.md`, line 1955 to the end), her summary `final9.txt`, her work notes and the coordinator's
question `reply9.txt`, and asks what Part IX changes in the record if it is correct. Everything marked
"conditional on Part IX" assumes the proofs of IX.3, IX.4, IX.5 and IX.7 pass their own referees. Two short proofs
were re-derived here because the report depends on their exact form: (IX.1) and the orientation reading of (IX.2).

Labels. **Proved**: checked line by line here. **Computed**: exact or reproducible computation with my own code
in `ref9_meaning/`. **Heuristic**: a count or model that no proof uses. **Gap**: a missing step. **False**: refuted.
**Conditional on Part IX**: a consequence for the record that holds if the named Part IX statement is correct.
**Assessment**: a judgement (probabilities, recommendations).

Notation (as in the record). A *cluster* is a set of `M` distinct points of `Z[i]` on `|z|^2 = N`, `R = sqrt N`, on an
arc of length `<= C sqrt R`. The *normalised* cluster is `z_i' = z_i/g`, `g` the Gaussian gcd, `N' = N/|g|^2`,
`R' = sqrt N'`, `W = log N'`, `s = omega(N')`, `P_+ = P_+(N')`. `a_ip` is the exponent row. `A_0`, `A_1(Q_0)`, `A_{2M}`
are as in `round31/README_VIII.md` (pair/cut/gap/character inputs; plus round-5 axiom (R) at odd prime powers
`<= Q_0`). `F^char_X`, `F^all_X` (P3-char, P3-all) are the classes of `round5/round6_sharp/sharp.md` Definition 1.1,
quoted in section 1.1 below. The *low-height window* is `P_+ <= M^{2-delta}`.

## 0. Results in brief

1. **Conditional on IX.3-IX.4, the residue axioms that are tested only through characters are capped at every
   height at once.** One object satisfies Definition 1.1 (a)-(d), (e-char), including the same-level clause at all
   parent primes, for **every** cutoff `X`, at `W <= (22/3 + o(1)) m log m`, `P_+ <= m^{11/6+o(1)}`. Hence P3-char and
   the point-residue class `A_X` are capped uniformly in `A` (the caveat of round 9 R1 is removed for them), in
   regime (1.4), and inside the low-height window for every `delta < 1/6`. The record scope limit "all primes of `N`
   must exceed the modulus bound `X`" is removed for these classes. Not covered: 2-adic collisions (Definition 1.1
   has odd `q` only, the round-5 referee asked for `q = 2` in axiom (R)); Haar-majority genericity;
   Lang-Waldschmidt rows 11-12; the analytic grafts of rows 7 and 17 (unchanged, still for fixed `A`).
2. **Conditional on IX.5, the single-label max-scale cross-gap is capped at constant one, outside the window.**
   `|sin(phi_j - arg varpi)| >= 1/max(p_j, N varpi)` for all Gaussian primes `varpi`, `>= c_*/max(p_j, N w)` for all
   Gaussian integers `w`, and the referee's finite test `>= 1/p_j` for `N w <= p_j` (on an unbounded sequence),
   together with all-height P3-char, at `W = O_C(m log m)` with `P_+ = m^{d}`, `d` a large constant. Still uncapped:
   the cross-gap inside the window, the composite and product-scale cross-gaps of round 20 (items (a), (b) of
   overview §5.5), and the union large sieve at the genuine constant as an operator-norm statement (equivalent to
   the cross-gap only up to absolute constants; Astra did not evaluate it).
3. **Conditional on IX.4, the window is capped for `0 < delta < 1/6` at coefficient `22/3`** (limsup statement,
   diagonal sequence), for all inputs of `A_{2M}` and for all-height P3-char. Open: `delta >= 1/6`, i.e.
   `(2 - o(1)) M log M <= P_+ <= M^{11/6}`.
4. **The one new unconditional theorem is an orientation theorem** (IX.2, (IX.1); **Proved**, re-derived here): for
   every nonzero `b in Z[i]`, `M^{3/2} << 2 + |b| (C sqrt R |sin(theta - arg b)| + C^2/2)`, with an absolute implied
   constant, where `theta` is the arc's midpoint direction. It holds in original coordinates; normalisation only
   sharpens it. It is a theorem about the original problem in a family of arc positions of proportion
   `~ (6/pi^2) T^2/sqrt R` for threshold `T` (**Computed**), and it extends the record's near-rational-direction
   bound (round 4 `hardness.md` Lemma 4.2(b), `2|b|(kappa + C)^2 + 2`) from `kappa = o((log R/loglog R)^{1/2})` to
   `kappa = o((log R/loglog R)^{3/2})`. All of IX.2's certificates (sum-two row combinations, conjugate pairs, marks of
   small primitive length) are certificates of orientation. It gives nothing for generic positions: with the best
   Gaussian direction the bound is `M << R'^{1/6}` (section 2.5; **Computed**: the optimal `|b|(1 + sqrt R' |sin|)` is
   `(0.1-1.4) R'^{1/4}` for 98% of directions at every `R'` from `10^8` to `10^40`).
5. **The residue column of the map is not complete.** It is closed at its two ends: character-tested residues at
   any height are capped (item 1), and all-monomial residues (e-all) at all heights with one rational lift per label
   are the problem itself (IX.7). The middle is open and is exactly the record's P3-all\*: e-all at polynomial
   cutoffs `X = 3M^A`. IX.7 says nothing there, and the e-all failure of IX.3's witnesses is at an unquantified
   cutoff. "Rational lift" is a uniformity-in-the-cutoff condition only: at each finite cutoff every compatible local
   system is the reduction of some `r_j in Q(i)^*` with `r_j conj(r_j) = p_j` (weak approximation on the norm-one
   torus; **Proved**, section 3.2).
6. **Her closing question is a restatement, not a route to the rate.** A positive answer would make the all-height
   all-monomial class equal to the set of actual clusters, which gives no rate by itself. The route through residues
   to the rate is a lower bound in P3-all\* at a polynomial cutoff, the record's question since round 6.
7. **Recommendation.** Ask for P3-all\* in regime (1.4), using her rational-lift machinery quantitatively (section 6).
   The exchange is converging on a wall: after Part IX every residue-type input is capped, is the problem, or is
   P3-all\*; every analytic input is capped or is a square-root-quality or sharp-constant statement; the remaining
   candidates (Gaussian cocycle, additive data) are, in the record's own words, the problem itself or open in both
   directions since round 31 Part IV. Probability that the next pass changes the growth rate for general `R`: about
   2-3% (Assessment).


## 1. The updated map of inputs

### 1.1 The record's wording, quoted

* **Definition 1.1 of `round5/round6_sharp/sharp.md`** (factored residue data). For every odd prime `q <= X` and
  every `j` with `p_j != q`, a unit `r_j(q)` of `Z[i]/q^{A_q}` with `Norm r_j(q) = p_j mod q^{A_q}`; point residues
  `rho_x(q) = i^{-k_x} prod_j r_j(q)^{a_xj} conj(r_j(q))^{e_j - a_xj}` for `q` not dividing `N`; level residues, the
  same product over `j != l`, at `q = p_l`.
  (e-char): for every `c in Z^M` with `sum c = 0`, `v(c) != 0`, and every unit `eta`,
  `|sin((c.delta - arg eta)/2)| >= 2^{-1/2} m_(c,eta) e^{-w(v(c))/2}`, `m_(c,eta) = prod_{q <= X, q not| N} q^{a_(q,eta)(c)}`,
  plus "the same-level pair clause at the primes `p_l <= X`".
  (e-all): in addition, for every `v in Z^r \ {0}` and `s in Z/4`,
  `|sin(<v,phi> - s pi/4)| >= nu_s m_(v,s) e^{-w(v)/2}`, `a_(q,s)(v) = max{a <= A_q : prod_j (r_j/conj r_j)^{v_j} = i^s mod q^a}`.
  `F^all_X subset F^char_X subset A_X` ("Relations" paragraph there).
* **What the record proved and called uncapped.** P3-char: `W_min <= (33 max(A,2) + o(1)) M log M` at moduli
  `<= 3M^A`, **for each fixed `A`**, with every prime of `N` above the cutoff (sharp.md Thm 5.1; round-6 referee:
  "Lemma 5.1b and lower.md Lemma 5.2 need every deleted prime to exceed `X'`"). P3-all: only
  `W <= (66A + o(1)) M log^2 M/loglog M`; **P3-all\*** ("Is `W_min^all_{3M}(M,C) = O(M log M)`?") open. Overview §5.4
  caveat (round 9 R1): "Rows 2-4 and 7 hold for residue moduli up to `M^A` with `A` **fixed in advance** ... A growth
  argument may choose `A` after the cluster's size parameter `c_0`, and **caps uniform in `A` are not proved** ...
  **no fake is known in regime (1.4)**, where `p_j <= M^{h_0}`, `X = 3M^A` and `A` is chosen after `c_0`." Overview
  §5.5: "**The sharp `(X+K)` large sieve**: up to constants equivalent to the max-scale cross-gap
  `|sin(phi_p - arg varpi)| >= c/max(p, Norm varpi)` ... open and unclassified", and round 20's composite form
  "`|sin(<v,phi> - arg alpha)| >= c/max(Nd, N alpha)`" (item (a)) and product-scale form (item (b)). Expert note §1:
  "P3-all\* and the sharp-constant large sieves ... remain uncapped, and the residue caps are not uniform in `A`".
  `round31/README_VIII.md` §3(e): uncapped "in every window: prime-level residues (P3-char, P3-all\*); point residues
  at moduli `M^A`, `A > 1`; the union large sieve at the genuine constant (max-scale cross-gap) ...;
  square-root-quality (GRH) prime sums and box counts; Haar-majority and full-torus Haar-almost-everywhere
  properties; Lang-Waldschmidt rows 11-12; (VIII.6) ...; the association of each prime norm with its genuine
  Gaussian argument", and every input for `delta >= 1/12`.
* **Twisted point collisions** (round-5 referee, `round5/referee_lower.md` §1.3(b); overview §5.2 item 2). For
  `sum c = 2k != 0`, the input `prod rho^c = u N^k mod q^a` implies
  `|sin(<v,phi> - psi_u)| >= D' e^{-w(v)/2}/sqrt2` with `v = sum c_x a_x - k e`. Round-5 witnesses meet it at `X = 2M`;
  "the twisted extension at `X = M^A`, `A > 1`, is open". Axiom (R) itself (`construct.md` Def 1.2) is stated for
  zero-sum `c` only.

### 1.2 Astra's wording, quoted, and the scope differences

* IX.3: "growth-scale fakes with `m asymp omega(N')`, satisfying P3-char for **every** finite cutoff `X` with the
  same prime labels, marks, real angles, windings and compatible prime residues. This includes both original
  Liouville gaps on all integer monomials, the strengthened gaps on every rational row character and every Gaussian
  unit, and same-level pairs at all odd parent primes." The witnesses are reductions of fixed
  `r_j = i^{u_j} pi_j prod_l gamma_l^{V_jl} in Q(i)^*` with `r_j conj(r_j) = p_j` exactly (IX.6).
  *Scope relative to Definition 1.1.* (i) Stronger in the cutoff: one object for all `X` (the record's caps were one
  object per fixed `A`). (ii) Stronger in the parent primes: Definition 1.1 allows `p_l <= X`, and the record's
  constructions had all `p_j > X`; IX.3 has `p_j <= m^{11/6+o(1)}` far below `X` and handles the level residues by
  deleting the singular factor. (iii) Same in the unit set (all four units) and in the norm constraint (exact, even
  globally). (iv) Odd prime powers only, as in Definition 1.1. The 2-adic clause that the round-5 referee added to
  axiom (R) is not claimed. (v) "Rational row character" means zero-sum `c` (IX.3 text: "a rational zero-sum row
  combination"). Twisted collisions (`sum c != 0`) are not claimed at any height. (vi) The weight bound `22/3` and the
  window exponent `11/6` come from IX.4; the all-height residue part itself only needs the physical-prime slack (IX.7).
* IX.5: (IX.8), quoted in item 2 of §0. *Scope relative to the record.* The record's cross-gap allows any absolute
  `c`; IX.5 gives `c = 1` against Gaussian primes at every norm, `c_*` (absolute, not one) against all Gaussian
  integers, and `c = 1` against all Gaussian integers of norm `<= p_j` (the referee's finite test of
  `ref8_spread.md`), on an unbounded sequence. It tests **single labels** `phi_j`; round 20's composite form
  (`<v,phi>` against `arg alpha` at scale `max(Nd, N alpha)`) and the product-scale form are not addressed. The large
  sieve itself is not evaluated (IX.5: "We have not evaluated the optimal operator norm of a large sieve"). The
  primes are of size `q^d` with `d > max(2, 2/beta)` (or `d > 2` with twelve copies), so the object is **outside** the
  low-height window.
* IX.4: "for fixed `0 < xi < 1/6` and `.99 < rho < 1`, take `q` sufficiently large depending on both. Then
  `P_+ <= 2q^{2-xi}`, `W/(m log m) <= 4(2 - xi)/rho + o_{rho,xi}(1)`. Consequently (IX.3) caps every fixed height
  window `P_+ <= m^{2-delta}` for `0 < delta < 1/6`. The endpoint in (IX.3) is a limsup statement."

### 1.3 The map (conditional on Part IX where marked)

| input class | before Part IX | after Part IX (conditional on the named block) | what remains |
|---|---|---|---|
| `A_0` (pairs, cuts, monomial gaps, rational characters) | capped; in the window for `delta < 1/12` at `23/3` | capped in the window for `delta < 1/6` at `22/3` (IX.4) | window `delta >= 1/6`; Haar-majority genericity |
| point residues, axiom (R) / `A_X` (zero-sum characters, four units, odd prime powers, level pairs) | capped for `X <= M^A`, `A` fixed, primes `> X`; at `X = 2M` inside the window for `delta < 1/12` | **capped at every `X` simultaneously**, primes `<= m^{11/6+o(1)}`, uniform in `A`, in regime (1.4), in the window for `delta < 1/6` (IX.3-IX.4; via `F^char_X subset A_X`) | 2-adic clause; twisted collisions; window `delta >= 1/6` |
| **twisted point collisions** (`sum c != 0`) | met at `X = 2M` by round-5 witnesses; open at `X = M^A`, `A > 1` | **not addressed** (IX.3 treats zero-sum `c` only; IX.6 treats `n != 0` words with the plain gaps only) | uncapped at `X = M^A`; in factored form these are e-all instances on row-span words `v = (1/2) S^T c` |
| P3-char (Def 1.1 (d), (e-char)) | capped for fixed `A`, primes `> X`; uncapped in the window | **capped at every cutoff** (IX.3), in the window for `delta < 1/6` (IX.4) | window `delta >= 1/6`; 2-adic |
| P3-all at `X = 3M^A` (P3-all\*) | open; bracket `[3, (66A + o(1)) log M/loglog M]` | **unchanged.** IX.3's witnesses fail e-all at some cutoff that is not quantified; nothing is claimed at `X = 3M^A` | the last open local class (section 3) |
| e-all at **all** cutoffs, one rational lift per label (finite exceptional sets) | not a record class | equals actual clusters (IX.7): the problem itself | — |
| e-all at all cutoffs, no rational lift (independent local witnesses) | not a record class | open (her closing question); contains no more than P3-all\* as a lower-bound route, since `F^all_infinity subset F^all_{3M^A}` | section 4 |
| max-scale cross-gap, single labels, prime directions, constant one | uncapped | **capped outside the window** (IX.5, `P_+ = m^{d}`, `d` large) | inside the window: uncapped |
| max-scale cross-gap against all lattice directions | uncapped | capped with an absolute `c_*` at all norms, and with `c = 1` for norms `<= p_j` (IX.5) | `c = 1` at all norms against composite `w`; inside the window |
| composite and product-scale cross-gaps (overview §5.5 (a), (b)) | uncapped | **unchanged** | — |
| union large sieve at the genuine constant (full operator form) | uncapped | **unchanged** (equivalent to the cross-gap only up to absolute constants) | — |
| square-root-quality prime sums and box counts (GRH quality) | uncapped in the window | **unchanged** | — |
| Lang-Waldschmidt rows 11-12, full LW with `eps -> 0` | uncapped | **unchanged** (not checked on any Part IX object) | — |
| (VIII.6) height-order inequality | untested on fakes | unchanged | — |
| Gram identity (G3) with real cofactor magnitudes `B_ij` | not in the record | satisfied by IX.3 and IX.5 objects with `max log B_ij = O_C(log m)` (IX.6 remark). This is close to automatic: `B_ij` is *defined* from `d_ij` and the chords, so only the size bound has content | — |
| Gaussian-integer cocycle `g_ij b_ij + g_jk b_jk + g_ki b_ki = 0` with `b_ij in Z[i]` | not a record class; integral models without the norm condition open both ways (Part IV) | **unchanged** (IX.6: "The complex Gaussian cofactor cocycle in (G4) is an additional requirement") | section 4 |
| additive / archimedean data (`Log((z_x+z_y)/(2z_u)) = O(Delta)`) | uncapped; "used in full ... the problem itself" | unchanged; IX.5's Schottky measure and IX.7's irrational rotations are archimedean in the angles only | section 4 |
| integer projection (IX.1) | not in the record | **new genuine input; violated by `A_1(2m)` fakes in rational-direction windows** (IX.6) | gives a rate only in orientation windows (section 2) |

**Summary of item (1).** Conditional on IX.3-IX.4: prime-level residues through characters at **all** moduli,
yes; point residues at moduli `M^A` with `A` not fixed in advance, yes for zero-sum characters (not for twisted
collisions); the window for `delta < 1/6`, yes, at coefficient `22/3`. Conditional on IX.5: the max-scale cross-gap
at constant one, yes for single labels against prime directions, outside the window only. What remains: e-all
(P3-all\*) at polynomial cutoffs; twisted point collisions at `X = M^A`; the union large sieve at the genuine
constant in operator form, the composite cross-gaps and every cross-gap inside the window; square-root-quality
prime sums; Lang-Waldschmidt; the Gaussian-integer cocycle; additive (archimedean) data; the window for
`delta >= 1/6`.

## 2. IX.2: the restricted rate theorem, stated for the record

### 2.1 The projection lemma

**Lemma P (Proved; Astra's (IX.1), re-derived here).** Let `M >= 2` distinct lattice points lie on `x^2 + y^2 = N`
(no arc needed). For `0 != b in Z[i]` let `X_b` be the diameter of the integers `Re(conj(b) z_i)`. Then, with an
absolute implied constant,
`2 + X_b >> M^2 exp(-sum_{p <= M, p = 1 (4), p | N|b|^2} log p/p) >> M^{3/2}`.

*Proof (as checked).* `x_i = Re(conj(b) z_i)` satisfies `x_i^2 + y_i^2 = D := N|b|^2`, so each value is taken at most
twice and there are `n >= M/2` distinct values. For odd `p` not dividing `D` the first coordinates of the conic
`x^2 + y^2 = D` mod `p` occupy `nu_p <= (p+3)/2` classes; for inert `p | D`, one class; for split `p | D`, at most `p`.
By Cauchy-Schwarz at least `n^2/nu_p - n` ordered pairs are congruent mod `p`. Comparing with
`prod_{i != j} |x_i - x_j| <= X_b^{n(n-1)}` and summing `log p` over `p <= M` gives
`log(2 + X_b) >= sum_{p <= M} log p/nu_p - O(1)`, and Mertens (with the progression `1 mod 4` for the bad primes)
gives the two bounds. The loss is at most `(1/2) log M + O(1)`, attained when every split `p <= M` divides `N`.
**Computed** (`ref9_meaning/orient_check.py`): on all 93,924 windows with `M >= 3` and chord constant `<= 6` on circles
`N <= 2*10^5`, and all `b` with `|b|^2 <= 50`, `M^{3/2}/(2 + X_b) <= 3.78` and `M^2 e^{-loss}/(2 + X_b) <= 8.70`.

### 2.2 Orientation theorem

**Theorem O (Proved, given Lemma P).** Let `M` lattice points of `|z| = R` lie on an arc of length `<= C sqrt R`
with midpoint direction `theta`. For every `0 != b in Z[i]`,
`M^{3/2} <<  2 + |b| (C sqrt R |sin(theta - arg b)| + C^2/2)`,
with an absolute implied constant. Equivalently, if `|sin(theta - arg b)| <= kappa R^{-1/2}`, then
`M << (|b|(C kappa + C^2))^{2/3} + 1`, and `M << (|b|(C kappa + C^2))^{1/2} + 1` when no split prime `p <= M` divides
`N|b|^2`.

*Proof.* `x = |b| R cos(phi - arg b)` on the arc `|phi - theta| <= Delta/2`, `Delta = C R^{-1/2}`; the derivative is at
most `|b| R (|sin(theta - arg b)| + Delta/2)`, so `X_b <= |b| R Delta (|sin(theta - arg b)| + Delta/2)`. Apply Lemma P.
The same holds for the normalised cluster with `R'`, `theta' = theta - arg g`; that version is sharper by a power of
`|g|` but needs the gcd.

**Corollary O1 (rate in orientation windows; Proved).** If `|b|(1 + kappa) = o((log R/loglog R)^{3/2})` then
`M = o(log R/loglog R)`, uniformly in `N`. For fixed `|b|` and `kappa`, `M` is bounded. No input from Theorem 1.1 or 1.2
is needed: Theorem O is a direct bound.

**Certificates of orientation (Proved; IX.2 reread).** Every hypothesis of IX.2 is a certificate that the arc points
at a Gaussian direction of small height.
* *Sum-two row combinations.* For `c in Z^M` with `sum c = 2`, `prod z_i'^{c_i}/N' = u alpha/conj(alpha)` with
  `|alpha|^2 = T(c)`, and `b = eta alpha` has `arg b = (1/2) sum c_i arg z_i' mod pi`, within `||c||_1 Delta/4` of
  `theta'`. So (IX.2) is Theorem O with `|b| <= sqrt(2T(c))`, `kappa <= ||c||_1 C/4`.
* *Conjugate pairs.* `c = e_i + e_j`: `z_i'/conj(z_j') = u alpha/conj(alpha)`, `T = N'/|gcd(z_i', conj z_j')|^2`. In words:
  **if two points of the arc have a Gaussian gcd with each other's conjugate of norm `>= N'/T`, then `M << T^{1/3}`.**
  `T = 1` means `z_j' = u conj(z_i')`: the arc straddles an axis or a diagonal, the classical trivial case.
  With `T <= P_+^2` (rows complementary up to two unit changes): `M << P_+^{2/3}`; with (VIII.4),
  `P_+ >= (2 - o(1)) M log M`, this is the statement that a near-complementary pair forces `P_+ >> M^{3/2}`.
* *A mark of small primitive length.* `c = 2e_i`: `T = |z_i'|^2/gcd(Re z_i', Im z_i')^2`. In original coordinates:
  **if the arc contains a lattice point `k w` with `k in Z`, `w` primitive, then `M << (C^2 |w|)^{2/3} + 1`.**

**Computed** (`ref9_meaning/tenpoint.py`). The record's ten-point cluster (`N = 1176852625 = 5^3 13^2 17 29 113`,
chord constant 8.2386) is symmetric under `x <-> y`: its five reflection pairs `(z, i conj z)` are the five `T = 1`
pairs, its midpoint is exactly the diagonal direction, and the projections `x + y` are `48503, 48505, 48509, 48511,
48515` (diameter 12, so `M^{3/2}/(2+X) = 2.26`). The record's largest explicit cluster is therefore an orientation
cluster of height `sqrt 2`.

### 2.3 Comparison with what was known

* The record's round 4 `hardness.md` Lemma 4.2(b): "an arc of length `C sqrt R` whose midpoint direction is within
  `kappa R^{-1/2}` of `arg w` holds at most `2|w|(kappa + C)^2 + 2` points" (whole-cap count). The trivial
  improvement (count integers in the projected interval of the arc itself) gives `2 + 2|w|(C kappa + C^2/2)`, linear in
  `kappa`. Theorem O improves the exponent to `2/3` (and `1/2` without small split primes). The range where the bound
  is `o(log R/loglog R)` grows from `kappa = o((log R/loglog R)^{1/2})` to `kappa = o((log R/loglog R)^{3/2})`.
* CG07 (as quoted in overview §1.1): arcs touching an axis point are trivial. Chan (IJNT 2015, as summarised in
  `lit18/direct.md` §0, §2.3; I read only the record's summary): at most 10 points with `|b| < n^{1/4}(log n)^{1/7}`
  for `n = N^2`, i.e. on the cap of half-length `sqrt R (log R)^{1/7}` about an axis point when `R` is an integer.
  Chan is better than Theorem O inside his range (bounded count; Theorem O gives `<< (log R)^{4/21}` there via the
  cap version of Lemma P). Theorem O covers every Gaussian direction, every `N`, and `kappa` up to
  `(log R)^{3/2-o(1)}`, with a count that grows.
* Novelty against the literature: the method is Gallagher's larger sieve on a conic fibre; I did not search for this
  application and cannot certify novelty.

### 2.4 For which parents, and which families

* **No parent-only consequence.** Theorem O depends on the arc's position, not on `N'`. For a fixed parent it covers
  the arcs within `kappa R^{-1/2}` of the directions of height `<= T/(1+kappa)`; the proportion of directions `theta`
  with `min_b |b|(1 + sqrt R' |sin(theta - arg b)|) <= T` is `~ (6/pi^2) T^2/sqrt R'` (**Computed**: at `R' = 10^12`,
  `T = 51`, sampled proportion `0.00149` over `10^5` directions, first-order prediction `0.00158`). At
  `T = (log R)^{3/2}` this is a vanishing proportion, about `(log R)^3` arc positions out of `sqrt R`.
* **Families where it is a new theorem about the original problem** (Proved, restricted): (F1) arcs within angle
  `kappa R^{-1/2}` of a fixed Gaussian direction, `1 << kappa = o((log R/loglog R)^{3/2})`: `M << kappa^{2/3}`;
  (F2) arcs containing a point `k w`, `w` primitive, `|w| = o((log R/loglog R)^{3/2})`: `M << |w|^{2/3}`;
  (F3) clusters with a pair of large conjugate gcd. Each is an orientation family. "Arcs within angle
  `R^{-1/2} M^{3/2-eta}` of an axis or diagonal have bounded `M`" is correct and is (F1).
* **Combined with Theorem 1.2 and (VIII.4).** Theorem O gives `M = o(s)` exactly when `F = o(s^{3/2})`; Theorem 1.2
  adds nothing to it. With (VIII.4) the pair certificate gives: deep-window clusters (`P_+ <= M^{3/2-eta}`) with `M`
  large have no near-complementary pair of rows. That is a structural restriction, not a rate.
* **Status of the input.** The integer projection is the first genuine input in Part IX that fakes violate in a
  regime where it bites: IX.6 places `A_1(2m)` fakes (`W = O(m log m)`, primes `~ m^A`, `A > 8(K+1)`) at projected
  diameter `O(m^{5/4})` in direction `2 + i`, where Lemma P forbids `m` genuine points. It carries the rate only on
  orientation windows.

### 2.5 Item (4c): the Dirichlet choice is vacuous

For an arc midpoint `theta'` and `Q >= 1`, Dirichlet gives `b`, `1 <= |b| <= Q`, with
`|sin(theta' - arg b)| <= 1/(|b| Q)`. Then `X_b <<_C |b|(1 + sqrt R'/(|b|Q)) <= Q + sqrt R'/Q`, minimised at
`Q = R'^{1/4}`: `X_b << R'^{1/4} = N'^{1/8}`. Lemma P then gives `M << N'^{1/12}`, far weaker than `r_2(N') = N'^{o(1)}`.
The condition `X_b = o(s^{3/2})` needs both `Q = o(s^{3/2})` and `sqrt R'/Q = o(s^{3/2})`, hence
`N'^{1/4} = o(s^3) = o((log N')^3)`, which holds for finitely many `N'` only. **Vacuous** (Proved). The Dirichlet
exponent is also the truth for generic directions: **Computed** (`ref9_meaning/dirichlet.py`, continued fractions
checked against brute force): the optimum `F(theta') = min_b |b|(1 + sqrt R' |sin(theta' - arg b)|)` satisfies
`0.10 R'^{1/4} <= F <= 1.39 R'^{1/4}` on 98% of 400 random directions at each of `R' = 10^8, 10^12, 10^16, 10^24,
10^40`, against `s^{3/2} = 33, 51, 71, 114, 210`. A projection argument for general positions would need a gain of
`R'^{1/4}`, i.e. Jarník-quality information; this is the record's "near rational direction of bounded height" class
(lit18 `direct.md` §2.3: "By Dirichlet, `H` can be as large as `R^{1/4}`, which is Jarník-quality").

## 3. The residue column after IX.3 and IX.7

### 3.1 What is proved at the two ends (conditional on IX.3 and IX.7)

* **Characters, any height (capped).** Every residue input that the record tests through zero-sum characters
  (axiom (R), `A_X`, P3-char with the same-level clause) admits a fake at `W = O(m log m)` that passes it at every
  cutoff at once, with prime residues that are reductions of fixed `r_j in Q(i)^*`, `r_j conj(r_j) = p_j`.
* **All monomials, all heights, rational lifts (the problem).** If each label's exterior residues come from one
  `r_j in Q(i)^*` of norm `p_j` (finitely many exceptional primes allowed), and e-all holds at every cutoff for
  `v = n e_j`, `n >= 1`, at unit `-1` (plus `v = n(e_i + e_j)` on a spanning tree), then `r_j` is a Gaussian prime and
  `phi_j = +- arg r_j mod pi` with synchronised signs: the model is an actual cluster (Lemma T of round 21 then
  realises it).

### 3.2 Why the column is not complete

1. **The middle is P3-all\*.** The record's open local class is e-all at `X = 3M^A` (regime (1.4): `A` chosen after
   `c_0`). Neither end reaches it. IX.3 proves nothing about e-all at a polynomial cutoff; its negative paragraph
   ("which stronger axiom fails") proves failure "at a sufficiently large finite cutoff" and says explicitly that the
   bad cutoff is not shown to be polynomial in `m`. IX.7 needs every cutoff: its proof runs `n` through arithmetic
   progressions to infinity and uses the full real part of `alpha^n`, whose prime factors exceed any fixed bound.
2. **Rational lifts are free at each finite cutoff (Proved).** Let `(r_j(q))_{q <= X}` be compatible local witnesses
   with `Norm r_j(q) = p_j mod q^{A_q}`. Write `r_j(q) = pi_j t_q` with `t_q` of norm one in `(Z[i]/q^{A_q})^*`. The
   norm-one torus `x^2 + y^2 = 1` is a rational curve over `Q` (it has the point `(1,0)`), so it has weak approximation:
   there is `t in Q(i)^*`, `t conj(t) = 1`, with `t = t_q mod q^{A_q}` for all `q <= X` (and units there). Then
   `r_j = pi_j t` is a rational lift at cutoff `X`. Its height grows with `X`. So IX.7's hypothesis is a
   **uniformity in the cutoff**, not an arithmetic restriction at any finite level, and a finite-cutoff version of
   IX.7 must bound the lift's height.
3. **Heights are what a polynomial cutoff sees (Proved, elementary).** For a rational witness with
   `r_j/conj(r_j) = u alpha/conj(alpha)`, `|alpha|^2 = J_j`, the singleton collision moduli at the four units are the
   odd exterior `X`-smooth parts of `Im alpha`, `Re alpha`, `Re alpha - Im alpha`, `Re alpha + Im alpha` (for `n = 1`;
   `alpha^n` for `n e_j`). e-all needs each to be `<= nu_s^{-1} sqrt(p_j)`. For a genuine witness `|Re pi_j|, |Im pi_j| <
   sqrt(p_j)`, so the test is automatic; for `J_j > p_j` it is a smoothness event for an integer of size `sqrt(J_j)`.
4. **Computed (`ref9_meaning/witness_height.py`).** IX.3-type witnesses (Paley `q`, all `(q-1)/2` marks of the class,
   integral right inverse from LLL, columns of `V` reduced against the integral kernel of `D` by Babai rounding;
   physical split primes just above `q^{1.9}`): at `q = 11` (5 marks, 10 labels, `P ~ 95`): row norms `||F_j||_1 = log J_j/log p` from 5 to 13 (median 9); 6 of 10 labels fail singleton e-all for every angle at `X = 3m`, 7 of 10 at `X = m^2`. At `q = 19` (9 marks, 18 labels, `P ~ 268`): `||F_j||_1` from 11 to 29 (median 23, about `1.2q`); 16 of 18 labels fail at `X = 3m` and at `X = m^2`. The run for `q = 23, 31, 43, 47, 59` did not finish within this session; any output is in `ref9_meaning/witness_height.out`. On these data `log J_j` grows like `q log P`, a constant fraction of `W`.
   So IX.3's witnesses, with this choice of `V`, are nowhere near P3-all at polynomial cutoffs. A different `V` might
   lower the heights; the minimum over `V` is a closest-vector problem that I did not solve exactly.
5. **Twisted point collisions sit on the e-all side.** In factored form, `prod rho_x^{c_x} = u N^k` with
   `sum c = 2k != 0` is the e-all instance on the row-span word `v = sum c_x a_x - k e`. These are point-level inputs
   (valid for actual circles, round-5 referee), they are open at `X = M^A` since round 5, and Part IX does not treat
   them. IX.7's singletons `n e_j` are, in the Paley profiles, outside the rational row span, so IX.7 uses
   prime-level information that no point residue carries.

### 3.3 The careful statement

**Statement (conditional on IX.3, IX.4 and IX.7).** Among residue-type inputs valid for actual circles:
(a) every input tested through zero-sum characters, at any cutoff and uniformly in the cutoff, is capped by a fake
at `W <= (22/3 + o(1)) M log M` (odd prime powers; 2-adic clause not covered);
(b) e-all at every cutoff, with one rational lift per label (finitely many exceptional primes per instance), holds
only on actual clusters.
**Not covered:** (c) e-all at a polynomial cutoff (P3-all\*), with or without rational lifts; (d) twisted point
collisions at `X = M^A`; (e) e-all at all cutoffs with local witnesses that have no common rational lift; (f) e-all
on rational lifts at a finite cutoff, where the lift may have any height (by 3.2(2) this is no restriction).
Caveats: the negative result of IX.3 (failure of e-all) has an instance-dependent cutoff; the rigidity of IX.7 is
ineffective in the cutoff (irrational rotation, Cesàro means); (b) is the problem itself and gives no rate.

So the sentence "every residue-type input is either capped by a fake at `W = O(M log M)` or is the problem itself"
is **False as stated**. The correct form is the statement above: the dichotomy holds at the two ends, and the
middle, P3-all\* together with twisted collisions at polynomial cutoffs, is exactly where the record's residue
question has stood since round 6.

## 4. The closing question, and what a proof of growth could still use

### 4.1 "Can the strengthened tests force independently chosen local prime witnesses to share a rational lift?"

* **At a finite cutoff the answer is trivially yes** (3.2(2)): any finite family of compatible local witnesses has a
  rational lift. So the question only has content at all cutoffs at once.
* **At all cutoffs a positive answer gives no rate.** It would show `F^all_infinity` = actual clusters (via IX.7),
  i.e. the all-height all-monomial class collapses onto the problem. A rate needs a class strictly larger than the
  actual clusters on which `W_min/(M log M) -> infinity`. A collapse result is a rigidity theorem of the type of
  round 21's Theorem C and Lemma T; the record already has several (round 14 Rig(K); round 21 Thm C; Part VIII's
  rigidity on the Question B thread), and none gave a count.
* **At all cutoffs the class of independent local witnesses is no easier for lower bounds than P3-all\*.**
  `F^all_infinity subset F^all_{3M^A}`, so a lower bound on the former is weaker than one on the latter only if it
  uses the infinitely many cutoffs, and the only known way to use them is IX.7's ineffective rotation argument.
* **Assessment:** a restatement. The route through residues to the rate is a lower bound in P3-all\* at a
  polynomial cutoff, or a fake there. Her machinery (rational lifts, heights `J_j`, clean real parts) is the right
  language for that question; her question as posed is not.

### 4.2 Candidates left after Part IX (Assessment; probabilities are of a growth-rate theorem for general `R`
within one further exchange, and in brackets of a useful restricted or barrier result)

1. **P3-all\* in regime (1.4)** (prime-level residues on every monomial at `X = 3M^A`, `A` after `c_0`). Neither capped
   nor the problem. Part IX isolates it: characters are capped at all heights, so a lower bound must use monomials
   outside the character lattice, i.e. heights of prime witnesses (3.2(3)). Probability 1.5% (15% for a fake or a
   quantitative rigidity statement). *First step:* the finite-cutoff, height-bounded form of IX.7: for
   `X >= p_j^K`, show that e-all on `n e_j`, `n <= n_0(K)`, at the four units forces either `J_j <= p_j^{1+o(1)}` or a
   collision modulus exceeding `sqrt(p_j)`; then decide whether Legendre-profile designs exist with all witness
   heights `J_j <= p_j^{1+o(1)}` (i.e. `F` close to `I`), which the integral-right-inverse identity `DF = H` controls.
2. **Twisted point collisions at `X = M^A`.** Point-level, valid, open since round 5, untouched by Part IX, and the
   residue counterpart of IX.2's sum-two words. Probability 0.5% (10% for a cap or a restricted theorem). *First step:*
   check whether IX.3's object can be extended to twisted words at `X = 2M` and then at all `X`; the words are
   `v = K_S^T c` with `sum c != 0`, whose angles involve the common phase `theta_0`, which IX.3 leaves free.
3. **The deep window `delta >= 1/6`** (`(2 - o(1)) M log M <= P_+ <= M^{11/6}`). Restricted. Probability 0.5% for the
   general rate, 5% for either a restricted theorem at `P_+ <= M(log M)^A` or a fake at `P_+ = M^{1+eps}`.
   *First step:* the digit sup norm `E = q^{1.5+a}` is the obstruction in IX.4; lifts of sup norm `q^{1+o(1)}`.
4. **The Gaussian-integer cocycle (a).** `g_ij b_ij + g_jk b_jk + g_ki b_ki = 0`, `b_ij in Z[i]`,
   mean `log|b_ij| <= W/(4(M-1)) + log C'`. This is an integral model without the norm condition; the record says
   (Part IV) "Whether integer-matrix models without the norm condition exist at `W = O(m log m)` is **open** in both
   directions", and abc on single triangles is vacuous on squarefree profiles (§5 item 6). Probability 1% (5%).
   *First step:* build integral models (`z_i in Z[i]`, all divisibilities `g_ij | z_i - z_j` with the rows' exact
   `g_ij`, norms equal to `N` only approximately) on the Legendre profile at `W = O(m log m)`; if they exist, the
   cocycle is capped and only the norm condition is left.
5. **The archimedean order (b).** Angle models already have an archimedean order of their marks on the arc; what they
   lack is the additive structure (`z_x + z_y` related to other marks). The record: "Used in full they are the problem
   itself: the additive gap is violated by almost every fake, but its universal form is strictly stronger than
   growth" (overview §5.5). Probability 0.5%. *First step:* state the weakest additive inequality that genuine
   clusters satisfy and IX.3's object violates (for instance a three-point second-difference bound), and test it on the
   object.
6. **Projection through simultaneous approximation (c).** Vacuous for generic positions (§2.5). Probability 0%.
   A two-dimensional larger sieve on the points themselves gives `log Q <= (1/2) log R + O(1)`, i.e. only what
   Theorem 1.1's pigeonhole already uses; the projection gains only through curvature, which needs a direction of
   height `<= R^{1/4}`.
7. **IX.7 at a finite cutoff (d).** By 3.2(2) it needs a height bound on the lift; with the bound it is candidate 1.
   Probability: included in 1.
8. **Square-root-quality prime sums, Lang-Waldschmidt rows 11-12, the composite cross-gap, the cross-gap inside the
   window.** Conditional or unclassified; no mechanism proposed in the record. Probability together 0.5%.

Total for the growth rate in one exchange: about 2-3%.

## 5. Proposed wording for the expert note

### 5.1 Section 1, replace the paragraph "Actual clusters satisfy `W >= 3M log M - O(M)`. Not everything is capped. ..."
up to "... the union large sieve at the genuine constant and square-root-quality prime sums and box counts
(`round31/README_VIII.md` §3(e))." by:

"Actual clusters satisfy $W\ge3M\log M-O(M)$. Not everything is capped. Residues that factor through prime
residues of the genuine norms and are tested through characters are now capped **at every height at once**: round
31, Part IX.3 gives fakes at $W\le(22/3+o(1))M\log M$, $P_+(N')\le M^{11/6+o(1)}$, satisfying
`round5/round6_sharp/sharp.md` Definition 1.1 (d), (e-char) at every odd prime power, with the same-level clause at
every parent prime, with witnesses that are fixed elements of $\mathbb Q(i)^*$ of norm $p_j$ (**proved**, refereed,
`round31/README_IX.md`). Hence point residues at moduli $M^A$ with $A$ chosen after the cluster are capped for
zero-sum characters (2-adic and twisted collisions excepted). The single-label max-scale cross-gap
$|\sin(\phi_j-\arg\varpi)|\ge1/\max(p_j,\mathrm{Norm}\,\varpi)$ is capped at constant one outside the low-height
window (Part IX.5). Still uncapped: P3-all$^*$ (prime-level residues on every monomial at moduli $3M^A$), twisted
point collisions at moduli $M^A$, the union large sieve at the genuine constant and the composite cross-gaps,
square-root-quality prime sums and box counts, Lang–Waldschmidt, and additive data. At all moduli at once, with one
rational lift per prime, the all-monomial test holds only on actual clusters (Part IX.7), so that form is the
problem itself. The low-height window $P_+(N')\le M^{2-\delta}$ ($N'$ the normalised norm, $P_+(N')$ its largest
prime factor) is **capped for $0<\delta<1/6$**, at coefficient $22/3$ along a diagonal sequence, for pair, cut,
monomial-gap, rational-character, point-residue and character-tested prime-residue inputs (Parts VIII.3 and IX.4).
It stays **open** for $\delta\ge1/6$, and in every window for the uncapped inputs just listed. Actual clusters have
$P_+(N')\ge(2-o(1))M\log M$ (Theorem 1.2)."

### 5.2 Section 1, add after Theorem 1.2 and its bullets

"> **Theorem 1.3** (**proved**; round 31, Part IX.2, `round31/README_IX.md`). For $M$ distinct lattice points on
> $x^2+y^2=N$ and $0\ne b\in\mathbb Z[i]$, the integers $\mathrm{Re}(\bar bz_i)$ have diameter $X_b$ with
> $2+X_b\gg M^{3/2}$ (absolute constant; $\gg M^2$ if no split $p\le M$ divides $N|b|^2$). Hence an arc of length
> $\le C\sqrt R$ whose midpoint direction $\theta$ satisfies $|\sin(\theta-\arg b)|\le\kappa R^{-1/2}$ holds
> $\ll(|b|(C\kappa+C^2))^{2/3}+1$ points.

This is a restricted (position-dependent) theorem: it gives $M=o(\log R/\log\log R)$ when
$|b|(1+\kappa)=o((\log R/\log\log R)^{3/2})$, against $o((\log R/\log\log R)^{1/2})$ from the sagitta count
(`round4/hardness.md` Lemma 4.2(b)). It says nothing about generic positions: the best direction gives only
$M\ll R^{1/6}$. Sum-two combinations of exponent rows, pairs of points with a large Gaussian gcd with each other's
conjugate, and points $kw$ with $w$ primitive and small are certificates of such an orientation. The record's
ten-point cluster is symmetric about the diagonal."

### 5.3 Section 5, item 1 (replace the caveat)

"Caveat: caps uniform in $A$ are proved for zero-sum characters (round 31, Part IX.3: one fake at every cutoff,
primes below the cutoff, factored witnesses), not for twisted collisions ($\sum c\ne0$) or for P3-all$^*$, which is
uncapped. With one rational lift per prime and all cutoffs, the all-monomial test is the problem itself (Part IX.7)."

### 5.4 Section 5, item 2 (replace the last sentence)

"Inside the low-height window the earlier exceptions (`round21/README.md`) stand only for $\delta\ge1/6$, or for the
inputs listed as uncapped in §1 (round 31, Parts VIII.3 and IX.4). The single-label max-scale cross-gap at constant
one is capped outside the window (Part IX.5); the union large sieve at the genuine constant (operator form) and the
composite cross-gaps are not."

### 5.5 Section 5, item 12, new bullet

"* ***Part IX: projecting genuine marks; factored residues at all heights*** (`round31/README_IX.md`).
  * *Orientation theorem* (Theorem 1.3; **proved**): Gaussian projections of $M$ genuine points have diameter
    $\gg M^{3/2}$; restricted rate near rational directions of height $o((\log R/\log\log R)^{3/2})$; vacuous for
    generic arcs. A round-5 fake can sit in such a window (IX.6), so this input separates circles from fakes, but only
    there.
  * *All-height P3-char fakes* (IX.3; **proved**): $m/\omega(N')\to1/2$, $W\le(22/3+o(1))m\log m$,
    $P_+\le m^{11/6+o(1)}$ (IX.4), rational witnesses; they fail the all-monomial test at an unquantified cutoff.
  * *Cross-gap fakes* (IX.5; **proved**): constant one against all Gaussian primes and against lattice directions of
    norm $\le p_j$, absolute constant against all lattice directions; primes of size $m^d$, $d$ large (outside the
    window); all-height P3-char.
  * *Rigidity* (IX.7; **proved**): with one rational lift per prime, the all-height all-monomial test at unit $-1$ on
    powers of single labels forces genuine angles. Ineffective in the cutoff. At any finite cutoff rational lifts
    exist for every compatible local system (weak approximation), so the hypothesis is a uniformity in the cutoff.
  * *Unchanged:* the growth rate for general $R$, the constant $2/3$, P3-all$^*$, twisted collisions at $M^A$, the
    window for $\delta\ge1/6$."

### 5.6 Pointer table

Replace the window row by "round 21, round 31 Parts VIII and IX: low-height window (capped for $\delta<1/6$ for
inputs up to character-tested prime residues at all heights; open otherwise) | `round21/README.md`,
`round31/README_VIII.md`, `round31/README_IX.md`", and add "Theorem 1.3, orientation (projection) theorem |
`round31/README_IX.md`".

### 5.7 Sentences the note must not contain

* "Point residues at moduli $M^A$ are capped uniformly in $A$" (only zero-sum characters; twisted and 2-adic
  excepted).
* "P3-all is capped" / "P3-all$^*$ is answered" / "every residue-type input is capped or is the problem itself".
* "Part IX gives a growth theorem" or "a rate in the deep window" (it gives an orientation theorem only).
* "IX.2 gives $M=o(\log R/\log\log R)$ for a class of radii" (it is position-dependent; no parent-only consequence).
* "The max-scale cross-gap is capped" without "single-label, outside the window"; "the sharp large sieve is
  capped".
* "The window is capped for $\delta\le1/6$" (strict inequality; limsup endpoint not reached) or "for all $\delta$".
* "Rational witnesses are a restriction at a finite cutoff".
* "IX.3's witnesses fail P3-all$^*$" as a theorem (computed for small Paley orders with one choice of $V$ only).
* "The Gram identity (G3) is capped" without "with real cofactor magnitudes".

## 6. The next question, and whether the thread is converging

**Recommendation.** Ask Astra for P3-all\* in regime (1.4), in the language of IX.7, aimed at the rate:

> Fix the Legendre profile of VIII.3 (or any profile with `W = O(M log M)`). For cutoffs `X = 3M^A` with `A` chosen
> after the weight constant, either (i) construct factored residue data passing e-all (every monomial, four units,
> all odd prime powers `<= X`, same-level clauses) at `W = O(M log M)`, which closes the last local class and would
> be the strongest barrier statement of the thread; or (ii) prove that e-all at such `X` forces witness heights
> `J_j <= p_j^{1+o(1)}` on a positive proportion of labels (a finite-cutoff, height-quantified IX.7), and then that
> near-genuine witnesses force `W >= f(M) M log M` with `f -> infinity`, which is a growth theorem for actual
> circles (sharp.md §6.7: "a lower bound ... is a theorem about actual clusters"). State the twisted collisions
> (`sum c != 0`) at `X = 2M` and at `X = M^A` on IX.3's object as a side item.

Secondary: lifts of sup norm `q^{1+o(1)}` in IX.4 (deep window), only if they come cheaply.

**Is the exchange converging on a wall? Yes, and Part IX makes the wall's shape clear.** Parts VIII and IX are the
two passes about integer circles. Each capped one more input class with a fake at `W = O(M log M)`, and each proved
one restricted theorem (Theorem 1.2; Theorem 1.3) that improves nothing for general `R`. After Part IX the inputs
fall into three groups: capped (all pair, gap, character and character-tested residue inputs at every height; the
single-label cross-gap; the window for `delta < 1/6`), the problem itself (all-height rational all-monomial residues;
additive data in full; genuine angles), and open since round 6 or round 31 Part IV (P3-all\*, twisted collisions at
polynomial height, the Gaussian cocycle without the norm condition, GRH-quality prime sums). The colleague's
constructions are strong, and Part IX's fakes are probably correct; but each further fake shrinks the map without
moving the rate. One more exchange is justified only on P3-all\* with the explicit alternative (fake or height
rigidity); if it returns another fake, the residue route should be closed for internal work.

## Files

* `ref9_meaning/orient_check.py`: Lemma P on all windows of chord constant `<= 6`, `N <= 2*10^5`, `|b|^2 <= 50`.
* `ref9_meaning/tenpoint.py`: the ten-point cluster is diagonal-symmetric; projections `x + y`.
* `ref9_meaning/dirichlet.py`: optimal Gaussian direction for generic arcs; proportion of orientation positions.
* `ref9_meaning/witness_height.py`, `witness_height.out`: heights of IX.3-type rational witnesses and the singleton
  e-all test at `X = 3m` and `X = m^2`.
