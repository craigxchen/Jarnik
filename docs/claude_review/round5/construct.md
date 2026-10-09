# Round 5: cheap angle models. `W_min = Theta(M log M)`, so the local + character + residue class is sharply capped

Target: angle models with `W = O(M log M)`, or at least `W = o(M^2 log M)`, for an infinite family of `M`.
Context: `round4/outside.md` §§1-2 (angle models, Lemma 0.1, Theorem A.1, Theorem A.3), its referee report,
`growth/walsh-extraction.md` (Theorem F, Lemma F.1, Prop P, Theorem C) and `two_thirds/proof.md`.
Scripts are in `round5/checks/`. Each runs from that directory with `python3 <script> [args]`, and its output is
stored next to it as `*_out.txt`.

Conventions, as in `outside.md`:
* *Theorem*, *Proposition* and *Lemma* mean proved here with all quantifiers, or cited with the exact source.
* *Exact check* means integer or `Fraction` arithmetic, or 50-digit `Decimal` logarithms compared with margins.
* *Illustration* means floating point or a single random sample. It is never used as a proof.
* `W = log N = 2 log R`. An arc of length `<= C sqrt R` has angular width `Delta <= C e^(-W/4)`.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The target is reached. Angle models with `W = O(M log M)` exist for every `M`, including residue-collision
strengthening. So the method class is sharply capped at order `log R / loglog R`, and no growth improvement is
possible inside it.**

**Attribution (read first).**
* *Without residue data* the result was already in the research notes, and `outside.md` and its referee missed it.
  * Item 392 (`paley_canonical_all_character_half_height.md`) proves the character margins of Theorem 1 below for
    the canonical diagonal-flip Paley profile. Its mechanism is the same one used here: a skew-Hadamard
    tournament, at most one beneficial full-magnitude flip, and a Parseval defect.
  * Item 394 (`paley_compatible_phase_coordinate_gap.md`) builds simultaneously admissible phases with actual prime
    weights at `W = 20 M log M + O_C(M)`. Combined with `outside.md` Theorem A.1 these are fake circles.
* *A concurrent write-up in this round*, `round5/lower.md`, reaches the same conclusion independently, with residue
  data, at `(48 + o(1)) M log M`. Its route is an `l^1`-stability lemma for the Paley transform ("Lemma S", built
  on item 393) plus a residue design of the same type as Lemma 4.1 here. It also makes the observation about
  items 392 and 394, and the correction of Prop 1.1.
* *What this write-up adds*:
  * the residue version at the better constant `27`, via the short rank-one-rectangle form of item 392, which needs
    no `l^1`-stability;
  * the constant `18` without residues;
  * the general skew-pairing statement, and the fact that Walsh cores of order `>= 16` have no skew pairing, which
    links the dichotomy to Theorem C;
  * an independent exact validation of the residue axiom on actual points, and exact concrete instances.

1. **Theorem 1 (character margins; item 392 of the notes in rank-one form).** Let `H` be a Hadamard matrix with a
   *skew pairing* (§2). Paley-I matrices of order `q+1`, `q = 3 mod 4` prime, have one. Take the fully flipped
   profile: `b` unflipped copies of every label, and one flipped copy per row, of the label paired with that row.
   Then every nonzero zero-sum
   integral character `c`, with `n = ||c||_1`, satisfies
   * `G(c) >= b - 4` for pairs;
   * `G(c) >= 2n + b - 8` for every other `+-1`-valued `c`;
   * `G(c) >= (b-3)n/2 + b` if `||c||_inf >= 2`.

   For `b = 8` this gives `G >= 4` on pairs and `G >= 2n` on every other character. Here
   `G(c) = sum_j (|(S^T c)_j| - 1)`, so `V_c - W/2 = (tau/2) G(c)` for equal weights `tau`.

   No `l^1`-stability (inverse-uncertainty) theorem is needed. Take a `+-1` character and the labels where it is
   exactly flat (`|T_a| = ||c||_1`); together with the support they form a rank-one rectangle of `H`. A skew
   pairing meets every rank-one rectangle in at most one row. So the obstacle (P1) of `outside.md` §2.2 was not
   a real obstacle; item 392 had already bypassed it.
2. **Lemma 4.1 (a residue design).** The point residues `rho_x = s * omega_x / conj(omega_x) mod q^a` are taken from
   small auxiliary Gaussian primes `omega_x` with `N(omega_x) = l_x > Q_0`. Every collision modulus of a
   character `c` then divides the nonzero Gaussian integer `Lambda_c - u conj(Lambda_c)`, whose norm is
   `<= 4 L^||c||_1`, with `L = max l_x`. So residue strengthening at every odd prime power `<= Q_0` costs only the
   factor `L^(n/2)`.
3. **Theorem 2 (fake circles at `W = O(M log M)`).** Take `C in (0,1]`, `Q_0 >= 3`, `q = 3 mod 4` prime and
   `M = q+1`. There is an angle model with the following properties:
   * genuine primes `= 1 mod 4` and a genuine exponent profile, with `b = 8`;
   * real angles satisfying **every** Lemma 0.1 gap, for every monomial;
   * residue data satisfying the residue-collision axiom (R) for every character at every odd prime power
     `<= Q_0`;
   * `M` distinct points on an arc of length `<= C sqrt R`.

   Its size is `W <= (27 + o(1)) M log M + 9 M log(1/C)` for `Q_0 = 2M`. Without residue data it is
   `(18 + o(1)) M log M + 9 M log(1/C)`. For `Q_0 = M^K` it is `(9(2+K) + o(1)) M log M + 9 M log(1/C)`.
4. **Corollary 6.1.** Let `W_min` be taken over the valid class (angle models with residue data, §1). For every
   `C > 0`,

   ```text
   (3 - o(1)) M log M  <=  W_min(M, C)  <=  (27 + o(1)) M log M + O_C(M)     (M = q+1),
   W_min(M, C) = O_C(M log M)                                               (all M).
   ```

   The lower bound is the `two_thirds` inequality on this class, as stated in the task (§1.2).

   **Corollary 6.2 (sharp barrier).** Consider any argument that uses only:
   * exponent identities;
   * Lemma 0.1 gaps for all monomials;
   * the residue axiom (R) at prime powers `<= 2M`;
   * Haar-generic properties of the angles.

   No such argument can prove `M <= (2/27 - eps) log R/loglog R`. The best constant for this class therefore
   lies in `[2/27, 2/3]`, and `two_thirds/proof.md` attains `2/3`.
5. **Correction to the framing of the question.** Read literally, "residue strengthening as in A.3" means
   strengthening every character gap by `Q_M = lcm(1..2M)` for every residue assignment (A.3, property 3).
   * That is not a valid input for actual circles: it fails already for `2+i, 1+2i` on `x^2+y^2 = 5`.
   * By itself it forces `W >= (8 - o(1)) M^2` for every model (Prop 1.1).

   The class for which "`W_min/(M log M) -> infinity` would imply a growth improvement" is the class with
   residue data and axiom (R) (§1.2). That is the class for which the task's lower bound `3 M log M - O(M)`
   holds, and it is the class used here.
6. **Consistency with refereed results.**
   * Walsh cores of order `>= 16` admit no skew pairing (exhaustive search plus a restriction argument,
     Remark 2.3). This is consistent with Theorem C of `walsh-extraction.md`, which excludes flipped Walsh
     profiles by aligned flats.
   * Prop P (fully flipped Paley passes each character **individually**) is upgraded here to **simultaneous**
     admissibility, which is what Theorem A.1 requires.
   * The fully flipped Hadamard phase lemma of `walsh-extraction.md` §10 cannot be proved by inputs (I)-(IV)
     plus (R). Problem (P2) of `outside.md` has a negative answer for skew-paired Paley profiles.

What remains outside the class is listed in §7. It consists of:
* residue data at the level of the individual Gaussian primes, which carries quadratic-reciprocity (parity)
  constraints that the present design does not satisfy;
* the actual angles: transcendence, e.g. Theorem L of `round4/flipped.md` under Lang-Waldschmidt;
* exact algebraic identities;
* multi-place approximation.

---------------------------------------------------------------------------------------------

## 1. Definitions: angle models with residue data

### 1.1 Angle models (recalled from `outside.md` §1)

Fix distinct primes `p_j = 1 mod 4` (`j <= r`), exponents `e_j`, rows `a_x in prod [0, e_j]` (`x <= M`), arc
positions `delta_x in [0, Delta]` and integers `k_x`. Put `N = prod p_j^(e_j)`, `W = log N`, `w(v) = sum |v_j| log p_j`.
* `phi in R^r` *realises* the profile if `<2a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x` for all `x`.
* It is *gap-admissible* if Lemma 0.1 holds for all `v != 0`:
  `dist(<v,phi>, (pi/2)Z) >= arcsin(e^(-w(v)/2))` and `dist(<v,phi>, pi/4 + (pi/2)Z) >= arcsin(e^(-w(v)/2)/sqrt 2)`.
* The model's points are `R i^(-k_x) exp(i <2a_x - e, phi>)`, with point units `eps_x := i^(-k_x)`.
* For an integral zero-sum `c`, the character monomial is `v(c) = sum_x c_x a_x`. Its value is
  `<v(c), phi> = (1/2) sum_x c_x delta_x + (pi/4) c.k`. The characters of the profile are the `v(c)`, together
  with the rational `v = sum lambda_x a_x in Z^r`.

### 1.2 Residue data and the residue axiom (R)

Fix `Q_0 >= 3`. Let `Q` be the set of odd prime powers `m = q^a <= Q_0` with `q` not dividing `N`.

**Definition 1.2.** *Residue data* is a family `rho = (rho_x^(m))`, `m in Q`, with:
* `rho_x^(m) in (Z[i]/m)^*`;
* `rho_x^(m) conj(rho_x^(m)) = N (mod m)`, i.e. the point lies on the conic;
* `rho_x^(q^a) = rho_x^(q^(a-1)) mod q^(a-1)`.

For a nonzero zero-sum `c in Z^M` and a unit `u in {+-1, +-i}`, put

```text
D_u(c) = prod_q q^(a_q),   a_q = max{ a : q^a in Q,  prod_x (rho_x^(q^a))^(c_x) = u  (mod q^a) }   (0 if none).
```

The **residue axiom (R)** is

```text
|sin(<v(c), phi> - psi_u(c))|  >=  D_u(c) e^(-w(v(c))/2) / sqrt 2,     psi_u(c) = arg(u)/2 + (pi/4) c.k,     (R)
```

for all such `c` and `u`. For pairs `c = e_x - e_y` and `u = 1` it is the collision-strengthened pair gap of
`two_thirds/proof.md` (Lemmas 3.1-3.3), stated for angles. If a small prime `q <= Q_0` divides `N`, one adds the
level-wise cofactor data of `two_thirds` Lemma 3.3 (S|). In all models built below every prime of `N` exceeds
`Q_0`, so that part is vacuous and we do not spell it out.

**Lemma 1.2 (validity).** Let `z_1, ..., z_M` be a normalised actual cluster (`two_thirds` Lemma 1.1), with
`rho_x^(m) = z_x mod m`. Then (R) holds.

*Proof.*
1. Write `z_x = eps_x prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj)`. Let `Z_+ = prod_{c_x>0} z_x^(c_x)`,
   `Z_- = prod_{c_x<0} z_x^(-c_x)`, `E = prod eps_x^(c_x)` and `A = A_v(c)`.
2. Since `sum c = 0`, we get `Z_+ / Z_- = E A / conj(A)`. Now `A` and `conj A` are coprime, because each `p_j`
   contributes only `pi_j` or only `conj pi_j`. Hence `conj(A) | Z_-`. Write `Z_- = conj(A) K`; then
   `Z_+ = E A K`.
3. If `prod rho_x^(c_x) = u mod q^a`, then `q^a | Z_+ - u Z_- = K (E A - u conj A)`. Since `Norm K` divides
   `N^(n/2)` and `q` does not divide `N`, we get `q^a | E A - u conj A`.
4. Products over distinct `q` divide as well, so `D_u(c) | E A - u conj A`.
5. Both `E A` and `u conj A` are `= 1 mod (1+i)`, because they have odd norm. So `(1+i) D_u(c)` divides
   `E A - u conj A`.
6. `E A - u conj A` is nonzero, since `A` is a conjugate-primitive nonunit. Hence
   `sqrt2 D_u(c) <= |E A - u conj A| = 2|A| |sin(theta + (arg E - arg u)/2)|`, with `theta = arg A`.
7. With `arg E = -(pi/2) c.k` (from `eps_x = i^(-k_x)`) this is (R). QED

`check_axiomR.py` checks the integer form `|Z_+ - u Z_-|^2 Norm(A_v) >= 2 D_u^2 N^(n/2)` exactly on 6,336
(character, unit) instances on four circles. Equality occurs, so the constant `sqrt 2` cannot be improved.

**Classes.**
* `A_0`: angle models.
* `A_1(Q_0)`: angle models with residue data satisfying (R) at all `m in Q`.
* `A_1 = A_1(2M)`, which matches `outside.md` Corollary A.4 (III).

By Lemma 1.2 every actual normalised cluster lies in `A_1(Q_0)` for every `Q_0`. Hence lower bounds on `W` proved
for a class transfer to actual circles. The task asserts that `two_thirds` gives
`3 M log M <= W + K_C M` on `A_1`. Every input of `two_thirds` Theorem 5.1 is:
* a pair instance of (R) or of its level-wise version;
* the cut identity;
* or pigeonhole on residue data in the conic fibres (`two_thirds` Lemma 3.2 counts the fibres).

### 1.3 The literal A.3 class is not a modelling class

`outside.md` Theorem A.3, property 3, requires `dist(<v,phi>, (pi/4)Z) >= arcsin(Q_M e^(-w(v)/2))` for **every**
character `v in L \ {0}`, with `Q_M = lcm(1, ..., 2M)`. Call the class of angle models with this property `A_U`.
It is contained in `A_1`, because `D_u(c) <= Q_M`.

**Proposition 1.1.**
(i) Every model in `A_U` has `W >= 4(M-1)(log Q_M + log(2/C)) = (8 - o(1)) M^2`.
(ii) Actual circles need not satisfy property 3.

*Proof.*
(i) For a pair `c = e_x - e_y`, `dist(<v(c),phi>, (pi/4)Z) <= |delta_x - delta_y|/2 <= Delta/2 = (C/2) e^(-W/4)`.
* Property 3 therefore forces `d_xy - W/2 >= 2 log Q_M + 2 log(2/C)`, where `d_xy = w(a_x - a_y)`.
* The cut identity (`two_thirds` Lemma 2.2) gives `sum_{x<y} d_xy <= W M^2/4`, hence
  `sum_{x<y} (d_xy - W/2) <= W M/4`.
* Comparing the two gives (i). The asymptotic form uses `log Q_M = psi(2M)`.

(ii) Take `z = 2+i` and `z' = 1+2i` on `x^2 + y^2 = 5`. The chord is `sqrt 2 <= 1.0 sqrt R`. With `M = 2`,
`Q_2 = 12`, and property 3 asks for `arcsin(12/sqrt 5)`, which does not exist. In the valid axiom (R), `z - z' = 1 - i`
collides modulo no odd prime, so `D = 1`. QED

So with the literal reading, `W_min >= 8M^2` would hold trivially. It would not transfer to actual circles,
and the implication in the task fails for that reading. The task's stated lower bound `3 M log M - O(M)` is the
`A_1` bound. `A.3`'s fakes lie in `A_U`, a subset of `A_1`, so they give the old upper bound
`W_min^(A_1) <= (30+o(1)) M^2 log M`.

**Lemma 1.3 (monotonicity).** `W_min^(A_i)(M', C) <= W_min^(A_i)(M, C)` for `M' <= M`, and
`W_min(M, C) <= W_min(M, min(C, 1))`.

*Proof.* Deleting rows keeps every inequality:
* the characters of the sub-cluster are characters of the cluster;
* the monomial gaps concern `phi` alone;
* the residue data restricts.

An arc of length `<= sqrt R` is an arc of length `<= C sqrt R` for `C >= 1`. QED

---------------------------------------------------------------------------------------------

## 2. The profile: skew-paired, fully flipped Hadamard cores

Let `H` be a Hadamard matrix of order `M >= 4`, normalised so that column `0` is all ones. Its rows are the
*rows*, and columns `1..M-1` are the *labels*.

**Definition 2.1 (skew pairing).** A skew pairing is a row `x_inf` together with a bijection
`pi : rows \ {x_inf} -> {1, ..., M-1}` such that for all distinct `x, y != x_inf`

```text
H(x, pi(x)) H(y, pi(y))  !=  H(x, pi(y)) H(y, pi(x)).
```

That is, every `2 x 2` "principal" minor on the paired rows and labels is nonzero. The condition is invariant
under row and column sign changes.

**Paley-I.** Let `q = 3 mod 4` be prime, `chi` the Legendre symbol, and index rows and columns by `{inf} cup F_q`.
* `H(inf, .) = 1`, `H(i, inf) = -1`, `H(i, j) = delta_ij + chi(j - i)`.
* After multiplying the rows `i in F_q` by `-1`, `H_n(i, j) = -delta_ij - chi(j - i)`.
* The pairing `pi(i) = i` is skew: `H_n(i,i) H_n(j,j) - H_n(i,j) H_n(j,i) = 1 - chi(-1) chi(j-i)^2 = 2`.

More generally, every Hadamard matrix equivalent to a skew-Hadamard matrix has a skew pairing.
`check_G_bound.py` (L) checks the minor value `2` for all 23 primes `q = 3 mod 4` with `7 <= q <= 199`.

**The profile `P(H, pi, b, a*)`.** Fix `b >= 3` and a label `a*`. Put `a(x) = pi(x)` for `x != x_inf` and
`a(x_inf) = a*`. The columns are:
* `b` unflipped copies `H(., a)` of every label `a = 1..M-1`;
* for every row `x`, one *flipped* copy of label `a(x)`, namely `H(., a(x))` with the entry in row `x` negated.

So `r = b(M-1) + M`, `S in {+-1}^(M x r)`, rows `a_x = (1 + S_x)/2`, `e = 1`. Every label carries `b + 1` copies,
except `a*`, which carries `b + 2`. This is the "fully flipped capacity-two" class of Prop P.
* *Twins.* The flipped copy of row `x` and any unflipped copy of `a(x)` differ exactly in row `x`. So
  `rank S = M`, the profile is realisable for every `delta`, and every rational character is integral:
  `L = {v(c)}`. This is the mechanism of `walsh-extraction.md` Lemma F.1 and `outside.md` Lemma A.7(ii).
* Every column is nonconstant, so the profile is normalised.

**Lemma 2.2 (exact decomposition).** For an integral zero-sum `c`, put
* `T_a = sum_x c_x H(x, a)`, so that `T_0 = 0`;
* `F = sum_{a>=1} |T_a|`;
* `h_x = H(x, a(x))`;
* `Delta_x = |T_a(x) - 2 c_x h_x| - |T_a(x)|`.

Then

```text
G(c) := sum_j (|(S^T c)_j| - 1)  =  (b+1)(F - M) + b + |T_a*| + sum_x Delta_x.
```

*Proof.*
1. The unflipped columns give `b sum_a (|T_a| - 1) = b(F - M + 1)`.
2. The flipped column of row `x` has `S^T c = T_a(x) - 2 c_x h_x`. These columns give
   `sum_x (|T_a(x)| - 1) + sum_x Delta_x`.
3. Since `a` hits every label once and `a*` twice, `sum_x |T_a(x)| = F + |T_a*|`. QED

`check_G_bound.py` (E) checks this identity exactly on up to 10,000 random integer characters (2,000 draws each for
`(q, b) = (11,5), (19,5), (23,7), (43,5), (67,12)`, zero vectors skipped).

**Remark 2.3 (Walsh cores have no skew pairing).** For the Walsh/Sylvester core `H(x, a) = (-1)^(x.a)` the
condition reads `(x + y).(pi(x) + pi(y)) = 1 (mod 2)` for all distinct `x, y`. Translations and the choice of
constant column reduce to `x_inf = 0` and `pi : F_2^t \ 0 -> F_2^t \ 0`.

Exhaustive backtracking:
* `sylvester_core_skew.py`: pairings exist for `t = 2, 3` and do **not** exist for `t = 4, 5`;
* `sylvester_skew.py`: full skew-equivalence of the Sylvester matrix exists for `t = 2, 3` and does **not** exist
  for `t = 4`;
* `sylvester_core_skew_relaxed.py`: there is no injective `pi : F_2^4 \ 0 -> F_2^4` with the property, even
  with the value `0` allowed.

Restricting a pairing for `t >= 4` to a 4-dimensional subspace `U` gives such an injection. Let
`pibar(x) = pi(x)|_U in U^*`. Then `(x+y).(pi(x)+pi(y))` depends only on the restrictions, and the condition forces
`pibar` to be injective. Hence **no Walsh core of order `>= 16` has a skew pairing**. This agrees with
`walsh-extraction.md` Theorem C. For flipped Walsh profiles, the Ramsey step of item 333 finds, for every flip
assignment, flats aligned up to boundedly many bad rows. By §3 an aligned flat `P` with functional `l` is a
rank-one rectangle `P x A_l` with `a(P) subset A_l`, i.e. one meeting the flips `|P|` times, while a skew pairing
allows at most twice.

---------------------------------------------------------------------------------------------

## 3. Theorem 1: character margins without an `l^1`-stability theorem

This is research-notes item 392, `paley_canonical_all_character_half_height.md`, equation (5), with `h = n`.
* The notes count all `b` physical copies of a label, the flipped one included. With that convention their
  `B - r/2 >= b/2 + h - 4` is `G >= b + 2h - 8`.
* Here the flipped copy comes on top of `b` unflipped ones, so the constants shift by one copy.

The proof below is the same argument in rank-one-rectangle form. It is stated for an arbitrary skew pairing,
because that is what separates Paley-I from Walsh (Remark 2.3).

**Theorem 1.** Let `H` have a skew pairing, `b >= 3`, and let `S = P(H, pi, b, a*)`. For every nonzero
integral `c` with `sum c = 0` and `n = ||c||_1`:

```text
(i)   G(c) >= b - 4                     if n = 2;
(ii)  G(c) >= 2n + b - 8                if ||c||_inf = 1 and n >= 4;
(iii) G(c) >= (b - 3) n / 2 + b         if ||c||_inf >= 2.
```

In particular, for `b = 8`: `G(c) >= 4` on pairs and `G(c) >= 2n` on all other characters.

*Proof.* Throughout, `F >= M`, because `sum_a T_a^2 = M ||c||_2^2` (`H H^T = M I`), `|T_a| <= n` and
`||c||_2^2 >= n` (`outside.md` Lemma A.7(i)). Also `Delta_x = 0` for `c_x = 0`, and `Delta_x >= -2|c_x|` always.

(iii) Let `K = ||c||_inf >= 2`.
1. From `c = H T / M` we get `|c_x| <= F/M`, so `F >= KM`.
2. With `sum Delta >= -2n`, Lemma 2.2 gives `G >= (b+1)(K-1)M + b - 2n`.
3. Since `n <= KM` and `(K-1)/K >= 1/2`, we have `(b+1)(K-1)M >= (b+1) n/2`.

(i) A pair `c = e_x - e_y` has `T_a in {0, +-2}`, nonzero on exactly `M/2` labels, so `F = M`. Lemma 2.2 then
gives `G = b + |T_a*| + Delta_x + Delta_y >= b - 4`.

(ii) Let `X = supp c`, so `|X| = n` and `c_x = +-1` on `X`. Parseval gives `sum_a T_a^2 = M n`, hence

```text
sum_{a>=1} |T_a| (n - |T_a|) = n (F - M),        every term >= 0.                       (3.1)
```

All `T_a` are even, and so is `n`. Split the labels with `T_a != 0` into two sets:
* the *big* set `A = {a : |T_a| = n}`;
* the *medium* set `L' = {a : 2 <= |T_a| <= n-2}`.

Each medium label contributes at least `2(n-2)` to (3.1), so `|L'| <= n(F - M)/(2(n-2))`.

* *Big labels form a rank-one rectangle.* `|T_a| = n` means that `c_x H(x,a)` has constant sign on `X`. Hence
  `H[X, A] = c|_X sigma^T` has rank one.
* *The pairing meets it at most once.* If `x != y` are paired rows in `X` with `pi(x), pi(y) in A`, then the
  `2 x 2` submatrix on rows `{x, y}` and labels `{pi(x), pi(y)}` is a submatrix of `H[X, A]`, so it has rank one.
  That contradicts Definition 2.1. Together with the unpaired row, `#{x in X : a(x) in A} <= 2`.
* *Flip terms.* For `x in X` with `T_a(x) = 0`, `Delta_x = |2 c_x h_x| = 2`. For the other `x in X`,
  `Delta_x >= -2`. Let `N_X = #{x in X : T_a(x) != 0}`. Then `sum_x Delta_x >= 2n - 4 N_X` and
  `N_X <= 2 + |L'|`. Here `a^(-1)(L')` has at most `|L'|` paired rows; the unpaired row is already counted in
  the `2`.
* Therefore

  ```text
  G >= (b+1)(F-M) + b + 2n - 8 - 4|L'| >= (b + 1 - 2n/(n-2)) (F - M) + b + 2n - 8 >= 2n + b - 8,
  ```

  because `2n/(n-2) <= 4 <= b+1` for `n >= 4`. QED

**Exact checks** (`check_G_bound.py`, all integer arithmetic). The bound in (i)-(ii), with the slack reported:

| `q` | `M` | `b` | characters checked exhaustively | min slack over bound |
|---|---|---|---|---|
| 11 | 12 | 5, 8 | all `+-1` characters, every `n <= 12` (73,788 per `b`) | 2 at pairs, 4 at half-sums and columns |
| 19 | 20 | 5 | all `+-1`, `n <= 8` (9,622,550) | 2 (pairs) |
| 23 | 24 | 5 | all `+-1`, `n <= 6` (2,756,228) | 2 (pairs) |
| 31 | 32 | 5 | all `+-1`, `n <= 6` (18,340,592) | 2 (pairs) |
| 43 | 44 | 5 | all `+-1`, `n <= 4` (816,398) | 2 (pairs) |

The general bound (iii) was checked exhaustively on every `c` with entries in `[-3,3]` and `||c||_inf >= 2`:
`||c||_1 <= 8` for `M = 12` (636,834 vectors, min slack 62) and `||c||_1 <= 6` for `M = 20` (796,480 vectors,
min slack 110). The script also lists `M = 24` for (iii). That part was stopped at the session time limit, is not
in the output file, and is not used.

`check_G_large.py` covers `M = 44, 104, 200` (`b = 8`) and `M = 44, 200` (`b = 5`). It runs over all pairs, all
`+-` columns, all two-column half-sums, 20,000 random `+-1` characters and 40 adversarial hill-climbs from random and
column starts.
* With the skew pairing the minimum slack is `2`, attained at pairs, in every case.
* With a random row-to-label assignment it is `0`, also at pairs, so the bound `b - 4` is attained.
* No other character came within the bound.

**Remarks.**
1. *What happened to (P1).* Items 392-394 of the notes already contained the answer, and `outside.md` §2.2 did not
   use them. That section bounded the flip terms by `-2n` and then needed
   `F - M >= gamma n` away from pairs, columns and half-sums: an `l^1`-stability theorem for the Paley transform,
   which looked as hard as Paley clique bounds.
   * Lemma 2.2 shows the flip term is `+2` on every support row whose paired label has `T = 0`.
   * Medium labels are paid for by `F - M` through (3.1).
   * The only dangerous labels are the exactly flat ones, which form a rank-one rectangle with the support.
   * The skew pairing lets such a rectangle hit the flips at most once, whatever the rectangle is.
   * So no structure theorem for flat or near-flat vectors of `H` is needed.
2. *The pairing matters.*
   * For a random flip assignment the same proof only gives `N_X <= #{x in X : a(x) in A} + |L'| + 1`. Large
     rank-one rectangles `X x A` with `|X cap a^(-1)(A)| >= n/2` would break it. For Paley such rectangles,
     other than pairs, columns and half-sums, are not known to be absent, which is again a Paley clique-type
     question.
   * For Walsh cores, the Ramsey step of item 333 (used in Theorem C) finds, for **every** assignment, flats that
     are aligned up to boundedly many bad rows.
3. *Sharpness.* The pair bound `b - 4` is not attained with the skew pairing (observed minimum `b - 2`). The
   column and half-sum values `G = 2n + b - 8 + 4` are close to (ii).

---------------------------------------------------------------------------------------------

## 4. A residue design whose collisions are controlled by a small Gaussian integer

**Lemma 4.1.** Let `Q_0 >= 3`, and let `l_1 < ... < l_M` be distinct primes `= 1 mod 4`, all `> Q_0`. Put `L = l_M`,
let `omega_x` be a Gaussian prime above `l_x`, and let `N` have no prime factor `<= Q_0`. Then:
* There is `s in Z[i]` with `s conj(s) = N (mod m)` for every `m in Q`.
* `rho_x^(m) := s omega_x conj(omega_x)^(-1) mod m` is residue data in the sense of Definition 1.2.
* For every nonzero zero-sum `c` and every unit `u`, with `Lambda_c := prod_{c_x>0} omega_x^(c_x) prod_{c_x<0} conj(omega_x)^(-c_x)`,

  ```text
  sqrt2 D_u(c)  <=  |Lambda_c - u conj(Lambda_c)|  <=  2 L^(n/2),   n = ||c||_1.
  ```

  In particular `D_1(c)` divides the odd part of `Im Lambda_c`.

*Proof.*
1. *Existence of `s`.* By `two_thirds` Lemma 3.2, the norm `(Z[i]/q^a)^* -> (Z/q^a)^*` has nonempty fibres
   for odd `q^a`. Take `s` modulo the largest power of each `q <= Q_0` and combine by CRT.
2. *Norms and compatibility.* Both are clear, since `omega_x conj(omega_x)^(-1)` has norm `1` modulo `m`, and
   `conj omega_x` is a unit modulo `m` because `q != l_x`.
3. *Collisions.* Since `sum c = 0`, `prod rho_x^(c_x) = Lambda_c / conj(Lambda_c)`. So a collision modulo `q^a`
   with unit `u` means `q^a | Lambda_c - u conj(Lambda_c)`.
4. *Nonvanishing.* `Lambda_c - u conj(Lambda_c) != 0`, because `Lambda_c / conj(Lambda_c)` has `omega_x`-valuation
   `+-|c_x| != 0` for some `x`.
5. *The factor `(1+i)`.* It divides the difference, since both terms have odd norm.
6. Combining 3-5 over the distinct `q` gives `sqrt2 D_u <= |Lambda_c - u conj Lambda_c| <= 2|Lambda_c| = 2 prod l_x^(|c_x|/2)`. QED

So the residues are those of the small numbers `omega_x / conj(omega_x)` (Cayley-type, norm one). A collision
modulus is a divisor of a fixed nonzero Gaussian integer of norm `<= 4 L^n`. This is what keeps residue
strengthening at `O(n log L)`.

A random assignment behaves differently: its worst pair collides at about `(log M)/(loglog M)` primes, each of
weight about `log M`. That heuristic is in `walsh-extraction.md` §7.2. It would cost a factor
`log M / loglog M` in `W`.

`instance.py` checks Lemma 4.1 exactly for `M = 12, 20, 24, 44`. Over all characters with `||c||_1 <= 6`
(`M = 12`, 66,594 of them) or `<= 4` (the others):
* `D_u(c)` is computed from the residue data;
* `D_1(c) | Im Lambda_c` holds;
* `D_u(c) <= 2 L^(n/2)` holds;
* every `rho_x` lies in the conic fibre.

---------------------------------------------------------------------------------------------

## 5. Theorem 2: fake circles with `W = O(M log M)`

**Theorem 2.** Let `C in (0, 1]`, `Q_0 >= 3`, `q >= 7` a prime `= 3 mod 4`, `M = q + 1`, `b = 8` and
`r = 9M - 8`. Let `L` be the `M`-th prime `= 1 mod 4` exceeding `Q_0`, and suppose

```text
P >= P_* := max{ 160 M^2 L / C,  21 r^2 },        p_1 < ... < p_r  primes = 1 mod 4 in [P, P e^(1/(4r))].
```

Then there is an angle model in `A_1(Q_0)` with:
* primes `p_j`, `e = 1`, the profile `P(H_q, pi, 8, a*)` of §2 (Paley-I, identity pairing), and `k = 0`;
* `M` distinct points on an arc of length `<= C sqrt R`;
* all Lemma 0.1 gaps for every `v in Z^r \ {0}`;
* the residue data of Lemma 4.1, satisfying (R) for every character and unit;
* `W <= r (log P + 1/(4r))`.

The same holds in `A_0` with `L` replaced by `1`.

**Window lemma.** Put `P_** = max{P_*, x_1}`, where `x_1` is the least `x >= x_0` with `x >= 8.4 r^2 log x`, and `x_0`
is such that `pi(2x;4,1) - pi(x;4,1) >= x/(3 log x)` for `x >= x_0`. Such an `x_0` exists by the prime number
theorem for the progression `1 mod 4`, and explicit versions of that theorem (e.g. Bennett-Martin-O'Bryant-Rechnitzer,
*Illinois J. Math.* 2018) make it explicit. Then a window as required exists with `P_** <= P <= 2P_**`.

*Proof.* `[log P_**, log 2P_**]` splits into `ceil(4r log 2) <= 2.8r` windows of log-width `<= 1/(4r)`, using
`r >= 34`. One of them holds `>= P_**/(8.4 r log P_**) >= r` primes `= 1 mod 4`. QED

Consequently, as `M -> infinity`:
* for `Q_0 = 2M`, `L ~ 2M log M` and `W <= (27 + o(1)) M log M + 9 M log(1/C)`;
* for `Q_0 = M^K`, `W <= (9(2 + K) + o(1)) M log M + 9 M log(1/C)`;
* in `A_0`, `W <= (18 + o(1)) M log M + 9 M log(1/C)`, where the window condition `P ~ 8.4 r^2 log P` dominates.

*Proof.*
1. *Weights.* Let `tau = log P`, so that `w_j in [tau, tau + eta]` with `eta <= 1/(4r)`. For a character `c` with
   `y = S^T c`, and since the terms with `y_j = 0` lose at most `eta/2` each,

   ```text
   V_c - W/2 = (1/2) sum_j w_j (|y_j| - 1) >= (tau/2) G(c) - eta r/2 >= (tau/2) G(c) - 1/8.
   ```

2. *Target.* We find `delta in [0, Delta]^M`, `Delta = C e^(-W/4)`, such that for every nonzero zero-sum `c`

   ```text
   dist(X_c, (pi/4) Z)  >=  arcsin(L^(n/2) e^(-V_c/2)),       X_c = (1/2) c.delta = <v(c), phi>.     (5.1)
   ```

   By Lemma 4.1, `D_u(c) e^(-w/2)/sqrt2 <= L^(n/2) e^(-w/2)`. Every `psi_u(c)` lies in `(pi/4)Z`, and
   `|sin t| >= sin dist(t, pi Z)`. So (5.1) implies (R) for all units. It also implies both Lemma 0.1 inequalities
   for every `v in L \ {0}`, since `L^(n/2) >= 1`.

   The argument of `arcsin` is tiny. We have `V_c >= (tau/2) sum_j |y_j| >= 4 tau F >= 4 tau n`, because there are
   `b = 8` unflipped copies, `F >= M ||c||_inf` and `n <= M ||c||_inf`. So `arcsin t <= 1.0368 t` applies.

3. *One character.* Draw `delta` uniformly. `X_c` has density `<= 2/Delta`: condition on all `delta_y` except one
   with `c_x != 0`. Its range has length `<= Delta n/2`, and it meets `<= 2 Delta n/pi + 2 <= n + 2` windows
   around points of `(pi/4)Z`; here `Delta <= 1`, because `W >= 4` and `C <= 1`. Using step 1,

   ```text
   P(5.1 fails for c) <= (4 * 1.0368 / Delta) L^(n/2) e^(-V_c/2) (n + 2)
                      <= (4.415 / C) L^(n/2) e^(-tau G(c)/4) (n + 2).
   ```

4. *Union bound.* There are at most `(2M)^n` characters with `||c||_1 = n`, since each is a sum of `n` signed unit
   vectors. Theorem 1 with `b = 8` gives `G >= 4` for `n = 2` and `G >= 2n` for `n >= 4`. Put
   `x := 4 M^2 L e^(-tau) <= C/40`. Then

   ```text
   P(some c fails) <= (4.415/C) [ 4M^2 L e^(-tau) * 4 + sum_{n>=4 even} (n+2) x^(n/2) ]
                   <= (4.415/C) (4x + 6.3 x^2) <= (4.415/C)(4.16)(C/40) < 0.46.
   ```

   So the set of good `delta` has measure `> 1/2`.

5. *Completion.* Fix a good `delta`. The twins give rank `M`, so the realisations form a coset `phi_0 + U`.
   * On it, every `v in L \ {0}` has value `X_c` and satisfies (5.1). The profile is therefore
     character-admissible.
   * Since `p_j >= 21 r^2`, `outside.md` Theorem A.1(ii) (with the referee's constant `1.3201`) gives a
     Haar-positive set of realisations satisfying the `(pi/4)Z` form of Lemma 0.1 for all `v` outside `L`.
   * The referee's Fubini argument (§3.1) upgrades "Haar-null per torus" to "Haar-null in `(R/2 pi Z)^r`".
6. *Points.* The angles of the points are `theta_0 + delta_x in theta_0 + [0, Delta]`, an arc of length
   `R Delta = C sqrt R`. They are distinct by the pair case of (5.1).
7. *Size.* `W = sum w_j <= r (tau + eta)`.
   * The window lemma puts `P <= 2P_**`. For `Q_0 = 2M`, `P_** = P_*` once `M` is large, since
     `160 M^2 L >> 8.4 r^2 log P`.
   * For `Q_0 = 2M`, the prime number theorem in the progression gives `L = (2 + o(1)) M log M`, so
     `log P_* = 3 log M + loglog M + O(1) + log(1/C)` and `r = 9M - 8`.
   * The other cases are the same computation.
   * In `A_0`, `P_** = x_1` for large `M`, so `log P_** = log(8.4 r^2 log x_1) = 2 log M + loglog M + O(1)`.

   QED

**Exact instances** (`instance.py`, the arguments `q b C Q0/M P0exp nmax`).

| `M` | `C` | `r` | primes | `W` | `W/(M log M)` | union bound `P(fail) <=` | `D_1` max (`n <= nmax`) | illustration: min log-margin of (R) |
|---|---|---|---|---|---|---|---|---|
| 12 | 1 | 100 | `10000121..10003529` | 1611.8 | 54.1 | `0.067` | 177905 (`n <= 6`) | 21.1 |
| 20 | 1 | 172 | `100000037..100005877` | 3168.4 | 52.9 | `0.034` | 25415 (`n <= 4`) | 21.0 |
| 24 | 1/2 | 208 | `100000037..100006981` | 3831.5 | 50.2 | `0.121` | 34017 (`n <= 4`) | 17.8 |
| 44 | 1 | 388 | `1000000009..1000016273` | 8040.6 | 48.3 | `0.040` | 247 (`n <= 2`) | 23.6 |
| 12 | 1/100 | 100 | `1000000009..1000004381` | 2072.3 | 69.5 | `0.065` | 3795 (`n <= 4`) | 23.4 |

The rows for `M = 44` and for `C = 1/100` are in `instance_q43_n2_out.txt` and `instance_q11_C001_out.txt`; the
others are in `instance_out.txt`. The last block of `instance_out.txt` is a longer `M = 44` run (`nmax = 4`). It was
stopped for time after the residue-data check and is not used.

For these instances the following is **exact**:
* the prime list, window and `p_j >= 21 r^2`;
* rank and twins;
* the conic-fibre condition of the residue data;
* every `D_u(c)` for `||c||_1 <= nmax`, the divisibility `D_1 | Im Lambda_c`, and `D_u <= 2L^(n/2)`;
* the union-bound sum, evaluated with 50-digit logarithms of the actual primes and the proved `G`-bounds.

The union-bound sum runs to `n = 4000`. Its tail is bounded by a geometric series with ratio `<= 0.008`.

A union bound `< 1` proves that a Haar-positive set of `delta` works for that concrete `(M, primes, residues)`.
Together with `1.33(Pi - 1) < 1`, which is printed, Theorem A.1(ii) applies, so these five fake circles exist.

Two things are illustrations only:
* the last column: one random `delta`, all characters with `||c||_1 <= nmax`;
* the check of 100,000 random monomials outside `L` on a generic realisation, with minimum ratio
  `dist/required` between 3.2 and 36.

**Growth of the constant** (`rates.py`, from the proved bounds, `b = 8`, `C = 1`; floating point in logs):

| `M` | `tau/log M` (A_1) | `W/(M log M)`, `A_1(2M)` | `W/(M log M)`, `A_0` |
|---|---|---|---|
| 12 | 5.72 | 50.0 | 49.2 |
| 104 | 4.53 | 41.7 | 37.3 |
| 10,004 | 3.82 | 35.1 | 28.3 |
| 1,000,004 | 3.59 | 32.8 | 25.1 |

The limits are `27` and `18`, approached at the rate `loglog M / log M`.

---------------------------------------------------------------------------------------------

## 6. Consequences

**Corollary 6.1.** Fix `C > 0`. Write `W_min^(A_1)(M, C)` for the least `W` over models in `A_1 = A_1(2M)` with an
`M`-point cluster on an arc of length `<= C sqrt R`, and `W_min^(A_0)` for the same without residue data.
* Along `M = q + 1`, `q = 3 mod 4` prime:
  `(3 - o(1)) M log M <= W_min^(A_1) <= (27 + o(1)) M log M + 9M log^+(1/C)`.
* `(kappa_C - o(1)) M log M <= W_min^(A_0) <= (18 + o(1)) M log M + 9M log^+(1/C)`, with
  `kappa_C = 1/ceil(C/sqrt2)` from `outside.md` Prop A.4.
* For all `M`, `W_min = O_C(M log M)`. By Breusch's theorem there is a prime `q = 3 mod 4` with
  `M <= q + 1 <= 2M`, and Lemma 1.3 applies.

So `W_min / (M log M)` stays bounded, and the alternative "`W_min/(M log M) -> infinity`" in the task is false.

**Corollary 6.2 (the class is sharply capped).** Consider an argument that proves `M <= f(R)` for all
clusters on arcs of length `<= C sqrt R`, and uses about the configuration only:
* (I) its exponent profile and exact exponent identities;
* (II) the Lemma 0.1 gaps for all monomials;
* (III) the residue axiom (R) at odd prime powers `<= 2M`, or `<= M^K`;
* (IV) properties of the angles valid off Haar-null sets.

Then `f(R) >= (2/27 - o(1)) log R/loglog R` for infinitely many `R`, or `(2/(9(2+K)) - o(1))` with `M^K`.

*Proof.* The fakes of Theorem 2 have `log R = W/2 <= (13.5 + o(1)) M log M`, and `x/log x` is increasing. QED

The class contains:
* every input of `two_thirds/proof.md`;
* pair, integral and rational character certificates;
* Theorem F (residue pigeonhole on windings);
* the aggregate collision and capacity inequalities of Prop P;
* every Fourier/Bohr/inverse-LO/energy argument on the angle GAP whose arithmetic input is (I)-(IV) (Theorem A.1).

So the best constant for the class lies in `[2/27, 2/3]`. **No growth improvement `M = o(log R/loglog R)` can come
from it.** In particular:
* (P1) of `outside.md` is settled affirmatively, with `b = 8`.
  * Without residues this was already implicit in items 392 and 394 together with Theorem A.1.
  * With residues it is new here, and independently in `lower.md`.
  * No `l^1`-stability is needed, but the flip assignment must be a skew pairing; `lower.md` instead proves an
    `l^1`-stability lemma valid for every capacity-two assignment.
* (P2) has a negative answer for skew-paired Paley profiles: a set of positions `delta` of positive measure,
  with `k = 0`, makes every character simultaneously admissible at `W = O(M log M)`.
* The "fully flipped Hadamard phase lemma" (`walsh-extraction.md` §10), if true for actual primes, needs an input
  outside (I)-(IV) and (R).
* Corollary A.4's ceiling `(log R/loglog R)^(1/2)` improves to the sharp order `log R/loglog R`, now also with
  residue strengthening.

The referee's per-character product-pigeonhole moduli `<= M^k/k!` (referee §3.2) are covered for every fixed `k`
by taking `Q_0 = M^k`.

**Consistency with refereed results** (no contradiction found).
* *`two_thirds`*: `3 M log M <= W + K_C M` holds, with `W ~ 27 M log M`.
* *Theorems B and B''* (`walsh-extraction.md`) need pure cores, or modifications confined to
  `<= M/(40 log M)` rows. Here every row is flipped.
* *Theorem C* excludes flipped Walsh profiles. Walsh cores of order `>= 16` have no skew pairing (Remark 2.3).
* *Prop P* is the individual version of Theorem 2, for the same profile class.
* *Theorem L of `round4/flipped.md`* is conditional on Lang-Waldschmidt. It excludes Hadamard-core profiles with
  flip mass `< (1/2 - delta) W` **for the actual angles of Gaussian primes**. Here the flip mass is `W/9`.
  So a two-logarithm transcendence input would exclude genuine versions of these fakes. That is consistent:
  the fakes' angles are not the actual ones.
* *`linear_allocation_affine_rigidity.md`*: `M <= r + 1 = 9M - 7`.
* *`outside.md` Prop A.4*: `W >= (1 - o(1)) M log M` for `C <= sqrt 2`, and here `W ~ 18-27 M log M`.

---------------------------------------------------------------------------------------------

## 7. What is not covered

1. **Prime-level residue data.** An actual circle also carries the residues `r_j = pi_j mod q^a` of the
   individual Gaussian primes, with `r_j conj(r_j) = p_j`. The point residues are the induced products.
   * With twins the induced map is onto, up to one parity per `q`: `rho_x / rho_y` must lie in a prescribed
     coset of the squares in the norm-one group `K_(q^a)`. For `l = N(omega)`, the coset of
     `omega/conj(omega)` is the Legendre symbol `(l/q)`, for inert and split `q` alike.
   * The design of Lemma 4.1 respects these cosets only if the `l_x` have prescribed quadratic characters modulo
     every `q <= Q_0`. That condition has density `2^(-pi(Q_0))`, so heuristically it puts `l_x` near
     `2^(pi(Q_0))`, which is useless here.
   * Whether `W = O(M log M)` holds with prime-level residue data is **open**. Two routes look possible, and
     neither is proved here:
     * a random prime-level assignment, which heuristically gives `O(M log^2 M/loglog M)`;
     * restricting every `p_j` to be a quadratic residue modulo every odd `q <= Q_0` (pseudosquare-type primes),
       with point residues `s (omega_x/conj omega_x)^2`. This makes all parities trivial at the cost of
       `log p_j >= pi(Q_0) log 2`, and so would give `W = O(M^2/log M)` if such primes exist at the conjectured
       size.
   * This parity (quadratic-reciprocity) information is the first arithmetic input beyond the class.
2. **Residue strengthening for monomials outside `L`.** `outside.md` A.3' has this, at `(log R)^(1/3)`. The point
   design does not define residues of single primes; see item 1.
3. **The actual angles.** Transcendence measures (Theorem L above), and residues at all moduli (`outside.md`
   §2.4).
4. **Exact algebraic identities** (Ptolemy, auxiliary polynomials) and **multi-place approximation**. Both are
   barred by earlier rounds.

---------------------------------------------------------------------------------------------

## 8. Files (`round5/checks/`)

| script | what it checks | type |
|---|---|---|
| `lib.py`, `gauss5.py` | Paley/Sylvester matrices, profile builder, Gaussian integers | helpers |
| `check_G_bound.py` | (L) skew minors, `q <= 199`; (E) Lemma 2.2; (T1) Theorem 1 (i)-(ii), exhaustive on `+-1` characters for `M = 12..44`; (T2) Theorem 1 (iii), exhaustive small norms | exact |
| `check_G_large.py` | Theorem 1 on `M = 44, 104, 200`: all pairs, columns, half-sums, random and hill-climbed characters; skew pairing against random assignment | exact (evidence) |
| `sylvester_skew.py`, `sylvester_core_skew.py`, `sylvester_core_skew_relaxed.py` | Remark 2.3: no skew pairing for Walsh cores of order `>= 16` | exact, exhaustive |
| `check_axiomR.py` | Lemma 1.2: validity of (R) on actual points; Prop 1.1 example | exact |
| `instance.py` | Theorem 2 instances `M = 12, 20, 24, 44`: primes, window, rank and twins, residue data, `D_u(c)` for small `c`, union bound | exact; one-`delta` check and monomial check are illustrations |
| `rates.py` | `W/(M log M)` from the proved bounds, `M <= 10^6` | arithmetic in floating-point logs |
| `milp_minG.py` | early exploration (random assignment, MILP). Not used in any proof | evidence |

The files `check_canonical.py`, `check_fake.py`, `check_lemmaS.py`, `check_medium.py`, `paley.py` and
`unionbound.py`, with their outputs, belong to the concurrent `round5/lower.md`, not to this write-up. So does
`round5/exact_checks/`.

The repository `/home/user/Jarnik` was not modified.
