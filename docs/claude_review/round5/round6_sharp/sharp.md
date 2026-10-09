# Round 6: fake circles for residue data that factor through the primes (P3)

Conventions are those of `round4/outside.md` and `round5/lower.md`:
* `W = log N`, `R = N^(1/2)`, and the arc has angular width `Delta = C e^(-W/4)`;
* `w(v) = sum_j |v_j| log p_j`, `L` is the character lattice, `v(c) = sum_x c_x a_x`, `n = ||c||_1`;
* the slack of a character is `s_c = w(v(c)) - W/2`.

*Theorem*, *Lemma* and *Proposition* mean proved here with all quantifiers, or cited exactly and
re-derived here. Results taken from the draft `round5/lower.md` (Lemma S, Prop. 3.3, Lemma 5.3,
Lemma 5.2, Prop. 1.2/1.3) are used only after the re-derivation in §1 and §4. *Evidence* means a
computation that no proof uses. Scripts are in `round6/sharp_checks/`; each runs from that
directory with `python3 <script>`, and its output is stored next to it as `*_out.txt`.

The task changed twice while this was written. The original target (the class III-small, built on
a Sylvester multi-block profile with Green-Sanders) is answered by Theorem 5.1 with `A = 1` and is
summarised in Appendix A. The main body treats the refined question (P3 of `lower.md`): residue data
that **factor through residues of the individual primes**, with the genuine norm constraint,
separately for characters (P3-char) and for all monomials (P3-all), and the question of moduli up
to `N^A`.

---------------------------------------------------------------------------------------------

## 0. Verdict

| class | residue strengthening | proved here | lower bound |
|---|---|---|---|
| III-small (task) `⊇` P3-char at `X = 3M` | characters, prime powers `< 3M` | `W_min <= (66 + o(1)) M log M` (Thm 5.1, `A = 1`) | `3 M log M - O(M)` (`two_thirds`) |
| **P3-char**, moduli `<= 3M^A` | characters, factored residues with genuine norms | `W_min <= (33 max(A,2) + o(1)) M log M` (Thm 5.1) | same |
| **P3-all**, moduli `<= 3M^A` | every monomial `A_v`, `v in Z^r` | `W_min <= (66 A + o(1)) M log^2 M / loglog M` (Thm 6.4) | same |
| any of these, moduli up to `N^A` | | not proved; where it breaks: §7 | same |

1. **P3-char is sharp at `Theta(M log M)`.** For every fixed `A >= 1` and `C > 0` there are angle
   models whose residue data factor through prime residues `r_j(q)` with `Norm r_j(q) = p_j`, at
   every odd prime power `<= 3M^A`, with every character gap strengthened by its actual collision
   modulus, and with `W <= (33 max(A,2) + o(1)) M log M` (Theorem 5.1). So no argument with these
   inputs proves `M = o(log R/loglog R)` (Corollary 5.2). This contains the task's class III-small
   and answers its parts (a) and (b), with `c_0 = 66 + o(1)`.
2. **What factorisation forces is exactly one bit per prime and modulus** (Prop. 2.1, Cor. 2.2).
   The norm-one group `T_Q` of `Z[i]/q^a` is cyclic of order `m = (q - chi(q)) q^(a-1)`. The class
   `ell_j` of `r_j/conj(r_j)` satisfies `ell_j = nu_j(q) mod 2`, where `nu_j(q) = 1` iff `(p_j/q) = -1`,
   and it is free otherwise:
   * `ell_j mod 4` is free;
   * prime powers add nothing;
   * the unit `i` flips parity iff `q = +-3 mod 8`, and `-1` never does.

   With twins, a point class map is realisable iff it respects the parity colouring
   `sigma_q(x) = sum_j a_xj nu_j(q) + [q = +-3 mod 8] k_x mod 2`. A parity-blind design is
   unrealisable at every one of the 535 moduli of the exact test (evidence).
3. **The design** (Def. 3.2) puts classes inside parity classes: `x -> (sigma_q(x), x mod p'(q))`.
   Here `p'(q) <= m_q/2` comes from a dyadic prime assignment of multiplicity `<= 9` (Lemma 3.1).
   * Max pair demand: `<= (10 + o(1)) log M` (measured `<= 6.8 log M` up to `M = 16384`).
   * Above `4M`, the maps are random parity-respecting injections.
   * Realisability through the exponent matrix is explicit (Cor. 2.2), including prime powers,
     units and level compatibility. It is checked exactly end to end on a Paley `M = 44`, `b = 33`
     profile with 1419 actual primes and all 535 odd primes `q <= 3872` together with their powers.
4. **The parity costs a constant, not the rate.** The coherent auxiliary-prime design of
   `lower.md` Lemma 4.1 (constant 48) is not available for factored data. Coherent objects would
   have to reproduce the Legendre pattern of every `p_j` at every `q <= X`, and only the genuine
   `pi_j` (up to squares) are known to do so. We therefore need uniform slack for every non-pair
   character, which `lower.md` Lemma S provides for `b >= 33` copies. The constant is
   `33 max(A,2)` (66 for `A <= 2`).
5. **P3-all.** All monomials strengthened at moduli `<= 3M^A`.
   * *Proved:* `W <= (66 A + o(1)) M log^2 M / loglog M` (Theorem 6.4). This gives the ceiling
     `M >= (1/(33A) - o(1)) log R logloglog R / (loglog R)^2`.
   * *Not `O(M log M)`.* With per-modulus independent residue choices, the necessary **height
     condition** `m_v <= sqrt2 e^(w(v)/2)` (Prop. 6.1) fails w.h.p. below
     `log P = (2 - eps) log^2 M/loglog M` (Prop. 6.5). So Theorem 6.4 is sharp for that method.
   * *Coherent choices.* Parity forces them to use the genuine `pi_j`, which makes pair
     collisions genuine-random; heuristically this gives the same rate (§6.6).
   * *Correction to the coordinator's estimate.* "`E m_v ~ 2^(pi(2M))`, hence `log P >~ M/log M`"
     is too pessimistic. That expectation is dominated by events of probability `e^(-Theta(M))`.
     Theta-moments give `log^2 M/loglog M`.
   * *Open (P3-all\*).* Is `W_min = O(M log M)` in P3-all? A negative answer would be a growth
     improvement from purely local inputs. No lower-bound mechanism beyond `two_thirds` is known;
     the question is a combinatorial design problem (§6.7).
6. **Moduli up to `N^A`: not proved, for free or factored residues.** It breaks at moduli of
   size `e^(Theta(M))` and beyond.
   * A fixed labelling collides everywhere (Siegel trap, Prop. 7.1, rigorous).
   * Coherent designs need auxiliary objects that are units modulo every `q <= X`. Their demand
     bounds are `>= 2n log X`. For `X >= N^(1/4)` that exceeds every possible slack
     `s_c <= (||c||_1 - 1) W/2` (§7.2).
   * Independent designs have Poisson tails: some characters with `||c||_1 ~ M` collide at about
     `M/log M` moduli of size `e^(Theta(M))` (expected-count computation, heuristic, §7.3).
   * A "balanced residue design" would suffice (Prop. 7.4), but none is known.
7. **The Sylvester multi-block route** of the original lead also works for III-small and for
   P3-char, with astronomically large constants (Green-Sanders). Appendix A gives its key lemmas
   and status. It is superseded by the Paley route.

---------------------------------------------------------------------------------------------

## 1. The classes

### 1.1 Definition

Fix `M`, `C > 0` and `X >= 3`. For an odd prime `q`, put `A_q = A_q(X) = max{a >= 0 : q^a <= X}`.
Write `T_Q` for the norm-one subgroup of `G_Q = (Z[i]/Q)^*`, and `m_Q = |T_Q|` (Prop. 2.1).

**Definition 1.1 (factored residue data; classes `F^char_X(M, C)` and `F^all_X(M, C)`).** A member
is a tuple `(p, e, (a_x), phi, theta_0, delta, k, r)` with the following properties.

* (a) *Profile.* As in `lower.md` Def. 1.1(a): distinct primes `p_j = 1 mod 4`; exponents `e_j >= 1`;
  pairwise distinct rows `a_x in prod_j [0, e_j]`; and `min_x a_xj = 0`, `max_x a_xj = e_j`.
  Put `N = prod p_j^(e_j)` and `W = log N`.
* (b) *Angles.* `phi in R^r` and `delta in [0, Delta]^M` with `Delta = C e^(-W/4)`, and
  `k in Z^M`, satisfying (1.1) of `outside.md`:
  `<2a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x`.
* (c) *(G).* Every `v != 0` satisfies both gap inequalities of Lemma 0.1.
* (d) *Factored residue data.* For every odd prime `q <= X` and every `j` with `p_j != q`: a unit
  `r_j(q)` of `Z[i]/q^(A_q)` with `Norm r_j(q) = p_j mod q^(A_q)`. Lower levels are reductions,
  so they are automatically compatible.
  * For `q` not dividing `N`, the *point residues* are
    `rho_x(q) = i^(-k_x) prod_j r_j(q)^(a_xj) conj(r_j(q))^(e_j - a_xj)`.
  * For `q = p_l`, the *level residues* are the same product over `j != l`.
* (e-char) For every `c in Z^M` with `sum c = 0` and `v(c) != 0`, and every unit `eta in {1, i, -1, -i}`:

  ```text
  |sin((c.delta - arg eta)/2)|  >=  2^(-1/2) m_(c,eta) e^(-w(v(c))/2),                        (1.1)
  m_(c,eta) = prod_(q <= X, q not dividing N) q^(a_(q,eta)(c)),
  a_(q,eta)(c) = max{a <= A_q : prod_x rho_x(q)^(c_x) = eta mod q^a}.
  ```

  Add the same-level pair clause at the primes `p_l <= X`, exactly as `lower.md` Def. 1.1(e).
* (e-all) In addition, for every `v in Z^r \ {0}` and every `s in Z/4`:

  ```text
  |sin(<v, phi> - s pi/4)|  >=  nu_s m_(v,s) e^(-w(v)/2),     nu_s = 1 (s even),  nu_s = 2^(-1/2) (s odd),   (1.2)
  m_(v,s) = prod_(q <= X, q not dividing N) q^(a_(q,s)(v)),
  a_(q,s)(v) = max{a <= A_q : prod_j (r_j/conj r_j)^(v_j) = i^s mod q^a}.
  ```

`F^char_X` consists of (a)-(d) with (e-char). `F^all_X` is `F^char_X` with (e-all) added. We write
**P3-char** and **P3-all** for these classes, and `W_min^char`, `W_min^all` for the infima of `W`.

**Relations.**
* *To `lower.md`.* `rho_x(q)` has norm `N mod q^(A_q)`. So a member of `F^char_X` induces a member
  of `lower.md`'s class `A_X` (point residues, free up to the norm). Condition (e) there is (1.1)
  with `eta = 1`. Hence `F^all_X ⊂ F^char_X ⊂ A_X`.
* *To the task's III-small.* That class has residue maps realisable through the exponent matrix,
  strengthened by the actual collision moduli at prime powers with `(q -+ 1) q^(a-1) < 2M`.
  Those prime powers satisfy `q^a < 3M`. So `F^char_(3M)` is a subclass of III-small, and a fake
  circle in `F^char_(3M)` is one in III-small.

### 1.2 Validity and the 2/3 bound inside the classes

**Proposition 1.2 (actual clusters are members).** Let `z_1, ..., z_M` be a normalised actual
cluster on an arc of length `<= C sqrt R` (`two_thirds` Lemma 1.1), with `phi_j = arg pi_j`. Define
`delta`, `k` by (1.1), and put `r_j(q) = pi_j mod q^(A_q)`. Then the tuple lies in `F^all_X(M, C)`
for every `X`.

*Proof.*
* (a)-(c) are `outside.md` §1 and Lemma 0.1. (d) holds since `Norm pi_j = p_j`.
* *Units.* From (1.1), `z_x = eps_0 i^(-k_x) prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj)` for a
  common unit `eps_0`. So `rho_x(q) = conj(eps_0) z_x mod q^(A_q)`, and a common unit cancels in
  `prod rho_x^(c_x)`, because `sum c = 0`.
* *(e-char).* This is `lower.md` Prop. 1.2, with `U = eta V` in place of `U = V`.
  * Put `U = prod_(c_x>0) z_x^(c_x)` and `V = prod_(c_x<0) z_x^(-c_x)`, and `G = gcd(U, V)`. Then
    `U/G = u_1 A_v` and `V/G = u_2 conj(A_v)`.
  * A collision `prod z_x^(c_x) = eta mod q^a` means `U = eta V mod q^a`. Then
    `q^a | D_eta := U/G - eta V/G`.
  * `D_eta != 0`, since `A_v` is conjugate-primitive. `(1+i) | D_eta`, since odd-norm elements and
    units are `1 mod (1+i)`. Distinct `q` are coprime. Hence `|D_eta| >= sqrt2 m_(c,eta)`.
  * `|U - eta V| = |G| |D_eta|`, with `|G| = N^(n/4) e^(-w/2)`. Also
    `|U - eta V| = N^(n/4) |e^(i c.delta) - eta| = 2 N^(n/4) |sin((c.delta - arg eta)/2)|`.
  * The same-level clause is `two_thirds` Lemma 3.3 (S|).
* *(e-all).* Let `A_v = prod pi_j^(v_j^+) conj(pi_j)^(v_j^-) = X_v + i Y_v`.
  * A collision `prod (r_j/conj r_j)^(v_j) = i^s mod q^a` means `A_v = i^s conj(A_v) mod q^a`
    (`A_v` is a unit mod `q`).
  * `A - i^s conj A` equals `2iY`, `(X - Y)(1 - i)`, `2X` or `(X + Y)(1 + i)` for `s = 0, 1, 2, 3`.
    So `q^a` divides the rational integer `D_s = Y, X - Y, X, X + Y`.
  * `D_s != 0` (conjugate-primitivity), and coprime moduli multiply, so `|D_s| >= m_(v,s)`.
  * Finally `|D_s| = nu_s^(-1) |A_v| |sin(arg A_v - s pi/4)|`, with `arg A_v = <v, phi>` and
    `|A_v| = e^(w(v)/2)`. QED

**Proposition 1.3 (two_thirds inside).** For `X >= 2M`, every member of `F^char_X(M, C)` satisfies
`3 M log M <= W + (14 + 4 log^+ C) M`.

*Proof.* `F^char_X ⊂ A_X`, and this is `lower.md` Prop. 1.3. That proposition was re-checked here
against `two_thirds/proof.md`, which uses the configuration only in two places:
* (2.5), i.e. the pair inequality (1.1) with `eta = 1` summed against the cut identity;
* (3.6), the pigeonhole with the class counts of Lemma 3.2. Those counts depend only on the norm
  of the point residues, `N mod q^a`, and on collisions at `q^a < 1.5 M`.

Factored data have point residues of norm `N`, so both hold. Units are irrelevant (`two_thirds`
§3, end). QED

**Inputs of `two_thirds` and where they sit.**

| input | where it sits |
|---|---|
| normalisation | (a) |
| pair chord identity and cut identity | exponent identities of the profile |
| collisions at inert `q`, at split `q` not dividing `N`, and same-level collisions at `p | N` | (1.1) for pairs, with the class counts |
| ramified factor `sqrt 2` | the `2^(-1/2)` in (1.1) |
| chain lemma, prime sums | pure combinatorics |

Theorem F of `walsh-extraction.md` is (G) on `L` together with the integral winding data. On
saturated profiles (twins) every rational character is integral (`lower.md` Lemma 3.1), so
Theorem F adds nothing.

---------------------------------------------------------------------------------------------

## 2. What factored residue data force

**Proposition 2.1.** Let `q` be an odd prime, `a >= 1`, `Q = q^a`, and `m = m_Q = (q - chi(q)) q^(a-1)`.
1. `T_Q` is cyclic of order `m`, and `4 | m`.
2. `r -> r/conj(r)` is a homomorphism of `G_Q` onto `T_Q`. Its fibres are the cosets
   `r (Z/Q)^*`.
3. There is a unique homomorphism `nu_Q : T_Q -> {+-1}` with `nu_Q(r/conj r) = (Norm r / q)`.
   It is onto, and its kernel is `T_Q^2`.
4. For `t in T_Q` and `n in (Z/Q)^*`: some `r` has `r/conj(r) = t` and `Norm r = n` iff
   `nu_Q(t) = (n/q)`.
5. `i in T_Q^2` iff `q = +-1 mod 8`. Always `-1 in T_Q^2`.
6. Reduction `G_(q^(a+1)) -> G_(q^a)` maps `T` onto `T` and generators to generators, and
   `nu_(q^(a+1)) = nu_(q^a) o red`.

*Proof.*
* *Item 2.* `r = conj(r) mod Q` iff `2i Im r = 0` iff `r in Z/Q`. So the kernel of the
  homomorphism `r -> r/conj(r)` is `(Z/Q)^*`, and the image lies in `T_Q`. `Norm : G_Q -> (Z/Q)^*`
  is onto, with fibres of size `m` (`two_thirds` Lemma 3.2). So
  `|G_Q| = m |(Z/Q)^*|`, and the image has `m = |T_Q|` elements.
* *Item 1, split `q`.* `Z[i]/q^a ≅ Z/q^a x Z/q^a`, with conjugation swapping the factors. Then
  `T ≅ {(u, u^(-1))} ≅ (Z/q^a)^*`, which is cyclic.
* *Item 1, inert `q`.* `T` is the product of the norm-one part of `F_(q^2)^*` (cyclic of order
  `q+1`) and the norm-one part of `(1 + qO)/(1 + q^a O)`, `O = Z_q[i]`. The latter is cyclic of
  order `q^(a-1)`: the `q`-adic logarithm identifies it with the trace-zero part of `qO/q^a O`.
  The orders are coprime, so `T` is cyclic.
* *Item 1, `4 | m`.* `q - chi(q) = 0 mod 4`.
* *Item 3.* `nu` is well defined: `r' = ru` with `u` rational changes the norm by `u^2`. It is
  multiplicative, and onto because `Norm` is onto. A cyclic group of even order has exactly one
  subgroup of index 2, namely the squares.
* *Item 4.* Take `r_0` with `r_0/conj(r_0) = t`. Then `n/Norm(r_0)` is a unit square mod `q`,
  hence mod `q^a` (Hensel, `q` odd): `n = s^2 Norm(r_0)`. Put `r = s r_0`.
* *Item 5.* `i = (1+i)/conj(1+i)`, and `Norm(1+i) = 2`, so `nu(i) = (2/q)`. Also `-1 = i^2`.
* *Item 6.* The reduction is a surjective homomorphism of cyclic groups (item 2 at both levels).
  Legendre symbols only see the residue mod `q`. QED

**Corollary 2.2 (factored data = parity-respecting class maps).** Fix `X`, an odd prime `q <= X`
not dividing `N`, its top power `Q = q^(A_q)`, and a generator `g` of `T_Q`. Write
`r_j/conj(r_j) = g^(ell_j)`, `i = g^iota`, and `nu_j(q) = [(p_j/q) = -1] = [(q/p_j) = -1]`
(quadratic reciprocity, `p_j = 1 mod 4`).
1. Data as in (d) correspond exactly to vectors `ell in (Z/m_Q)^r` with `ell_j = nu_j(q) mod 2`.
   Nothing else is constrained.
2. The point residues are `rho_x = B g^(kappa_x)` with `B` common and
   `kappa_x = <a_x, ell> - iota k_x`.
3. `c` collides at `q^a` with unit `i^s` iff `sum_x c_x kappa_x = iota s mod m_(q^a)`.
4. Every `kappa` arising this way satisfies `kappa_x - kappa_y = sigma_q(x) - sigma_q(y) mod 2`, with

   ```text
   sigma_q(x) = sum_j a_xj nu_j(q) + lambda_q k_x  mod 2,     lambda_q = [q = +-3 mod 8].       (2.1)
   ```

5. *Twins.* Suppose the profile has, for each row `x`, columns `f(x)`, `s(x)` with
   `a_(.,f(x)) - a_(.,s(x)) = eps_x e_x`, `eps_x = +-1`, all `2M` columns distinct. Then
   conversely every `kappa` with `kappa - sigma_q` constant mod 2 arises. An explicit preimage is

   ```text
   ell = nu + 2 ell',      ell' = sum_x y_x eps_x (e_(f(x)) - e_(s(x))),      y = (kappa - kappa_0 + iota k - A nu)/2.   (2.2)
   ```

*Proof.*
* Item 1 is Prop. 2.1(4) per `j`.
* Item 2: `rho_x = i^(-k_x) prod_j conj(r_j)^(e_j) prod_j (r_j/conj r_j)^(a_xj)`.
* Item 3 follows from item 2, using `sum c = 0`.
* Item 4: `iota` is odd iff `i` is a non-square (Prop. 2.1(5)), i.e. iff `lambda_q = 1`.
* Item 5: `A ell' = sum_x y_x eps_x^2 e_x = y`. So `A ell = A nu + 2y = kappa - kappa_0 + iota k`.
  Here `y` is defined since `kappa - kappa_0 + iota k - A nu` is even.
* The `r_j` then exist by Prop. 2.1(4), and lower levels follow by reduction (Prop. 2.1(6)). QED

**Remarks 2.3.**
1. *Units.* `-1` never changes parity. `i` changes it iff `q = +-3 mod 8`. Windings with odd `k_x`
   therefore move row `x` between parity classes at exactly those `q`. The construction below uses
   `k = 0`.
2. *The quartic part is free.* `ell_j mod 4` is not constrained by the norm. For actual Gaussian
   primes and `8 | m`, `ell_j mod 4` is an invariant of `pi_j` (the quartic character). It is not a
   function of `(p_j/q)`: both values occur among primes with `(p_j/q) = 1` (evidence,
   `check_parity.py` (f)).
3. *Prime powers* add nothing beyond level compatibility (Prop. 2.1(6)).
4. *Norms versus full residues.* Suppose a class fixed the residues `pi_j mod q^a` themselves.
   * Once `prod_(q <= X) q > 2 sqrt(p_j)`, these residues determine `pi_j`. Two Gaussian integers
     of norm `p_j` congruent modulo an integer `Q > 2 sqrt(p_j)` are equal, since a nonzero
     Gaussian multiple of `Q` has modulus `>= Q`.
   * For `p_j <= M^(O(1))` and `X >= 3M` this always happens. So "full prime residues" means
     "genuine Gaussian primes `pi_j`, fake angles `phi`". That is a different, essentially global,
     class (§6.6).
   * The local content of the genuine primes is their norms, i.e. exactly the parity of item 1.
     That is the class studied here. It is P3 of `lower.md` §8, whose Prop. 8.1 is the case
     `a = 1`, `q = +-1 mod 8` of Prop. 2.1(3).

*Exact checks* (`check_parity.py`, all 89 odd prime powers `<= 400`, exhaustive over `Z[i]/Q`, 211
Gaussian primes `p <= 3000`): cyclicity; the fibres in item 2; `nu` well defined and equal to
square-ness; `i` square iff `q = +-1 mod 8`; the realisability criterion of item 4 for all `t` and
40 values of `n` per modulus; `pi/conj(pi)` square iff `(p/q) = 1`; parity compatibility across
levels. ALL PARITY CHECKS PASSED.

---------------------------------------------------------------------------------------------

## 3. The parity-respecting residue design

### 3.1 A prime assignment

**Lemma 3.1.** Define a map `q -> p'(q)` from odd primes to primes.
* `p'(3) = p'(5) = 2` and `p'(7) = 3`.
* For `k >= 3`, list the odd primes `q > 7` in `[2^k, 2^(k+1))` as `q_0 < ... < q_(n-1)`, and the
  primes in `[2^(k-2), 2^(k-1))` as `t_0 < ... < t_(t-1)`. Put `p'(q_i) = t_(floor(i t/n))`.

Then:
1. `2 p'(q) <= m_q` and `q < 8 p'(q)` for every odd prime `q`, and `p'(q) <= (q-1)/2` for `q >= 5`;
2. every prime has at most 9 preimages;
3. for every integer `d >= 1`,
   `Lambda'(d) := sum_(q : p'(q) | d) log q <= 9 log d + 18.72 omega(d)`.

*Proof.*
1. `p'(q) < 2^(k-1) <= q/2` with `q` odd, so `p'(q) <= (q-1)/2 <= m_q/2`; and `q < 2^(k+1) <= 8 p'(q)`.
   The three small cases are checked directly: `m_3 = m_5 = 4 = 2 * 2` and `m_7 = 8 >= 2 * 3`.
2. A target in block `k` receives at most `ceil(n/t)` primes. Target blocks of different `k` are
   disjoint dyadic intervals, and the special cases add at most 2 to the primes 2 and 3.
   * *`k <= 22` (`q < 2^23`).* Exact count (`check_design.py`, Part 1): the multiplicity is at
     most 4.
   * *`k >= 20`.* Use the Rosser-Schoenfeld bounds `x/ln x < pi(x)` (`x >= 17`) and
     `pi(x) < 1.25506 x/ln x` (`x > 1`). With `n <= pi(2^(k+1)) - pi(2^k) + 1` and
     `t >= pi(2^(k-1)) - pi(2^(k-2)) - 1`, dividing by `2^k/ln 2` gives

     ```text
     n/t <= R(k) := (2.51012/(k+1) - 1/k + ln2 2^(-k)) / (0.5/(k-1) - 0.313765/(k-2) - ln2 2^(-k)).
     ```

     Direct evaluation for `20 <= k <= 2000` gives `max R(k) = 8.1047`. For `k > 2000`, the
     numerator is `<= 1.51013/k` and the denominator is `>= (0.186235 - 1/k - 10^(-500))/k`. So
     `R(k) <= 8.14` there. Hence the multiplicity is at most 9.
3. `Lambda'(d) = sum_(p | d) sum_(p'(q) = p) log q <= sum_(p | d) 9 log(8p)`. QED

By Robin's bound `omega(d) <= 1.3841 log d / loglog d` (`d >= 3`), item 3 gives
`Lambda'(d) <= (9 + 26/loglog d) log d`.

### 3.2 The design

**Definition 3.2 (design `D_A`).** The inputs are the following.
* A profile with twins as in Cor. 2.2(5).
* Rows labelled `0, ..., M-1`.
* Primes `p_j > X = 3 M^A` with `A >= 1`, and windings `k = 0`.

For each odd prime `q <= X`, with top power `Q = q^(A_q)`, put
`sigma_q(x) = sum_j a_xj nu_j(q) mod 2` (this is (2.1) with `k = 0`) and define
`kappa^(q^a) in Z/m_(q^a)` for `1 <= a <= A_q` as follows.

* **Small `q`** (`q < 4M`). The level-1 class is

  ```text
  kappa^(q)_x = sigma_q(x) + 2 (x mod p'(q))   in Z/m_q.
  ```

  Put `a*(q) = min{a >= 1 : p'(q) q^(a-1) >= M}`.
  * For `2 <= a <= a*(q)`: `kappa^(q^a)_x` is the element of
    `Z/m_(q^a) ≅ Z/m_q x Z/q^(a-1)` with components `(kappa^(q)_x, x mod q^(a-1))`.
  * For `a > a*(q)`: a random lift `kappa^(q^a) = kappa^(q^(a-1)) + m_(q^(a-1)) u^(q,a)`, with
    `u^(q,a)` uniform in `(Z/q)^M`.
* **Large `q`** (`4M < q <= X`). `kappa^(q)_x = sigma_q(x) + 2 y^(q)_x`, where `y^(q)` is a
  uniformly random injection of the rows into `Z/(m_q/2)` (possible, since `m_q/2 >= 2M`). Higher
  levels are random lifts as above.

All random choices are independent. By construction `kappa^(Q) = sigma_q mod 2`, and every level is
the reduction of the next one. Cor. 2.2(5) (with `kappa_0 = 0`), applied at the top power `Q`, gives
residues `r_j(q) mod Q` with `Norm r_j(q) = p_j`. Their point classes are `kappa^(q^a)` at every level.

**Lemma 3.3 (properties of `D_A`).** Let `d = |x - y|`.
1. *Pairs, unit 1.* `rho_x = rho_y mod q^a` iff `q` is small, `a <= a*(q)`,
   `sigma_q(x) = sigma_q(y)`, `p'(q) | d` and `q^(a-1) | d`. Hence

   ```text
   log m_(xy,1) <= Lambda'(d) + log d <= mu_M log M,     mu_M := 10 + 26/loglog M.
   ```

2. *Structured part.* For every `c`, `prod_(q<4M) q^(min(a_q(c), a*(q))) <= m_struct := prod_(q<4M) q^(a*(q))`,
   and `log m_struct <= pi(4M) log(32 M^2) <= 11 M` for `M >= M_0`.
3. *Random levels.* Let `c != 0` with `sum c = 0`. Pick `x_0 in supp c` with `|c_(x_0)| = c_min`,
   the least nonzero `|c_x|`.
   * At a large `q`, level 1, each unit: probability `<= 4 c_min/m_q`; some unit: `<= 32 c_min/q`.
   * At a random lift `a` (large `q` with `a >= 2`, or small `q` with `a > a*(q)`), given a
     collision with unit `i^s` at level `a-1` and everything else: probability `<= gcd(c_(x_0), q)/q`.
   * Different `q` are independent.
4. *Pairs at random levels* never collide with unit 1. With a unit `!= 1`, at a large `q`,
   level 1: probability `<= 12/m_q`.

*Proof.*
1. *Small `q`, `a = 1`.* As an integer, `kappa^(q)_x - kappa^(q)_y = (sigma_x - sigma_y) + 2((x mod p') - (y mod p'))`
   has absolute value `<= 2p' - 1 < m_q`. So it is `0 mod m_q` iff both parts vanish.
   * *Levels `a <= a*`.* Use the CRT components.
   * *Level `a*` is injective,* since `p' q^(a*-1) >= M`. Lifts preserve injectivity, so no pair
     collides above `a*`.
   * *Large `q`.* `y` is injective, and a difference `sigma_x - sigma_y + 2(y_x - y_y)` with
     `y_x != y_y mod m/2` is nonzero mod `m`.
   * *Demand.* `q` contributes at most `(1 + v_q(d)) log q`, and only when `p'(q) | d`. Then use
     Lemma 3.1(3).
2. `p'(q) q^(a*-2) < M`, so `q^(a*) < M q^2/p' < 8 M q <= 32 M^2`. Then use
   `pi(x) < 1.26 x/ln x`.
3. *Large `q`, level 1.* Condition on `y_x` for `x != x_0`. Then `y_(x_0)` is uniform on at least
   `m/2 - M >= m/4` values. The congruence `c_(x_0) y = t mod m/2` has at most
   `gcd(c_(x_0), m/2) <= c_min` solutions, and the parity condition can only lower the count.
   *Lifts.* `sum_x c_x kappa^(q^a)_x = sum_x c_x kappa^(q^(a-1))_x + m_(q^(a-1)) sum_x c_x u_x`.
   Given the collision at level `a-1`, a collision at level `a` is one congruence mod `q` for
   `c_(x_0) u_(x_0)`. (The unit `i` reduces compatibly: Prop. 2.1(6).)
4. Same computations with `c = e_x - e_y`. QED

*Exact check of realisability* (`check_realise.py`).
* *Set-up.* Paley `q_0 = 43` (`M = 44`), `b = 33`, canonical profile, `r = 1419` actual primes
  `= 1 mod 4` above `10^7`, and `X = 2 M^2 = 3872`.
* *Scope.* All 535 odd primes `q <= 3872` and their 564 levels.
* *Construction.* Build `D_A`. Solve (2.2). Construct `r_j(q)` with `Norm r_j = p_j mod q^(A_q)`
  by Hilbert 90 and a square root.
* *Checks, at every level:*
  * `Norm rho_x = N`;
  * `rho_x/rho_y = eta` iff `kappa_x - kappa_y = iota s`, for all pairs and all four units;
  * the line-1 pair pattern is exactly that of Lemma 3.3(1);
  * for 300 random characters per level, `prod rho^c = eta` iff `sum c kappa = iota s`.
* *Parity.* `sigma_q` is nonconstant at all 535 primes, so a parity-blind design is unrealisable
  at every one of them.
* ALL REALISATION CHECKS PASSED.

*Pair demand* (evidence, `check_design.py` Part 3, worst case over parity colourings): the maximum
of `Lambda'(d) + sum (levels)` is `5.73, 5.94, 6.61, 6.77, 6.54` times `log M` for
`M = 64, 256, 1024, 4096, 16384`.

---------------------------------------------------------------------------------------------

## 4. Profile and slack (from `lower.md`, re-derived)

**The profile.**
* Let `q_0 = 3 mod 4` be a prime with `q_0 >= 43`, `M = q_0 + 1`, and `H` the normalised Paley-I
  matrix of `lower.md` §2.
* Take `b >= 33` physical columns per nonconstant label, and the canonical assignment:
  * `alpha(x) = x` for finite `x`, with flip column `f(x) = (x, 0)` and sibling `s(x) = (x, 1)`;
  * `alpha(infinity) = a*`, with flip column `(a*, 2)` and sibling `(a*, 3)`.
* `S` is the sign matrix, `a_x = (1 + S_x)/2`, `e = 1`, and `r = b(M-1)`.
* The twins satisfy `a_(.,f(x)) - a_(.,s(x)) = eps_x e_x` with `eps_x = -H(x, alpha(x))`, and the
  `2M` twin columns are distinct.

**Lemma 4.1** (`lower.md` Lemma 3.1).
1. `L = {v(c) : c in Z^M, sum c = 0}`, and `v(c) != 0` for `c != 0`. So every rational character
   is integral.
2. Every column is nonconstant.
3. `rank S = M`.

**Lemma 4.2** (`lower.md` Lemma S). Let `T(c)_a = sum_x c_x H(x,a)` and `F(c) = sum_(a in F_q0) |T_a|`.
Every nonzero zero-sum `c in Z^M` other than a pair, a signed column `+-H_a` or a signed half-sum
`(s H_a + s' H_b)/2` has `F(c) >= (17/16) M`. If `max|c_x| >= 2`, then `F(c) >= 2M`.

*Re-derivation (checked line by line).* Put `u_a = |T_a|/n` and `kappa = #{u_a > 1/2}`.

* *Defect identity.* `sum_a |T_a| (n - |T_a|) = nF - M ||c||_2^2`, and `||c||_2^2 >= n`.
* *Window inequality.* `|F - kappa n| <= 2(F - M)`. This holds because
  `|sum u - kappa| <= 2 sum u(1-u) <= 2(F - M)/n`.
* *Inversion.* `c = H T / M` gives `F >= M max|c_x|`.

For ternary `c` of support `h >= 4` with `F < (17/16) M`, the window confines
`kappa h` to `((7/8) M, (19/16) M)`. Three cases remain.

* *`kappa >= 3`.* Then `h < 19M/48 <= (q_0 - 1)/2`. The transposed prime-field uncertainty
  (Lemma 2.2 there) applies: `T_a = -G(a) mod q_0` off the support, where
  `G(Y) = sum_x c_x (Y - x)^((q_0-1)/2) - c_infinity` is nonconstant by a Vandermonde argument.
  * So at least `M/2 - h` labels have `T_a != 0`, and at most `M/h` have `|T_a| = h`.
  * Each remaining nonzero `T_a` has `2 <= |T_a| <= h - 2`, so `|T_a|(h - |T_a|) >= 2(h - 2)`.
  * This forces `F - M >= (2(h-2)/h)(M/2 - h - M/h) >= M/16` for `M >= 44`.
* *`kappa = 1`.* `1 - u_a < 2(F-M)/h` gives `|T_a| > (3/4) M`. Column peeling then gives
  `F >= min(2|T_a|, 2M) > (3/2) M`.
* *`kappa = 2`.* Half-sum peeling gives
  `F >= M + phi(|T_a|) + phi(|T_b|) > (5/4) M`, with `phi(t) = t - |t - M/2|`.

*Evidence* (`check_lemmaS6.py`, independent code):
* exhaustive support 4: minimum `F = 60, 64, 84` for `q_0 = 43, 47, 59` (against `17M/16`);
* random supports up to `M`: minimum `F = 68, 76, 96`;
* amplitude `>= 2`: minimum `F = 212, 244, 356` (against `2M`).

**Proposition 4.3 (slack).** Let the weights satisfy `w_j in [tau, tau + 1/r]`, put
`tau_j = c^T S_j` and `2B(c) = sum_j |tau_j|`. Then `s_c >= (tau/2)(2B(c) - r) - 1/2`, and:

```text
pairs:                            2B - r >= b - 4;
every other c != 0:               2B - r >= kappa_b M,   kappa_b := min(1, b/16 - 2)  (>= 1/16 for b >= 33);
amplitude h = ||c||_inf >= 2:     2B - r >= h M (b - 4)/2 + b.
```

*Proof.*
* *Slack.* `s_c = (1/2) sum_j w_j (|tau_j| - 1)`. Every term with `tau_j != 0` has
  `|tau_j| - 1 >= 1`, and there are at most `r` terms with `tau_j = 0`, each at most `1/r` above
  `tau`.
* *Identity.* `lower.md` (3.1): `2B = bF + sum_x (|T_(alpha(x)) - 2 c_x H(x, alpha(x))| - |T_(alpha(x))|)`.
  Each flip term is `>= -2|c_x|`, so `2B - r >= b(F - M + 1) - 2n`.
* *Cases.*
  * Pairs: `F = M`, and only two flip terms move, each by `+-2`.
  * Columns: `b + 2M - 8`; half-sums: `b + M - 16` (`lower.md` Prop. 3.3). Both are `>= M`.
  * Other ternary `c`: Lemma 4.2 and `n <= M` give `b(M/16 + 1) - 2M`.
  * Amplitude `h >= 2`: `F >= hM` and `n <= hM` give `b(hM - M + 1) - 2hM >= hM(b-4)/2 + b`.

QED. *Evidence* (`check_lemmaS6.py`, `b = 33`):
* pairs: `2B - r = 29 = b - 4`, attained;
* every other tested character: `2B - r >= 69, 73, 85`, against the proved `35.8, 36.0, 36.8`.

---------------------------------------------------------------------------------------------

## 5. The theorem for P3-char

**Lemma 5.1a (one character; `lower.md` Lemma 5.3).** Let `delta` be uniform on `[0, Delta]^M`,
`c != 0` with `sum c = 0` and `n = ||c||_1`, and `eps > 0`. Then

```text
Pr[ dist(c.delta/2, (pi/4) Z) < arcsin(min(1, eps)) ]  <=  4 eps (n + pi/Delta).
```

*Proof.*
* `c.delta/2` has density `<= 2/Delta`, and takes values in an interval of length `n Delta/2`.
* At most `2n Delta/pi + 2` bad intervals meet that range, each of length `<= pi eps`. QED

**Lemma 5.1b (restriction; `lower.md` Lemma 5.2 for factored data).** Suppose `F^char_X(M', C)`
(respectively `F^all_X`) has a member with `e = 1`. Then for every `M <= M'` and `X' <= X`,
`F^char_(X')(M, C)` (respectively `F^all_(X')`) has a member with `W <= W'`.

*Proof.* Keep `M` rows, and delete the columns that are constant on them.
* *Residues.* Keep the `r_j` of the kept columns. The point residues change by a factor common
  to all kept rows, which cancels in every `prod rho^c`.
* *Inequalities.* Every character and monomial of the restriction is one of the original model,
  with the same angle and weight. Its collision modulus can only drop, since fewer moduli are
  used. So (c), (1.1) and (1.2) are inherited.
* *Normalisation.* `Delta` grows, and the profile remains normalised. QED

**Theorem 5.1 (P3-char).** Fix `A >= 1`, `b >= 33` and `C > 0`. There is `M_0 = M_0(A, b, C)` such
that for every `M >= M_0`, `F^char_(3M^A)(M, C)` has a member with `e = 1` and

```text
W  <=  (b max(A, 2) + o(1)) M log M      (M -> infinity).
```

At a prime Paley order `M = q_0 + 1` one has, explicitly,
`W <= b(M-1)(max(A,2) log M + 2 loglog M + log(42 b^2)) + 1`.

*Proof.* First let `M = q_0 + 1` be a prime Paley order, and use the profile of §4.

1. *Weights.*
   * Put `Z = max(3 M^A, 21 r^2) (log M)^2`.
   * For large `M`, `[Z, 2Z]` contains more than `r^2` primes `= 1 mod 4` (PNT in progressions).
     One of the `r` intervals of length `Z/r` that partition `[Z, 2Z]` contains `r` of them. Call
     them `p_1, ..., p_r`, all in `[A_0, A_0 + Z/r]`.
   * Put `tau = log A_0`. Then `w_j in [tau, tau + 1/r]` and `W <= r tau + 1`.
   * `p_j > X = 3 M^A`, so `N` is prime to every modulus of the class.
   * `p_j >= 21 r^2`, so `outside.md` Theorem A.1(ii) applies.
   * `tau <= max(A, 2) log M + 2 loglog M + log(42 b^2)`.
2. *Residues.* Use the design `D_A` of Def. 3.2 with `k = 0`. By Cor. 2.2 it gives data as in (d).
3. *Events.* Draw `delta` uniform on `[0, Delta]^M`, independent of the design. For every `c`
   with `sum c = 0` and `v(c) != 0` (by Lemma 4.1, these are all characters), define
   `Fail_c = {dist(c.delta/2, (pi/4)Z) < arcsin(min(1, eps_c))}`, where:
   * for non-pairs, `eps_c = m^gen_c e^(-w(v(c))/2)` and `m^gen_c = prod_q q^(max_eta a_(q,eta)(c))`;
   * for pairs, `eps_xy = m_(xy,1) e^(-w/2)`, plus the event `Off_xy = {log m^off_xy > W/4 - 1}`,
     where `m^off` is the collision modulus with units `!= 1`.

   *Sufficiency.* If no event occurs, then (1.1) holds for every `c` and every `eta`.
   * *Non-pairs.* `|sin((c.delta - arg eta)/2)| >= sin dist(c.delta/2, (pi/4)Z) >= m^gen e^(-w/2)`.
   * *Pairs, `eta = 1`.* `|c.delta/2| <= Delta/2 < pi/8`, so the distance is `|c.delta/2|`.
   * *Pairs, `eta != 1`.* The left side is `>= sin(pi/8) > 0.38`. The right side is
     `<= 2^(-1/2) e^(log m^off - w/2) <= 2^(-1/2) e^(-1)`, since `w >= W/2`.

   Since `m >= 1`, the non-occurrence also gives (G) on `L`.
4. *Union bound.* Put `theta = 1/(A+1)` and `K_c = 4(n + pi/C)`. By Lemma 5.1a and
   `min(1, x) <= x^theta`:

   ```text
   Pr[Fail_c] <= E[min(1, K_c m^gen_c e^(-s_c/2))] <= K_c^theta e^(-theta s_c/2) m_struct^theta E[m_rand(c)^theta],
   ```

   with `m^gen <= m_struct m_rand` (Lemma 3.3(2)).
   * *The random factor.* `X^theta <= 3 M^(A/(A+1))` and
     `sum_(q<=X) q^(theta-1) <= X^theta/theta =: Sigma_theta = O_A(M^(A/(A+1))) = o(M)`.
     Lemma 3.3(3) and independence across `q` give, for `q` not dividing `c_(x_0)`:
     * `E[q^(theta a_q)] <= 1 + 64 c_min q^(theta - 1)` (large `q`);
     * `E[q^(theta (a_q - a*)^+)] <= 1 + 2 q^(theta - 1)` (small `q`).

     The at most `log_2 c_min` primes dividing `c_(x_0)` give factors `<= 1 + 100 c_min A_q X^theta`.
     Hence

     ```text
     log E[m_rand(c)^theta] <= S(c) := (64 c_min + 2) Sigma_theta + log_2(c_min) log(1 + 100 c_min (A+1) M log M) = c_min o(M).
     ```

   * *Ternary non-pairs.* There are at most `3^M` of them, with `s_c >= tau kappa_b M/2 - 1/2`
     (Prop. 4.3) and `c_min = 1`. Their total is at most

     ```text
     3^M (4(M + pi/C))^theta exp(-theta(tau kappa_b M/4 - 1/4) + 11 theta M + 66 Sigma_theta)
       = exp(-M [theta(tau kappa_b/4 - 11) - log 3 - o(1)]),
     ```

     which tends to 0 since `tau -> infinity`.
   * *Amplitude `h >= 2`.* There are at most `(2h+1)^M` of them, with `n <= hM`, `c_min <= h` and
     `s_c >= 7 tau h M - 1/2`. The sum over `h >= 2` is at most
     `sum_h exp(-hM [3.5 theta tau - log(2h+1)/h - 11 theta - o(1)])`, which tends to 0.
   * *Pairs, `eta = 1`.* `s >= (b-4) tau/2 - 1/2` and `m_(xy,1) <= M^(mu_M)` (Lemma 3.3(1)). The
     total is at most `(M^2/2) 4(2 + pi/C) M^(mu_M) e^(-(b-4)tau/4 + 1/4)`. With `tau >= 2 log M`
     and `b >= 33`, this is `O_C(M^(2 + mu_M - 14.5)) -> 0`, because `mu_M -> 10`.
   * *Pairs, units `!= 1`.* `m^off_xy <= m_struct m_rand(xy)`, so
     `Pr[Off_xy] <= exp(-theta(W/4 - 1 - 11M) + 66 Sigma_theta)`. Here `W >= r tau >= 64 M log M`,
     so this is `<= M^(-3)` for large `M`.

   Hence for `M >= M_0(A, b, C)` the total probability is `< 1`. The dependence on `C` is only
   through `log(1/C)`. Fix a design and a `delta` for which no event occurs.
5. *Completion.* `rank S = M`, so `delta` is realisable. Theorem A.1(ii) gives a Haar-positive set
   of realisations `phi` satisfying (G) for `v` outside `L`. On `L`, (G) holds by step 3. The
   points are distinct, by the pair gaps. All of (a)-(d) and (e-char) hold.
6. *All `M`.* Let `M'` be the least prime Paley order `>= max(M, 44)`. Then `M' = (1 + o(1)) M`
   (PNT for primes `= 3 mod 4`). Apply Lemma 5.1b with `X' = 3M^A <= 3M'^A`. QED

**Corollary 5.2 (sharp barrier for factored residues).** Fix `C > 0` and `A >= 1`. Let `I` be any
set of inputs about a cluster consisting of:
* (I) the exponent profile and exact exponent identities;
* (II) the Lemma 0.1 gaps for all monomials;
* (III_F) the residue-strengthened character gaps (1.1), for **some** residue assignment that
  factors through prime residues of the genuine norms, at odd prime powers `<= 3M^A`. This
  includes all pigeonhole collisions used by `two_thirds`, Theorem F with every denominator, and
  every parity (quadratic reciprocity) consequence of the factorisation;
* (IV) any property of the angles valid off a Haar-null set (of each realisation torus or of the
  full torus).

Every actual cluster supplies such inputs (Prop. 1.2). An argument deriving `M <= f(R)` from `I`
alone has `f(R) >= (2/K - o(1)) log R/loglog R` along a sequence `R -> infinity`, with
`K = 33 max(A, 2)`. In particular no such argument proves `M = o(log R/loglog R)`. For `A <= 2`,
the best constant attainable from `I` lies in `[1/33, 2/3]`.

*Proof.* Theorem 5.1. Haar-positivity, also in the full torus, follows as in `referee_outside.md`
§3.1, since the good `delta` have positive measure. Then
`log R = W/2 <= (K/2 + o(1)) M log M`. QED

---------------------------------------------------------------------------------------------

## 6. P3-all: every monomial strengthened

### 6.1 A necessary condition on the data

**Proposition 6.1 (height condition).** In `F^all_X`, every `v != 0` and every `s` satisfy
`m_(v,s) <= nu_s^(-1) e^(w(v)/2)`.

*Proof.* Take (1.2) with `|sin| <= 1`. For actual circles this is `m_(v,s) | D_s(A_v)` together
with `|D_s(A_v)| <= nu_s^(-1) |A_v|`. QED

The condition concerns only residues and weights, not angles. For `v` outside `L`, the angle
`<v, phi>` is Haar-uniform on the realisation torus (`outside.md` Thm A.1(i)). So the completion
step needs the slab measures `sum_(v not in L) min(1, 4 m_v e^(-w(v)/2))` to be small, with
`m_v = prod_q q^(max_s a_(q,s)(v))`.

### 6.2 Degeneracy and the random fibre

**Lemma 6.2 (degeneracy).** For the profile of §4 (or any profile with `2M` distinct twin columns
as in Cor. 2.2(5)): if `v in Z^r \ L` and `v in L + G Z^r` for an integer `G >= 1`, then
`G <= 2 ||v||_1`.

*Proof.* Suppose `||v||_1 < G/2`, and let `v = c^T a + G u` with `sum c = 0`.
* Twin columns give `v_(f(x)) - v_(s(x)) = eps_x c_x mod G`.
* Put `c'_x = eps_x (v_(f(x)) - v_(s(x)))`. Then `c' = c mod G`, and
  `||c'||_1 <= ||v||_1 < G/2`, because the twin columns are distinct. So `sum c' = 0`, since it is
  `0 mod G` and of absolute value `< G`.
* Hence `v = c'^T a mod G`, and both sides have entries of absolute value `< G/2`. So
  `v = c'^T a in L`, a contradiction. QED

**Definition 6.3 (random fibre).** In the design `D_A`, replace the particular solution (2.2) by
`ell = nu + 2(ell'_0 + k)`. Here `ell'_0` is (2.2), and `k` is uniform in

```text
K'_Q = {k in (Z/(m_Q/2))^r : A k in (Z/(m_Q/2)) 1},
```

independently for each `q` and of everything else.

**Lemma 6.3.**
1. The point classes change only by a common shift `2t`, where `A k = t 1`. So all pair and
   character collisions are those of `D_A`.
2. A uniform `k` is obtained by drawing `k_j` (`j` not a twin column), `u in (Z/(m/2))^M` and `t`
   uniformly, and putting `k_(f(x)) = eps_x (t - (A k)_x^(non-twin) - sum_y a_(x,s(y)) u_y)` and
   `k_(s(x)) = u_x - k_(f(x))`.
3. For `v` outside `L`, every level `a <= A_q` and every `s`:

   ```text
   Pr[<v, ell^(q^a)> = iota s mod m_(q^a)] <= 4 ||v||_1 / m_(q^a)      (and <= 16 ||v||_1/m_(q^a) for some s),
   ```

   independently across `q`. Labels of non-twin columns are independent and uniform in their
   parity class.

*Proof.*
1. `A ell = A nu + 2A ell'_0 + 2t 1`.
2. Substitute `k_(s(y)) = u_y - k_(f(y))` into `(Ak)_x`. This gives
   `(Ak)_x = (non-twin part) + sum_y a_(x,s(y)) u_y + eps_x k_(f(x))`. The map from the free
   variables to `K'_Q` is a bijection.
3. *Reduction.* Reduction mod `m_(q^a)/2` maps `K'_Q` onto `K'_(q^a)`. To see this, lift
   `kbar in K'_(q^a)` arbitrarily to `k~`. Then `A k~ = t~ 1 + (m_(q^a)/2) z`; choose `w` with
   `A w = z` (twins) and take `k~ - (m_(q^a)/2) w`. So `ell' mod m_(q^a)/2` is uniform on a coset
   of `K'_(q^a)`.
   * *The event.* It requires the right parity and `<v, ell'> = const mod m_(q^a)/2`. Its
     probability is 0 or `1/|v(K'_(q^a))| = d/(m_(q^a)/2)`, where `d` is the index of the image.
   * *If `d > 1`.* `<v, .> mod d` vanishes on `K'_d = A^(-1)((Z/d) 1)`, by the same lifting.
     Since `A` is onto, any functional vanishing on `ker A` has the form `psi o A`. Vanishing on
     `A^(-1)((Z/d)1)` means `psi = lambda^T` with `sum lambda = 0 mod d`. Lifting `lambda` to an
     integral zero-sum vector shows `v in L + d Z^r`. So `d <= 2||v||_1` by Lemma 6.2.
   * The last sentence of item 3 follows from item 2. QED

### 6.3 The theorem

**Theorem 6.4 (P3-all).** Fix `A >= 1`, `b >= 33` and `C > 0`. For all large `M`,
`F^all_(3M^A)(M, C)` has a member with `e = 1` and

```text
W  <=  (2 A b + o(1)) M log^2 M / loglog M          (b = 33:  (66 A + o(1)) M log^2 M / loglog M).
```

Consequently, take inputs (I), (II), (IV) and the factored residue strengthening (1.2) of
**every** monomial at prime powers `<= 3M^A`. An argument using only these has, along a sequence
of `R`, `f(R) >= M >= (1/(Ab) - o(1)) log R logloglog R / (loglog R)^2`.

*Proof.* Use the profile of §4 at a prime Paley order, and put `X = 3M^A`.

1. *Parameters.*
   * Put `e^beta = log(4er) / (58 (loglog X + 1) loglog M)`. Then
     `beta = loglog M - O(logloglog M)`.
   * Put `theta = beta/log X`, which is `<= 1/2` for large `M`.
   * Put `tau_* = (2 log X/beta)(log(4er) + 14 + 20 beta + 10 loglog X + 58 e^beta (loglog X + 1))`.
     The bracket is `(1 + o(1)) log(4er)`, so `tau_* = (2A + o(1)) log^2 M/loglog M`.
   * Choose `r` primes `= 1 mod 4` within an interval of length `Z/r` inside `[Z, 2Z]`,
     `Z = e^(tau_*)`, as in Thm 5.1. Then `tau in [tau_*, tau_* + 1]`, `W <= r tau + 1`, and
     `p_j > X`.
2. *Residues and positions.* Use the design `D_A` with the random fibre (Def. 6.3), `k = 0`, and
   `delta` uniform. By Lemma 6.3(1) the character events of Thm 5.1, step 3, are unchanged. Their
   union bound (step 4) holds a fortiori, since `tau` is larger. So with probability `>= 2/3`,
   (1.1) and (G) hold on `L`, together with (1.2) for `v in L`. The last follows because
   `|sin| >= sin dist >= m^gen e^(-w/2) >= nu_s m_(v,s) e^(-w/2)`, and with unit `i^s` the
   monomial collisions of `v(c)` are the point collisions of `c` (`k = 0`).
3. *Monomials outside `L`.* For `v` outside `L`, the realisation torus maps onto the circle by
   `phi -> <v, phi>`, with Haar measure going to uniform measure. For each `s`, the set where
   (1.2) fails has measure `(2/pi) arcsin(min(1, nu_s m_(v,s) e^(-w/2))) <= min(1, m_v e^(-w/2))`.
   (1.2) also implies (c) for `v`. So

   ```text
   E_design[ Haar(bad) ] <= sum_(v not in L, v != 0) E[min(1, 4 m_v e^(-w(v)/2))] <= sum_v 4 e^(-theta w(v)/2) E[m_v^theta].
   ```

4. *Moments.* Let `kk = ||v||_1 <= X`.
   * By Lemma 6.3(3) and `m_(q^a) >= (2/3) q^a`: `Pr[a_q(v) >= a] <= min(1, 24 kk q^(-a))`.
   * *Primes `q > 24 kk`:* `E[q^(theta a_q)] <= 1 + 58 kk q^(theta - 1)`.
   * *Primes `q <= 24 kk`:* `E <= 1 + (log_3(24X) + 3.4) e^(2 beta)`.
   * With `pi(y) <= 1.26 y/log y` and `sum_(q<=X) 1/q <= loglog X + 1`,

     ```text
     log E[m_v^theta] <= kk [ 20 beta + 10 loglog X + 10 + 58 e^beta (loglog X + 1) ].
     ```

   * For `kk > X` use `m_v <= prod_(q^a <= X) q <= e^(1.04 X) <= e^(1.04 kk)` (Rosser-Schoenfeld)
     and `theta = 1`.
   * There are at most `(4er)^kk` vectors with `||v||_1 = kk`, and `w(v) >= kk tau`. By the
     choice of `tau_*`,

     ```text
     E[Haar(bad)] <= 4 sum_(kk>=1) e^(-4 kk) + 4 sum_(kk > X) (4 e r e^(1.04 - tau/2))^kk < 1/4.
     ```

5. *Conclusion.* By Markov, `Haar(bad) <= 1/2` with probability `>= 1/2`. Together with step 2,
   some design and some `delta` satisfy all character conditions and have `Haar(bad) < 1`. A
   realisation `phi` outside the bad set satisfies (c), (1.1) and (1.2). All `M` follow by
   Lemma 5.1b.
6. *The corollary.* `log R = W/2 <= (Ab + o(1)) M log^2 M/loglog M`, with
   `log M = (1 + o(1)) loglog R` and `loglog M = (1 + o(1)) logloglog R`. QED

### 6.4 Sharpness of the method

**Proposition 6.5 (independent fibres cannot reach `O(M log M)`).** Fix `eps > 0`. Use the
construction of Theorem 6.4 (profile of §4, equal weights `w_j in [tau, tau + 1/r]`, design `D_A`
with random fibre, any `X >= M`), but with `tau <= (2 - eps) log^2 M/loglog M`. Then with
probability `-> 1` as `M -> infinity`, some single-prime monomial `e_j` has
`m_(e_j) > sqrt2 e^(w_j/2)`. So the height condition fails, and no choice of `delta` or `phi`
completes the data to a member of `F^all`.

*Proof.* Fix `eta_0 > 0`, and let `Q_0` be the set of primes `q = +-3 mod 8` with
`M^(1-eta_0) <= q <= M`.
* *Hit probability.* For `q in Q_0`, `m_q = 4 mod 8`, so the 4-torsion `{0, m/4, m/2, 3m/4}`
  contains exactly two elements of each parity. By Lemma 6.3(2), for the at least `(b-2) M`
  non-twin columns `j` the labels `ell_j^(q)` are independent and uniform on the `m_q/2` elements
  of parity `nu_j(q)`. So `e_j` collides at `q` (level 1) with probability exactly `4/m_q`,
  independently over `j` and `q`.
* *Expected hits.* `lambda := sum_(q in Q_0) 4/m_q -> 2 log(1/(1 - eta_0))` (Mertens in
  progressions mod 8).
* *Many hits.* For `kk = floor((1 - eta_0) log M/loglog M)`, a fixed column has at least `kk`
  hits with probability `>= (1 - o(1)) e^(-2 lambda) lambda^kk/kk!`. This is
  `exp(-(1 - eta_0 + o(1)) log M)`.
* *Some column.* The expected number of such columns is `>= (b-2) M exp(-(1 - eta_0 + o(1)) log M) -> infinity`.
  By independence over columns, one exists with probability `-> 1`.
* *Height violated.* That column has `log m_(e_j) >= kk (1 - eta_0) log M = (1 - eta_0)^2 (1 - o(1)) log^2 M/loglog M`,
  whereas the height condition needs `log m_(e_j) <= tau/2 + 1`. Choose `eta_0` with
  `(1 - eta_0)^2 > 1 - eps/2`. QED

So for this method the threshold is `tau = (2 + o(1)) log^2 M/loglog M` at `X = 3M`. Theorem 6.4
attains it. The obstruction is the Poisson tail of independent per-modulus choices. It is a
property of the method, not of the class.

### 6.5 Correction to the expectation estimate

Suppose the `r_j(q)` are solved independently for each modulus. Then `E[m_v] = prod_q (1 + Pr[hit](q - 1)) ~ 2^(pi(2M))`
or more. Demanding `sum_v E[m_v] e^(-w/2) < 1` would force `log P >~ M/log M`. That expectation is
dominated by events of probability `e^(-Theta(M))`, such as one column colliding at every
`q <= 2M`. Among `r = O(M)` columns such events do not occur. Theta-moments (Theorem 6.4) give the
true rate `log P ~ log^2 M/loglog M` for independent designs, and Prop. 6.5 shows this is not an
artefact of the union bound.

### 6.6 Coherent designs and the parity

A design is *coherent* if `r_j(q) = zeta_j mod q^(A_q)` for fixed Gaussian integers `zeta_j` and
every `q`, after rescaling by rationals (which leaves `zeta/conj(zeta)` unchanged).

* *Advantage.* The collision moduli of a monomial then divide a fixed integer:
  `m_(v,s) <= |D_s(prod zeta_j^(v_j))| <= nu_s^(-1) prod |zeta_j|^(|v_j|)`. If
  `|zeta_j|^2 <= M^(O(1))`, both the height condition and the completion sum are automatic at
  `log P = O(log M)`. This is how `lower.md` Lemma 4.1 handles characters for free point residues.
* *Parity obstruction.* By Cor. 2.2(1), coherence with genuine norms needs
  `(Norm zeta_j / q) = (p_j / q)` for every prime `q <= X`. That is, `p_j Norm(zeta_j)` must be a
  quadratic residue modulo every prime `<= X`.
  * The only known way to guarantee this for every `j` is `Norm zeta_j = p_j` times a square,
    i.e. `zeta_j` = `pi_j` up to unit, rational and square factors.
  * A non-square `p_j Norm(zeta_j)` that is a QR modulo all `pi(X)` primes `<= X` is a
    pseudosquare-type integer. Heuristically the least one has `log ~ pi(X) log 2 ~ M log 2/log M`.
    Unconditionally, Burgess-type bounds for the least quadratic non-residue only give a
    polynomial lower bound in `X`, so polynomial size is not excluded rigorously. (Heuristic
    statement.)
* *Coherent with genuine `pi_j`.* The point residues are those of genuine Gaussian integers
  `zhat_x = prod zeta_j^(a_xj) conj(zeta_j)^(1-a_xj)`, and pair collisions are divisibilities of
  the genuine reduced chords `(zhat_x - zhat_y)/gcd`. For random choices of the primes, the
  maximal pair demand over `M^2` pairs has the same Poisson tail as in Prop. 6.5, now with
  `2 log M/loglog M` hits. That forces `(b-4) tau/2 >~ log^2 M/loglog M`. (Heuristic: making it
  rigorous would need equidistribution of Gaussian primes of polynomial size modulo products of
  many primes `<= 3M`. At `M = 44`, `check_p3all.py` Part C shows nothing asymptotic.)

So both natural routes stall at the same rate. One is independent and structured for pairs, the
other coherent and structured for monomials. The tension is that point classes `kappa = A ell` are
a Hadamard-type transform of the column labels. Structure for pairs (`x mod p'`) and structure for
small monomials (coherent or Sidon-like labels) do not survive that transform simultaneously. With
twins, `M` labels are always solved from the point design, and solved labels are pseudo-random.

### 6.7 The open problem

**(P3-all\*).** Is `W_min^all_(3M)(M, C) = O(M log M)`?

* *Sufficient for "yes" (Prop. 6.6).* For the profile of §4 with `log P = O(log M)`, factored
  residue data at the prime powers `<= 3M` with:
  * every pair demand `log m_xy = O(log M)`;
  * `m_v <= M^(beta ||v||_1)` for all `v` outside `L`, for a constant `beta`.

  Then Theorem 5.1 plus the completion sum `sum_kk (4er M^beta e^(-tau/2))^kk` give
  `W = O(M log M)`. The proof is steps 3-4 of Theorem 6.4 with `theta = 1`.
* *"No" would be a growth improvement.* Every actual circle lies in `F^all` (Prop. 1.2). So a
  lower bound `W_min^all(M) >= g(M) M log M` with `g -> infinity` is a theorem about actual
  clusters: `M = o(log R/loglog R)`.
* *What such a proof must use.* Angles outside `L` are free modulo the slab conditions. So it must
  come from pair pigeonhole (`two_thirds`), the height condition (Prop. 6.1), or a genuine
  covering statement for the slabs. Pigeonhole forces only average demands. The height condition
  is satisfiable for each monomial separately.
* *Assessment.* We found neither a construction nor a lower-bound mechanism. P3-all is the one
  local input class in this programme whose growth rate is still undecided. The proved bracket is
  `3 M log M - O(M) <= W_min^all <= (66 + o(1)) M log^2 M/loglog M` at `X = 3M`.

---------------------------------------------------------------------------------------------

## 7. Moduli up to `N^A`: where it breaks

Throughout, `X = N^A = e^(AW)` with `W = O(M log M)`, so `log X = Theta(M log M)`. Two new features
appear.
* Moduli `q = p_j <= X` divide `N`, so level data enter (Def. 1.1(d)). For factored data, the
  parity at `q = p_l` is given by the symbols `(p_j/p_l)`, i.e. reciprocity among the varying
  primes.
* Products of `k` points force collisions of `2k`-term characters at moduli up to about
  `binom(M+k-1, k)`. Reaching `N^A` needs `k >= M^(cA)`.

Neither forces anything by itself. What fails is every design we can analyse.

**Proposition 7.1 (fixed labellings collide everywhere).** Let `𝒬` be a set of prime powers `Q`
(`q` not dividing `N`), and suppose the point classes at every `Q in 𝒬` are induced by one
integer labelling: `kappa^Q_x = lambda_x mod m_Q`. Then there is `c != 0` with `sum c = 0` and
`||c||_inf <= (M max(1, max|lambda_x|))^(2/(M-2))` that collides (unit 1) at every `Q in 𝒬`.
Condition (1.1) then forces `sum_(Q in 𝒬) log q <= w(v(c))/2 + 1 <= ||c||_1 W/4 + 1`.

*Proof.*
* *Siegel's lemma* (2 equations, `M` unknowns, coefficients `<= max(1, max|lambda|)`) gives `c`
  with `sum c = 0` and `sum c_x lambda_x = 0`. Then `sum c_x kappa_x = 0 mod m_Q` for every `Q`.
* *The inequality.* `|sin| <= 1` in (1.1), and `w(v(c)) = (1/2) sum_j w_j |c^T S_j| <= ||c||_1 W/2`. QED

If `𝒬` contains all primes of an interval `(Y, Y']` with `sum_(Y < q <= Y') log q > ||c||_1 W/4 + 1`
(for instance `Y' = N^A`, `Y = Y'/2`), this is violated. This is the coordinator's trap. A
labelling such as `kappa_x = x` at all large moduli is excluded, as is any labelling with
`max|lambda| <= e^(o(M))`.

**7.2 Coherent designs.** These need auxiliary objects that are units modulo every `q <= X`. In
`lower.md` Lemma 4.1 that means auxiliary primes `ell_x > X`.
* Its bound `m_c^gen <= Y^(2n)/4` with `Y > X` gives, for pairs, only `log m <= 4 log Y`, and
  `4 log Y > 4 log X = 4AW`.
* The average pair slack is `<= W/(2(M-1))` (cut identity, `lower.md` §6). So once `4 log X` exceeds
  `W/(2(M-1))`, i.e. for every `X` beyond a fixed power of `M` when `W = O(M log M)`,
  the bound certifies nothing.
* For a general character the bound is `2n log Y > 2nAW`. For `A >= 1/4` this exceeds every
  possible slack, since `s_c <= (1/2) sum_j w_j |c^T S_j| - W/2 <= (n - 1) W/2`.
* For factored data, coherence also meets the parity obstruction of §6.6.

**7.3 Independent designs (heuristic).** Suppose the class maps at different moduli are
independent, with collision probability about `1/m_Q` per character per modulus.
* The characters with `||c||_1 <= 4M` number about `e^(5.3 M)`.
* The expected number of them with `kk` collisions at moduli in `[e^y, e^(2y)]` is about
  `e^(5.3M) (log 2/y)^kk / kk!`. This is `>= 1` for `kk ~ 5.3 M/log(kk y)`.
* At `y ~ M` this gives characters of demand `kk y ~ 5 M^2/log M`.
* The slack of such (random-looking) characters is about `(W/2) ||c||_2 = O(M^(3/2) log M)`
  (CLT heuristic). So the design fails at moduli around `e^(Theta(M))`, inside `[M^A, N^A]`.
* Avoiding short characters at large moduli (`B_h`-type maps, or "avoid all `c` with
  `#{||c'||_1 <= ||c||_1} < m_Q/2`") removes only characters with `||c||_1 <~ log m_Q/log M`. It
  does not help at `||c||_1 ~ M`. No second-moment computation was done.

**Proposition 7.4 (what would suffice).** Suppose that for the profile of §4 there are factored
residue data at every odd prime power `<= X(M)`, including level data at the `p_j <= X`, such
that:
* every pair has `log m_xy <= C_1 log M`;
* every non-pair character has `log m^gen_c <= C_1 ||c||_1 log M`.

Then `W_min^char_(X(M))(M, C) = O_(C_1)(M log M)`.

*Proof.* Step 4 of Theorem 5.1 with these bounds in place of `m_struct m_rand`:
* *Ternary non-pairs:* `tau kappa_b M/4 > C_1 M log M + M log 3` once `tau >= (4C_1/kappa_b + 1) log M`.
* *Amplitude `h`:* `n <= hM` and the slack `7 tau h M` absorb `C_1 h M log M`.
* *Pairs:* as before.
* *Completion:* plain gaps outside `L`, Thm A.1(ii). QED

So the extension to `N^A` reduces exactly to the existence of such a **balanced residue design**.
Independent designs have Poisson tails, fixed labellings fail by Prop. 7.1, and coherent ones fail
by §7.2. We found no other candidate. Whether such designs exist is open, for free and for
factored residues alike.

---------------------------------------------------------------------------------------------

## 8. Evidence (no proof depends on it)

| script | what it does | result |
|---|---|---|
| `gres.py` | Gaussian residues mod `q^a`, `T_Q`, generators, Hilbert 90, square roots | helper |
| `profile.py` | independent builder of the canonical one-flip Paley profile | helper |
| `check_parity.py` | Prop. 2.1 (1)-(6) exhaustively for all 89 odd prime powers `<= 400`; 211 Gaussian primes | ALL PARITY CHECKS PASSED; quartic part of `ell` not a function of `(p/q)` whenever `8 | m` |
| `check_design.py` | Lemma 3.1 (multiplicity: exact `<= 4` for `q < 2^23`; R-S bound `8.11` for `k >= 20`), pair demand of `D_A` | max demand `5.7 - 6.8 log M` for `M = 64 ... 16384`; ALL DESIGN CHECKS PASSED |
| `check_realise.py` | end-to-end exact realisation of `D_A` with forced parity: `M = 44`, `b = 33`, 1419 actual primes, all 535 odd `q <= 3872` and 564 levels | ALL REALISATION CHECKS PASSED; parity-blind design unrealisable at all 535 primes |
| `check_lemmaS6.py` | Lemma 4.2 and Prop. 4.3 for `q_0 = 43, 47, 59` (exhaustive support 4, random, amplitude `>= 2`, columns, half-sums; `b = 33` profile) | min `F/M = 1.36, 1.33, 1.40` at support 4; pair bound `b - 4` attained; others `>= 69` vs proved `35.8` |
| `check_p3all.py` | random fibre at `M = 44`: single-column and 2-column collision moduli; abstract independent model up to `M = 3*10^4`; genuine-coherent pair demand | single columns `log m/log M` mean 3.3, max 8.0; abstract model max `9.9, 10.6, 12.5, 12.7` for `M = 10^2 ... 3*10^4` (slow growth, consistent with `log M/loglog M` tails; small `M` is not asymptotic); genuine-coherent pair demand `3.8 log M` vs structured worst case `5.6 log M` at `M = 44` (no asymptotic information) |

---------------------------------------------------------------------------------------------

## Appendix A. The Sylvester multi-block route (original lead; superseded)

`lead6/notes.md` §7-9 proposed:
* `t` blocks of the Sylvester matrix (nonconstant columns), with independent uniformly random row
  permutations;
* one twin column per row: a copy of a block-`t` column, flipped at that row, with the twin labels
  covering every label (one repeat).

It works, with an absolute but astronomically large constant (proportional to the Green-Sanders
bound), and it is subsumed by §5. The parity-respecting design of §3 applies verbatim, since only
twins are used. These are the statements that were checked. They are written at the level of
detail below; a referee pass is still owed.

* **Margin identity.** Put `E_i = F_i - (M-1)`, with `F_i` the transform `l^1` norm in block `i`.
  For equal weights, `s_c >= (tau/2) G(c) - 1/2` with `G(c) = sum_i E_i(c) + Tw(c)`, where the
  twin term satisfies
  `Tw >= max(-M, E_t - 1 - 2n)` (onto labels give `sum_x |T_(a(x))| >= F_t`).
  * *Pairs:* `E_i = 1` (two Sylvester rows differ in exactly `M/2` nonconstant columns) and
    `Tw >= -4`, so `G >= t - 4`.
  * *`n >= 2M`:* `F >= n` (inversion) gives `G >= (t-1) n/2`.
* **Lemma A.1 (Hoelder alignment).** Let `c` be `+-1` with support `s >= 4` in one block, and
  `E < 1 + M/(s-2)`. Then the permuted support lies in an affine subspace of `F_2^m` of dimension
  `<= ceil(log2 s)`. If some `|c_x| >= 2`, then `E >= 1 + 2M/n`.

  *Proof.* `sum_a |T_a|(n - |T_a|) = n(E-1)` for `+-1` vectors. This gives
  `#{|T_a| = n} >= M/n - (n-2)(E-1)/(2n) > M/(2n)`. The aligned set is a coset of
  `span(S - S)^perp`, so `|span(S - S)| < 2s`. The non-`+-1` case is the defect identity with
  `||c||_2^2 - n >= 2`.
* **Lemma A.2 (level sets; Boolean version suffices).** Let `f : F_2^m -> Z` with
  `||f||_A := sum_a |fhat(a)| <= 3`, where `fhat(a) = E_x f(x)(-1)^(a.x)`. Then `f` is a signed sum
  of at most `L = 12 L_GS(225)` coset indicators.

  *Proof.* The Wiener algebra is a Banach algebra with `||1||_A = 1`, so
  `||1_(f=v)||_A <= prod_(u != v, |u| <= 3) (3 + |u|)/|v - u| <= 225` (the maximum, at `v = +-1`), and
  `f = sum_(v != 0) v 1_(f=v)`. Then apply the idempotent theorem:
  * Green-Sanders, *A quantitative version of the idempotent theorem in harmonic analysis*, Ann.
    of Math. 168 (2008) 1025-1054, with `L_GS(K) <= exp exp(C K^4)` for indicator functions; the
    `F_2^n` case is Green-Sanders, *Boolean functions with small spectral norm*, GAFA 18 (2008).
  * Normalisation: `||1_H||_A = 1` for cosets `H`.
  * In `F_2^n` a coset is a difference of two subgroups, so "cosets" versus "subgroups" costs a
    factor 2.
  * For integer-valued `f` directly: Sanders, *Bounds in Cohen's idempotent theorem*
    (arXiv:1610.07092), `exp(K^(4+o(1)))` coefficients. The statement is taken from a search
    summary and was not re-read.
* **Counting.**
  * `#{affine subspaces of F_2^m} <= 3.47 (m+1) 2^(m^2/4 + m)`.
  * `Pr[pi(S) lies in a d-dimensional affine subspace] <= 3.47 2^(d(s-d-1)) M^(d+1-s)`. Support 4
    gives exactly `1/(M-3)`.
  * `Pr[f_i structured] <= |S_L| / Mult(mu_c)`, with `log |S_L| = O(L log^2 M)` and
    `Mult(mu) >= binom(M, ceil(s/4))` for entries `<= 3`.
* **Assembly** (`t = 5`, `P = M^(A')`, `A' = O(L)`). W.h.p. over the permutations:
  * every character of support `>= 4 L log M` is unstructured (`||f||_A > 3`, hence `E > 2M`) in
    some block;
  * every `+-1` character of support `4 <= s < 4 L log M` is non-aligned in some block (support 4
    needs `t >= 5`);
  * non-`+-1` small supports have `E >= 1 + 2M/n` in every block.

  So every non-pair has `G >= M/(4L log M)` (or more). The bound `m_c <= e^(2.6M)` (all
  accessible moduli) is then absorbed once `A' >~ L`, and pairs use the design of §3. This gives
  `W <= (6 A' + o(1)) M log M`, with `A'` proportional to `L_GS(225)`.

Paley profiles avoid all of this. Lemma S replaces Green-Sanders by an elementary prime-field
uncertainty argument, and a single block with `b = 33` copies replaces the random permutations.

---------------------------------------------------------------------------------------------

## Appendix B. The original task's checklist

| item | answer |
|---|---|
| (a) actual clusters in the class; `two_thirds` inputs valid | Prop. 1.2 and Prop. 1.3. The class contains III-small's realisability (twins make `A` onto) and its strengthening. Factored data add only the parity constraint |
| (b) fakes with `log N <= c_0 M log M` for every `C` | Thm 5.1 with `A = 1`: `c_0 = 66 + o(1)` (Paley, `b = 33`), uniform in `C` (only `M_0` depends on `log(1/C)`) |
| twin-column penalty | Paley: flip terms `>= -2|c_x|` (Prop. 4.3); pairs `b - 4`. Sylvester: `Tw >= max(-M, E_t - 1 - 2n)`, pairs `t - 4` |
| pair margins vs residue demand vs positional term | Thm 5.1 step 4: `(b-4) tau/4 > (2 + mu_M + o(1)) log M`, with `mu_M -> 10` and the factor `4(2 + pi/C)` |
| characters with large entries | amplitude classes: `2B - r >= hM(b-4)/2 + b` (Prop. 4.3), count `(2h+1)^M` |
| rational characters | saturation (Lemma 4.1): every rational character is integral, so windings drop out |
| realisability of residue maps (prime powers, units, CRT) | Cor. 2.2. One residue mod the top power gives all levels. Units: `k = 0`; `i` flips parity iff `q = +-3 mod 8`. CRT across `q` is free. Exact end-to-end check (§3) |
| union bounds converge with `t` absolute | Paley: no `t`; `b = 33` fixed. Sylvester: `t = 5` (Appendix A) |
| short-interval primes | not needed: one interval of length `Z/r` in `[Z, 2Z]` (PNT in progressions mod 4 plus pigeonhole) |
| Green-Sanders statement | Appendix A (Boolean version suffices via level sets) |
| exact small-`M` computations | §8 (residue design demand, parity, realisability, Lemma S). The Sylvester support-4 alignment data are in `lead6/notes.md` §8 and were not recomputed |

---------------------------------------------------------------------------------------------

## Files (`round6/sharp_checks/`)

`gres.py`, `profile.py` (helpers); `check_parity.py`, `check_design.py`, `check_realise.py`,
`check_lemmaS6.py`, `check_p3all.py` with outputs `check_*_out.txt`. The repository
`/home/user/Jarnik` and `round5/` were not modified.
