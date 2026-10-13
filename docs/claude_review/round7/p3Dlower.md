# Round 7, P3-all\*: can designed prime residues be forced to cost `>> M log M`? (lower-bound attempt)

Status: DRAFT, unrefereed. Companion to `round7/p3lower.md`, `round7/p3G.md` and
`round5/round6_sharp/sharp.md`. Scripts: `round7/p3Dlower_checks/` (each runs from that directory
with `python3 <script> [args]`; outputs are stored next to them as `*_out.txt`). The repository and
the round-5/6/8 folders were not modified.

Conventions are those of `sharp.md` (Def. 1.1, Prop. 2.1, Cor. 2.2):
* `W = log N`; `w(v) = sum_j |v_j| log p_j`; `L` is the character lattice; `n = ||c||_1`.
* For an odd prime `q` not dividing `N`, `Q = q^(A_q)`, `T_Q` is cyclic of order `m = m_Q`, `4 | m`.
  A residue datum at `q` is a *label vector* `ell = ell(q) in (Z/m)^r` (`r_j/conj r_j = g^(ell_j)`),
  and `iota = log_g(i)`, so `{0, iota, 2 iota, 3 iota} = (m/4) Z/m` is the log of `mu_4 = {1, i, -1, -i}`.
* `v` *collides at `q^a` with unit `i^s`* if `<v, ell> = s iota mod m_(q^a)`; `a_(q,s)(v)` is the
  largest such `a <= A_q`; `m_(v,s) = prod_q q^(a_(q,s)(v))`, and `m_v = prod_q q^(max_s a_(q,s)(v))`
  is the all-unit modulus used in the completion sum (`sharp.md` §6.1).
* **Parity** (`sharp.md` Prop. 2.1(3), Cor. 2.2(1)): factored data with genuine norms are exactly the
  label vectors with `ell_j = nu_j(q) mod 2`, `nu_j(q) = [(p_j/q) = -1]`.
* **P3-all\*** = `F^all_X(M, C)` of `sharp.md` with the per-unit height condition
  `m_(v,s) <= nu_s^(-1) e^(w(v)/2)` (necessary, `sharp.md` Prop. 6.1 as corrected) and the
  completion requirement for `v` outside `L` (`sharp.md` §6.1). Actual clusters are members
  (`sharp.md` Prop. 1.2).
* `M log M` is the target scale; *growth* means `W_min >= f(M) M log M` with `f -> infinity`.

*Theorem/Proposition/Lemma* = proved here with all quantifiers, or a cited classical result named as
such. *Conditional* results name their hypothesis. *Evidence* = a computation no proof uses.
*Heuristic* = an argument with a named unproved independence assumption.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The target was not reached.** No lower bound `W >= f(M) M log M` with `f -> infinity` is proved for
P3-all\*, nor for its near-uniform-margin subclass, nor for any natural subclass. No growth
improvement for actual circles is claimed. The proposed mechanism fails at two identified points,
proved here, and the remaining gap is isolated as a single combinatorial "inverse" step.

| # | statement | class | status |
|---|---|---|---|
| 1.1 | **parity localisation:** at `q = +-3 mod 8` parity-respecting and unconstrained labels have the same collision patterns and the same all-unit moduli `m_v` (units are permuted); at `q = +-1 mod 8` the parity only vetoes collisions of `v` with `<v, nu(q)>` odd | factored data (P3-char, P3-all, P3-all\*), actual | proved |
| 2.1 | **twin decomposition:** on twin profiles, `(ell_J, U, Lambda)` is a parity-compatible change of basis; point classes are `kappa = F(ell_J, U) + Lambda`; every monomial is the same-sparsity combination (2.2) of the basis labels; every `kappa` (parity colouring respected) combines with **every** parity-respecting free label vector | twin profiles in P3-all\* (incl. the round-6 fakes) | proved |
| 2.2 | consequence: low pair demand never forces incoherence of the `r - M` free labels; it acts on monomials only through `Lambda = kappa - F` (twin-asymmetric monomials) | same | proved |
| 3.1 | additive labels (`ell = g_q lambda`, `|lambda| <= M^h`) force `W >= (1/4 - o(1)) M X/(h+2)` (Siegel trap for monomials) | P3-all\*, twin profiles, `p_j > X` | proved |
| 3.2 | exactly Gaussian-coherent labels with the right parity at all `q = +-1 mod 8 <= X` have `Norm zeta_j >= p_j`; so coherent certificates have zero slack and can never certify completion | P3-all\* | **proved under GRH** (effective Chebotarev) |
| 4.3 | twisted coherence (Def. 4.1, parity repaired by a coherent non-square) is admissible and has explicit **quadratic** caps `log m_(v,s) <= (T(2k+1)+4) log(sqrt2|Z_v|) + T k(k+1)/2 log ell_T`; linear part only from `q = +-3 mod 8` | P3-all\* | proved |
| 4.4 | twisted coherence certifies the height condition only up to norm `k_1 = O(tau/(T log M))` | P3-all\* | proved |
| 4.5 | twisted coherence lands at the random rate `tau ~ log^2 M/loglog M` | model | heuristic (independence of divisibility events) |
| 5.1 | per modulus, parity forces no collision among fewer than `m_q/8` monomials; counting arguments cannot produce concentration | P3-all\* | proved |
| 5.2 | every analysed design (additive, coherent, genuine, random fibre, twisted, `D_A` + coherent free labels) needs `tau >~ log^2 M/loglog M` or more | -- | table: proved for additive and random fibre, GRH for exact coherence, heuristic otherwise |

**Answer to the task.**
1. *The proposed mechanism does not work as stated, for two proved reasons.*
   * (M1) "low pair demand forces incoherent residues" is **false** for the `r - M` free labels on
     twin profiles (Cor. 2.2). Only the `M` twin labels `Lambda = kappa - F` are constrained.
   * (M2) "pigeonhole/entropy over residue configurations" **cannot** force a height violation
     (Lemma 5.1 with `p3lower.md` B.3/B.5). The parity is invisible at half the moduli
     (Thm 1.1) and is only a veto at the rest. Counting yields collisions on average, which fit in
     the budget at `tau = Theta(log M)`; growth needs concentration, which counting never forces.
2. *What is proved instead about designed residues.* Certified designs cannot reach `O(M log M)`:
   * additive certificates collapse (Prop. 3.1);
   * exact multiplicative certificates with the forced parity are of genuine size, with zero slack
     (Thm 3.2, GRH);
   * the best parity repair costs a factor `k` (Thm 4.3).
   So any fake at `W = O(M log M)` in P3-all\* would have to use **uncertified, non-random** labels
   whose sparse combinations nevertheless have bounded agreement with `mu_4` at every scale.
3. *What a growth theorem needs* (§5.3): (I1) an inverse theorem "bounded monomial demand at all
   scales implies certified labels" (open, and of discrepancy type); (I2) certified implies genuine
   size (proved here under GRH for exact coherence); (I3) genuine size implies Poisson tails
   (Hypothesis `H` of `p3lower.md`; open). The bracket for P3-all\* is unchanged:
   `3 M log M - O(M) <= W_min <= (66 A + o(1)) M log^2 M/loglog M` (`X = 3M^A`).

**Validity for actual circles.** Theorems 1.1 and 2.1 are identities about factored residue data
and hold for actual clusters (genuine labels `log(pi_j/conj pi_j)`). Proposition 3.1, Theorem 3.2
and Theorem 4.3 are statements about *designs*; for actual clusters the labels are genuine, so
Theorem 3.2 is consistent (`zeta_j = pi_j`, `Norm zeta_j = p_j`) and Proposition 3.1 does not apply
(genuine labels are not additive). Nothing here bounds actual clusters beyond `two_thirds`.

---------------------------------------------------------------------------------------------

## 1. Where the forced parity acts

**Theorem 1.1 (parity localisation).** Fix an odd prime `q <= X` not dividing `N`.

1. *(`q = +-3 mod 8`: parity is invisible.)* For every label vector `ell' in (Z/m)^r` (no parity
   constraint) and every `nu in {0,1}^r` there is `u in (Z/4)^r` such that `ell := ell' + iota u`
   satisfies `ell_j = nu_j mod 2` for all `j`, and for every `v in Z^r`, every level `a <= A_q` and
   every `s in Z/4`:

   ```text
   v collides at q^a with unit i^s under ell   <=>   v collides at q^a with unit i^(s - <v,u>) under ell'.   (1.1)
   ```

   In particular `max_s a_(q,s)(v)` (hence the all-unit modulus `m_v`, the completion sum and the set
   of `(v, a)` that collide with *some* unit) is the same for `ell` and `ell'`.
2. *(`q = +-1 mod 8`: one bit, and only as a veto.)* For parity-respecting `ell`, if `v` collides at
   `q` with any unit, then `<v, nu(q)> = 0 mod 2`. Conversely, for `ell` uniform on the parity
   class and `v` with `gcd(v_1, ..., v_r, m_q/2) = 1`, the probability of a level-1 collision with some
   unit is exactly `8/m_q` if `<v, nu(q)>` is even and `0` if it is odd (it is `4/m_q` without the
   parity constraint).

*Proof.*
1. By `sharp.md` Prop. 2.1(5), `i` is a non-square in `T_Q` iff `q = +-3 mod 8`, so `iota` is odd.
   Put `u_j = [ell'_j != nu_j mod 2]`; then `ell_j = ell'_j + iota u_j` has parity `nu_j`.
   `iota u` is the log of `i^(u)`, a vector of units, so its reduction to level `a` is
   `iota_a u` with `iota_a = log i` at level `a` (Prop. 2.1(6)). Hence
   `<v, ell> = <v, ell'> + iota <v,u> mod m_(q^a)`, which is (1.1). For fixed `(v, a)` the four
   events `<v, ell> = s iota mod m_(q^a)` (`s in Z/4`) are disjoint because `4 | m_(q^a)`, so
   `max_s a_(q,s)` is unchanged.
2. If `q = +-1 mod 8`, `i` is a square, `iota` is even, so a collision forces `<v, ell>` even, and
   `<v, ell> = <v, nu> mod 2`. For the probability: write `ell = nu + 2 ell''` with `ell''` uniform
   in `(Z/(m/2))^r`. Since `gcd(v, m/2) = 1`, `<v, ell''>` is uniform on `Z/(m/2)`, so `<v, ell>`
   is uniform on the `m/2` residues of parity `<v, nu>`. Here `8 | m` (`m = q - chi(q)` and
   `q = +-1 mod 8`), so all four elements of `(m/4) Z/m` are even. Without the constraint `<v,ell>`
   is uniform on `Z/m`. QED

**Remarks 1.2.**
* Item 1 is the global form of `sharp.md` Remark 2.3(1) ("the unit `i` flips parity iff
  `q = +-3 mod 8`"): at half of the moduli (Dirichlet density 1/2), *every* collision pattern
  available to unconstrained factored data is available to parity-respecting data, with the same
  all-unit moduli. A parity-based lower bound can therefore only use moduli `q = +-1 mod 8`.
* Item 2 is the single-modulus content of `p3lower.md` Prop. B.3: the parity is a *veto* (odd
  `<v, nu>` cannot collide) and a factor 2 on the even side. It forces no collision.
* Neither item depends on the profile, on pairs or on angles.

---------------------------------------------------------------------------------------------

## 2. Twin profiles: the pair channel touches only the twin labels

The round-6 constructions (the profile of `sharp.md` §4, used in its Theorems 5.1 and 6.4) are
*twin profiles* (`sharp.md` Cor. 2.2(5), Lemma 6.2): `e = 1`, and for every row `x` there are columns
`f(x), s(x)` with `a_(., f(x)) - a_(., s(x)) = eps_x e_x`, `eps_x = +-1`, all `2M` columns distinct.
Let `J` be the other `r - 2M` columns and `sigma_(yx) := a_(y, s(x))`. A lower-bound mechanism for
P3-all\* must in particular work on twin profiles (the designer may choose one), so this section
determines exactly what pair data can say there. Units `k = 0` for simplicity (windings add the
fixed shift `-iota k_y` to `kappa_y` and change nothing below).

**Theorem 2.1 (twin decomposition).** Fix a twin profile and an odd prime `q` not dividing `N`
(top modulus `Q`, `m = m_Q`). Define the *basis labels*

```text
ell_J = (ell_j)_(j in J),     U_x := ell_(f(x)) + ell_(s(x)),     Lambda_x := eps_x ell_(f(x)),                      (2.1)
F_y := sum_(j in J) a_(yj) ell_j + sum_x sigma_(yx) U_x .
```

1. `ell -> (ell_J, U, Lambda)` is a bijection of `(Z/m)^r` (inverse:
   `ell_(f(x)) = eps_x Lambda_x`, `ell_(s(x)) = U_x - eps_x Lambda_x`), compatible with reduction to
   every level `q^a`. `ell` is parity-respecting iff `ell_J = nu_J`, `U_x = nu_(f(x)) + nu_(s(x))` and
   `Lambda_x = nu_(f(x))` mod 2.
2. The point classes are `kappa_y = <a_y, ell> = F_y + Lambda_y`. Every point-level datum (point
   residues up to a common factor, pair and character collisions at every level and unit) is a
   function of `kappa` alone.
3. For every `v in Z^r` put `vt := (v_J, (v_(s(x)))_x, c'(v))` with `c'(v)_x := eps_x (v_(f(x)) - v_(s(x)))`. Then

   ```text
   <v, ell> = <v_J, ell_J> + <v_s, U> + <c'(v), Lambda>          (all levels),                                     (2.2)
   ```

   and `v -> vt` is a unimodular bijection of `Z^r` with `||vt||_1 <= 2||v||_1` and `||v||_1 <= 2||vt||_1`.
4. *(Decoupling.)* For every `kappa in (Z/m)^M` and every parity-respecting `(ell_J, U)` there is
   exactly one `Lambda` with `F + Lambda = kappa`, namely `Lambda = kappa - F`, and the resulting
   `ell` is parity-respecting iff `kappa_y = F_y + nu_(f(y)) mod 2` for all `y`. That condition
   reads `kappa_y = sigma_q(y) + const mod 2`, with `sigma_q` the parity colouring of `sharp.md`
   (2.1), and does not involve `ell_J` or `U` beyond their (fixed) parities.

*Proof.*
1. The inverse formula is immediate, and all maps are `Z`-linear with integer coefficients, so they
   commute with reduction mod `m_(q^a)`. The parity statement follows from the inverse formula.
2. `<a_y, ell> = sum_J a_yj ell_j + sum_x (a_(y,f(x)) ell_(f(x)) + a_(y,s(x)) ell_(s(x)))` and
   `a_(y,f(x)) = sigma_(yx) + eps_x delta_(xy)`, `a_(y,s(x)) = sigma_(yx)`, so the twin sum equals
   `sum_x sigma_(yx) U_x + eps_y ell_(f(y))`. Point residues are `rho_y = B g^(kappa_y)` with `B`
   common (`sharp.md` Cor. 2.2(2)); collisions of a character `c` are conditions on
   `sum c_y kappa_y` (Cor. 2.2(3)).
3. `v_(f(x)) ell_(f(x)) + v_(s(x)) ell_(s(x)) = v_(s(x)) U_x + (v_(f(x)) - v_(s(x))) ell_(f(x))`, and
   `ell_(f(x)) = eps_x Lambda_x`. The inverse is `v_(f(x)) = v_(s(x)) + eps_x c'_x`. The norm bounds
   follow from `|c'_x| <= |v_(f(x))| + |v_(s(x))|` and `|v_(f(x))| <= |v_(s(x))| + |c'_x|`.
4. Immediate from items 1 and 2: the parity of `Lambda_y = kappa_y - F_y` must be `nu_(f(y))`, and
   `F_y = sum_J a_yj nu_j + sum_x sigma_(yx)(nu_(f(x)) + nu_(s(x))) mod 2` depends only on the
   parities. Comparing with `sharp.md` (2.1): `sigma_q(y) = sum_j a_yj nu_j = F_y + nu_(f(y)) mod 2`. QED

**Corollary 2.2 (what pair data can and cannot constrain).** On twin profiles:
* *Free labels are unconstrained by pairs.* Whatever the pair demands are (in particular for every
  pair design with `max D_xy = O(log M)`, such as `D_A` of `sharp.md` §3), the free labels
  `(ell_J, U)` at each modulus may be **any** parity-respecting vector: `(r - M)` coordinates per
  modulus, completely decoupled from pairs. So "low pair demand forces incoherent residues" is
  false for the free labels. The twin-symmetric monomials (`c'(v) = 0`, a sublattice of rank
  `r - M`) see only free labels, so pair data say nothing about them.
* *Everything goes through `Lambda = kappa - F`.* A monomial `v` interacts with the pair data only
  through `<c'(v), Lambda>`, i.e. only if it is asymmetric on some twin pair. By (2.2) such `v` with
  `||c'(v)||_1 = 1` and free part `w := (v_J, v_s)` collides at `q` iff

  ```text
  Lambda_y  in  -+ <w, (ell_J, U)> + (m/4) Z/m         (y = the asymmetric twin, sign of c'_y).                (2.3)
  ```

  So the only route from "pairs are cheap" to "some monomial is expensive" is: **the twin labels
  `Lambda_y = kappa_y - F_y` must lie close (in the sense (2.3), at many moduli) to sparse
  combinations of the free labels.**

*Remark 2.3.* The decoupling is exact and holds for every `X`, every level and every unit, with the
genuine norm constraint in force. It also shows that the task's proposed mechanism has to be stated
for the twin labels: incoherence of the free labels is never forced.

---------------------------------------------------------------------------------------------

## 3. Certified designs: additive labels die, coherent labels are forced to genuine size

A monomial bound is *certified* if it comes from a global object: the collisions of `v` at all
moduli are divisors of one integer whose size is controlled by `||v||_1`. Two kinds are available:
additive (labels are reductions of fixed integers) and multiplicative (labels are logs of fixed
Gaussian integers, as in `lower.md` Lemma 4.1 and `sharp.md` §6.6). This section shows what each
costs in P3-all\*.

**Proposition 3.1 (additive labels: the Siegel trap for monomials).** Let a member of P3-all\* be
given, a set `𝒬` of odd primes `q <= X` not dividing `N`, integers `lambda_1, ..., lambda_r` with
`|lambda_j| <= H`, and for each `q in 𝒬` elements `g_q in Z/m_Q` and `u(q) in (Z/4)^r` with

```text
ell(q) = g_q lambda + iota u(q)    in (Z/m_Q)^r         (q in 𝒬).                                              (3.1)
```

Let `r' = ceil(r/2)` and let `k >= 1` satisfy `binom(r', k) > 2kH + 1`. Then, with
`Pi_𝒬 = prod_(q in 𝒬) q^(A_q)`,

```text
W  >=  (r/2) (log Pi_𝒬 - log 2) / (4k).                                                                        (3.2)
```

*Proof.* Let `J'` be `r'` columns of smallest weight, and `w_* = max_(j in J') w_j`. At least
`r - r' + 1 >= r/2` columns have weight `>= w_*`, so `W >= (r/2) w_*`. The `binom(r', k)` sums
`sum_(j in S) lambda_j` (`S` a `k`-subset of `J'`) lie in `[-kH, kH]`, so two distinct subsets
`S, S'` give equal sums. Put `v = 1_S - 1_(S')`: `v != 0`, `||v||_1 <= 2k`, `<v, lambda> = 0`. By
(3.1), `<v, ell(q)> = iota <v, u(q)>` at every level, so `v` collides at `q^(A_q)` with the unit
`i^(<v,u(q)>)` for every `q in 𝒬`. Hence `prod_s m_(v,s) >= Pi_𝒬`. The per-unit height condition
gives `m_(v,s) <= nu_s^(-1) e^(w(v)/2)`, so `Pi_𝒬 <= 2 e^(2 w(v)) <= 2 e^(4k w_*)`. Combine. QED

*Consequence.* Take `𝒬` = all odd primes `<= X = 3M^A` not dividing `N`, and suppose (as in every
construction of rounds 5-7) that all `p_j > X`. Then `log Pi_𝒬 >= theta(X) - log 2 = (1 - o(1)) X`
(prime number theorem). For twin profiles `r >= 2M`, and if `H <= M^h` (`h` fixed), then
`k = h + 2` satisfies `binom(r',k) > 2kH + 1` for large `M`, so `W >= (1/4 - o(1)) M X/(h+2)`,
i.e. at least of order `M^2/h`. Additive certification is therefore useless unless, roughly,
`log H >~ X loglog M/log M` (so that the relation `v` has norm `>~ X/log M`), and then the certificate `|<v, lambda>| <= ||v||_1 H` no longer bounds anything at the scale `log M`.
(This is `sharp.md` Prop. 7.1 transposed from characters to monomials; the point is that `D_A`
is additive on *points* only, and Theorem 2.1 shows its labels are then `Lambda = kappa - F`.)

**Theorem 3.2 (multiplicative labels with the right parity are of genuine size; conditional on
GRH).** Assume the Generalized Riemann Hypothesis for Dedekind zeta functions of biquadratic fields.
There is an absolute constant `c_1` such that the following holds. Let `p = 1 mod 4` be a prime and
`zeta` a Gaussian integer with `n := p Norm(zeta)`, and let `X >= c_1 (log 2n)^2 (loglog 3n)^2`.
Assume (as in every construction, and in `sharp.md` Thms 5.1/6.4) that every prime of `N` exceeds
`X`. Suppose a parity-respecting label of column `j` (with `p_j = p`) satisfies

```text
ell_j(q) in log(zeta/conj zeta) + (m_q/4) Z/m_q       for every prime q <= X, q = +-1 mod 8, q not dividing n.   (3.3)
```

Then `n` is a square or twice a square. In particular `Norm(zeta) >= p`, i.e.
`|zeta| >= p^(1/2) = e^(w(e_j)/2)`.

*Proof.*
* *Parity.* At `q = +-1 mod 8`, `m_q/4` is even (`8 | m_q`), so (3.3) and the parity constraint give
  `log(zeta/conj zeta) = nu_j(q) mod 2`. By `sharp.md` Prop. 2.1(3) the left side is even iff
  `(Norm zeta/q) = 1`, and `nu_j(q) = 0` iff `(p/q) = 1`. Hence `(n/q) = 1` for every prime
  `q <= X`, `q = +-1 mod 8`, `q` not dividing `n`.
* *Chebotarev.* Suppose `n` is neither a square nor twice a square. Then `sqrt n` is not in
  `Q(sqrt 2)`, so `L = Q(sqrt n, sqrt 2)` is biquadratic, with Galois group `{1, sigma, tau, sigma tau}`
  where `sigma` fixes `sqrt 2` and negates `sqrt n`. An odd prime `q` not dividing `n` has Frobenius
  `sigma` iff `(2/q) = 1` and `(n/q) = -1`, i.e. iff `q = +-1 mod 8` and `(n/q) = -1`. The
  discriminant satisfies `log d_L <= 3 log(32 n)`. By the effective Chebotarev theorem under GRH
  (Lagarias-Odlyzko 1977, in Serre's form), the number of such `q <= x` is
  `x/(4 log x) + O(x^(1/2) (log d_L + log x))`. For `x >= c_1 (log 2n)^2 (loglog 3n)^2` with `c_1`
  large this exceeds `omega(n) <= log_2 n`, so some such `q <= X` does not divide `n` (nor `N`, whose
  primes exceed `X`). That contradicts the first bullet.
* *Size.* If `p Norm(zeta) = t^2` or `2t^2`, then `p | t` (as `p` is odd), so `Norm(zeta) >= p`. QED

**Corollary 3.3 (coherent certificates have no slack).** Suppose every column's labels satisfy
(3.3) with Gaussian integers `zeta_j` and `X >= c_1 (log 2 n_max)^2 (loglog 3 n_max)^2`. Under GRH:
* The multiplicative certificate for a monomial `v` is (as in `sharp.md` §6.6)
  `m_(v,s)^(+-1 mod 8) <= prod_(s' in Z/4) |D_(s')(Z_v)| <= 2 |Z_v|^4`, and `<= nu_s^(-1) |Z_v|` if the
  `mu_4`-shifts in (3.3) vanish; here `Z_v = prod zeta_j^(v_j^+) conj(zeta_j)^(v_j^-)` and the left
  side is the part of `m_(v,s)` from `q = +-1 mod 8` (by Theorem 1.1(1), the moduli `q = +-3 mod 8`
  can be made coherent with arbitrary `zeta_j`, so only this part matters). The proof is Lemma 4.2
  below with `delta = 0`.
* By Theorem 3.2, `|Z_v| >= e^(w(v)/2)`. So the certificate is never better than the trivial genuine
  bound `m_(v,s) <= nu_s^(-1) e^(w(v)/2)` (`sharp.md` Prop. 1.2), which is the height condition with
  zero slack. The completion sum `sum_v m_v e^(-w(v)/2)` (`sharp.md` §6.1) needs slack
  `e^(-||v||_1 log(4er))`, so it can **never** be certified by coherence. Whether it holds then
  depends on the actual `X`-smooth parts of the integers `D_s(Z_v)`, i.e. on the genuine-type
  arithmetic studied in `p3G.md` (where `zeta_j = pi_j` collapses the class to actual clusters).

*Remarks 3.4.*
* Theorem 3.2 is the rigorous form of `sharp.md` §6.6's "pseudosquare" heuristic, restricted to the
  moduli `q = +-1 mod 8`, which by Theorem 1.1 are the only ones where it can bite.
* Unconditionally, the least prime in the class `sigma` is only known to be `<= d_L^(c)`
  (Lagarias-Montgomery-Odlyzko 1979), i.e. polynomial in `n`. That gives a threshold `X >= n^(c')`;
  since `p_j > X` in every construction, the unconditional version says nothing useful, and we do
  not state it.
* The theorem concerns *exact* coherence at all `q = +-1 mod 8` below `X`. Section 4 relaxes it.

---------------------------------------------------------------------------------------------

## 4. Twisted coherence: the best certified repair of the parity, and its quadratic cost

Theorem 3.2 forbids exact coherence with small `zeta_j`. The natural repair keeps the coherent
labels and corrects the parity by a *fixed* non-square, chosen coherently from a short list.

**Definition 4.1 (twisted coherent labels).** Let `ell_0 = 2 < ell_1 < ell_2 < ...` be `2` followed
by the primes `= 1 mod 4`, with Gaussian primes `rho_0 = 1 + i`, `rho_t | ell_t`. Let `zeta_j` be
Gaussian primes over distinct primes `p'_j = 1 mod 4` with `p'_j > X` and `p'_j` not in the list
`{ell_t : t <= T(X)}` below. For an odd prime `q <= X` not dividing `N`:
* `t(q)` is the least `t` with `(ell_t/q) = -1`. It exists (Dirichlet: some prime `= 1 mod 4` lies
  in a non-residue class mod `q`), and `t(q) = 0` iff `q = +-3 mod 8`. Put `T(X) = max_(q<=X) t(q)`.
* `delta_j(q) := [log(zeta_j/conj zeta_j) != nu_j(q) mod 2]` (the parity mismatch, `sharp.md`
  Prop. 2.1(3)).
* `ell_j(q) := log(zeta_j/conj zeta_j) + delta_j(q) log(rho_(t(q))/conj rho_(t(q)))` in `Z/m_Q`.

Since `(ell_t/q) = -1`, `log(rho_t/conj rho_t)` is odd, so `ell(q)` is parity-respecting: these are
factored residue data with genuine norms, i.e. admissible in P3-all\*.

**Lemma 4.2 (collision criterion).** For `v in Z^r \ {0}`, put `Z_v = prod_j zeta_j^(v_j^+) conj(zeta_j)^(v_j^-)`,
`d_q(v) = <v, delta(q)> in Z`, and `W_(t,d)(v) = Z_v rho_t^d` (with `rho^d := conj(rho)^(|d|)` for
`d < 0`). Then for every level `a <= A_q` and unit `s`:

```text
v collides at q^a with i^s     <=>     q^a  |  D_s( W_(t(q), d_q(v))(v) ),                                         (4.1)
```

where `D_s(X + iY) = Y, X - Y, X, X + Y` for `s = 0, 1, 2, 3` (`sharp.md` Prop. 1.2). Moreover
`D_s(W_(t,d)(v)) != 0` for all `t, d, s`.

*Proof.* `<v, ell(q)>` is the log of `(Z_v/conj Z_v)(rho/conj rho)^d = W/conj W`, `W = W_(t(q), d_q(v))`.
All `zeta_j`, `rho_t` are units mod `q` (`q <= X < p'_j`, and `ell_t != q`), so `W` is a unit mod
`q`, and `W/conj W = i^s mod q^a` iff `q^a | W - i^s conj W`. As in `sharp.md` Prop. 1.2,
`W - i^s conj W` is `D_s(W)` times `2i`, `1 - i`, `2` or `1 + i`, which are units mod odd `q`, and
a rational integer divisible by `q^a` in `Z[i]` is divisible by `q^a` in `Z`. If `D_s(W) = 0` then
`W/conj W = i^s`, a unit; but `zeta_j/conj zeta_j` has valuation `+1` at `zeta_j` and `-1` at
`conj zeta_j` (`p'_j` split), the `p'_j` are distinct and differ from `ell_t`, so the valuation of
`W/conj W` at `zeta_j` is `v_j`, forcing `v = 0`. QED

**Theorem 4.3 (quadratic caps).** In the setting of Definition 4.1, let `v != 0`, `k = ||v||_1`.
Then for every unit `s`:

```text
log m_(v,s)  <=  4 log|Z_v| + log 2                                       (moduli q = +-3 mod 8)
               + sum_(1 <= t <= T(X)) sum_(|d| <= k) [ log(sqrt2 |Z_v|) + (|d|/2) log ell_t ]      (q = +-1 mod 8)
             <=  (T(X)(2k + 1) + 4) log(sqrt2 |Z_v|) + (T(X) k(k+1)/2) log ell_(T(X)),                              (4.2)
```

with `log|Z_v| <= (k/2) max_j log p'_j`. The first line alone is linear in `k`; the quadratic term
comes only from `q = +-1 mod 8`.

*Proof.* Group the moduli `q` by `(t(q), d_q(v))`; `|d_q(v)| <= sum_j |v_j| delta_j(q) <= k`.
* For `t >= 1`, the prime powers `q^(a_(q,s)(v))` with `(t(q), d_q(v)) = (t, d)` are pairwise
  coprime divisors of the nonzero integer `D_s(W_(t,d)(v))` (Lemma 4.2), whose absolute value is at
  most `nu_s^(-1) |W_(t,d)| <= sqrt2 |Z_v| ell_t^(|d|/2)`.
* For `t = 0` (`q = +-3 mod 8`), `rho_0/conj rho_0 = i`, so `W/conj W = (Z_v/conj Z_v) i^d` and
  (4.1) reads `q^a | D_(s - d)(Z_v)`. Grouping by `s - d mod 4`, the `t = 0` part of `m_(v,s)` divides
  `prod_(s' in Z/4) D_(s')(Z_v)` (each factor nonzero by Lemma 4.2 with `d = 0`), of absolute
  value `<= 2 |Z_v|^4`.
Multiply over groups, and sum `sum_(|d|<=k) |d|/2 = k(k+1)/2`. QED

**Corollary 4.4 (range certified by twisted coherence).** Let all column weights be `>= tau`, let
`max_j log p'_j <= h log M`, and write `T = T(X)`, `lam = log ell_T`. Then the per-unit height
condition `m_(v,s) <= nu_s^(-1) e^(w(v)/2)` holds for every `v` with

```text
1 <= ||v||_1 = k  <=  k_1 := ( tau - (T + 4)(h log M + 0.7) - T lam ) / ( T (2 h log M + 1.4 + lam) ).            (4.3)
```

Nothing is certified for `k > k_1`. (Completion: since (4.2) does not depend on `s`,
`m_v <= prod_s m_(v,s) <= e^(4 cap(k))`, `cap(k)` the right side of (4.2). The `(4er)^k` monomials of
norm `k` then contribute `<= 4 e^(-k)` to `sum_v min(1, 4 m_v e^(-w(v)/2))` whenever
`4 cap(k) <= k tau/2 - k(log(4er) + 1)`, i.e. for `k <= k_1'` := (4.3) with `tau` replaced by
`tau/4 - (log(4er) + 1)/2`.)

*Proof.* By (4.2) with `log(sqrt2 |Z_v|) <= (k h/2) log M + 0.35` and `w(v) >= k tau`, it suffices that
`(T(2k+1) + 4)(k h log M + 0.7)/2 + T k(k+1) lam/2 <= k tau/2`. Dividing by `k/2` and using
`0.7/k <= 0.7`, it suffices that `(2Tk + T + 4)(h log M + 0.7) + T(k+1) lam <= tau`, which is (4.3)
after collecting the terms in `k`. QED

*Size of `T(X)`.* `T(X)` is the number of list primes needed to supply a non-residue to every
`q <= X`, `q = +-1 mod 8`. Under GRH (Lagarias-Montgomery-Odlyzko, least prime in a Chebotarev class
of `Q(i, sqrt q)`), `ell_(T(X)) <= c (log X)^2`, so `T(X) = O((log X)^2/loglog X)`; unconditionally
the Lagarias-Montgomery-Odlyzko bound only gives `ell_(T(X)) <= X^(c)` for an absolute `c`. Even with `T = O(1)` (true for a positive
proportion of `q`, false for all `q`), (4.3) gives `k_1 = O(tau/log M)`: with `tau = C log M`,
twisted coherence certifies only monomials of norm `O(C)`, while the completion requires every
norm up to `~ 2.08 X/tau`.

**Remark 4.5 (heuristic: twisted coherence lands at the random rate).** For `k > k_1` the cap (4.2)
is not the binding bound; the actual demand is decided by which primes `q` divide the integers
`D_s(W_(t,d)(v))` *and* have `d_q(v) = d`. If these events behave like independent events of
probability `~ 4/q` (the model of `sharp.md` Prop. 6.5), the hits of the worst `v` of norm `k`
are the Poisson tail `~ k log M/log(k log M)` truncated at the cap `~ T h k^2`. The two meet at
`k_* ~ log M/(T h loglog M)`, where the per-unit demand is `~ log^2 M/loglog M`. So the
per-unit budget `tau/2` must be `>~ log^2 M/loglog M`: the same rate as the uniform random fibre
(Prop. 6.5) and as genuine labels (`p3G.md`). (Heuristic: the independence of the divisibility events
across the `2k + 1` values of `d` is assumed.) Evidence (§6) agrees at small sizes: twisted labels
behave like coherent ones for `k = 1` and like random ones for `k = 2`.

---------------------------------------------------------------------------------------------

## 5. The proposed mechanism, step by step, and what is left

The target mechanism was: *(M1)* low pair demand (needed by near-uniform margins) forces incoherent
prime residues; *(M2)* incoherent residues violate the per-unit height condition of some low-weight
monomial, by a pigeonhole/entropy count of the `e^(Theta(M))` residue configurations against the
number of low-weight monomials.

### 5.1 (M2) as a counting argument cannot work

**Lemma 5.1 (no forced collisions).** Fix an odd prime `q <= X` not dividing `N`, a parity vector
`nu`, and a finite set `V` of vectors `v` with `gcd(v_1, ..., v_r, m_q/2) = 1`. If `|V| < m_q/8`,
some parity-respecting `ell(q)` has no level-1 collision (with any unit) for any `v in V`.

*Proof.* For `ell` uniform on the parity class, `<v, ell>` is uniform on the `m_q/2` residues of
parity `<v, nu>` (proof of Theorem 1.1(2)). The four elements of `(m_q/4) Z/m_q` are all even if
`q = +-1 mod 8` and split two/two if `q = +-3 mod 8` (`m_q/4` odd). So each `v` collides with
probability `<= 8/m_q`. Union bound. QED

So at each modulus the parity forces *nothing*: the only forced collisions are pigeonhole ones,
which already exist without parity, and parity changes their count by at most a factor 2
(`p3lower.md` Prop. B.3). Counting over all moduli at once, the forced collisions are the ones of
`p3lower.md` Prop. B.5, which give only `W >= r log min(r, X) - O(r)`, a consequence of the prime
number theorem. The reason is structural: counting shows that collisions exist **on average**
(about `4 log X` of demand per monomial, within the budget `||v||_1 tau/2` at `tau = Theta(log M)`). Growth needs
**concentration**: some `v` of norm `k` hit with demand `>> k log M`. Random designs have it (the
Poisson tail, `sharp.md` Prop. 6.5), but no counting argument forces it, because the per-modulus
choices can be correlated across moduli at will (Lemma 5.1 is per modulus; Theorem 2.1 shows the
pair data leave `r - M` coordinates per modulus completely free).

### 5.2 (M1) is false for the free labels, and what survives

By Corollary 2.2, on twin profiles low pair demand constrains only the twin labels
`Lambda = kappa - F`. What survives of (M1) is:

> *(M1') On a twin profile, the twin labels `Lambda_y(q) = kappa_y(q) - F_y(q)` must avoid, at all but
> `O(||w||_1 + 1)` moduli, the values `-+ <w, (ell_J, U)(q)> + (m_q/4)Z` for every sparse free vector
> `w` (by (2.3)), while `kappa = F + Lambda` must be a point code with maximal pair agreement
> `O(log M)` (near-uniform margins, `p3lower.md` Prop. B.1).*

**The fate of every design we can analyse.** In each line `tau` is the per-column weight needed; growth
would require every design to need `tau >> log M`.

| design (labels per modulus) | in P3-all\*? | pairs (near-uniform) | monomials | `tau` needed | status |
|---|---|---|---|---|---|
| additive: `ell = g_q lambda`, `|lambda| <= M^h` | yes if parities fit | -- | Siegel relation collides everywhere | `W >= (1/4 - o(1)) M X/(h+2)` | **proved** (Prop. 3.1) |
| coherent `log(zeta_j/conj zeta_j)`, small `zeta_j` | **no** (parity) | dense coherent pair integers | linear caps | `O(log M)` for monomials | proved (caps); excluded by Thm 3.2 (GRH) |
| coherent with right parity | yes | as above | caps `>= e^(w/2)`: no slack | `-` | **proved under GRH** (Thm 3.2, Cor. 3.3) |
| genuine `pi_j` | yes | genuine pair integers | height automatic, completion fails at poly sizes | `~ log^2 M/loglog M` | `p3G.md` (verified sets, model theorem) |
| random fibre `D_A` | yes | `O(log M)` | Poisson tail | `(2 + o(1)) log^2 M/loglog M` | **proved** (`sharp.md` Thm 6.4, Prop. 6.5) |
| twisted coherent (Def. 4.1) | yes | dense coherent | quadratic caps (Thm 4.3) | `~ log^2 M/loglog M` | caps proved; rate heuristic (Rem. 4.5) |
| `D_A` points + coherent or twisted free labels | yes | `O(log M)` | free part: linear/quadratic caps; twin part `Lambda = kappa - F` random-like | `~ log^2 M/loglog M` | heuristic (below) |

For the last line: by (2.3) a twin monomial `+-e_(f(y)) + w` collides at `q` iff `Lambda_y(q)` hits
one coherent value. With `kappa` from `D_A` (additive in `y`) and `F` coherent and dense,
`Lambda_y(q)` has no relation to the coherent values, and if the hits are independent with
probability `~ 4/q` the `4M (2r)^(k-1)` such monomials of norm `k` have the Poisson tail
`~ k log M/loglog M` hits. Their budget is `w(v)/2 <= k (tau + 1)/2`, so again
`tau >~ 2 log^2 M/loglog M` for **all** columns (the free part of the monomial carries `k - 1` of
the `k` units of norm). Making the twin columns alone heavy does not help.

### 5.3 What a growth proof would have to contain

Every analysed design fails for one of three reasons, and a growth theorem for (near-uniform
members of) P3-all\* would have to show that **every** design fails for one of them:

* **(I1) Inverse step.** If every sparse monomial has per-unit demand `<= beta ||v||_1 log M` for
  all norms up to the completion range `~ 2.08 X/tau`, then the labels are (close to)
  *certified*: induced by a global object of controlled height. No such inverse theorem is known,
  and Lemma 5.1/Theorem 2.1 show that it cannot come from counting. It is a statement about
  concentration over `e^(Theta(k log M))` monomials at once, of discrepancy type.
* **(I2) Certified implies genuine size.** Proved for exact Gaussian coherence under GRH
  (Theorem 3.2), and for additive certification unconditionally (Proposition 3.1). Approximate
  coherence (Definition 4.1) costs a factor `k` (Theorem 4.3).
* **(I3) Genuine size implies Poisson tails.** This is the arithmetic of dense coherent integers
  (pair integers `B_xy`, twin values): Hypothesis `H` of `p3lower.md` Theorem C, or the smooth-part
  statistics of `p3G.md`. It is open even for genuine primes.

(I2) is the only step proved here. (I1) is the missing combinatorial core. A design defeating it
would be a family of parity-respecting labels, not induced by any global object, whose sparse
combinations have bounded agreement with `mu_4` at all scales: a "multiplicative code" with
prescribed (Legendre) parities. We found none, and we found no argument excluding it.

---------------------------------------------------------------------------------------------

## 6. Checks and evidence

**Exact checks** (`check_theorems.py`, output `check_theorems_out.txt`: `ALL OK`).
* Theorem 1.1(1): for 12 primes `q = +-3 mod 8` up to 67 and 200 random label/parity/monomial
  triples each, adding `iota u` permutes the collision unit by `<v, u>` exactly as in (1.1).
* Theorem 2.1: on a random twin profile (`M = 5`, `r = 14`, random signs `eps_x`), the identities
  `kappa_y = F_y + Lambda_y` and (2.2) hold for 300 random label vectors and monomials.
* Lemma 4.2 / Theorem 4.3: for all odd `q <= 200`, five twisted-coherent columns (Definition 4.1,
  `zeta_j` over the first split primes above 200, genuine `p_j` of size about `10^6`), and all
  `5^5 - 1` monomials with entries in `[-2, 2]`, every level-1 collision with unit `i^s` holds iff
  `q | D_s(Z_v rho_(t(q))^(d_q(v)))`, and `D_s != 0` in every case (18140 collisions tested). The
  parity of every twisted label equals `[(p_j/q) = -1]`.
* `tgroup.py` self-test: discrete logs in `T_q` for 14 primes; `log i` is `m/4` or `3m/4`; and
  `log(z/conj z)` is even iff `(Norm z/q) = 1` for all `z = a + bi`, `0 <= a < 12`, `1 <= b < 12`.

**Evidence** (`designs.py`, output `designs_out.txt`; no proof uses it). Columns `r = M`, moduli
`q in [M^(1/2), 3M]`, level 1, genuine `p_j ~ M^3`. Entry: `max_v max_s` (per-unit demand)`/log M`
over monomials of norm `k` (`k = 3` only for `M <= 200`).

| `M` | RAND `k=1,2,3` | COH (not in class) | TWIST (in class) | GEN (in class) |
|---|---|---|---|---|
| 50  | 4.73, 5.92, 6.17 | 1.63, 2.93, 4.89 | 2.37, 5.04, 6.77 | 2.51, 4.56, 6.38 |
| 100 | 3.37, 5.58, 7.32 | 1.55, 3.09, 4.68 | 2.36, 5.39, 6.87 | 2.60, 5.28, 6.31 |
| 200 | 4.68, 6.92, 7.31 | 1.51, 3.07, 4.58 | 1.77, 5.53, 7.32 | 2.78, 5.11, 7.07 |
| 400 | 4.71, 6.37, - | 1.48, 2.97, - | 2.08, 5.30, - | 3.03, 5.48, - |
| 800 | 4.44, 6.59, - | 1.45, 2.90, - | 2.20, 6.19, - | 3.00, 6.17, - |

Reading (evidence only): coherent labels have demand `~ 1.5 k log M`, linear in `k`, as certified
(their `zeta_j` have norm `~ X`). Twisted labels match coherent ones at `k = 1` (the cap of
Theorem 4.3 has only `d in {-1, 0, 1}`) but are already at the random level at `k = 2`, in line with
Remark 4.5. Genuine labels sit with the random ones. The sizes are far from asymptotic
(`log M <= 6.7`), so no rate is read off.

---------------------------------------------------------------------------------------------

## 7. Open points

1. **(I1), the inverse step.** Is there a parity-respecting label system (any profile, `X = 3M`)
   whose sparse monomials all have per-unit demand `<= beta ||v||_1 log M` up to norm `~ X/log M`,
   with `beta` fixed, that is not certified? A "no" (an inverse theorem) plus (I2), (I3) would give
   growth for near-uniform members. A "yes" combined with Theorem 2.1 and a pair design for `kappa`
   would be the first step towards `W_min^all = O(M log M)`; the twin channel (M1') would remain.
2. **Twin channel (M1').** Even with ideal free labels, the `M` twin labels `Lambda = kappa - F`
   must avoid all sparse free combinations while `kappa` is a pair code. Is there a statement of the
   form "a point code with pair agreement `O(log M)`, translated by a dense combination of labels
   with bounded sparse agreement, has a coordinate with unbounded sparse agreement"?
3. **Unconditional Theorem 3.2.** Only GRH gives the threshold `X >~ log^2 n`. Unconditionally the
   least prime in the needed Chebotarev class is only known to be `<= n^(c)`, which is useless here
   because `p_j > X`.
4. **Non-twin profiles.** Theorem 2.1 uses twins. Without twins the point classes `kappa = A ell`
   still determine all pair data, but the free/twin split is replaced by a choice of complement to
   the row lattice. The decoupling of a rank-`(r - M)` block of free labels persists whenever `A`
   has a unimodular `M x M` minor (same proof). Profiles without such a minor were not examined.

## Files (`round7/p3Dlower_checks/`)

| file | content |
|---|---|
| `tgroup.py` | exact arithmetic and discrete logs in `T_q`; self-tests in `check_theorems.py` |
| `check_theorems.py`, `check_theorems_out.txt` | exact checks of Theorem 1.1(1), Theorem 2.1, Lemma 4.2/Theorem 4.3 (`ALL OK`) |
| `designs.py`, `designs_out.txt` | evidence table of §6 (`python3 designs.py M r 3 1 kmax 0.5`) |
| `notes.md` | working notes (not part of the write-up) |
