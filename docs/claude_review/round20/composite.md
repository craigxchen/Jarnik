# Round 20 — composite-supported analytic inputs: classification

*Lens: the one unclassified class of the barrier map (`round10/README.md` §7.2; `lit18/README.md` §2.2 item 3;
`overview.md` §5.5): composite-supported analytic inputs. Scripts: `composite_checks/`. Everything here is the author's
and **unrefereed**.*

## 0. Summary and verdict

**Verdict: the class is capped at the rate `log R/loglog R`, with the same constants and the same caveat as
`round10/analytic.md` Thm B. Composite support adds no separating power. No growth, no uniform bound.**

1. **Two composite grafts** (Def. 1.1). The *minimal graft* changes the angles of the S-smooth ideals only. The
   *multiplicative graft* rotates each S-fibre `{d beta : beta coprime to S}` by the constant `rho(d)` and keeps complete
   multiplicativity. Both agree with the fake on all its points, primes and divisors of its norm.
2. **The all-ideal large sieve has no sharp form** (Prop. 1.5, proved): its optimal constant is
   `>= max(#I(X), (2K+1) floor(X^(1/2)))`, because rational multiples share an angle. So every true all-ideal Hecke large
   sieve is automatically robust; the round-10 "sharp large sieve" problem has no all-ideal analogue.
3. **Theorem M (proved).** The minimal graft satisfies every `(c, gamma)`-robust composite input (large sieve over any
   support; twisted sums with any weights of divisor-bounded size on S-smooth ideals, which covers Type I/II sums,
   divisor sums, Selberg/GPY sieve weights; weighted polar-box counts with any support) up to `1 + O(M^(1-A/2) + 2^B M^(1-A gamma))`.
   The reason: there are at most `Psi_S(X) <= X^((1+o(1))/A)` S-smooth ideals in `I(X)`.
4. **Theorem P (proved).** If the argument also uses multiplicativity across S-smooth and S-coprime factors, the
   multiplicative graft satisfies: the all-ideal large sieve over any support (P1, unconditional; this is the "dual
   large-sieve argument" the round-10 referee asked for); every multilinear (Type II) bound with free coefficients,
   exactly (P2); Halasz-type bounds (P2'); sharp polar-box counts with error `>= c x^(1/2)` (P3, unconditional);
   divisor-weighted boxes with error `>= c x^(3/4)` (P4); and every family uniform in sector position and height (P5).
   The only quantities that enter are the fibre masses `tau_sigma = O(s P_0^(-sigma))`: each fibre is a rotated genuine
   sub-ensemble, which no position-uniform input can detect.
5. **Corollary 2.2 (the cap).** Grafts of the `round8/lw.md` 6.1 fakes (`W <= (2A + o(1)) M log M`) and of the P3-char
   fakes (`W <= (bA + o(1)) M log M`) satisfy all of these inputs when `A > max(2, 1/gamma)`, `A >= 4`. Proposal: move
   the class from `overview.md` §5.5 to the capped table §5.4 (kind F), with row 7's R1 caveat.
6. **Non-robust composite inputs collapse.** Smooth all-ideal Hecke sums are `O(X^(-B))` (true; Poisson) and are
   violated by both grafts. But Proposition 3.1 (proved) shows that, with multiplicativity, smooth decay at a single
   `k = +-1` forces `phi_p = +-theta_p (mod pi/2)` for every `p`: the fake **is** an actual cluster. This is the composite
   Theorem D: a separating statement that is the problem itself. The only other non-robust composite input, the sharp
   large sieve on primitive support, is (transferred from the round-10 referee) equivalent to a composite max-scale
   cross-gap. That is a badly-approximable condition and joins the existing unclassified sharp-large-sieve item.
7. **Exact strength requirement** (§5). To defeat every fake family, a composite input must have error
   `X^(gamma(c_0))` with `gamma(c_0) < 2/c_0`, uniformly in position, at heights `>= M^(c_0/2)`. For sharp counts this is
   false (rational edges give `X^(1/2)`); for smooth sums it is true but collapses (Prop. 3.1). **No known theorem
   reaches a fake-violated, non-collapsing composite statement.**
8. **Computations** (§4; evidence/exact, no proof uses them). Large-sieve constants change by `<= 3.2%` under the
   multiplicative graft and `<= 0.09%` under the minimal graft. Smooth sums are `~1e-12` genuine against `0.03-3` grafted,
   and linear in `epsilon` for near-genuine angles. Polar-box differences are `<= 0.0008` of the proved P3 bound.
   `Psi_S` matches `binom(2s+L, L)`.

## 1. Setting, notation, and four lemmas about S-smooth ideals and lattice points

All statements in §§1-3 are **proved here (unrefereed)** unless labelled otherwise.

### 1.1 The genuine angle function and the fake data

* `I` is the monoid of nonzero ideals of `Z[i]`; `I(X) = {alpha in I : N alpha <= X}`. Each ideal has exactly one
  generator `x + iy` with `x > 0, y >= 0`.
* The **genuine angle function** is `vartheta : I -> R/(pi/2)Z`, `vartheta(alpha) = arg(generator) mod pi/2`. It is
  completely additive (`vartheta(alpha beta) = vartheta(alpha) + vartheta(beta)`), `vartheta(conj alpha) = -vartheta(alpha)`,
  and `vartheta(n) = 0` for `n in Z`. The Hecke Grossencharacters of `Q(i)` are `lambda_k(alpha) = e^{4ik vartheta(alpha)}`.
* **Fake data** (as in `round10/analytic.md` §1): a finite set `S` of rational primes `p = 1 mod 4`, `s = |S|`,
  `P_0 = min S`; for each `p in S` the genuine `pi_p` with `theta_p = vartheta(pi_p)`, and a real **fake angle** `phi_p`.
  An ideal is **S-smooth** if all its prime ideals lie above `S`. For S-smooth
  `d = prod_p pi_p^(u_p) conj(pi_p)^(w_p)` put
  `vartheta^S(d) = sum_p (u_p - w_p) phi_p` and `rho(d) = vartheta^S(d) - vartheta(d) = sum_p (u_p - w_p)(phi_p - theta_p)`.
  For an angle model (fake circle) with primes `S` and angles `phi`, the fake points
  `z_x = prod pi^(a_x) conj(pi)^(e - a_x)` are S-smooth and `vartheta^S(z_x) = <2a_x - e, phi>` is exactly the fake
  angle of `round4/outside.md` (5.1); likewise for every divisor of `N`.
* **Fibre decomposition.** Every `alpha in I` factors uniquely as `alpha = alpha_S alpha'` with `alpha_S` S-smooth and
  `alpha'` coprime to `S`. The **S-fibre** of an S-smooth `d` is `F_d = {d beta : beta coprime to S}`.

**Definition 1.1 (two composite grafts).**
* The **multiplicative graft** is `vartheta^x(alpha) = vartheta^S(alpha_S) + vartheta(alpha')`. Equivalently, on the fibre
  `F_d` it is the genuine angle function rotated by the constant `rho(d)`; the fibre `F_1` (ideals coprime to `S`) is
  unchanged. It is completely additive, conjugation-equivariant and zero on `Z` (a rational `n` has
  `n_S = prod (pi_p conj pi_p)^v`, so `rho(n_S) = 0`). Its characters are `lambda_k^x = lambda_k * u_k`, where
  `u_k(alpha) = e^{4ik rho(alpha_S)}` is a completely multiplicative unimodular function that is `1` off the S-smooth part.
* The **minimal graft** is `vartheta^min(alpha) = vartheta^S(alpha)` if `alpha` is S-smooth and `vartheta(alpha)` otherwise.
  It changes only the S-smooth ideals. It is conjugation-equivariant, zero on `Z`, and additive on S-smooth ideals and
  on ideals coprime to `S` separately, but not across mixed products.

Both agree with the prime graft of `round10/analytic.md` Def. 3.1 on prime ideals (mod `pi/2`,
`vartheta^S(conj pi_p) = -phi_p = pi/2 - phi_p`), and both agree with the fake on every S-smooth ideal, in particular on
the fake's points, primes and all divisors of its norm. So a statement about the ensemble that an argument applies to
"the angle function whose values at the cluster's Gaussian integers are the cluster angles" makes sense for either graft.
The multiplicative graft is the referee's "natural extension" (`round10/referee_analytic.md` Slip 3).

### 1.2 Counting S-smooth ideals

**Lemma 1.2.** Let all primes of `S` be `>= P_0`.
1. For `sigma > 0` with `4 s P_0^(-sigma) <= 1`:
   `tau_sigma := sum_{d != 1 S-smooth} N d^(-sigma) = prod_{p in S}(1 - p^(-sigma))^(-2) - 1 <= 8 s P_0^(-sigma)`.
   More generally, for `B >= 0` and `D_B = sum_{a >= 1}(a+1)^B 2^(1-a)`, if `2 s D_B P_0^(-sigma) <= 1` then
   `tau^(B)_sigma := sum_{d != 1} d(d)^B N d^(-sigma) <= 4 s D_B P_0^(-sigma)` (`d(d)` = number of ideal divisors).
2. `Psi_S(X) := #{d S-smooth : N d <= X} <= binom(2s + L, L)` with `L = floor(log X / log P_0)`, hence
   `Psi_S(X) <= max(1, X^kappa)`, `kappa := log(4 e s)/log P_0`. Every S-smooth `d` with `N d <= X` has `Omega(d) <= L`
   and `d(d) <= 2^L <= X^(log 2/log P_0)`. For `X < P_0` the only S-smooth ideal in `I(X)` is `1`.

*Proof.* (1) The Euler product is over the `2s` prime ideals above `S`, each of norm `p >= P_0`. With
`t = p^(-sigma) <= 1/4`, `-log(1 - t) <= t/(1 - t) <= 2t`, so `log prod <= 4 s P_0^(-sigma) =: u <= 1` and
`e^u - 1 <= (e - 1)u <= 2u`. The weighted case: `sum_{a>=0}(a+1)^B t^a <= 1 + t D_B` for `t <= 1/2`, then the same.
(2) An S-smooth `d = prod pi^(v_pi)` with `N d <= X` has `sum v_pi log p <= log X`, so `sum v_pi <= L`; the number of
`v in N^(2s)` with `sum v <= L` is `binom(2s + L, L)`. If `L <= 2s`, `binom(2s+L, L) <= (e(2s+L)/L)^L <= (4es)^L`; if
`L > 2s`, `binom(2s+L, L) <= 2^(2s+L) <= 4^L`. Either way `<= (4es)^L <= X^kappa` since `L <= log X/log P_0`.
`Omega(d) = sum v <= L` and `d(d) <= 2^Omega`. QED

For the fakes used below `s <= 2M` (or `s = b(M-1)`) and `P_0 >= M^A`, so `kappa = (1 + O(1/log M))/A` and
`tau_sigma = O(M^(1 - A sigma))`.

### 1.3 Lattice points in sectors and the all-ideal large sieve

**Lemma 1.3 (polar boxes).** Let `I` be an arc of `R/(pi/2)Z` of length `omega <= pi/2` (open, closed or half open)
and `y >= 1`. Then `#{alpha in I(y) : vartheta(alpha) in I} = omega y/2 + E`, `|E| <= 2.53 y^(1/2) + 1.58`. For
`1 <= y_2` and `0 <= y_1 <= y_2`: `|#{y_1 <= N alpha <= y_2, vartheta(alpha) in I} - omega(y_2 - y_1)/2| <= 8.3 y_2^(1/2)`.

*Proof.* Lift `I` to `[c, c + omega]` with `c in [0, pi/2)`. Choosing for each ideal its generator with argument in
`[c, c + pi/2)` identifies the ideals counted with the nonzero lattice points of the closed sector
`K = {|z|^2 <= y, arg z in [c, c+omega]}` (up to the arc's endpoints, which the same argument handles). `K` is convex
(`omega < pi`), of area `omega y/2` and perimeter `P = (2 + omega) y^(1/2) <= 3.571 y^(1/2)`. With `r = 2^(-1/2)` and the
unit squares `z + [-1/2, 1/2]^2`: these are disjoint, contained in the outer parallel body `K_(+r)` for `z in K`, and
they cover the inner parallel body `K_(-r)`, each point of which is within `r` of a lattice point lying in `K`. So
`A - rP <= area(K_(-r)) <= #(Z^2 cap K) <= area(K_(+r)) = A + rP + pi r^2`. Subtract `1` for the origin:
`|E| <= rP + pi/2 <= 2.53 y^(1/2) + 1.58`. The annular count is the difference of two such counts (the inner one is `0`
with area `< pi/4` when `y_1 < 1`), so the error is `<= 2(2.53 + 1.58) y_2^(1/2) <= 8.3 y_2^(1/2)`. QED

**Lemma 1.4 (genuine all-ideal large sieve, upper bound).** For `Y >= 1`, `K >= 1` and complex `b_k`, `|k| <= K`:
`sum_{alpha in I(Y)} |sum_k b_k lambda_k(alpha)|^2 <= C_1 (Y + K Y^(1/2)) ||b||^2`, with `C_1 = 56`.

*Proof.* Put `t_alpha = 4 vartheta(alpha)/(2 pi) in R/Z`, `T(t) = sum_k b_k e(kt)`, `h = |T|^2`. For any `t` and the arc
`J_t` of length `delta` centred at `t`, `h(t) <= delta^(-1) int_(J_t) h + int_(J_t) |h'|` (average the identity
`h(t) = h(u) + int_u^t h'` over `u in J_t`). Summing over `alpha`, every `u` lies in at most `Q(delta)` arcs `J_(t_alpha)`,
where `Q(delta)` is the maximal number of points `t_alpha` in an arc of length `delta`. With
`int h = ||b||^2` and `int |h'| <= 2 ||T||_2 ||T'||_2 <= 4 pi K ||b||^2`,
`sum_alpha h(t_alpha) <= Q(delta)(delta^(-1) + 4 pi K) ||b||^2`. Take `delta = 1/K`; an arc of length `1/K` of `t`
is an arc of length `pi/(2K) <= pi/2` of `vartheta`, so by Lemma 1.3, `Q <= pi Y/(4K) + 2.53 Y^(1/2) + 1.58`. Hence the
sum is `<= (1 + 4 pi)(pi Y/4 + 4.11 K Y^(1/2)) ||b||^2 <= 56 (Y + K Y^(1/2)) ||b||^2` (using `1.58 K <= 1.58 K Y^(1/2)`). QED

**Proposition 1.5 (the all-ideal large sieve has no sharp form).** Let `Lambda*(X, K)` be the optimal constant in
`sum_{|k| <= K} |sum_{alpha in I(X)} a_alpha lambda_k(alpha)|^2 <= Lambda ||a||^2`. Then for `X >= 1`, `K >= 0`:
`max(#I(X), (2K + 1) floor(X^(1/2))) <= Lambda*(X, K) <= C_1 (X + K X^(1/2))` (upper bound for `K >= 1`).
The same lower bound holds for every support `T subset I(X)` that contains `I(X)` (first term) or contains the rational
integers `n <= X^(1/2)` (second term). Consequently every **true** all-ideal Hecke large-sieve inequality has
`Lambda >= c_* (X + K X^(1/2))` with an absolute `c_* > 0` (`c_* = 1/4` for `X >= 100`): it is automatically
`(c_*, gamma)`-robust in the sense of `round10/analytic.md` Def. 3.3 (R1).

*Proof.* Upper bound: Lemma 1.4 and duality (the operator `a -> (sum_alpha a_alpha lambda_k(alpha))_k` and its adjoint
have the same norm). Lower bounds: `a = 1` on `I(X)` and `k = 0` alone gives `#I(X)^2 = #I(X) ||a||^2`;
`a = ` indicator of `{n in Z : 1 <= n <= X^(1/2)}` has `lambda_k(n) = 1` for all `k`, so the left side is
`(2K+1) floor(X^(1/2))^2 = (2K+1) floor(X^(1/2)) ||a||^2`. For `X >= 100`, `#I(X) >= pi X/4 - 2.53 X^(1/2) - 1.58 >= X/2`
and `floor(X^(1/2)) >= X^(1/2)/2`, so `Lambda* >= (X/2 + (2K+1) X^(1/2)/2)/2 >= (X + K X^(1/2))/4` (max `>=` mean).
(Checked numerically at `X = 100, 200, 10^3, 10^4`, `K = 0, 1, 10, 1000`: ratio `>= 0.59`.) QED

*Meaning.* Round 10 left the **sharp** `(X + K)` large sieve over primes unclassified (`round10/README.md` §3 item 6).
Over all ideals there is no such statement: rational multiples share an angle, which forces the `K X^(1/2)` term. The
sharp form survives only for supports with at most one ideal per angle, such as primitive ideals (§3.3).

## 2. The cap: composite grafts satisfy every robust composite-supported input

### 2.1 The input classes

Fix constants `c > 0`, `gamma in (0, 1]`, `B >= 0`, `C_c >= 1`. A **`(c, gamma)`-robust composite input** is a family of
inequalities about an angle function `vartheta'` on `I`, each of one of the following forms (the support `T`, the
weights and all parameters are arbitrary unless stated):

* **(C1) all-ideal large sieve.** `sum_{|k|<=K} |sum_{alpha in T} a_alpha e^{4ik vartheta'(alpha)}|^2 <= Lambda ||a||^2`
  for all `a` supported on a set `T subset I(X)`, with `Lambda >= c (X + K X^(1/2))`. By Prop. 1.5 this is **every true**
  all-ideal large sieve (`T = I(X)`), with `c = c_*`.
* **(C2) weighted twisted sums** (Type I/II sums, divisor sums, sieve-weighted sums, Hecke sums over all ideals):
  `|sum_{alpha in T} w(alpha) f(N alpha) e^{4ik vartheta'(alpha)} - Main| <= E`, with `|f| <= 1` supported in `[1, X_1]`,
  and `E >= c X_1^gamma`.
* **(C3) weighted polar boxes** (lattice points of `Z[i]` in sectors and annuli, with composite support, divisor weights,
  sieve weights `(sum_{e | alpha} lambda_e)^2`, almost-prime or rough supports, congruence conditions):
  `L <= sum_{alpha in T, N alpha in J, vartheta'(alpha) in I} w(alpha) <= U`, with `L = 0` or `L >= c x^gamma`, and `U = infinity`
  or `U >= c x^gamma`, `x = sup J`.

In (C2)-(C3) the weights may be arbitrary, except that on **S-smooth** ideals `|w(alpha)| <= C_c d(alpha)^B`. (This covers
`1, mu, lambda, Lambda, d_k` with `k <= 2^B`, and Selberg/GPY weights `(sum_{e|alpha} lambda_e)^2` with `|lambda_e| <= 1`,
which are `<= d(alpha)^2`.) `Main`, `L`, `U` may be anything. The inputs used in Selberg/GPY-type sieves (level of
distribution in sectors, upper and lower bound sieves in polar boxes) are of the forms (C2)-(C3) for the sets and weights
they involve. Prime-supported inputs are the special case covered by `round10/analytic.md` Thm B.

### 2.2 Inputs used without multiplicativity: the minimal graft

**Theorem M (minimal-graft barrier).** Let `S`, `phi` be fake data, `kappa = log(4es)/log P_0`,
`kappa_B = (log(4es) + B log 2)/log P_0`. The minimal graft `vartheta^min` satisfies:
1. every (C1) inequality with `Lambda` replaced by `(1 + 3 t^(1/2)) Lambda`, `t = 131 s P_0^(-1/2)/c`, provided
   `kappa < 1/2` and `t <= 1`;
2. every (C2) inequality with `E` replaced by `(1 + eta) E`, and every (C3) inequality with `L, U` replaced by
   `(1 - eta) L, (1 + eta) U`, where `eta = (2 C_c/c)(4es) 2^B P_0^(-gamma)`, provided `kappa_B < gamma`.

*Proof.* Only S-smooth ideals other than `1` change (`vartheta^S(1) = 0`), and in `I(Y)` there are at most
`Psi_S(Y) - 1` of them; none if `Y < P_0` (Lemma 1.2(2)).
(1) Put `D_k = sum_{alpha in T} a_alpha (e^{4ik vartheta^min(alpha)} - e^{4ik vartheta(alpha)})`. Then
`|D_k| <= 2 sum_{S-smooth alpha != 1} |a_alpha| <= 2 Psi_S(X)^(1/2) ||a||` (Cauchy-Schwarz), so
`sum_k |D_k|^2 <= 12 K Psi_S(X) ||a||^2`. By Minkowski in `l^2(k)`, the grafted left side is
`<= (Lambda^(1/2) + (12 K Psi_S(X))^(1/2))^2 ||a||^2 = (1 + t_X^(1/2))^2 Lambda ||a||^2` with
`t_X = 12 K Psi_S(X)/Lambda <= 12 Psi_S(X)/(c X^(1/2))`. For `X >= P_0`, `Psi_S(X) X^(-1/2) <= X^(kappa - 1/2) <= P_0^(kappa - 1/2) = 4es P_0^(-1/2)`,
so `t_X <= 48 e s P_0^(-1/2)/c <= t`, and `(1 + t^(1/2))^2 <= 1 + 3 t^(1/2)` for `t <= 1`.
(2) The change of the sum in (C2) is at most `2 sum_{S-smooth alpha != 1, N alpha <= X_1} |w(alpha)|
<= 2 C_c 2^(BL) Psi_S(X_1) <= 2 C_c X_1^(kappa_B)` (Lemma 1.2(2)), and it is `0` if `X_1 < P_0`. Divided by
`E >= c X_1^gamma` this is `<= (2C_c/c) X_1^(kappa_B - gamma) <= (2 C_c/c) P_0^(kappa_B - gamma) = eta`. (C3) is the same
with `|w|` summed over the S-smooth ideals in the norm range; robustness gives `L, U >= c x^gamma` whenever they are finite
and nonzero. QED

*Remark 2.1.* Theorem M uses nothing about the form of the inequality except that it changes by at most the total
weight of the changed ideals. So it covers any functional of the angle function that is Lipschitz in the Hamming
metric with the stated weights. That includes statements with additive shifts, for example GPY-type counts of `alpha`
with all `alpha + h_j` in prescribed sets and sectors (a summand changes only if some `alpha + h_j` is S-smooth, at most
`(#shifts) Psi_S(x)` summands), and Kubilius-type statements about the joint law of prime factors and angles of a random
ideal.

### 2.3 Inputs used together with multiplicativity: the multiplicative graft

An argument may use that `lambda_k(alpha beta) = lambda_k(alpha) lambda_k(beta)` across S-smooth and S-coprime factors
(for instance to factor a sum into a bilinear form). The minimal graft violates that identity; the multiplicative graft
keeps it. Write `B(J, I) = {alpha : N alpha in J, vartheta(alpha) in I}`.

**Theorem P (multiplicative-graft barrier).** Let `S`, `phi` be fake data with `4 s P_0^(-1/2) <= 1`. The
multiplicative graft `vartheta^x` satisfies:

* **(P1) all-ideal large sieve, any support.** Every (C1) inequality with `Lambda` replaced by `(1 + 56 tau_(1/2)/c) Lambda`,
  `tau_(1/2) <= 8 s P_0^(-1/2)`. *No other hypothesis.*
* **(P2) multilinear forms with free coefficients (Type II).** Every inequality of the form
  `|sum_{m_1,...,m_r} a^(1)_(m_1) ... a^(r)_(m_r) W(N m_1, ..., N m_r) lambda_k(m_1 ... m_r)| <= B` that holds for all
  coefficient vectors in a class closed under multiplication by unimodular sequences (any `l^p` or pointwise bounds)
  holds for `vartheta^x` with the **same** `B`.
* **(P2') Halasz-type inputs.** Every inequality "for every completely multiplicative `g : I -> {|z| <= 1}`,
  `|sum_{alpha in I(X)} g(alpha) lambda_k(alpha)| <= H(X, g lambda_k)`" holds for `vartheta^x` with the same `H`. If `H` depends
  on `g` alone through a pretentious distance `D(g, .)`, the twist changes `D^2` by at most `4s/P_0`.
* **(P3) sharp polar boxes, all ideals or ideals coprime to a fixed `q` with `(q, S) = 1`, weight 1.** For every arc `I`
  with `|I| <= pi/2` and every `J subset [0, x]`,
  `|#(T cap B^x(J, I)) - #(T cap B(J, I))| <= 16.6 Pi_q (1 + tau'_S) tau_(1/2) x^(1/2)`,
  where `B^x` is the box for `vartheta^x`, `Pi_q = prod_{frak p | q}(1 + N frak p^(-1/2))` and
  `tau'_S = prod_{p in S}(1 + p^(-1/2))^2 - 1`. So every (C3) count of this kind with `gamma >= 1/2` survives with
  `eta = 16.6 Pi_q (1 + tau'_S) tau_(1/2)/c`.
* **(P4) divisor-weighted polar boxes.** For `w = d` (number of ideal divisors), `T = I`:
  `|sum_{B^x} d - sum_B d| <= 420 (1 + tau'_S) tau^(1)_(3/4) x^(3/4)`. So (C3) with `w = d` and `gamma >= 3/4` survives.
* **(P5) uniform families.** Let `w` be multiplicative with `w(1) = 1`, and let `h` be the Dirichlet inverse of `w`
  restricted to S-smooth ideals, so that `w 1_{coprime to S} = w * h`. Suppose that the genuine ensemble satisfies, **for
  all** boxes with `x >= 1` (all positions `I`, all `J subset [0, x]`),
  `|sum_{alpha in B(J,I)} w(alpha) - M_w(J, |I|)| <= C_E x^gamma` (main term depending on `I` only through `|I|`), or, for
  twisted sums, `|sum_alpha w(alpha) f(N alpha) lambda_k(alpha)| <= C_E X_1^gamma` for all `f` in a class closed under
  dilations `f -> f(n .)`. Then the same statements hold for `vartheta^x` with an extra error at most
  `2 C_E x^gamma tau^(w)_gamma (1 + tau^(h)_gamma)`, where `tau^(g)_gamma = sum_{d != 1 S-smooth} |g(d)| N d^(-gamma)`.
  If `|w(frak p^a)|, |h(frak p^a)| <= C_0 (a+1)^B`, this is `O_(C_0,B)(s P_0^(-gamma)) x^gamma`.

*Proof.*
(P1) By duality it suffices to bound `sum_{alpha in T} |sum_k b_k e^{4ik vartheta^x(alpha)}|^2`. Group by fibres:
`alpha = d beta`, `beta` coprime to `S`, and `e^{4ik vartheta^x(d beta)} = e^{4ik vartheta^S(d)} lambda_k(beta)`. So the sum is
`sum_d sum_{beta coprime to S, d beta in T} |sum_k b^(d)_k lambda_k(beta)|^2` with `b^(d)_k = b_k e^{4ik vartheta^S(d)}`,
`||b^(d)|| = ||b||`. The term `d = 1` is a sub-sum of the genuine dual sum over `T`, so it is `<= Lambda ||b||^2`. For
`d != 1`, drop the conditions on `beta` other than `N beta <= X/N d`, and apply Lemma 1.4 (the sum is empty if
`X/N d < 1`): it is `<= 56 (X/Nd + K (X/Nd)^(1/2)) ||b||^2`. Summing over `d != 1`:
`<= 56 (X tau_1 + K X^(1/2) tau_(1/2)) ||b||^2 <= 56 tau_(1/2) c^(-1) Lambda ||b||^2`.
(P2) `lambda^x_k(m_1 ... m_r) = prod_j lambda_k(m_j) u_k(m_j)` since `u_k` is completely multiplicative. So the grafted form
equals the genuine form with coefficients `a^(j)_m u_k(m)`, which lie in the same class.
(P2') `g lambda^x_k = (g u_k) lambda_k` and `g u_k` is again completely multiplicative with values in the unit disc. And
`D^2` changes only at the prime ideals above `S`, each by at most `2/p` (since `|1 - Re z| ` changes by at most `2` for
`|z| <= 1`, with weight `1/N frak p`), in total by `<= 4s/P_0`.
(P3) Write the count as `sum_d N_d(I - vartheta^S(d))` versus `sum_d N_d(I - vartheta(d))` with
`N_d(I') = #{beta coprime to qS : N beta in J/N d, vartheta(beta) in I'}` (`d` is coprime to `q`). The term `d = 1` is the
same on both sides. Inclusion-exclusion over the squarefree products `e` of prime ideals dividing `qS`:
`N_d(I') = sum_e mu(e) #{gamma : N gamma in J/(N d N e), vartheta(gamma) in I' - vartheta(e)}`. By Lemma 1.3 each count is
`|I'| |J|/(2 N d N e) + O(8.3 (x/(N d N e))^(1/2))`, and it is `0` when `x < N d N e`. The main terms do not depend on the
rotation, so `|N_d(I_1) - N_d(I_2)| <= 16.6 (x/N d)^(1/2) sum_e N e^(-1/2) = 16.6 (x/N d)^(1/2) Pi_q (1 + tau'_S)` for
`|I_1| = |I_2|`. Sum over `d != 1`.
(P4) The fibre weights split, `d(d beta) = d(d) d(beta)` for `beta` coprime to `S`. For the fibre sum
`sum_{beta coprime to S in a box} d(beta) = #{(gamma, delta) coprime to S : gamma delta in box}`, use the hyperbola
split at `y^(1/2)` (`y` the top of the fibre's norm range): the pairs with `N gamma <= y^(1/2)`, plus those with
`N delta <= y^(1/2)`, minus those with both. Each piece is a sum over `N gamma <= y^(1/2)` of a polar-box count of `delta`
(coprime to `S`, rotated by `vartheta(gamma)`), whose main term is rotation-invariant and whose error is
`<= 8.3 (1 + tau'_S)(y/N gamma)^(1/2)` (inclusion-exclusion as in P3). With
`sum_{N gamma <= Y} N gamma^(-1/2) <= 8.4 Y^(1/2)` (partial summation from `#I(t) <= 4.2 t`), the rotation-dependent part of
the fibre sum is `<= 3 * 8.3 * 8.4 (1 + tau'_S) y^(3/4) <= 210 (1 + tau'_S) y^(3/4)`. Two rotations differ by at most twice
that; weight the fibre `d` by `d(d)`, put `y = x/N d`, and sum over `d != 1`.
(P5) As in P3: the grafted minus genuine weighted count is
`sum_{d != 1} w(d) [F_d(I - vartheta^S(d)) - F_d(I - vartheta(d))]`, with
`F_d(I') = sum_{beta coprime to S, N beta in J/N d, vartheta(beta) in I'} w(beta)`. Since `w 1_{coprime} = w * h`,
`F_d(I') = sum_e h(e) sum_{gamma in B(J/(N d N e), I' - vartheta(e))} w(gamma)`. By the uniform hypothesis, each inner sum is
`M_w(J/(N d N e), |I'|)` up to `C_E (x/(N d N e))^gamma`, and both rotations have the same main term (the inner sums
vanish when `x < N d N e`). So the difference is `<= 2 C_E x^gamma sum_{d != 1} |w(d)| N d^(-gamma) sum_e |h(e)| N e^(-gamma)`.
For twisted sums the grafted sum is `sum_d w(d) e^{4ik vartheta^S(d)} G_d` with
`G_d = sum_e h(e) lambda_k(e) sum_gamma w(gamma) f(N d N e N gamma) lambda_k(gamma)`, each inner sum a member of the family
at height `X_1/(N d N e)`. (If a family has a main term `Main(f, k)` for `k != 0`, closure under rotations of
`vartheta` forces `|G| <= E`, because `|e^{4ikt} G - Main| <= E` for every `t`.) The Euler-product bounds are Lemma 1.2(1). QED

### 2.4 The cap

**Corollary 2.2 (composite-supported robust inputs are capped at the rate).** Fix `c, gamma, B, C_c` (and, for P5, `C_0`)
**before** the target constant. Let `A > max(2, 1/gamma)`, `A >= 4`. For every Paley order `M = q + 1 >= M_4(A, c, gamma, B, C_c)`
the angle model `F_M` of `round8/lw.md` Thm 6.1 with `P = M^A` (as in `round10/analytic.md` Thm B: `M` distinct points on
an arc of length `<= C_0 R^(1/2)`, every monomial gap of `round4/outside.md` Lemma 0.1, Haar-generic angles,
`W <= (2A + o(1)) M log M`, `s = 2q < 2M` primes, all `>= M^A`) has:
* a minimal graft satisfying every `(c, gamma)`-robust composite input (C1)-(C3) up to `1 + eta_M`;
* a multiplicative graft satisfying every input of the kinds P1-P5 (P3 needs `gamma >= 1/2`, P4 `gamma >= 3/4`, P5 a uniform
  family at exponent `gamma`) up to `1 + eta_M`,

with `eta_M = O(M^(1 - A/2) + M^(1 - A gamma) 2^B) -> 0`. The same holds for the P3-char fakes of
`round5/round6_sharp/sharp.md` Thm 5.1 (`s = b(M-1)`, primes `>= 3M^A (log M)^2`, `W <= (bA + o(1)) M log M`), so residue data
through characters at moduli `<= 3M^A` may be added. Hence no argument that uses the exponent profile, the monomial gaps
(or P3-char residue data), Haar-generic properties, and composite-supported analytic inputs that are robust with fixed
`(c, gamma)` (used either without cross multiplicativity, or with it in the forms P1-P5) can prove `M = o(log R/loglog R)`.

*Proof.* `P_0 >= M^A`, `s < 2M` (resp. `b(M-1)`), so `kappa, kappa_B = (1 + o(1))/A < min(1/2, gamma)`, `t = O(M^(1-A/2))`,
`eta = O(2^B M^(1 - A gamma))`, `tau_sigma = O(M^(1 - A sigma))` (Lemma 1.2). Apply Theorems M and P. The fake's local
properties are unchanged. QED

*Remarks.*
1. **Same constants as Theorem B.** The composite inputs cost nothing beyond the prime-supported ones: `A > max(2, 1/gamma)`.
   This settles `round10/referee_analytic.md` Slip 3: the "dual-large-sieve argument" the referee asked for is P1, and the
   Lemma-3.2 obstruction (`K X s/P_0`) disappears because each S-fibre is a **rotation** of a genuine sub-ensemble, which
   the large sieve cannot detect (only the fibre masses `tau_(1/2)` enter).
2. **The round-9 R1 caveat stands, in sharper form.** The cap bars arguments with `gamma` fixed in advance. Unlike the
   prime-supported case, composite inputs with **arbitrarily small** `gamma` do exist and are true: smooth all-ideal Hecke
   sums are `O(X^(-B))` (Poisson summation). They are the subject of §3.1, where they turn out to collapse the fake to the
   actual cluster.
3. **Prop. C already covers composite counting.** `round10/analytic.md` Prop. C is stated for arbitrary finite candidate
   sets `Pi_j`, so absolute-value Fourier counting over composite candidates is rotation-invariant and never certifies
   emptiness, with no change.
4. **What P5 needs.** Published theorems are uniform in the position of the sector and in the height, which is the P5
   hypothesis. Theorem P does not need uniformity in a coprimality modulus (which would cost `d(prod S)^O(1) = 2^O(M)`): the
   coprimality to `S` is removed by inclusion-exclusion at Euler-product cost `1 + O(s P_0^(-gamma))`.
5. **Twisted sums at large frequency.** For a twisted sum at frequency `k` with `|k| <= X^theta`, a genuine uniform family
   at height `Y < k^2` holds only with `C_E >= Y^(1-gamma)`. So P5 at frequencies up to `X^theta` caps only for
   `gamma > 2 theta + 1/A` unless the input's own constant already absorbs this (Theorem M has no such restriction).

## 3. Non-robust composite inputs: what fakes violate, and why it does not help

### 3.1 Smooth all-ideal Hecke sums: true with error `X^(-B)`, violated by grafts, and they collapse the fake

For `k != 0` and `f in C_c^infinity((0, infinity))`, the genuine sum `F_k(X) = sum_{alpha in I} f(N alpha/X) lambda_k(alpha)`
equals `(2 pi i)^(-1) int_(2) f~(s) X^s L(s, lambda_k) ds` (`f~` the Mellin transform). `L(s, lambda_k)` is entire of
polynomial growth in vertical strips and `f~` decays rapidly, so shifting to `Re s = -B` gives `F_k(X) = O_(f,k,B)(X^(-B))`
for every `B`. This is a **true composite-supported input with error below every power of `X`**; it is not robust, and
§2 does not apply to it.

**Proposition 3.1 (entire-L collapse).** Let `vartheta'` be a completely additive angle function on `I` that agrees with
`vartheta` on every prime ideal not above `S`, is conjugation-equivariant, and has `vartheta'(pi_p) = phi_p` (for example the
multiplicative graft). Fix `k != 0`. Suppose that for every `f in C_c^infinity((0, infinity), C)` there is `delta_f > 0`
with `sum_alpha f(N alpha/X) e^{4ik vartheta'(alpha)} = O_f(X^(-delta_f))` as `X -> infinity`. Then for every `p in S`,
`4k phi_p = +- 4k theta_p (mod 2 pi)`. For `k = +-1` this says `phi_p = +- theta_p (mod pi/2)` for all `p`: the fake's
"primes" are the genuine `pi_p` or `conj(pi_p)` up to units, so the fake **is an actual cluster** (on the same circle).

*Proof.* Let `a_p = lambda_k(pi_p) = e^{4ik theta_p}` and `b_p = e^{4ik phi_p}`. For `Re s > 1`, by complete multiplicativity
and agreement off `S`,
`L'(s) := sum_alpha e^{4ik vartheta'(alpha)} N alpha^(-s) = L(s, lambda_k) prod_{p in S} E_p(s)`,
`E_p(s) = (1 - a_p p^(-s))(1 - conj(a_p) p^(-s)) / ((1 - b_p p^(-s))(1 - conj(b_p) p^(-s)))`.
This continues `L'` meromorphically to `C`.
*Mellin step.* For `f` as in the hypothesis, `G(X) = sum_alpha f(N alpha/X) e^{4ik vartheta'(alpha)}` vanishes for `X` below
`1/sup supp f`, and `int_0^infinity G(X) X^(-s-1) dX = f~(s) L'(s)` for `Re s > 1`. The hypothesis makes the left side
absolutely convergent and analytic on `Re s > -delta_f`. So `f~ L'` has no pole on `Re s = 0`.
*Poles.* Suppose `b_p notin {a_p, conj(a_p)}` for some `p`. Let `s_m = i(arg b_p + 2 pi m)/log p`, `m in Z`. The factor
`1 - b_p p^(-s)` has a simple zero at each `s_m`, and the numerator of `E_p` is nonzero there (`1 - a_p/b_p != 0`,
`1 - conj(a_p)/b_p != 0`). A numerator factor of another `E_(p')` vanishes at `s_m` only if
`(arg b_p + 2 pi m)/log p = (+-arg a_(p') + 2 pi m')/log p'`; two such coincidences for the same `p'` and sign would give
`(m_1 - m_2) log p' = (m'_1 - m'_2) log p`, hence `m_1 = m_2` by unique factorisation. So at most `2(s-1)` of the `s_m`
are cancelled by other numerators. Other denominators only add poles. `L(s, lambda_k)` is entire and has no zero on
`Re s = 0`: by the functional equation `Lambda(s) = W Lambda(1 - s, conj lambda_k)` with
`Lambda(s) = (2 pi)^(-s) Gamma(s + 2|k|) L(s, lambda_k)` (conductor 1, `|W| = 1`), `L(it)` is a nonzero multiple of
`L(1 - it, conj lambda_k) != 0` (nonvanishing on `Re s = 1`), and `Gamma(it + 2|k|)` is finite and nonzero. Hence `L'` has a
pole at all but finitely many `s_m`. Pick one such `s_m = i t_m` and `f(u) = f_0(u) u^(-i t_m)` with `f_0 >= 0`,
`f_0 != 0` in `C_c^infinity`: `f~(s_m) = int f_0(u) u^(-1) du > 0`, so `f~ L'` has a pole at `s_m`. Contradiction.
So `{b_p, conj b_p} = {a_p, conj a_p}` for all `p`. QED

*Meaning.* This is the composite analogue of `round10/analytic.md` Thm D: a composite input precise enough to see the
graft pins the fake's prime angles to genuine Gaussian-prime angles, i.e. it **is** the original problem. It supplies
no mechanism: after it, one must still bound actual clusters.

*Quantitative form (evidence, §4 check 2; heuristic statement).* With `phi_p = theta_p + epsilon u_p`, `|u_p| <= 1`, the
grafted smooth sum at `X = 10^6` is proportional to `epsilon` (`9.3e-2, 8.8e-3, 8.5e-4, 8.5e-5` for
`epsilon = 10^-1, ..., 10^-4`). The residues at the poles `s_m` are proportional to `b_p - a_p`, so precision `E` at height
`X` constrains roughly `sum_p |b_p - a_p|^2 <~ E^2` (averaging the phases `X^(i t_m)` over `log X`). Since
`|b_p - a_p| <= 2`, this is vacuous unless `E^2 < 4s`, i.e. unless the error is `O(M^(1/2))`, far below `X^gamma` for any
fixed `gamma > 0` at heights `X >= M^A`. And a precision `E = O(1)` pins the angles only to `|b_p - a_p| = O(E)`.
Whether gap-admissible fakes exist with all `|phi_p - theta_p| <= epsilon` (they would satisfy every such input at
precision `>~ epsilon s^(1/2)`) is **not decided** here. A dimension count (the realising set is a thickened subtorus of
dimension `~ r - M + 1`, and its image under the pair map covers the torus once `epsilon r^(1/2) >~ 1`) suggests
`epsilon ~ M^(-1/2)` is attainable for the angle constraints alone, with gap admissibility unchecked (**heuristic**).

### 3.2 Ring structure

If an argument uses both multiplicativity and the archimedean additive structure of the ensemble on all composites
(`arg(alpha + beta)` determined by `alpha, beta` as complex numbers), the ensemble is the embedding of the ring `Z[i]` into
`C`, so the angle function is genuine (or conjugate) everywhere, and the fake collapses to an actual cluster. The
minimal graft keeps the additive structure exactly among ideals that are not S-smooth; relations that involve the fake's
own points (chords `z_x - z_y`, sums `z_x + z_y`) are the **additive data** of `overview.md` §5.5, which are not modelled by
angle models and which, used in full, are the problem itself (`round10/additive.md`). Composite-supported analytic inputs
add nothing to that item.

### 3.3 The sharp large sieve on primitive support

By Prop. 1.5 the only composite large sieve that can have the sharp form `Lambda = O(X + K)` is one whose support has at
most one ideal per angle, for example **primitive** ideals (no rational integer factor `> 1`). On that support the genuine
angles of `alpha != beta` with `N alpha, N beta <= X` satisfy `|sin(vartheta(alpha) - vartheta(beta))| = |Im(alpha conj beta)|/(|alpha||beta|) >= 1/X`,
so the sharp form holds (Gallagher with `delta = 1/X`). By the round-10 referee's argument for primes (`referee_analytic.md`,
sharp large sieve; transferred, not re-derived here), up to constants it is equivalent to a **composite max-scale
cross-gap**: for every fake S-smooth primitive `d` (angle `<v, phi>` with `e^(w(v)) = N d`) and every primitive genuine `alpha`,
`|sin(<v, phi> - arg alpha)| >= c/max(N d, N alpha)`. Actual clusters satisfy it with `c = 1`, since
`|Im(A_v conj alpha)| >= 1` gives `|sin| >= e^(-w(v)/2) |alpha|^(-1) >= 1/max(e^(w(v)), N alpha)`.
* It contains round 10's prime version (`v = e_j`), which is unclassified (`round10/README.md` §3 item 6, §7.2).
* For a single angle `psi` it is a **badly-approximable** condition on `tan psi` (writing `alpha = x + iy`,
  `|sin(psi - arg alpha)|` is comparable to `|tan psi - y/x| x^2/N alpha` near `psi`; by Hurwitz every irrational slope violates
  it for `c` above an absolute constant, and by Khintchine almost every slope violates it for every `c > 0`), not a collapse: unlike LCG (Thm D), it does not
  force `<v, phi> = arg beta`. Whether gap-admissible fakes at barrier scale can satisfy it for all `v` simultaneously
  is a Schmidt-game question that this round does **not** decide.
* No growth mechanism from it is known. The record's assessment of the prime version is that at best it
  reconstitutes the actual problem up to `O(M)` in `W` (`round10/README.md` §3 item 6).
So composite support adds one family (monomial cross-gaps) to an item that was already unclassified; it creates no new
open class.

### 3.4 The cluster's own polar box

A polar box that contains an actual cluster at the critical scale has angular width `Delta = C N^(-1/4)` and outer arc
length `Delta N^(1/2) = C N^(1/4)`. Sharp lattice counts there have intrinsic error of the order of the boundary
(Lemma 1.3: up to `2.53 N^(1/2)` for the radial edges, and the arc itself carries the representation numbers of the
integers in the norm range), so no robust count in that box can see `M = O(log N)` points; a count with error `o(M)` in
that box is a statement about `r_2(n)` on the arc, i.e. the problem. Smoothing the box radially to width `>= N^(1/2)` in
norm makes Poisson summation exact up to `O(N^(-B))`, which is the §3.1 situation (it sees the graft, and collapses).

## 4. Computations (scripts in `composite_checks/`; outputs `*_out.txt`)

Labels: **exact** for integer counts, **evidence** for floating-point linear algebra or sums. No proof uses them.

1. **`ls_check.py`** (evidence; double-precision eigenvalues). Optimal all-ideal large-sieve constants
   `Lambda(X, K) = lambda_max(T(k - l))`, `T(m) = sum_{N alpha <= X} e^{4im vartheta(alpha)}`, for `X = 2*10^4, 10^5, 3*10^5`,
   `K = 5, 20, 60, 150`, `S` = 4, 6, 8 primes `>= 101, 211, 401`, random `phi`.
   * Prop. 1.5 lower bound `max(#I(X), (2K+1) floor(X^(1/2)))` and Lemma 1.4 upper bound hold in all 12 cases.
   * Multiplicative graft: relative change `+0.05%` to `+3.2%`; minimal graft: `<= 0.09%`. (Here `tau_(1/2) ~ 1.2`, far
     outside the asymptotic regime `s P_0^(-1/2) -> 0`; the changes are small anyway.) This extends the round-10 referee's
     `composite_graft.py` (`<= 7.5%`, two primes).
2. **`smooth_check.py`** (evidence). Smooth sums `F_k(X)` with a `C^infinity` bump on `(1/2, 2)`, `X <= 2*10^6`, `k = 1, 2, 5`,
   `S = {101, 109, 113, 137, 149, 157}`.
   * Genuine: `|F_k| <= 2.2e-8` (at `X = 2*10^4`) and `<= 4.1e-12` for `X >= 10^5` (rounding level).
   * Multiplicative graft: `|F_k|` between `0.03` and `2.1`. Minimal graft: `0.03` to `3.2` when S-smooth norms lie in
     `(X/2, 2X)`, and equal to genuine when none does (`X = 10^5, 5*10^5`; S-smooth norms here are `~10^4` and `~10^6`).
   * Pinning (Prop. 3.1, quantitative): `phi = theta + epsilon u` gives `|F_1(10^6)| = 9.3e-2, 8.8e-3, 8.5e-4, 8.5e-5` at
     `epsilon = 10^-1..10^-4`, linear in `epsilon`.
3. **`box_check.py`**.
   * (a) exact counts: over 400 random sectors, `max |count - area|/(2.53 y^(1/2) + 1.58) = 0.058 <= 1` (Lemma 1.3).
   * (b) exact counts: 400 random polar boxes, `x <= 2*10^6`: `max |graft - genuine| = 107`, at most `0.0008` of the
     proved bound of P3.
   * (c) exact: `Psi_S(X) = binom(2s + L, L)` exactly in the three cases (`455, 1771, 861`), consistent with Lemma 1.2(2).

## 5. The positive direction: what a separating composite statement would need

**The exact strength requirement.** Take a fake at `W = c_0 M log M`; Corollary 2.2 needs `P_0 = M^A` with
`A ~ c_0/2`. A composite input sees the minimal graft only if its error at height `X >= P_0` is below the graft's
footprint, the total weight of the S-smooth ideals, `Psi_S(X) max |w| <= X^(kappa_B) = X^((2 + o(1))/c_0)`; it sees the
multiplicative graft (in the forms P1-P5) only if its error is below `O(s P_0^(-gamma)) X^gamma`. So a composite input
that defeats **every** fake family of the record, for every `c_0`, must have, for each `c_0`,

> error `<= X^(gamma(c_0))` with `gamma(c_0) < 2/c_0`, uniformly over the boxes or sums it uses, at heights `X >= M^(c_0/2)`,

so `gamma -> 0`: precision `X^(o(1))`. (This is the composite form of the round-9 R1 caveat: `gamma` must be chosen after
`c_0`.) What exists at that precision:
* **Sharp** counts (lattice points in sectors and annuli, with any of the supports of §2.1): precision `X^gamma`,
  `gamma < 1/2`, is **false** uniformly in position. Sectors with rational edges of small height carry `>> X^(1/2)` excess
  points on the edge, and the arc of the cluster's own box carries `r_2(n)` (§3.4). GRH-quality inputs give
  `gamma >= 1/2` and are capped (P3, P5).
* **Smooth** all-ideal sums: precision `X^(-B)` is **true** (Poisson). They see the graft, but in exactly one way: they pin
  the fake's prime angles to genuine Gaussian-prime angles (Prop. 3.1), which is the problem itself. The quantitative
  version pins `|b_p - a_p|` only to the precision of the input, and is vacuous above precision `O(M^(1/2))` (§3.1,
  evidence plus heuristic).
* **Sharp large sieve on primitive support**: true, non-robust, equivalent up to constants (transferred, not
  re-derived) to the composite max-scale cross-gap of §3.3, whose status is the round-10 unclassified status. No growth
  mechanism from it is known.

**Verdict on the positive direction.** No composite-supported statement was found that (i) actual clusters satisfy,
(ii) fakes violate, (iii) is not averaged away (Prop. C, Theorems M and P), and (iv) is weaker than the problem. Every
fake-violated composite statement identified is either a lattice-structure statement that collapses the fake to the
actual cluster (Prop. 3.1, Thm D, §3.2), or the sharp primitive large sieve, which only joins the existing unclassified
round-10 item. **No known theorem reaches a non-collapsing fake-violated composite statement**, and none of them gives
growth.

## 6. Barrier checklist, item by item

The cap (Corollary 2.2) is itself a barrier result, so the checklist is applied to the **positive** candidates of §§3, 5.

1. *Inputs that fake circles satisfy* (pairs, gcds, cut, residues, characters, GAP-Fourier, energy, `L^p`). The
   candidates of §3 are not of these kinds: they are ensemble statements. Robust composite ensemble statements are now
   added to this list by Corollary 2.2 (grafted fakes satisfy them).
2. *Robust or averaged analytic statements* (Thm B, Prop C, Thm D). Theorems M and P extend Thm B to composite support at
   the same constants; Prop. C already covers composite candidate sets; Prop. 3.1 is the composite Thm D. The candidates
   that survive item 2 (smooth sums at precision `X^(-B)`, sharp primitive large sieve) are pointwise, hence fall under
   the Thm D dichotomy: collapse (§3.1) or unclassified cross-gap (§3.3).
3. *Auxiliary polynomials* on `X_M` (`tau* < 2`) or on `Mbar_(0,k)`. Not used.
4. *Single linear forms in logarithms; SLW*. Not used. (Prop. 3.1 uses only the Euler product, the functional equation
   and nonvanishing of `L(s, lambda_k)` on `Re s = 1`.)
5. *S-integer arguments with thresholds depending on the primes only through `|S|`*. Not used. All bounds of §2 depend on
   `S` only through `s` and `P_0`, but they are **upper bounds on graft footprints** (used to build the cap), not thresholds of
   a growth argument; the cap holds for every fake of the record.
6. *Route 5's box statement and CIRCLE(q)*. Not used.
7. *`(L_5)`/`(L_6)`-type family mechanisms at `k = 8`; computational tests of `(L_8)`*. Not used.
8. *The abc unit-pencil lane*. Not used.
9. *The exact-Dirichlet-exponent fact* (`round4/hardness.md`). Consistent: robust inputs cannot reach it (they have a fixed
   power saving, i.e. they act at heights where the fake's footprint is negligible), and the inputs that do reach it
   (Prop. 3.1, §3.3) are pointwise lattice statements.

**Project constraints.** No uniform exceptional set from a fixed-`S` Subspace Theorem is assumed (none is used). No
Gaussian gcd or denominator cost is discarded: the grafts change only angles; profile, gcds and residues of the fake are
those of the refereed fakes, and inclusion-exclusion and Dirichlet inverses over `S` are carried with their full
Euler-product costs (`tau`, `tau'`, `Pi_q`). No drifting common phase is replaced by zero: rotations `rho(d)` are kept
exactly, and the cap uses only that inputs are blind to them. No template-dependent threshold is treated as uniform in
the radius: `M_4(A, c, gamma, B, C_c)` depends only on the input's constants, which are fixed before the target constant,
exactly as in Thm B; the caveat that arguments choosing `gamma` after `c_0` are not barred is stated (Cor. 2.2 Remark 2,
§5).

## 7. Status labels and open points

| item | status |
| --- | --- |
| Lemmas 1.2-1.4, Prop. 1.5 (all-ideal large sieve is automatically robust) | proved here, unrefereed |
| Theorem M (minimal graft: every robust composite input, without cross multiplicativity) | proved here, unrefereed |
| Theorem P (multiplicative graft: P1 large sieve any support; P2 Type II / multilinear; P2' Halasz; P3 sharp polar boxes `gamma >= 1/2`; P4 divisor-weighted `gamma >= 3/4`; P5 uniform families) | proved here, unrefereed |
| Corollary 2.2 (composite robust inputs capped at `W <= (2A + o(1)) M log M`, `A > max(2, 1/gamma)`, `A >= 4`) | proved here, unrefereed; R1 caveat applies |
| Prop. 3.1 (entire-L collapse) | proved here, unrefereed |
| quantitative pinning; near-genuine fakes | evidence + heuristic |
| sharp primitive large sieve / composite max-scale cross-gap | open (inherits round-10 unclassified status) |
| computations §4 | exact (counts) / evidence (floating point) |

**Proposed change to the barrier map** (`overview.md` §5.4-5.5), subject to refereeing: move "composite-supported analytic
inputs" from §5.5 (unclassified) to §5.4 as a capped row (kind F, witness: minimal and multiplicative grafts of the
`lw.md` 6.1 and P3-char fakes, `W <= (2A + o(1)) M log M` resp. `(bA + o(1)) M log M`, same R1 caveat as row 7), and list in
§5.5 only the sharp large sieve on primitive support as part of the existing sharp-large-sieve item.

**Open points.**
1. P5 at small `gamma` with frequencies up to `X^theta` needs `gamma > 2 theta + 1/A` (Remark 5 of §2.4); Theorem M has no such
   restriction but does not allow cross multiplicativity. An input that is non-uniform in position, has `gamma < 1/2`, and
   is used with cross multiplicativity is not covered by Theorem P. No published input of that kind is known to me.
2. Whether gap-admissible fakes can have all `|phi_p - theta_p| <= epsilon` with `epsilon -> 0` (near-genuine fakes). They
   would satisfy the quantitative smooth-sum inputs at every precision `>~ epsilon s^(1/2)`.
3. Whether a non-multiplicative ensemble consistent with a fake can satisfy the smooth-sum decay of §3.1 (a compensation
   question). Irrelevant for growth: with multiplicativity the decay collapses the fake.
4. The composite max-scale cross-gap (§3.3): satisfiable by barrier-scale fakes or not.
