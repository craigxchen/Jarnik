# The fully flipped Hadamard phase lemma: two-logarithm reduction, a conditional proof, and the exact transcendence gap

Round 4, target: the "fully flipped Hadamard phase lemma" of Section 10 of
`docs/claude_review/growth/walsh-extraction.md` (and its referee report). Checks are in
`round4/flipped_checks/` (paths below are relative to that directory; every script is run from it
with `python3 <script>`; outputs are stored next to the scripts as `*_out.txt`).

## Status (honest)

**The lemma is not proved unconditionally here, for Paley cores or for any other non-Walsh type.**
For Sylvester (Walsh) cores it is already known unconditionally: it is the case `F = all rows` of
Theorem C of walsh-extraction.md (item 333 of the notes), which gives `W/(M log M) -> infinity`.
No growth saving is obtained for Paley or general cores without a hypothesis.

What is new (all statements proved below with full quantifiers, unless labelled heuristic/data):

| # | result | type |
|---|---|---|
| Prop 2.1-2.2 | Every label pinning (and every label-pair pinning) of a Hadamard-core profile with arbitrary flips is a **nonzero linear form in two logarithms of two varying numbers of Q(i)** (a block and a Hadamard product of flipped primes) plus `log i`, with error `<= M Delta`, uniformly in `M` | exact reduction, new |
| Thm L | **Conditional on a two-logarithm lower bound of Lang–Waldschmidt type over Q(i)**, every Hadamard-core profile (any normalised `H`, any flips, any capacity, any `W`) whose flip mass is `F <= (1/(2 kappa) - delta) W` has `M < M_0(kappa, delta, C)`. Under the Lang–Waldschmidt conjecture (`kappa = 1+eps`) this covers `F < (1/2 - delta) W`, hence the whole Prop-P barrier class (`F/W ~ 1/b <= 1/5`) and every fully flipped profile whose flipped primes are not heavier on average than the others | conditional theorem, new |
| Prop 3.1 / Rem 3.2 / Cor 3.3 | Liouville's inequality on these forms is exactly the column character; it gives only `D <= 2F + O(W/M)` and otherwise misses a contradiction by exactly `F/2 - D/4 + O(W/M)` | unconditional, quantified |
| Prop B | Every bound of Baker type (`log 1/|Lambda| <= c h_1 h_2 L`) is useless for the lemma: in any configuration `h_1 h_2 >= (W+D) min((1/8 - o(1)) log M, 0.87)` for single-label forms (`1/16` for pairs), i.e. known bounds miss by a factor `~ c log M` | unconditional, the sharp obstruction |
| Prop U, Lemma 6.2 | **Flipped primes are not pinned**: no rational consequence of the phase system isolates flipped angles; isolating `t` of them drags in `>= s(t)` unflipped blocks, `s(t) >= (M/2 - t)/(B+1)` for prime-order Paley (integer prime-field uncertainty lemma) | unconditional, new |
| Lemma 6.3, Prop 6.4 | Sector/grid counting: if the flipped angles were pinned to `<= M-5` directions mod `pi/2` at precision `O(M Delta)`, the lemma would follow (via `prod q_x >= (1/(2 eps))^(M-G)`); by Prop U the system supplies no such pinning | unconditional |
| Prop 7.1 + data | For Paley cores in the **heavy-flip regime `F >= (1/2 + delta) W`**, every bounded-support form has support weight `> W/4`, so even the Lang–Waldschmidt route needs a different input there; this regime remains open even conditionally | partial obstruction |

**Sharpest obstruction (summary).** In the regime `W_F < W/2` (which contains the whole Prop-P
barrier class) the lemma is exactly one transcendence statement away: each label condition is a
two-logarithm form `M log u_a - 2 log v_a` with `h(u_a) ~ W/(2M)` and
`h(v_a) = W_F/2 >= (M/2) log(4M/e)`. A contradiction needs a lower bound `|Lambda| >= e^(-(W+D)/4 + O(log MC))`.
The naive heights `log H(u_a) + log H(v_a) ~ W/M + W_F` sum to less than `W/4` when `W_F < W/4`
(label pairs: `W_F < W/2`), which is enough under Lang–Waldschmidt, while every known unconditional
bound depends on the *product* `h(u) h(v) >= (1/8 - o(1)) (W+D) log M`, which exceeds `(W+D)/4` by a
factor of order `log M`. Liouville (one-number arguments, i.e. all integral/rational/fractional
characters) misses by exactly `F/2 - D/4 + O(W/M)`, and forces `D <= 2F + O(W/M)`.

---------------------------------------------------------------------------------------------

## 0. Setting and conventions

`H` is a Hadamard matrix of order `M >= 8` (so `4 | M`) with **first column all ones** (column `0`,
the constant label); `A = {1, ..., M-1}` are the nonconstant labels; `H^T H = M I`, and
`sum_x H(x,a) = 0` for `a in A`. Rows may not be renormalised (a row sign change is not a symmetry
of the problem); columns may (it reorients a block).

**Profiles with modifications.** Columns `j in J` are conjugate-primitive nonunit Gaussian integers
`kappa_j` (single split primes in the lemma; products with pairwise disjoint rational prime supports
are allowed throughout), `w_j = log Norm(kappa_j)`, with label `a(j) in {0} u A`, orientation
`sigma_j = +-1` and a *modification set* `F_j` of rows. With `kappa^(1) = kappa`,
`kappa^(-1) = conj(kappa)` (and `kappa^(e)` for `e in Z` accordingly, always a Gaussian integer),

```text
z_x = g eps_x prod_j kappa_j^(s_xj),     s_xj = sigma_j H(x,a(j)) (-1)^[x in F_j],     eps_x in {1,i,-1,-i}.
```

Unmodified label-0 columns are absorbed into `g`; `W = sum w_j` over the remaining columns,
`D = log Norm g`, `R^2 = e^(W+D)`. The **flip mass** is `F = sum_j |F_j| w_j`.
The lemma's class ("fully flipped, one flip per row at distinct primes, capacity `B`") is:
`F_j = {x}` for `j = j(x)`, `x -> j(x)` injective, `a(j(x)) != 0`, each label carrying `b >= 5`
columns of which at most `B` are flipped; then `F = W_F := sum_x w_(j(x))`, and we write
`q_x = Norm kappa_(j(x))`, `sigma_x = H(x, a(x))`.

Blocks and flip products (all Gaussian integers, conjugate-primitive by disjointness):

```text
G_a = prod_{a(j)=a} kappa_j^(sigma_j),        W_a = log Norm G_a,          sum_{a in A} W_a <= W,
c_j(a) = sum_{x in F_j} H(x,a) H(x,a(j)),     P_a = prod_j kappa_j^(sigma_j c_j(a)),     F_a = log Norm P_a <= F,
d_j(a,b) = (c_j(a) - c_j(b))/2 in Z,          Q_ab = prod_j kappa_j^(sigma_j d_j(a,b)),  F_ab = log Norm Q_ab.
```

(`d_j` is an integer because `H(x,a) - H(x,b) in {0, +-2}`.) With
`D_ab = {x : H(x,a) != H(x,b)}` (exactly `M/2` rows), `F_ab <= sum_j w_j |F_j cap D_ab|`, with
equality in the one-flip-per-column case. In the fully flipped class `c_(j(x))(a) = H(x,a) sigma_x`
for every `a`, so `F_a = W_F` for every label and `F_ab = sum_{x in D_ab} log q_x`.

The points lie on an arc of length `<= C sqrt R`, i.e. of angular width `Delta <= C R^(-1/2) = C e^(-(W+D)/4)`.
As in walsh-extraction.md (1.1), with real lifts `phi_j` of `arg kappa_j`:
`(S phi)_x = theta' + delta_x + (pi/2) k_x`, `|delta_x| <= Delta/2`, `k_x in Z`.

**Heights.** For `alpha in Q(i)`, `H(alpha)` is the naive height (largest absolute coefficient of the
primitive minimal polynomial over `Z`) and `h(alpha)` the absolute logarithmic height. If `gamma` is a
conjugate-primitive nonunit with odd norm `N`, then `u = gamma/conj(gamma)` has primitive minimal
polynomial `N X^2 - 2 Re(gamma^2) X + N`, so `H(u) <= 2N` and `h(u) = (1/2) log N`
(primitivity: a prime dividing `N` and `Re gamma^2 = x^2 - y^2` divides `x` and `y`; checked in
`check_two_log_identity.py`).

**Axis gap** (walsh-extraction Lemma 1.1, used repeatedly): a conjugate-primitive nonunit Gaussian
integer of odd norm is neither on an axis nor on a diagonal; if its argument is within `e` of
`(pi/4) Z` then `|gamma| sin e >= 1/sqrt 2`.

## 1. The modified pinning, exactly

**Lemma 1.1.** For every `a in A` and every pair `a != b` in `A`,

```text
M Phi_a - 2 sum_j c_j(a) sigma_j phi_j            = (pi/2) m_a        + sum_x H(x,a) delta_x,                 (1.2)
M (Phi_a - Phi_b) - 4 sum_j d_j(a,b) sigma_j phi_j = (pi/2)(m_a - m_b) + sum_x (H(x,a) - H(x,b)) delta_x,      (1.3)
```

with `Phi_a = sum_{a(j)=a} sigma_j phi_j`, `m_a = sum_x H(x,a) k_x in Z`, and both error terms of
absolute value `<= M Delta/2`. Equivalently, in `Q(i)` (exact identity, `zeta_x = z_x/conj(z_x)`):

```text
prod_x zeta_x^(H(x,a)) = +- (G_a/conj G_a)^M (conj P_a / P_a)^2.                                              (1.4)
```

*Proof.* `sum_x H(x,a) s_xj = sigma_j [M delta_(a,a(j)) - 2 c_j(a)]` because
`sum_x H(x,a) H(x,a(j)) = M delta_(a,a(j))` (also for `a(j) = 0`) and each modified row flips the
sign of one term. Multiply (1.1) by `H(x,a)` and sum: `theta'` cancels since `sum_x H(x,a) = 0`.
The error is `<= (Delta/2) sum_x |H(x,a)| = M Delta/2`; for (1.3), `H(x,a) - H(x,b)` is `+-2` on `M/2`
rows and `0` elsewhere, again `<= M Delta/2`. For (1.4): `zeta_x = eps_x^2 prod_j (kappa_j/conj kappa_j)^(s_xj)`,
`g/conj g` cancels, `prod_x eps_x^(2H(x,a)) = +-1`. QED

For the lemma's class this is Lemma 7.2 of walsh-extraction.md (with `psi_x = sigma_(j(x)) phi_(j(x))`).
As the referee noted, the pinnings are **necessary** conditions only (the row conditions are recovered
from them with error `(M-1) Delta/2`, and `m` must lie in `H^T Z^M`); everything below uses them only
as necessary conditions, so no equivalence is needed.

*Check.* `check_two_log_identity.py` verifies (1.4) and its pair version exactly (big Gaussian
integers, no floating point) on literal profiles: Walsh-8/16 and Paley-12/20, random capacity-two
flip assignments, random orientations and units, plus profiles with several flips per column and
flipped content columns (2,600 exact identities). It also checks that wrong exponents fail
(negative control) and the minimal-polynomial height facts.

## 2. Two-logarithm reduction

**Proposition 2.1 (single label).** Let `a in A`, `u_a = G_a/conj G_a`, `v_a = P_a/conj P_a`
(elements of `Q(i)` of modulus one), `Log` the principal logarithm, `Log i = i pi/2`. There is
`b_3 in Z`, `|b_3| <= 2M + 5`, such that

```text
Lambda_a := M Log u_a - 2 Log v_a + b_3 Log i     satisfies   |Lambda_a| <= M Delta,
H(u_a) <= 2 e^(W_a),     H(v_a) <= 2 e^(F_a)   (H(v_a) = 1 if P_a is a unit).
```

If some column `j` of label `a` has `|F_j| < M/2` (always true in the lemma's class), then `Lambda_a != 0`.

**Proposition 2.2 (label pair).** For `a != b` in `A`, with `U_ab = G_a conj(G_b)`,
`u_ab = U_ab/conj U_ab`, `w_ab = Q_ab/conj Q_ab`, there is `b_3 in Z`, `|b_3| <= 2M + 9`, with

```text
Lambda_ab := M Log u_ab - 4 Log w_ab + b_3 Log i,      |Lambda_ab| <= M Delta,
H(u_ab) <= 2 e^(W_a + W_b),      H(w_ab) <= 2 e^(F_ab),
```

and `Lambda_ab != 0` as soon as some column of label `a` or `b` has `|F_j| < M/4`.

*Proof.* `Log u_a = 2 i arg G_a + 2 pi i n_1` and `Log v_a = 2 i arg P_a + 2 pi i n_2` for integers
`n_1, n_2`, and `arg G_a = Phi_a`, `arg P_a = sum_j c_j(a) sigma_j phi_j` modulo `2 pi`. By (1.2),
`M Log u_a - 2 Log v_a = 2 i eta + i pi m_a + 2 pi i n` with `|eta| <= M Delta/2` and `n in Z`. Put
`b_3 = -2(m_a + 2n)`; then `Lambda_a = 2 i eta`. Since `|Log| <= pi`,
`|b_3| pi/2 <= M pi + 2 pi + M Delta`, giving the bound on `|b_3|` whenever `M Delta <= pi/2` (the only
case used below: if `M Delta > pi/2` then `(W+D)/4 < log(2MC/pi)` and all statements below are
trivial). The heights follow
from Section 0 (`G_a`, `P_a` are conjugate-primitive of odd norm).

Nonvanishing. If `Lambda_a = 0` then `M arg G_a - 2 arg P_a in (pi/4) Z`. This is the argument of the
Gaussian integer `Gamma_a = G_a^M conj(P_a)^2`. Combining `kappa^(s) kappa^(t) = Norm(kappa)^(min) kappa^(s+t)`
for opposite signs, `Gamma_a` is a positive rational integer times the conjugate-primitive
`Gamma'_a = prod_j kappa_j^(sigma_j e_j)`, where `e_j = M 1[a(j)=a] - 2 c_j(a)`. For the column `j` of
label `a` with `|F_j| < M/2`, `|c_j(a)| <= |F_j|` gives `e_j != 0`. So `Gamma'_a` is a nonunit of odd
norm with argument in `(pi/4) Z`, contradicting the axis gap.

The pair case is identical with (1.3), `Gamma_ab = U_ab^M conj(Q_ab)^4`, and `e_j = +-M - 4 d_j(a,b)`,
`|d_j| <= |F_j| < M/4`; the bound on `b_3` gains `4 pi` instead of `2 pi`. QED

*Remark 2.3 (what is new).* The notes use two-logarithm bounds (Matveev, items 276-277 and
`minimal_prime_rank_width_bound.md`) only for forms `N log alpha - m log(-1)`: **one** varying number
and a root of unity, i.e. a fixed target. In (2.1)-(2.2) **both** numbers vary with the configuration.
This is the precise meaning of the "moving targets" of Section 10: the target of `u_a` is the
Hadamard product `v_a` of all flipped primes. The reduction to two numbers, uniformly in `M`, is
special to cores with a square inverse of small denominator (`H^T/M`): the general `M`-row profile of
item 532 (all patterns on few rows) has no such inversion (Section 8).

## 3. One-number arguments: what Liouville gives

**Proposition 3.1 (column character, unconditional).** In every configuration as in Section 0 and
every `a in A` having a column with `|F_j| < M/2`,

```text
(M/2) W_a + F_a  >=  (W + D)/2 - 2 log( M C / (2 sqrt 2) ).                                                    (3.1)
```

*Proof.* Halve (1.2) (`M` is even): `(M/2) arg G_a - arg P_a` lies within `M Delta/4` of `(pi/4) Z`.
This is the argument of `gamma_a = G_a^(M/2) conj(P_a)`, whose conjugate-primitive part is a nonunit
of odd norm (exponent `M/2 - c_j(a) != 0` at the column of label `a`). The axis gap gives
`|gamma_a| sin(M Delta/4) >= 1/sqrt 2`, i.e. `log |gamma_a| >= log(2 sqrt 2/(M Delta))`.
Use `log|gamma_a| = ((M/2) W_a + F_a)/2` and `Delta <= C e^(-(W+D)/4)`. QED

(3.1) is the integral character `lambda = H(., a)/2` (half-integral entries, integral column image), so
it is one of the inequalities that Prop P of walsh-extraction.md already shows to be satisfied by the
barrier profiles. Lemma F.1 there shows that in the fully flipped class every rational character
(every Theorem-F consequence, every modulus `m`) is half-integral. Moreover, for any rational
`lambda` with fractional column image of denominator `e`, the only Gaussian integer the method sees is
`beta(lambda)^e = beta(e lambda)`, and `e lambda` is again half-integral. So **every one-number
(Liouville) argument extracted from the phase system is one of these integral character inequalities.**

**Remark 3.2 (exact size of the Liouville deficit).** In the language of Section 2, Liouville gives
`log(1/|Lambda_a|) <= (M/4) log H(u_a) + (1/2) log H(v_a) + O(1)`. A contradiction would need
`(M/4) W_a + F_a/2 < (W+D)/4 - O(log MC)` for some `a`. Since `min_a W_a <= W/(M-1)`, the left side
can be as small as `W/4 + F/2 + O(W/M)`: **Liouville misses by exactly `F/2 - D/4 + O(W/M)`.** For
primitive clusters (`D = 0`) this is the flipped mass `F/2`. For pure cores (`F = 0`) the deficit
vanishes, and Theorem B of walsh-extraction.md closes the gap with the residue pigeonhole. In the
fully flipped class, `F = W_F >= M log(4M/e)` (Section 5). The same accounting for any rational
`lambda` gives the deficits computed in items 391-393 (all positive for Prop-P profiles).

**Corollary 3.3 (content is controlled by the flips).** In the fully flipped class (where `F_a = W_F`
for all `a`), summing (3.1) over `a in A` and using `sum_a W_a <= W` gives
`D <= 2 W_F + W/(M-1) + 4 log(MC/(2 sqrt 2))`. So Liouville alone gives a contradiction exactly
when the common content exceeds twice the flipped mass. Otherwise it gives nothing.

## 4. A conditional proof

**Hypothesis `(H_kappa)` (two logarithms over Q(i)).** There are `kappa >= 1`, `kappa' >= 0` and
`c_* > 0` such that for all `alpha_1, alpha_2 in Q(i)` with `|alpha_1| = |alpha_2| = 1` and all
integers `b_1, b_2, b_3` with `Lambda = b_1 Log alpha_1 + b_2 Log alpha_2 + b_3 Log i != 0`:

```text
|Lambda| >= c_* B^(-kappa') ( H(alpha_1) H(alpha_2) )^(-kappa),          B = max(2, |b_1|, |b_2|, |b_3|).
```

*Relation to known conjectures.* This is the case "three logarithms, one of them `log i`, numbers in
`Q(i)`" of the **Lang–Waldschmidt conjecture** on linear forms in logarithms. That conjecture asserts
`|Lambda| >= C(eps)^m B / (|b_1|...|b_m| A_1...A_m)^(1+eps)`, with `A_j = max(H(alpha_j), e^(|log alpha_j|))`
(see M. Waldschmidt, *Diophantine Approximation on Linear Algebraic Groups*, Springer 2000, Chapter 1).
With `m <= 3`, `|Log| <= pi` and `H(i) = 1`, it gives `(H_kappa)` with `kappa = 1 + eps`,
`kappa' = 2 + 3 eps`, `c_* = C(eps)^3 e^(-3 pi (1+eps))`. I state `(H_kappa)` separately because the
**height normalisation matters**. With the absolute height `e^(h)` in place of the naive height, the
statement is false for these forms. `check_lw_scale.py` (`M = 12`, all conjugate-primitive
`G, P` with norms `<= X`) finds `min d sqrt(Norm G Norm P)` falling from `3.6*10^-2` (`X = 500`) to
`1.5*10^-5` (`X = 2*10^4`), where `d = dist(M arg G - 2 arg P, (pi/2)Z) ~ |Lambda|/2`. With naive
heights, `min d Norm G Norm P` stays between `0.3` and `0.7` up to `X = 8000`. It then drops to
`0.011` because of a Roth-type approximation of a fixed direction (`P = 6+i`, `G = 71+98i`). This is
compatible with exponent `1 + eps` and is evidence about the normalisation only, not about the truth
of the hypothesis. Exact `Lambda = 0` pairs (for example `P` an associate of `(2+i)^6` against
`G = 2+i`) are excluded, as the hypothesis requires. A refined form of the
Lang–Waldschmidt conjecture is known to imply a version of the abc conjecture (A. Baker,
*Logarithmic forms and the abc-conjecture*, 1998). It is open, and so is every `(H_kappa)`.

**Theorem L (conditional on `(H_kappa)`).** Consider a configuration as in Section 0 (any normalised
Hadamard matrix of order `M`, any orientations, units, content and modification sets) on an arc of
length `<= C sqrt R`, and assume every nonconstant label has a column with `|F_j| < M/4`.
Then, for every pair `a != b`,

```text
(W + D)/4  <=  kappa (W_a + W_b + F_ab) + 2 kappa log 2 + kappa' log(2M + 9) + log(M C / c_*),               (L.1)
```

and hence

```text
(W + D)/4  <=  2 kappa W/(M-1) + kappa F (M-1)/(2(M-2)) + E_M,      E_M := 2 kappa log 2 + kappa' log(2M+9) + log(MC/c_*).   (L.2)
```

Consequently, if `0 < delta < 1/(2 kappa)` and `F <= (1/(2 kappa) - delta) W`, then `M < M_0`, where
`M_0` is the least integer `>= 2 + (1 + 8 kappa)/(kappa delta)` such that
`(kappa delta/4)(M - 1) log 5 > E_M` holds for all `M >= M_0` (it exists because `E_M = O(log M)`).
No bound on `W` (no `K`), on the capacity, or on the prime sizes is needed.

*Proof.* We may assume `c_* <= 1` and `M Delta <= pi/2`; otherwise `(W+D)/4 < log(2MC/pi)` and (L.1)
is trivial. (L.1): by Prop 2.2, `0 < |Lambda_ab| <= M Delta <= M C e^(-(W+D)/4)`. Apply `(H_kappa)` with
`alpha_1 = u_ab`, `alpha_2 = w_ab`, `B <= 2M + 9`, `H(u_ab) <= 2 e^(W_a+W_b)`, `H(w_ab) <= 2 e^(F_ab)`,
and take logarithms.
(L.2): average (L.1) over the `binom(M-1, 2)` pairs. First, `sum_{a<b} (W_a + W_b) = (M-2) sum_a W_a <= (M-2) W`,
so the average of `W_a + W_b` is `<= 2W/(M-1)`. Second, `F_ab <= sum_j w_j |F_j cap D_ab|`, and a fixed row
`x` lies in `D_ab` for exactly `P_x N_x <= (M-1)^2/4` pairs, where `P_x, N_x` count the `+-1` entries
of row `x` on `A`. So the average of `F_ab` is `<= F (M-1)/(2(M-2))`. Some pair is at most the
average of `W_a + W_b + F_ab`.
Conclusion: with `F <= (1/(2 kappa) - delta) W`,
`kappa F (M-1)/(2(M-2)) <= (1/4 - kappa delta/2) W (1 + 1/(M-2))`. So (L.2) gives
`(kappa delta/2) W <= W (1/4 + 2 kappa)/(M-2) + E_M`. For `M - 2 >= (1 + 8 kappa)/(kappa delta)` the first
term is `<= (kappa delta/4) W`, leaving `(kappa delta/4) W <= E_M`. Every nonconstant label has a
nonunit column, so `W >= (M-1) log 5`, and this contradicts the choice of `M_0`. QED

**Corollaries.**

* **(L.3) The lemma under Lang–Waldschmidt.** Under `(H_(1+eps))`, for every `C` and `delta > 0`, fully
  flipped Hadamard-core clusters with `W_F <= (1/2 - delta) W` have bounded `M` (any `W`, any `K`).
  This is the fully flipped phase lemma of Section 10, in a stronger form, on the regime
  `W_F < W/2`.
* **(L.4) The Prop-P barrier class is excluded.** With `b` copies per label and nearby primes
  (`tau <= log p_j < tau + 1/r`), `W_F/W <= (tau + 1/r) M/(tau b (M-1)) <= 1/4` for `b >= 5`,
  `M >= 6`, `tau >= 1`. More generally, whenever the flipped primes are on average no heavier than
  the unflipped ones, `W_F/W <= M/(b(M-1)) < 1/2` for `b >= 3`. In both cases `(H_kappa)` with any
  `kappa < 2` (respectively `kappa < b/2` asymptotically for nearby primes) suffices. Exact
  bookkeeping on literal nearby-prime profiles (`check_margins.py`, Paley-12 to Paley-200 and
  Walsh-8 to Walsh-128, `b = 5`): the largest admissible exponent
  `kappa_2* = (W/4)/min_(a,b)(W_a + W_b + F_ab)` grows from `0.86` (M=12) through `1.19` (M=20) to
  `2.26` (M=200), approaching `b/2 = 2.5`. The single-label version
  `kappa_1* = (W/4)/min_a (W_a + W_F)` tends to `b/4 = 1.25`.
* **(L.5) The intermediate flip range** of the referee report (Paley with `M/(40 log M) < f < M` flipped
  rows) is also covered under `(H_kappa)` whenever the flipped mass is `< W/(2 kappa)`. The theorem is
  insensitive to how many rows are flipped.
* **(L.6) Variants.** The proof only needs a two-logarithm bound
  `log(1/|Lambda|) <= kappa_1 log H(alpha_1) + kappa_2 log H(alpha_2) + O(log B)` with
  `2 kappa_1/(M-1) + kappa_2 F/(2W) < 1/4 - delta'`. So the exponent on the **block side may grow
  like `eps M`**. Only the exponent on the flipped side matters. For example, a single-label bound
  improving Liouville's `(M/4, 1/2)` to `((1-eta) M/4, 1/2)` would exclude all profiles with
  `F < (eta/2 - o(1)) W`.

**Why the conditional statement is not "conjecture in, conjecture out".** `(H_kappa)` is a statement
about **two** numbers of `Q(i)`, uniform in everything. It does not imply the general
Cilleruelo–Granville bound by this route: for a general profile, every rational row combination has
support weight `>= ~W/2` (the pair conditions are critical; Section 8). It does settle the
structured class on which every rational-character, collision and capacity argument provably
fails (Prop P). This locates that barrier inside classical transcendence theory: it is a
"sum-of-heights versus product-of-heights" problem for two logarithms (Section 5). It is not a
problem of simultaneous approximation in growing dimension, contrary to the description in
Section 10 of walsh-extraction.md: one label pair suffices.

## 5. The gap to known transcendence bounds (the sharp obstruction)

All known unconditional lower bounds for linear forms in two (or a few) logarithms have the shape

```text
log(1/|Lambda|) <= c h(alpha_1) h(alpha_2) L,        L >= 1,                                                  (5.1)
```

up to lower-order terms. Here `h` is the absolute logarithmic height and `L` collects `log B`-type
factors. Examples are Baker–Wüstholz (1993), Matveev (2000; the notes use its degree-two
specialisation with `c <= 2^40`), Laurent–Mignotte–Nesterenko (1995) and Laurent (2008).
For our forms, `h(u_a) = W_a/2`, `h(v_a) = F_a/2`, `h(u_ab) = (W_a + W_b)/2`, `h(w_ab) = F_ab/2`.

**Lemma 5.1.** Any `n` distinct split primes have `sum log p >= sum_{k <= n} log(4k+1) > n log(4n/e)`.
(*Proof:* the `k`-th smallest prime `= 1 mod 4` is `>= 4k+1`; `n! >= (n/e)^n`.)

**Proposition B.** Let a fully flipped Hadamard-core configuration (one flip per row at distinct
split primes, every label with at least `b >= 5` columns at distinct split primes) lie on an arc of
length `<= C sqrt R`, and let `M >= 8`. Then for every label `a`,

```text
h(u_a) h(v_a) >= (W_F/2) max{ (W + D - 2 W_F)/(2M) - (2/M) log(MC/(2 sqrt 2)) ,  6.99 }.
```

Write `V = W + D`. In particular, as `M -> infinity` with `C` fixed, `h(u_a) h(v_a) >= V min((1/8 - o(1)) log M, 0.87)`.
For the pair forms, `h(u_ab) h(w_ab) >= V min((1/16 - o(1)) log M, 0.87)` for every pair. Hence a bound
of shape (5.1) yields a contradiction for a single label only if `c L <= (2 + o(1))/log M`, and for a
pair only if `c L <= (4 + o(1))/log M`, or else in either case only if `c L < 0.29`. All published
explicit two-logarithm constants have `c L >= 1` by a wide margin; the notes use Matveev's degree-two
specialisation with `c <= 2^40`. **So no known two-logarithm bound gives anything once `log M > 4`,
and the deficit grows like `c log M`.**

*Proof.* `F_a = W_F` (Section 0). By (3.1), `W_a >= (V - 2 W_F)/M - (4/M) log(MC/(2 sqrt 2))`. Also
`W_a >= log(5 * 13 * 17 * 29 * 37) > 13.98`, since label `a` has `b >= 5` columns at distinct split
primes. Multiply `h(u_a) = W_a/2` by `h(v_a) = W_F/2`. If `W_F <= V/4`, the first bound gives
`h_1 h_2 >= (W_F/4)(V/(2M))(1 - o(1))`, and `W_F >= M log(4M/e)` (Lemma 5.1), so `h_1 h_2 >= (1/8 - o(1)) V log M`.
If `W_F > V/4`, the second gives `h_1 h_2 > 6.99 V/8 > 0.87 V`. A contradiction from (5.1) requires
`c h_1 h_2 L < V/4 - log(MC)`, which fails under the stated conditions on `cL`.
For pairs: `D_ab` consists of `M/2` flipped rows with distinct primes, so `F_ab >= (M/2) log(2M/e)`.
Halving (1.3) and arguing as in Prop 3.1 (with `gamma = U_ab^(M/2) conj(Q_ab)^2`; the axis gap applies
because `|d_j| < M/4`) gives `(M/2)(W_a + W_b) + 2 F_ab >= V/2 - 2 log(MC/(2 sqrt 2))`. Hence, if
`F_ab <= V/8`, then `W_a + W_b >= V/(2M)(1 - o(1))` and `h_1 h_2 = (W_a + W_b) F_ab/4 >= (1/16 - o(1)) V log M`.
If `F_ab > V/8`, then `W_a + W_b > 27.9` gives `h_1 h_2 > 0.87 V`. QED

*Data* (`check_margins.py`, literal nearby-prime profiles, `b = 5`): `min_a W_a W_F/W` is `10.8`
(Paley-12), `15.5` (Paley-44), `21.3` (Paley-200), i.e. about `tau = 4 log M`. A Baker-type bound
needs it below `1/(4c)`.

**The three bounds side by side** (single-label form; need: right side `< (W+D)/4 - O(log MC)`):

| bound | `log(1/|Lambda_a|) <=` | at the best label | verdict |
|---|---|---|---|
| Liouville (= all characters) | `(M/4) W_a + F/2` | `W/4 + F/2 + O(W/M)` | fails by `F/2 - D/4` (`>= -O(W/M)` by Cor 3.3) |
| Baker type (proved) | `c (W_a/2)(F/2) L` | `>= (1/8 - o(1)) c (W+D) log M` | fails by factor `~ c log M` |
| Lang–Waldschmidt (conjectural) | `(1+eps)(W_a + F) + O(log M)` | `(1+eps)(W/M + F)` | succeeds iff `F < W/4` (pairs: `F < W/2`) |

## 6. Flipped primes are not pinned; what grid pinning would give

The suggested route was: pin the `M` flipped angles to a grid modulo `pi/(2M)` up to `O(Delta)`, then
use sector uniqueness and the distinctness of their norms. This section proves that the phase system
does not contain such a pinning. It also states exactly what pinning would suffice.

**Proposition U (no isolation).** Fully flipped class with capacity `<= B`, every label having an
unflipped column. Every zero-sum `lambda in Q^M` is `lambda = H_A n/M` for a unique `n in Q^A`
(`H_A` = nonconstant columns). Its column image `v = lambda^T S` equals, up to the orientation sign,
`n_a` on every unflipped column of label `a`, and `n_(a(x)) - 2 sigma_x (H_A n)_x / M` on the flipped
column of row `x`. Hence:

1. no nonzero `lambda` has column image supported on flipped columns (if `n_a = 0` for all `a` then `lambda = 0`);
2. if `v` is nonzero on at most `t` flipped columns, then `s = |supp n|` satisfies `s (t + B s) >= M`;
3. for a prime-order Paley core (`M = q + 1`, `q = 3 mod 4`), moreover `s >= (M/2 - t)/(B + 1)`.

So any rational relation that involves at most `t` flipped angles also involves the unflipped blocks
of at least `s(t)` labels: about `sqrt(M/B)` labels in general and `(M/2 - t)/(B+1)` for Paley.
Its Gaussian-integer height therefore contains the unflipped blocks of all these labels, which are
unknown. The targets
of the flipped angles are thus never fixed: they move with `Omega(sqrt M)` (Paley: `Omega(M)`)
unknown blocks.

*Proof.* `[1 | H_A]` is invertible and `1^T lambda = 0`, so `lambda = H_A n/M` with
`n = H_A^T lambda`. The formula for `v` is the computation of Section 0 with `H_A^T H_A = M I`.
(1) is immediate. (2): `supp(H_A n)` consists of rows `x` with `a(x) in supp n` (at most `B s` rows)
and rows counted in `t`, so `|supp H_A n| <= t + B s`. The Hadamard uncertainty principle gives
`|supp n| |supp H_A n| >= M`: indeed `||H_A n||_2^2 = M ||n||_2^2` and
`||H_A n||_inf <= ||n||_1 <= sqrt(s) ||n||_2`. (3): by Lemma 6.2 at least `M/2 - s` finite rows have
`(H_A n)_x != 0`; at most `B s` of them have `a(x) in supp n`; the rest are among the `t`. QED

**Lemma 6.2 (integer prime-field uncertainty for Paley).** Let `q = 3 mod 4` be prime and `H` the
normalised Paley-I matrix (`H(x,a) = -chi(a - x) - [a = x]` for `x, a in F_q`, row `infinity` all
ones). For every nonzero `n in Z^(F_q)` with `s = |supp n| <= (q+1)/2`,

```text
#{ x in F_q : (H n)_x != 0 }  >=  (q + 1)/2 - s.
```

*Proof.* Dividing by a power of `q` we may assume `n` is not `0 mod q`. Since
`chi(y) = y^((q-1)/2) mod q` (also for `y = 0`), `(H n)_x = -P(x) - n_x (mod q)` with
`P(X) = sum_a n_a (a - X)^d`, `d = (q-1)/2`. The coefficient of `X^k` in `P` is
`(-1)^k binom(d,k) sum_a n_a a^(d-k)`, and `binom(d,k)` is nonzero mod `q`. If `P = 0` mod `q`,
then `sum_a n_a a^m = 0` for `m = 0..d`. These are `d + 1 >= s` equations with an invertible
Vandermonde matrix on `supp n`, so `n = 0` mod `q` on its support, a contradiction. So `P` is a
nonzero polynomial of degree `<= d` with at most `d` roots. Every `x notin supp n` with
`(Hn)_x = 0` is a root, so at least `(q - s) - d` rows have `(Hn)_x != 0`. QED

This extends the ternary statement of item 393 to arbitrary integer vectors. *Check:*
`check_paley_uncertainty.py` exhausts 866,762 vectors (`q = 7, 11, 19, 23, 43`, supports up to 4,
entries up to 3) and spot-checks the congruence. The observed minimum of `#nonzero + s` is
`(q+1)/2 + 1` in every case.

**Lemma 6.3 (sector lemma).** Let `0 < 2 eps < pi/2`. A sector `{z : |arg z - t| <= eps, Norm z <= X}`
contains at most `1 + 2 eps X` primitive lattice points. Any two distinct ones satisfy
`Norm(u) Norm(v) >= sin(2 eps)^(-2) > (2 eps)^(-2)`.
*Proof.* Distinct primitive points in a sector of opening `< pi` are not collinear with `0`, so the
triangle `0 u v` is a nondegenerate lattice triangle: `|u||v| sin angle(u,v) = |det(u,v)| >= 1`.
Ordering the points by argument, the consecutive triangles are interior-disjoint, have area
`>= 1/2`, and lie in the sector, whose area is `eps X`. QED
(`check_sector_lemma.py`: 400 random rational sectors, exact enumeration.)

**Proposition 6.4 (grid pinning would suffice).** Let `omega_1, ..., omega_n` be pairwise
non-associate conjugate-primitive Gaussian integers. Suppose each `arg omega_x` is within `eps < pi/8`
of a set `T` of `G` directions mod `pi/2`. Then `sum_x log Norm omega_x >= (n - G) log(1/(2 eps))`.
In a fully flipped configuration, a pinning of the oriented flipped primes to `G <= M - 5`
directions mod `pi/2` with precision `eps <= 4 M C e^(-(W+D)/4)` is therefore impossible once
`W > 20 log(8MC)`. Indeed it would give `W >= W_F >= 5((W+D)/4 - log(8MC))`.
*Proof.* Rotate each `omega_x` by a unit into the sector of its grid point. Non-associate means the
rotated points are distinct primitive lattice points. By Lemma 6.3, each sector holds at most one
point of norm `< 1/(2 eps)`. QED

Pinning modulo `pi/(2M)` (`G = M` classes) would not suffice by counting alone. More importantly,
by Prop U the phase system pins no flipped angle to anything independent of the unflipped blocks.
**The suggested sector/grid mechanism is therefore unavailable, and Prop U gives the exact reason**
(Hadamard, respectively prime-field, uncertainty).

## 7. The heavy-flip regime `W_F >= W/2` (open even conditionally)

Theorem L uses label pairs, whose support weight is `W_a + W_b + F_ab ~ 2W/M + F/2`. Under
`(H_kappa)`, any rational row combination `lambda` whose column image takes boundedly many values
gives a form in boundedly many logarithms, with height sum equal to the **support weight**
`sum_{j in supp v} w_j`. The Lang–Waldschmidt route needs this weight to be `< W/(4 kappa)`.

**Proposition 7.1.** In a fully flipped prime-order Paley profile with capacity `B`, every rational
combination with `s = |supp n| <= (q+1)/2` has at least `M/2 - (B+1) s` flipped columns in its support.
If every flipped prime has `log q_x >= (1 - delta/4) W_F/M` and `W_F >= (1/2 + delta) W` with
`0 < delta <= 1/2`, then every combination with `s <= delta M/(4(B+1))` has support weight `> W/4`.
*Proof.* The first claim is Prop U(3). For the second, the weight is
`>= (1/2 - delta/4)(1 - delta/4)(1/2 + delta) W = (1/4 + 5 delta/16 - 11 delta^2/32 + delta^3/16) W > W/4`. QED

So in the balanced heavy regime only combinations involving a positive proportion of all labels
could work. *Data* (`check_support_minimum.py`; exhaustive over all `n in {-2..2}^11` for Paley-12, all
ternary `n` for Walsh-16, supports `<= 6` (ternary) and `<= 3` (`|n_a| <= 2`) for Paley-20; capacity
two, `b = 5`, unflipped weight 1, flipped weight `f`). For Paley-20 the whole Pareto frontier of
(unflipped, flipped) support counts consists of single labels `(3, 20)` and pairs `(7, 11)`. The
minimal support/`W` exceeds `1/4` as soon as `W_F/W >= 0.35`; the asymptotic threshold for pairs is
`1/2`. For Paley-12 (a small, highly symmetric order) one assignment admits support-6 vectors with
only 3 flipped primes, `(24, 3)`. These beat pairs only for `f >= 6` and never reach `1/4` in the
tested range: their support/`W` is `(24 + 3f)/(43 + 12f) > 1/4` for every `f`. For Walsh cores,
coset vectors `n = 1_(c+V)` (`c notin V`, `|V| = sqrt M`) have `H n` supported on `V^perp`. Their
column image takes at most five values, and their support weight is `~ sqrt M (W/M + (B+1) f) << W/4`
when the flipped primes have comparable sizes `~ f`. So under a five-logarithm Lang–Waldschmidt
hypothesis, flipped **Walsh** profiles with comparable flipped primes are excluded in every regime;
they are already excluded unconditionally by item 333. **Paley in the heavy regime is the precise
residual case, even conditionally.**

## 8. Where the obstruction sits, and what would finish

1. **Exact reduction.** Every label (or label-pair) condition of a Hadamard-core profile is a
   nonzero two-logarithm form over `Q(i)` with error `<= M Delta`. The block side has height
   `~ W/M`, the flipped side `F_a = W_F` (or `F_ab ~ W_F/2`). Prop 2.1-2.2.
2. **One-number methods are exhausted.** They are the integral characters and miss by exactly
   `F/2 - D/4 + O(W/M)` (Prop 3.1, Rem 3.2, Cor 3.3; Prop P of walsh-extraction.md). Fractional
   characters reduce to integral ones by taking the `e`-th power.
3. **Known transcendence bounds** (product of heights) miss by a factor `~ c log M` in every
   configuration (Prop B). The **sum-of-heights** (Lang–Waldschmidt) bound would prove the lemma
   whenever `W_F < W/2` (Thm L), in particular on the whole Prop-P barrier class. A constant-factor
   improvement of Liouville on the power side would prove it for `W_F < (eta/2) W` (L.6).
4. **Flipped primes are not individually pinned** (Prop U, Lemma 6.2). Grid pinning would suffice
   (Prop 6.4), but the system does not provide it.
5. **Residual case even under Lang–Waldschmidt:** Paley (quasi-random) cores whose flipped primes
   carry at least half of the weight. Bounded-support forms are blocked there (Prop 7.1), so one
   would need a form involving `Theta(M)` labels with few flipped primes. No such form is known, and
   none appears in the exhaustive small-order data except at order 12.
6. **Why this does not extend to general profiles.** The inversion `H^T/M` has denominator `M` and
   `l_1`-norm `1` per label, so the error grows only by `M`. For the general `k`-row profile of item
   532 (columns = all sign patterns `T` on the selected rows, with balanced block weights), take a
   zero-sum `lambda != 0` and a row `i` with `lambda_i != 0`. The column value at `T` is `2 lambda(T)`,
   and `lambda(T)`, `lambda(T xor {i})` cannot both vanish. So every form is nonzero on at least half
   of the patterns, i.e. it has support weight `>= (1/2 - o(1)) W`, and `(H_kappa)` gives only the
   (already critical) pair condition. The
   Hadamard-core class is special because its square structure produces forms whose support weight
   is far below `W/2`.

**Precise missing input (unconditional).** A lower bound for the two-logarithm form
`M Log(U/conj U) - 4 Log(Q/conj Q) + b Log i`, uniform over conjugate-primitive Gaussian integers
`U, Q`, of the form `log(1/|Lambda|) <= kappa_1 log Norm U + kappa_2 log Norm Q + O(log M)` with
`kappa_1 = o(M)` and `kappa_2 < W/(2F)`. For the single-label form `M Log(G/conj G) - 2 Log(P/conj P)`
the condition is `kappa_2 < W/(4F)`. Equivalently: a power of a small Gaussian integer cannot be
extremely close in direction to a square of a much larger one, with a polynomial (not
product-of-heights) dependence on the larger height. This is a special case of the Lang–Waldschmidt
conjecture and is not known.

## 9. Relation to the notes and to walsh-extraction.md (novelty check)

* Searched the 1,170 research notes for `Waldschmidt`, `Lang`, `two logarithm`/`two-logarithm`,
  `Laurent`, `moving target`, `Kummer`. Two-logarithm bounds (Matveev) appear in items 276-277,
  `minimal_prime_rank_width_bound.md` and `hadamard_augmentation_short_character.md`, always with one
  varying number and `log(-1)`, i.e. fixed targets. `Waldschmidt` appears only as a continued-fraction
  reference. `Lang` appears only for Bombieri–Lang/Lang's fixed-place theorem.
  No note reduces the flipped pinning to a two-number form, uses the Lang–Waldschmidt conjecture, or
  quantifies the Baker gap. Lemma 6.2 extends item 393 from ternary to integer vectors. Prop 3.1 is
  a column character, contained in Prop P / items 391-393.
* Corrections to Section 10 of walsh-extraction.md: (i) the pinning system is necessary, not
  equivalent (as the referee said); (ii) at least for `W_F < W/2`, the obstruction is not "`M - 1`
  simultaneous approximations that no tool controls jointly". A **single** label pair already gives a
  contradiction under a classical (open) two-logarithm conjecture. The real gap is product-of-heights versus
  sum-of-heights for two logarithms, a factor `~ log M` (Prop B).
* The Sylvester case of the lemma is item 333 / Theorem C of walsh-extraction.md (unconditional).

## 10. Files (`round4/flipped_checks/`)

* `fcommon.py`: exact Gaussian arithmetic, deterministic Miller–Rabin (validated against trial
  division below `2*10^5`), Paley-I and Sylvester matrices (Hadamard property asserted), profile builder.
* `check_two_log_identity.py`: identity (1.4) and its pair version, nonvanishing (axis gap),
  minimal-polynomial heights; exact; negative control. Output `check_two_log_identity_out.txt`. PASS.
* `check_margins.py`: height bookkeeping on literal nearby-prime and heavy-flip profiles
  (Cor L.4, Prop B data). Output `check_margins_out.txt`.
* `check_paley_uncertainty.py`: Lemma 6.2, exhaustive. Output `check_paley_uncertainty_out.txt`. PASS.
* `check_sector_lemma.py`: Lemma 6.3, exhaustive on random sectors. Output `check_sector_lemma_out.txt`. PASS.
* `check_support_minimum.py`: Section 7 data (Pareto frontiers of support counts). Output `check_support_minimum_out.txt`.
* `check_lw_scale.py`: height normalisation of `(H_kappa)` (Section 4). Output `check_lw_scale_out.txt`.

The repository `/home/user/Jarnik` was not modified.
