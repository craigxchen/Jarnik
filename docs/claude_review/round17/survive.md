# Round 17: survival attempt for CIRCLE(q) (lower bound / dual certificate)

## 0. Summary and status

**No growth improvement and no uniform bound.** Nothing here moves `M <= (2/3 + eps) log R/loglog R`. CIRCLE(q) is
**not decided**: no regime witness was constructed, so neither a quantitative "no" nor a cap is proved.

**Proved (new):**
1. **Liouville-free equivariant bound (Thm 3.1).** For `sdc_q` boxes and `m < q`, an integral certificate `c in V(A,m)`
   with `<c, M_i> != 0` (`M_i = hat k(1)^i hat k(-1)^{m-i}`) gives, for **every** `theta` and every `H >= |c|_inf`,
   `kappa_H(hat u_theta) >= 2^{-m} binom(m,i)/(|c|_1 q^m)`. The multiplier group realises the Galois group, so one evaluates
   at the largest conjugate (`>= 1` by the norm) instead of dividing by the others. Translation averaging (valid for
   `m < q`) makes the bound uniform on the circle with no `(m+1)^{-1/2}` loss. This removes the referee's
   `e^{-Theta(q^2)}` Liouville obstruction.
2. **Gap principle (Thm 3.4).** On the Fourier circles `kappa_H` is either identically `0` or at least
   `2^{-m}/(|A| H q^m)`. So CIRCLE(q) reduces to two **exact** questions: `C_1 subset TC_{H'}` (cap) versus
   `C_1 not subset TC_{H_*}` with `ln H_* = (1/2 - o(1)) q tau` (survival). Only a height window between them is left.
3. **The passage to `kappa_H` is free on the circles (§5).** The bound is a `kappa_H` bound from the start, and the
   record's `H' = e^{m'W/40}` exceeds every useful witness height.
4. **Survival criterion (Cor 3.2-3.3).** CIRCLE(q) survives at every barrier scale if there is an integral witness of
   degree `m'` with `ln|c|_1 < q tau_D/2 - m' ln 2q - ln(Cq/(m'+1))`. At `q = 3251` the budget is `4.63 * 10^5 = 142.5 q`,
   against `q ln 3 = 3572` for a witness with `O(1)` entries on all of `A` (`regime_budget_out.txt`).
5. **What the witness cannot be (Props 4.1, 4.2, Remark 4.3).**
   * It cannot be product-type (zonotopes, sub-cube tensors, or their convolutions): these have order
     `<= 8jq/(q+1) < m'`. This uses the exact Paley-Parseval identity `sum_a v(P'_a)^2 = (q+1)(K^2 + |v|_2^2)/4`.
   * It cannot be supported on a union of coordinate sub-cubes in `A`: these have order `<= 4j`.
   * Counting cannot produce it: `ln r_0 >= 1073` at `q = 3251`, so Bombieri-Vaaler heights are `exp(e^{1073})`.

**Exact checks:**
* the `AGL` and Galois identities (§6.1);
* the record's `sdc7` datum: the unique certificate is `AGL`-invariant with `<c, M_4> = -28224 in Z`, and it is
  LP-optimal (§6.2);
* toy witnesses in `Z[zeta]`, with the bound within a factor 7 (for `L <= 2`) to 33 (for `(7,3)`) of the LP truth (§6.3; the LP values are evidence).

**Open, stated sharply (§7):** (W_q), the existence at one regime `q` of a non-product, sign-correlated integral
certificate of degree `m'` that does not vanish on `C_1`, of height `e^{O(q)}`. Per the round's rule, route 5 stops
unless the cap side succeeds.

## 1. Setting, notation, the statement to prove

**Labels.** *Proved* = proved here with all quantifiers (or a one-line consequence of a refereed record result, named).
*Exact* = integer / exact-rational / exact finite-field computation, script named. *Evidence* = everything else
(floating point, finite search, toy sizes); nothing proved uses it. Record = `docs/claude_review/` (read-only).
Scripts: `survive_checks/` (`x.py` writes `x_out.txt`).

**Data.** `q` prime, `q = 3 mod 4`, `chi` the quadratic character, `zeta = e^{2 pi i/q}`. Relative `sdc_q` patterns (base row
`inf` omitted), `a in F_q`:

```text
P_a = {a} u {x : chi(a-x) = 1},      P'_a = {a} u {x : chi(a-x) = -1}        (2q patterns, |P| = (q+1)/2).
```

For an integer `L >= 0` put `A = A_q(L) = {k in {-1,0,1}^{F_q} : |k(P)| <= L for all 2q patterns}`; the regime box is
`L = j`. Here `n = q + 1`, `d = n - 1 = q` (coordinates = non-base rows), `r = 2q`, and as in the record
`w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`. Regime: `q >= 3251`; there `m' < q`
(R16-J2) and `m' <= 8j + 8`.

**Certificates.** For `m >= 1`: `V(A,m) = {c in R^A : sum_k c_k P(k) = 0 for all P in Poly_{<m}}`; `<c,f> := sum_k c_k f(k)`;
`Phi_c(u) = <c, l_u^m>`, `l_u(k) = <k,u>`; `kappa_c(u) = |Phi_c(u)|/(|c|_1 d^m)`; `kappa_V = sup_{c in V(A,m)} kappa_c`;
`V_H = span_R{c in V(A,m) cap Z^A : |c|_inf <= H}`, `kappa_H = sup_{c in V_H} kappa_c`. Duality (nearstrata Prop 2.1):
`d^m kappa_V(u) = Cheb_A(u) := min_{deg P < m} max_A |l_u^m - P|`. For a complex-valued `f` on `A` put
`Cheb^C_A(f) = min_{P in Poly_{<m} (x) C} max_A |f - P|`. For real `c in V(A,m)` and every complex `P` of degree `< m`,
`|<c,f>| = |<c, f - P>| <= |c|_1 max_A |f - P|`, so

```text
(1.1)        Cheb^C_A(f) >= |<c,f>| / |c|_1       for every real c in V(A,m), c != 0.
```

**The circle.** `u_theta(x) = cos(theta + 2 pi x/q)` (`x in F_q`, `u_inf = 0`), `C_1 = {u_theta}`. With
`l_s(k) := sum_x k_x zeta^{s x}` (`= hat k(s)`) one has `l_{u_theta} = (e^{i theta} l_1 + e^{-i theta} l_{-1})/2`, hence

```text
(1.2)   l_{u_theta}^m = 2^{-m} sum_{i=0}^m binom(m,i) e^{i theta (2i-m)} M_i,      M_i := l_1^i l_{-1}^{m-i},
        Phi_c(u_theta) = 2^{-m} sum_{i=0}^m binom(m,i) e^{i theta (2i-m)} <c, M_i>.
```

So `Phi_c` vanishes identically on `C_1` iff `<c, M_i> = 0` for all `i`, and `C_1 subset TC_V(A,m)` iff every `M_i`
lies in `Poly_{<m}|_A (x) C` (Cor 3.3 of barrier.md, here directly from (1.1)-(1.2) and duality).

**Barrier scale and thresholds** (barrier.md §3.3 with referee corrections). `W = r tau = 2 q tau`, `tau = K ln n >= tau_D`,
`eps_0 = C d e^{-W/4}/(m'+1) = C q e^{-q tau/2}/(m'+1)`, `ln H' = m' W/40 = m' q tau/20`. CIRCLE(q) asks whether
`inf_theta kappa_{H'}(hat u_theta) <= T_cap := eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)})` (cap) or not. Since MIN^pair needs
`kappa_H > eps_0` and `T_cap < eps_0`, the natural survival statement is the stronger
`inf_theta kappa_{H'}(hat u_theta) > eps_0`; this is what is targeted below. Normalisation:
`hat u = u/|u|_inf`, `|u_theta|_inf in [cos(pi/q), 1]`, and `kappa(hat u) = kappa(u)/|u|_inf^m >= kappa(u)`, so lower bounds
for `kappa(u_theta)` transfer to `hat u_theta` for free.

**Symmetry.** `AGL(1,q)` acts on `F_q`, on `k` by `(g k)_x = k_{g^{-1} x}` and on measures by `(g c)_k = c_{g^{-1} k}`. Write
`tau: x -> x + 1` and `mu_s : x -> s x` (`s in F_q^*`).

## 2. Dual form of kappa_V and kappa_H on the circle; isotypic reduction

**Lemma 2.1 (proved).** (i) Every `g in AGL(1,q)` maps the pattern set onto itself; hence `gA = A`, `g` preserves
`Poly_{<m}|_A`, `V(A,m)`, integrality, `|c|_1`, `|c|_inf`, `V_H`, and `Phi_{gc}(u) = Phi_c(g^{-1}u)`.
(ii) `l_s(tau k) = zeta^s l_s(k)` and `l_1(mu_s k) = l_s(k)`. Hence `M_i(tau k) = zeta^{2i-m} M_i(k)` and
`M_i(mu_s k) = sigma_s(M_i)(k)`, where `sigma_s in Gal(Q(zeta)/Q)`, `zeta -> zeta^s`, acts on the coefficients
(`k` is integral).
(iii) For every measure `c` on `A`: `<tau^a c, M_i> = zeta^{a(2i-m)} <c, M_i>`, and for rational `c`,
`<mu_s c, M_i> = sigma_s(<c, M_i>)`.

*Proof.* (i) `P_a + b = P_{a+b}`, `P'_a + b = P'_{a+b}`. For `s in F_q^*`: `s P_a = {sa} u {y : chi(a - y/s) = 1} =
{sa} u {y : chi(s)^{-1} chi(sa - y) = 1}`, which is `P_{sa}` if `chi(s) = 1` and `P'_{sa}` if `chi(s) = -1`; likewise for `P'_a`.
`(gk)(gP) = k(P)`, so `gA = A`; the rest is bookkeeping (coordinate permutations). (ii) `l_s(tau k) = sum_x k_{x-1} zeta^{sx} =
zeta^s l_s(k)`; `l_1(mu_s k) = sum_x k_{s^{-1}x} zeta^x = sum_y k_y zeta^{sy} = l_s(k)`, and `l_{-1}(mu_s k) = l_{-s}(k)`;
`sigma_s(M_i) = l_s^i l_{-s}^{m-i}`. (iii) `<g c, f> = sum_k c_{g^{-1}k} f(k) = sum_k c_k f(gk)`; apply (ii) and use that `sigma_s`
commutes with rational linear combinations. QED

(The multiplier part of (i) is the referee's AGL(1,q) invariance, referee_barrier §4, re-proved here because the proof
of Theorem 3.1 uses the explicit Galois compatibility (ii).) Exact check: `survive_checks/sym_check.py`.

**Isotypic decomposition (proved).** `V(A,m) (x) C` splits into the eigenspaces `V_t = {c : tau c = zeta^t c}`,
`t in F_q`. By (1.2) the circle sees `c` only through the `m+1` numbers `<c, M_i>`, and by (iii) `M_i` pairs only with
`V_t` for `t = 2i - m`. **Because `m < q` and `q` is an odd prime, the `m+1` exponents `2i - m`
(`0 <= i <= m`) are pairwise distinct mod `q`** (`q | 2(i - l)` forces `q | i - l`, and `|i - l| <= m < q`). So the
projection

```text
(2.1)      pi_i(c) := (1/q) sum_{a in F_q} zeta^{-a(2i-m)} tau^a c        satisfies   <pi_i(c), M_l> = [l = i] <c, M_i>,
           Phi_{pi_i(c)}(u_theta) = 2^{-m} binom(m,i) e^{i theta(2i-m)} <c, M_i>      (constant modulus in theta).
```

This is the "matching translation eigenspace" of the task: `pi_i(c)` lies in `V_{2i-m}` (the `zeta^{2i-m}`-eigenspace of
`tau` acting on measures), and its `Phi` has **constant modulus on the whole circle**. `pi_i(c)` is complex; it is used
only through (1.1) and through the real measures `tau^a c` it averages.

**Dual form of `kappa` on the circle (proved).** By (1.1), (2.1) and duality, for every `theta` and `i`:

```text
(2.2)   Cheb_A(u_theta) >= 2^{-m} binom(m,i) Cheb^C_A(M_i),
        Cheb^C_A(M_i) = sup { |<c, M_i>| / |c|_1 : c in V(A,m) real }   ( = the two-measure form: c = c_+ - c_-,
        c_+, c_- measures on A of equal mass with equal moments of order < m, and <c_+ - c_-, M_i> = "the m-th moment
        difference along the complex direction psi_1, isotypic part i").
```

*Proof of the inequality.* Let `P` be real of degree `< m` with `max_A |l_{u_theta}^m - P| = E := Cheb_A(u_theta)`. Since `tau` is
a box automorphism and `l_{u_theta}^m o tau^{-a} = l_{tau^a u_theta}^m = l_{u_{theta - 2 pi a/q}}^m`, each
`P o tau^{-a}` approximates `l_{u_{theta-2pi a/q}}^m` within `E` on `A`. By (1.2) and the distinctness of the exponents,
`(1/q) sum_a e^{-i(theta - 2 pi a/q)(2i-m)} l_{u_{theta-2pi a/q}}^m = 2^{-m} binom(m,i) M_i`. The same combination of the `P o tau^{-a}`
is a complex polynomial of degree `< m` within `E` of it on `A`. The equality is finite-dimensional LP duality over `C`
restricted to real certificates: `>=` is (1.1); `<=` is not needed and not claimed. QED

The `(m+1)^{-1/2}` loss of barrier Thm 3.2 does not occur in (2.2): the lower bound is uniform in `theta` directly.

## 3. AGL(1,q)-equivariant certificates: Liouville's loss removed

The referee's objection to a "no" answer (referee_barrier §4, Thm 3.5 item) was: exact non-membership of `C_1` in
`TC` plus Liouville at the algebraic point `u_0` gives only `kappa >= e^{-Theta(q^2)}`, because `Phi_c(u_0)` lies in a field
of degree `(q-1)/2` and Liouville divides by all the other conjugates. The point of this section is that **the Galois
conjugates are realised by box automorphisms**, so one may always evaluate at the *largest* conjugate, and Liouville's
loss disappears completely.

**Theorem 3.1 (equivariant witness bound; proved).** Let `q` be a prime, `q = 3 mod 4`, `L >= 0`, `A = A_q(L)`, and
`1 <= m < q`. Let `c in V(A,m) cap Z^A` and `0 <= i <= m` with `alpha := <c, M_i> != 0`. Then:
1. `alpha in Z[zeta]`, and there is `s in F_q^*` with `|sigma_s(alpha)| >= 1`.
2. For that `s` and **every** `theta in R` there is `a in F_q` such that the integral certificate `c_* := tau^a mu_s c`
   (which lies in `V(A,m) cap Z^A` and has `|c_*|_1 = |c|_1`, `|c_*|_inf = |c|_inf`) satisfies
   `|Phi_{c_*}(u_theta)| >= 2^{-m} binom(m,i)`.
3. Consequently, for every `theta` and every `H >= |c|_inf`,

```text
(3.1)      kappa_V(hat u_theta) >= kappa_H(hat u_theta) >= kappa_H(u_theta) >= 2^{-m} binom(m,i) / ( |c|_1 d^m ),      d = q.
```

*Proof.* (1) `M_i` has coefficients in `Z[zeta]` and `k in Z^q`, so `M_i(k) in Z[zeta]`; `c` is integral, so `alpha in Z[zeta]`.
`alpha` is a nonzero algebraic integer, so its norm `N(alpha) = prod_{s in F_q^*} sigma_s(alpha)` is a nonzero rational
integer, `|N(alpha)| >= 1`, and the largest factor has modulus `>= 1`.
(2) Put `c' = mu_s c`. By Lemma 2.1(i) `c' in V(A,m) cap Z^A` with the same norms, and by Lemma 2.1(iii)
`<c', M_i> = sigma_s(alpha)`. Let `f_a(theta) = Phi_{tau^a c'}(u_theta)`. By (1.2) and Lemma 2.1(iii),
`f_a(theta) = 2^{-m} sum_l binom(m,l) e^{i theta(2l-m)} zeta^{a(2l-m)} <c', M_l>`, so by (2.1)
`(1/q) sum_a zeta^{-a(2i-m)} f_a(theta) = 2^{-m} binom(m,i) e^{i theta(2i-m)} sigma_s(alpha)`. A mean of `q` numbers of modulus
`|f_a(theta)|` has modulus at most `max_a |f_a(theta)|`. Hence `max_a |f_a(theta)| >= 2^{-m} binom(m,i) |sigma_s(alpha)| >= 2^{-m} binom(m,i)`.
(3) `c_* in V_H` because it is integral with `|c_*|_inf <= H`; so `kappa_H(u_theta) >= kappa_{c_*}(u_theta)`, which is the right-hand
side of (3.1). `kappa_H <= kappa_V` since `V_H subset V(A,m)`; the normalisation step is §1. QED

*Remarks.* (a) Nothing about `theta` is assumed: the bound is uniform on the circle with no `(m+1)^{-1/2}` loss and no
Diophantine input. (b) The bound is insensitive to which isotypic index `i` carries the witness, up to `binom(m,i) >= 1`.
(c) The same proof works verbatim for every circle `C_b` (`u(x) = cos(theta + 2 pi b x/q)`, `b in F_q^*`), by
`mu_b`-equivariance, and for `kappa` computed with any `AGL(1,q)`-stable family of certificates.
(d) For `m` even and `c` translation-invariant, only `i = m/2` pairs, and if moreover `c` is `F_q^*`-invariant then
`<c, M_{m/2}>` is a rational integer: this explains the record's exact datum `sdc7`, `L = 2`, `m = 8`: `dim V = 1`,
`Phi_c = -7717.5 = -1975680/2^8` constant on every circle (referee_barrier §6), and (3.1) is checked there exactly in §6.

**Definition (witness height).** `lambda_q(L, m) := min { |c|_1 : c in V(A_q(L), m) cap Z^A, Phi_c not identically 0 on C_1 }`,
with `min(empty) = +infinity`. By (1.2), `lambda_q < infinity` iff `C_1 not subset TC_V(A,m)` (if `C_1` is not in `TC_V`, some real
`c in V` has `Phi_c != 0` somewhere on `C_1`; `V` is defined over `Q`, so a rational, hence an integral, such `c` exists).

**Corollary 3.2 (survival criterion for CIRCLE(q); proved).** In the regime (`L = j`, `m = m'`, `W = 2 q tau`), for every
`tau > 0` with

```text
(3.2)        ln lambda_q(j, m') < q tau/2 - m' ln(2q) - ln( C q/(m'+1) ),
```

one has `inf_theta kappa_{H'}(hat u_theta) > eps_0 > T_cap`, i.e. the circle `C_1` (and every `C_b`) carries **no**
`eps_0`-near-zero of `kappa_{H'}`, at every point and with the record's height `ln H' = m' W/40`. In particular the
circle test of barrier Thm 3.5 does not cap route 5 at that `tau`.

*Proof.* Take a minimising witness `c`. Then `|c|_inf <= |c|_1 = lambda_q < e^{q tau/2} <= e^{m' q tau/20} = H'` (`m' >= 10`), so
(3.1) applies with `H = H'`, `binom >= 1`: `kappa_{H'}(hat u_theta) >= 2^{-m'} q^{-m'}/lambda_q > C q e^{-q tau/2}/(m'+1) = eps_0` by (3.2). QED

**The passage to `kappa_H` is free on the circle.** (3.1) is a `kappa_H` statement from the start: the certificates
used are AGL-images of one integral witness, and the record's height `H' = e^{m' W/40}` exceeds the `|c|_1` allowed by
(3.2) by a factor `e^{(m'/10 - 1) q tau/2}`. So neither (P1) nor (P2) of reduction §4 is needed on the circle: the only
question left is the size of `lambda_q`.

**Corollary 3.3 (the quantitative test is the qualitative test plus a height; proved).** At the regime parameters:
1. If `C_1 subset TC_V(A_q(j), m')`, route 5 is capped (barrier Thm 3.5, H_light assumed) at every barrier `tau`.
2. If `C_1 not subset TC_V(A_q(j), m')`, the circle test is passed at every `tau >= tau_*(q) := (2/q)( ln lambda_q(j,m')
   + m' ln(2q) + ln(Cq/(m'+1)) )`.
3. So CIRCLE(q), at all barrier scales `tau >= tau_D` at once, is **decided in the survival direction** as soon as one
   exhibits an integral `c in V(A_q(j), m')` with `Phi_c|_{C_1} != 0` and `ln |c|_1 < q tau_D/2 - m' ln(2q) - ln(Cq/(m'+1))`.
   The only scales at which the test can still cap with `C_1 not subset TC_V` are `tau_D <= tau < tau_*(q)`, which is a
   non-empty range only if `ln lambda_q > q tau_D/2 - O(m' ln q)`, i.e. only if **every** witness has
   `ln |c|_1 >= (tau_D/2 - o(1)) q`.

*Proof.* (1) A point of `TC_V` has `kappa_V = 0 >= kappa_H`; Thm 3.5. (2), (3): Corollary 3.2. QED

**Theorem 3.4 (gap principle at every height; proved).** For `1 <= m < q`, every `L` and every `H >= 1`, exactly one of
the following holds for `A = A_q(L)`:
* `C_1 subset TC_H(A,m)`, i.e. `kappa_H = 0` on all of `C_1` (and then on every `C_b`);
* `inf_theta kappa_H(hat u_theta) >= 2^{-m} / (|A| H d^m)`.

*Proof.* If some `c in V_H` has `Phi_c` not identically zero on `C_1`, then, `V_H` being spanned by integral vectors of
height `<= H` and `c -> Phi_c|_{C_1}` being linear, one of those spanning vectors, `c_0`, is a witness with
`|c_0|_1 <= |A| H`; apply Theorem 3.1. Otherwise `kappa_H = 0` on `C_1`, and on `C_b` by `mu_b`-invariance. QED

So on the Fourier circles **`kappa_H` has no small positive values**: it is either identically zero or at least
`2^{-m}/(|A| H d^m)`. At the regime this turns CIRCLE(q) into two *exact* membership questions:

```text
  C_1 subset TC_{H'}(A, m')                ==>  route 5 capped (H_light assumed; barrier Thm 3.5),   ln H' = m' q tau/20;
  C_1 not subset TC_{H_*}(A, m')           ==>  circle test passed at tau,   ln H_* := q tau/2 - m' ln(2q) - ln(Cq/(m'+1)) - ln|A|.
```

Only `C_1 subset TC_{H_*}` together with `C_1 not subset TC_{H'}` is left undecided: a "height window"
`ln H in [ (1/2 - o(1)) q tau , (m'/20) q tau ]` (at `q = 3251`, `tau = tau_D`: `[4.6 * 10^5, 9.9 * 10^7]`).

This answers the referee's objection to "exact non-membership does not settle it": **it does, provided the witness of
non-membership has `ln |c|_1 <= (tau_D/2 - o(1)) q`.** A witness with `O(1)` entries supported on all of `A` has
`ln |c|_1 <= q ln 3 + O(1)`, so the budget allows entries up to `e^{(tau_D/2 - ln 3 - o(1)) q}`.

## 4. What the witness must look like: product-type, sign-free and counting witnesses all fail

By §3 the whole survival question is the existence of one integral witness of moderate height. This section proves
that none of the certificate types in the record can be that witness, and that counting (Siegel / Bombieri-Vaaler)
cannot produce it either. So the witness, if it exists, is the "first non-zonotopal certificate" of README-16 §7.2,
and nothing weaker will do.

**Lemma 4.0 (range of a pattern sum; proved).** For `v in Z^{F_q}` put `S_a(v) = v(P'_a)`, `S^-_a(v) = v(P_a)`, `K = sum_x v_x`.
Then `sum_a S_a(v)^2 = sum_a S^-_a(v)^2 = (q+1)(K^2 + |v|_2^2)/4`, and
`h'(v) := sum_{P} |v(P)| >= (q+1)|v|_2^2/(2|v|_1) >= (q+1)/2` for `v != 0`.

*Proof.* `P'_a = a + D'` with `D' = {0} u QR` (because `chi(a-x) = -1` iff `x - a in QR`, as `-1` is a non-residue), and
`P_a = a + D`, `D = {0} u NQR`. So `S = v` correlated with `1_{D'}`, and `|hat S(b)| = |hat v(b)| |hat 1_{D'}(b)|` with
`|hat 1_{D'}(0)| = (q+1)/2` and, for `b != 0`, `hat 1_{D'}(b) = 1 + (chi(b) g - 1)/2` (`g = i sqrt q`, Gauss sum), of modulus
`sqrt(q+1)/2`; the same for `D`. Parseval: `sum_a S_a^2 = (1/q)[((q+1)/2)^2 K^2 + ((q+1)/4)(q|v|_2^2 - K^2)] = (q+1)(K^2+|v|_2^2)/4`.
Then `h'(v) >= (sum_a S_a^2 + sum_a (S^-_a)^2)/max_a(|S_a|,|S^-_a|)` and `|S_a|, |S^-_a| <= |v|_1`; finally `|v|_2^2 >= |v|_1`
for integer `v`. QED (Exact check of the identity on 66990 vectors at `q = 3251`: `zono_cost.py`.)

**Proposition 4.1 (no product-type witness; proved).** Let `c = c_1 * c_2 * ... * c_t` (convolution of finitely
supported measures on `Z^q`) be supported in `A_q(j)`, where each factor is either
* (a) a first difference `delta_{k_l + v_l} - delta_{k_l}` (`v_l in Z^q \ 0`; order `m_l = 1`), or
* (b) a measure of order `m_l` whose support is exactly a product set `k_l + prod_{x in Y_l} T_x` (`T_x subset {-1,0,1}`,
  `|T_x| >= 2`, other coordinates fixed). This includes the tensor certificates `(x)_x mu^{alpha_x}` that span `V(cube, m)`.

If `Phi_c` is not identically zero, then `c` has order `M = sum_l m_l <= 8jq/(q+1) < 8j <= m'`. In particular no
zonotope certificate, no product of sub-cube certificates and no convolution of such is a witness in `V(A_q(j), m')`.

*Proof.* (i) *Orders and Phi.* If `c_l` annihilates `Poly_{<m_l}` then `c` annihilates `Poly_{<M}`, and expanding
`l_u(k_1 + ... + k_t)^M` gives `Phi_c(u) = (M!/prod_l m_l!) prod_l <c_l, l_u^{m_l}>` (a factor that receives degree `< m_l`
kills the term). So `c` has exact order `M` when `Phi_c != 0`, and `M = m'` for a witness.
(ii) *Extreme points add.* Let `lambda` be a linear form taking distinct values on each finite support `S_l = supp c_l`
(true for `lambda = sigma_P + eps rho` with `rho` generic and `eps > 0` small). The `lambda`-maximum of `S_1 + ... + S_t` is
attained only at the sum of the `lambda`-maximisers `k_l^+ in S_l`, so the coefficient of `c` there is
`prod_l c_l(k_l^+) != 0`: this point lies in `supp c`. Letting `eps -> 0`, `k_l^+` maximises `sigma_P := k -> k(P)` on
`S_l`. With the minimisers likewise, `range_P(supp c) >= sum_l range_P(S_l)`; and `range_P(supp c) <= 2j` on `A`.
(iii) *Counting.* Summing over the `2q` patterns: `sum_l sum_P range_P(S_l) <= 4jq`. Type (a): `range_P(S_l) = |v_l(P)|`,
so `sum_P range_P(S_l) = h'(v_l) >= (q+1)/2 = (q+1) m_l/2` (Lemma 4.0). Type (b): every function on a product grid is a
polynomial of degree `<= sum_x (|T_x| - 1)`, so `m_l <= sum_x (|T_x| - 1) <= sum_x range(T_x)`; and `range_P(S_l) =
sum_{x in P cap Y_l} range(T_x)`, so `sum_P range_P(S_l) = (q+1) sum_x range(T_x) >= (q+1) m_l` (each `x` lies in `q+1`
patterns). Hence `(q+1) M/2 <= 4jq`. Finally `m' = ceil(4w) >= 4w >= 8 floor(w/2) = 8j > 8jq/(q+1)`. QED

*Evidence that 4.1 is far from tight:* the true minimum of `h'` over small `v` at `q = 3251` is `q + 1` (support `<= 3`:
exhaustive up to `AGL`; support `4..8`: 15000 random vectors, minimum `>= 1.48(q+1)`; `zono_cost_out.txt`). With
`h' >= q+1` the bound becomes `M <= 4jq/(q+1)`: product-type witnesses stop at **half** the needed degree.

**Proposition 4.2 (no sign-free witness; proved).** Let `F` be a family of subsets `Y subset F_q` with
`{-1,0,1}^Y subset A_q(j)` (identifying `{-1,0,1}^Y` with vectors supported in `Y`), and `B_F = union_{Y in F} {-1,0,1}^Y`.
Then every `Y in F` has `|Y| <= 2j`, every function on `B_F` is the restriction of a polynomial of degree `<= 4j`, and
`V(B_F, m) = 0` for `m > 4j`. Since `m' > 4j`, no witness is supported on such a union (in particular none on the
ternary Hamming ball `{|supp k| <= j} subset A`).

*Proof.* `1_Y in A` gives `|Y cap P| <= j` for all `P`. If `Y = F_q` then `|P_a| + |P'_a| = q + 1 > 2j`, impossible; so pick
`a notin Y`: `|Y cap P_a| + |Y cap P'_a| = |Y|`, hence `|Y| <= 2j`. Let `F^` be the down-closure of `F` (still admissible), so `B_F = {k : supp k in F^}`. The monomials
`k^alpha` with `alpha in {1,2}^Z`, `Z in F^` (support exactly `Z`) are linearly independent on `B_F`: in a vanishing
combination take `Z` minimal with a nonzero coefficient and evaluate at the `2^{|Z|}` points with support `Z`; only the
monomials supported on `Z` survive, and on `{+-1}^Z` they are the distinct characters `prod_{alpha_x = 1} k_x`. Their number
is `sum_{Z in F^} 2^{|Z|} = |B_F|`, so they form a basis of `R^{B_F}`, of degrees `<= 2|Z| <= 4j`. A measure on `B_F`
annihilating `Poly_{<m}` with `m > 4j` annihilates this basis, hence is `0`. QED

**Remark 4.3 (counting cannot produce the witness; proved numbers).** `V(A,m')` is huge and has a basis of small
height in the Bombieri-Vaaler sense, but the witnesses form the complement of `V' = {c : <c, M_i> = 0 for all i}`,
whose codimension in `V` is at most `(m'+1)(q-1)`. Bombieri-Vaaler controls only the product of the heights of a
basis, `prod_l |c_l|_inf <= U^{r_0/2}` with `r_0 = dim Poly_{<m'}|_A`, so it puts a short vector outside `V'` only at
height `<= U^{r_0/(2 codim V')}`. By the Hamming-ball argument of Prop 4.2 (restriction to `{|supp k| <= j}` is onto),
`r_0 >= sum_{s <= j} binom(q,s) 2^s`, i.e. `ln r_0 >= 1073.3` at `q = 3251` (`r0_bounds_out.txt`; upper bound `3229.4`). So the
counting height is `exp(e^{1073})`, against the budget `ln|c|_1 < 4.6 * 10^5` of Corollary 3.2. Equivalently: the
witness question **is** the question whether `C_1 subset TC_{H_*}(A, m')` for `ln H_* ~ (tau_D/2) q`, and reduction Thm 4.1
(codimension of `V_H`) is `e^{1073}` orders of magnitude too weak to decide it. Siegel's lemma *outside a union of
sublattices or a hypersurface* (Fukshansky, Monatsh. Math. 2006, arXiv math/0409374 [L, abstract level]) bounds the
height of a vector avoiding the bad set in terms of the height (covolume) of the ambient lattice, i.e. of the same
`U^{r_0/2}`-type quantity; read at abstract level it does not change this verdict.

**Remark 4.4 (the `L^2` route does not supply a witness either; a check of one method, not a theorem about `kappa`).**
A witness-free route would be `Cheb_A(u) >= Psi_nu(u)^{1/2}` for a measure `nu` on `A`, with `nu` a product measure
conditioned on `A` and the tail handled by hypercontractivity. The budget of Cor 3.2 does not help, because that
chain closes only when the tail term is below the main term (all or nothing). For the `p`-biased measure
(`P(k_x = +-1) = p/2`) the chain needs `G(p) = m' ln 3 + (m'/2) ln(2/p) + (1/2) ln nu_p(A^c) < -ln 2`. With the exact lower
bound `nu_p(A^c) >= binom(|P|, j+1)(p/2)^{j+1}(1-p)^{|P|-j-1}`, `min_p G(p) = +1823` at `q = 3251` (`biased_l2_out.txt`; grid
`p in [10^{-20}, 1]`). The reason is structural: the tail decays like `p^{j}` and the hypercontractive cost grows like
`p^{-m'/2} = p^{-4j}`. This is the record's "regime slabs are thick" (the cut is only polynomially rare under the
uniform measure) in its biased form.

**Remark 4.5 (heuristic, not used).** A witness `c = c_+ - c_-` is a pair of measures on `A` with equal moments of all
orders `< m'`: a multidimensional Prouhet-Tarry-Escott configuration inside `A`. Prouhet-type solutions are exactly the
zonotopes of Prop 4.1. So (W_q) asks for a non-Prouhet PTE configuration, of degree `~ 8j`, inside a regime box, of
height `e^{O(q)}`, that is not orthogonal to the circle.

**What a witness must look like (summary).** Combining 4.1-4.3 with Lemma 2.1: a witness can be taken
`tau`-isotypic (in `V_{2i-m'}` for one `i`; then `mu_s` moves it to another eigenspace `V_{s^{+-1}(2i-m')}` and Theorem 3.1 evaluates the best
conjugate), integral with `ln |c|_1 <= 142 q` at `q = 3251`, and it must be **neither** product-type (4.1) **nor**
supported on a union of coordinate sub-cubes (4.2): it must use sign cancellations inside the patterns, at degree
between `8jq/(q+1)` and `m' <= 8j + 8`. For `m'` even the natural candidate is `F_q^*`-invariant in `V_0`, paired with
`M_{m'/2} = |hat k(1)|^{m'}`, where `<c, M_{m'/2}>` is a rational integer (Remark 3.1(d)); for `m'` odd no
`F_q^*`-invariant `tau`-isotypic witness exists (the stabiliser of `V_t`, `t != 0`, in `F_q^*` is trivial), and the
witness's `F_q^*`-orbit has `q - 1` members.

## 5. Passage to kappa_H

**5.1 On the circle the passage `kappa_V -> kappa_H` is free (proved).** Reduction §4 needs (P1) or (P2) to pass from
`kappa_V` (equivalently `Psi_A`) to `kappa_H`, and both are open in the regime because Bombieri-Vaaler controls only
codimensions (Remark 4.3). On `C_1` neither is needed: Theorem 3.1 produces its bound from the `AGL(1,q)`-images of
**one integral witness**, all of the same height, so the bound is a `kappa_H` bound for every `H >= |c|_inf`. The
record's height `ln H' = m' W/40 = (m'/20) q tau` exceeds the largest witness height that Corollary 3.2 can use,
`ln|c|_1 < q tau/2`, by the factor `m'/10 >= 205` in the exponent. So every witness that is good enough for the `kappa_V`
bound is automatically inside `V_{H'}`. The `H`-gap of the record (`ln H` between `2 ln U` and `(E/2) ln E`) is irrelevant on
the circles: what matters is the exact question of Theorem 3.4 at height `H_*`.

**5.2 What Theorem 3.4 does and does not say about `kappa_H` versus `kappa_V`.**
* It does **not** give `kappa_H >= e^{-o(n ln n)} kappa_V` on `C_1`. A circle could lie in `TC_{H_*}` and not in `TC_V`
  (all witnesses of height `> H_*`), and then Theorem 3.4 gives only `kappa_{H'} >= 2^{-m'}/(|A| H' q^{m'}) = e^{-(m'/20 + o(1)) q tau}`,
  far below `eps_0 = e^{-(1/2 + o(1)) q tau}`.
* It **does** say that a cap through the circles needs either an exact identity (`C_1 subset TC_{H'}`, e.g. `C_1 subset TC_V`,
  the cap agent's target) or the window case: all witnesses of height `<= H_*` vanish on `C_1`. In the window case
  fakes on `C_1` would have `kappa_{H'}` positive but below `eps_0`; Theorem 3.4 does not decide this case.

**5.3 Transfer to other symmetric boxes (remark).** The proof of Theorem 3.1 uses two facts: (1) a cyclic box
automorphism of odd order `N > m` acting on the circle, so that the exponents `2i - m` are distinct mod `N`; (2) every
Galois conjugate `sigma_s`, `s in (Z/N)^*`, is realised by a box automorphism. For `sdc_q` both hold (Lemma 2.1). For the
actual-type Sylvester boxes with a Singer cycle (`N = d = 2^k - 1`, odd), (1) holds for `m < d`, but the box
automorphisms normalising the Singer cycle realise only the Frobenius powers `s in <2>`, a subgroup of index `phi(d)/k`
[standard: the normaliser of a Singer cycle in `GL(k,2)` is `GammaL(1,2^k)`; quoted from memory, M]. Then
`|N(alpha)| >= 1` gives only `max_{s in <2>} |sigma_s alpha| >= B^{-(phi(d)/k - 1)}`, with `B = |c|_1 max_A |M_i| <= |c|_1 d^m`
bounding the unrealised conjugates: a Liouville loss with exponent `~ phi(d)/k`. `pal_q` fails (1) (no translation
symmetry, barrier §3.2). **So the Liouville-free bound is special to `sdc_q`-type data, whose symmetry contains the
whole multiplier group; it does not transplant to Sylvester as it stands.**

## 6. Calibration computations (evidence / exact)

All toys are far outside the regime (no exact model is a regime box; the first regime `q` is 3251). They test the
statements of §2-§3, not the open witness question.

**6.1 Symmetry (exact; `sym_check.py`).** `AGL(1,q)` preserves the pattern set for all `q = 3 mod 4`, `7 <= q <= 83`; the
identities `l_s(tau k) = zeta^s l_s(k)`, `l_1(mu_s k) = l_s(k)` hold exactly in `Z[x]/(x^q - 1)` on random integer vectors
(also at `q = 3251`), and `2i - m` are distinct mod `q` for every `m < q`.

**6.2 The record's sdc7 datum (exact up to the LP value; `sdc7_check.py`).** `q = 7`, `L = 2`, `|A| = 757`, `m = 8`: `dim V = 1`
(fmpz nullspace), the generator `c` has `|c|_1 = 1008`, `|c|_inf = 42`, is fixed by all 42 elements of `AGL(1,7)`, and pairs
only with `M_4`: `<c, M_4> = -28224`, a **rational integer**, as Remark 3.1(d) predicts. `Phi_c(u_theta) = 2^{-8} binom(8,4)(-28224)
= -7717.5` (the referee's value) on the whole circle, and `|Phi_c|/|c|_1 = 7.65625` **equals** the LP Chebyshev value
`Cheb_A(u_theta)` (`theta = 0, 0.7`; HiGHS, floating point): with `dim V = 1` the equivariant certificate is optimal.
(Here `m = 8 > q`; the translation-isotypic argument still applies because only `i = 4` satisfies `2i = m mod 7` in
`0 <= i <= 8`.) At `m = 7` (`dim V = 44`) the LLL basis vectors are sign-characters of `AGL(1,7)` and pair with all
`M_i`.

**6.3 Theorem 3.1 against the true Chebyshev value (witnesses exact in `Z[zeta]`; LP values evidence; `calib.py`).**
Witness = shortest (in `l^1`, among the 60 shortest LLL-reduced basis vectors of the exact integer kernel) integral
`c in V(A,m)` with some `<c, M_i> != 0`; "bound" = `2^{-m} binom(m,i) max_s |sigma_s <c,M_i>| / |c|_1` (the proof of Theorem
3.1 before using `max >= 1`); "LP" = `min` over `theta in {0, 0.3, pi/q, 1.1}` of `Cheb_A(u_theta)`.

| `(q, L)` | `|A|` | `m` with `dim V > 0` | circle in `TC_V`? | `|c|_1` of witness | LP / bound |
| --- | --- | --- | --- | --- | --- |
| (7, 1) | 57 | 1, 2, 3 (`V = 0` for `m >= 4`) | no (exact witnesses) | 2, 4, 768 | 3.90, 2.58, 2.51 |
| (7, 2) | 757 | 1..6 | no | 2, 6, 1024, 16, 32, 128 | 3.92, 5.69, 6.71, 5.58, 4.87, 4.88 |
| (7, 3) | 1765 | 4, 5, 6 (only these computed) | no | 12, 24, 36 | 19.9, 33.0, 20.1 |
| (11, 1) | 133 | 1, 2, 3 (`V = 0` for `m >= 4`) | no | 2, 4, 6 | 3.96, 2.63, 2.60 |
| (19, 1) | 381 | 1, 2, 3 (`V = 0` for `m = 4..6`) | no | 2, 4, 6 | 3.99, 2.65, 2.65 |

*Reading (evidence).* (i) The equivariant bound is within a factor `<= 7` of the true value for `L <= 2`, and within `<= 33` for
`(7, 3)`, where the LLL-selected witness is short but not LP-optimal. So the Liouville-free mechanism loses little
beyond the choice of witness. (ii) In every toy case with `V(A,m) != 0` the circle
is **not** in `TC_V` (rigorous for these toys: an exact witness was found). (iii) The toys never reach the regime's
defining feature: they either have `V = 0` (small `L`: the box is exhausted at order 3, a toy "cap" by emptiness that has
no regime analogue, since `|A| >= 3^q/2` there) or degree `m <= 8Lq/(q+1)` where product-type witnesses exist
(Prop 4.1); for every toy above the witnesses found have `|c|_inf <= 2` except two LLL artefacts (`|c|_inf = 128`). A toy
with `8Lq/(q+1) < m < q`, `dim V(A,m) > 0` and thin cuts needs roughly `q >= 64 c^2` with `L ~ c sqrt q`, i.e. `q` in the
hundreds: out of exact reach.

## 7. What is proved, what is open

| statement | status | where |
| --- | --- | --- |
| `AGL(1,q)` acts on `sdc_q` boxes; multipliers realise the Galois group of `Q(zeta_q)` on `M_i` | **proved** | Lemma 2.1 |
| Uniform circle bound `Cheb_A(u_theta) >= 2^{-m} binom(m,i) Cheb^C_A(M_i)`, all `theta`, `m < q` (no `(m+1)^{-1/2}` loss) | **proved** | (2.2) |
| Liouville-free witness bound: integral `c` with `<c,M_i> != 0` gives `kappa_H(hat u_theta) >= 2^{-m} binom(m,i)/(|c|_1 q^m)` for all `theta`, all `H >= |c|_inf` | **proved** | Thm 3.1 |
| Gap principle: on `C_1`, `kappa_H` is `0` or `>= 2^{-m}/(|A| H q^m)` | **proved** | Thm 3.4 |
| CIRCLE(q) survives at `tau` if `ln lambda_q(j,m') < q tau/2 - m' ln 2q - ln(Cq/(m'+1))` (budget `142.5 q` at `q = 3251`, `tau = tau_D(0.184)`) | **proved** implication | Cor 3.2, `regime_budget_out.txt` |
| Passage `kappa_V -> kappa_{H'}` on the circles | **not needed** (free) | §5.1 |
| No zonotope / product-type witness at `m'` (order `<= 8jq/(q+1) < m'`) | **proved** | Prop 4.1 |
| No witness on a union of coordinate sub-cubes inside `A` (order `<= 4j`) | **proved** | Prop 4.2 |
| Counting (Bombieri-Vaaler, Siegel variants) cannot give the witness height | **proved** numbers (`ln r_0 >= 1073` at `q = 3251`) | Remark 4.3 |
| Equivariant certificate optimal and rational on the record's sdc7 datum | **exact** (LP value floating) | §6.2 |
| Toys: circle never in `TC_V` when `V != 0`; bound within factor 7-33 of truth | **exact** witnesses / LP evidence | §6.3 |
| **Witness existence with height:** `C_1 not subset TC_{H_*}(A_q(j), m')`, `ln H_* ~ (tau_D/2) q`, at one regime `q` | **open** | — |
| `C_1 subset TC_{H'}(A_q(j), m')` (cap) | **open** (cap agent's side) | — |

**Main result (proved).** The quantitative CIRCLE(q) test is equivalent, up to a height window, to two *exact*
membership statements (Theorem 3.4). In particular the referee's obstacle "exact non-membership gives only
`e^{-Theta(q^2)}` by Liouville" is removed: on the `sdc_q` circles the Galois conjugates are box automorphisms, so exact
non-membership witnessed at height `e^{(tau_D/2 - o(1)) q}` gives the full survival bound, uniformly in `theta`, for
`kappa_{H'}` directly.

**What is not achieved.** No witness of degree `m'` was constructed, and no lower bound on `inf_theta kappa_V(u_theta)` at
the regime degree was proved. A quantitative "no" therefore remains open, and so does the cap. The open statement is
now sharp and exact:

> **(W_q)** There is an integral `c in V(A_q(j), m')` with `<c, M_i> != 0` for some `i` and
> `ln |c|_1 < q tau_D/2 - m' ln(2q) - ln(Cq/(m'+1))` (`~ 142 q` at `q = 3251`).

(W_q) at one regime `q` proves that the Fourier circles of `sdc_q` do not cap route 5 at any barrier scale (Cor 3.3).
By Props 4.1 and 4.2 its witness must be non-product and sign-correlated, of degree between `8jq/(q+1)` and `8j + 8`.
This is the "first non-zonotopal certificate at an irrational direction of an overdetermined regime box" of
README-16 §7.2, now with an explicit and generous height budget, and with the `H`-passage and the Diophantine and
Liouville losses all removed.

**Consequence for route 5 (per the round's rule).** CIRCLE(q) is **not decided** by this attempt. Survival is reduced to
(W_q), which no available method reaches: product witnesses stop at degree `8jq/(q+1)` (and, by the exact cost data,
in practice at about `4j`), and counting stops `e^{1073}` orders of magnitude short. Unless the cap attempt decides
the other side, route 5 stops here. The growth rate is unchanged: `M <= (2/3 + eps) log R/loglog R`.


## 8. Files (`survive_checks/`)

| script | output | content | label |
| --- | --- | --- | --- |
| `slib.py` | — | `sdc_q` patterns and box enumeration | library |
| `sym_check.py` | `sym_check_out.txt` | `AGL(1,q)` invariance (`q <= 83`), Galois identities, exponent distinctness | exact |
| `regime_budget.py` | `regime_budget_out.txt` | `tau_D`, witness-height budget of Cor 3.2 | 50-digit |
| `r0_bounds.py` | `r0_bounds_out.txt` | `ln r_0 in [1073.3, 3229.4]` at `q = 3251` | 50-digit |
| `zono_cost.py` | `zono_cost_out.txt` | Paley-Parseval identity; `min h'(v)` for small `v` at `q = 3251` | exact / sampled |
| `biased_l2.py` | `biased_l2_out.txt` | failure of the biased-measure `L^2` chain | 40-digit |
| `sdc7_check.py` | `sdc7_check_out.txt` | `sdc7`, `L = 2`, `m = 7, 8`: certificate, `AGL` character, `<c, M_i>`, LP | exact + LP evidence |
| `calib.py` | `calib_q*_L*_out.txt` | toy witnesses in `Z[zeta]` vs LP Chebyshev values | exact witnesses, LP evidence, `r_0` mod `p` |
