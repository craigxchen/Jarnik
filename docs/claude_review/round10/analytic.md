# Analytic mechanisms (Hecke characters, Gaussian primes in sectors) for short-arc clusters

Round 10, fresh-mechanism panel, target "analytic number theory". Scripts: `analytic_checks/` (each runs with
`python3 <script>` from that directory; output in `*_out.txt`).

Labels: **proved** (all quantifiers, here); **quoted** (a published theorem, cited from memory with its source;
the network was not used); **conjecture**; **exact** (integer arithmetic); **evidence** (a computation no proof
uses); **heuristic**.

## 0. Verdict

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` is still the best bound. The
round's question was whether the distribution of Gaussian primes in narrow sectors, Hecke characters and their
L-functions (large sieve, zero density, GRH, Kubilius/Ricci/Harman-Lewis sector theorems, sieve upper bounds in
sectors and short intervals), applied to the actual points of a cluster, can beat the fake-circle barrier. The
answer is a dichotomy, proved except where marked.

1. **Proposition A (proved).** Hecke characters restricted to the norm fibre of `N` (sector counts, `S_k(N)`,
   divisor sums) are functions of the angle vector, and therefore are the GAP Fourier analysis already capped by
   `outside.md` A.3 and `lw.md` 6.1.
2. **Theorem B, the graft barrier (proved).** Any *robust* ensemble input, meaning one whose error term is
   `>= c X^gamma` at height `X` (or a large sieve with constant `>= c(X + K X^(1/2))`), is satisfied up to a factor
   `1 + o(1)` by the genuine angle function with the fake's `O(M)` angles grafted in. This covers every published
   Hecke-type statement listed in §3. The grafted fakes have `W <= (2A + o(1)) M log M` for any `A > max(2, 1/gamma)`.
   It also holds together with residue data used through characters (P3-char). So these inputs, in any
   combination with the profile, every monomial gap and Haar-generic properties, cannot prove
   `M = o(log R/loglog R)`. The reason is selection: a cluster uses `O(M)` of the `~ M^A/log M` primes at its height,
   and robust statements cannot see `O(M)` points.
3. **Proposition C (proved).** Counting prime tuples with a Fourier majorant and absolute-value bounds on Hecke
   sums gives a bound invariant under column-wise rotation of the angles, hence `>= 1` always. It can never
   certify that a configuration does not exist, which is what an adversarial `N` needs. With amplification, the
   square-root heuristic still loses by `e^(Theta(M log M))`. Averaging statements of Kubilius type give exceptional
   sets, never empty ones (§4.2).
4. **Theorem D, collapse (proved).** The pointwise separation of a cluster prime from external lattice points,
   `|sin(theta_j - arg alpha)| >= (p_j Norm alpha)^(-1/2)` for `Norm alpha <= p_j + 1` (true for actual primes, as
   `Im(pi_j conj alpha)` is a nonzero integer), is **false for every fake that is not the shadow of an actual
   cluster of no larger size**: by Minkowski it forces `phi_j = arg beta_j`, `Norm beta_j <= p_j`, and then
   `prod beta_j^(a) conj(beta_j)^(e - a)` is an actual cluster. So the class it defines has exactly the actual
   minimum: fake-killing, but no reduction. Lemma 0.1 of `outside.md` is the case `Norm alpha <= 2`.
5. **The analytic version (conditional).** With external Gaussian *primes* in place of lattice points, the same
   collapse follows from Conjecture GPA (§5.1): every direction that is not a lattice direction of norm `<= h` is
   approximated by Gaussian prime directions with `|sin| < (h Norm varpi)^(-1/2)`. GPA is beyond GRH (GRH reaches
   sectors of angle `X^(-1/2)(log X)^(3+o(1))`; Harman-Lewis reach `X^(-0.381)`); at lattice directions it is of
   Bunyakovsky type. It is supported by the computations of §7. Under PCG the genuine angle is an isolated
   admissible point (Lemma D.2, proved; hypothesis checked for all `p <= 10^5`).

**Status: barrier.** Narrow-sector distribution of Gaussian primes either cannot see the cluster primes
(statistical use, Theorem B and Proposition C) or is the lattice itself (pointwise use, Theorem D and
Proposition D.1). Neither half gives growth. The input that is true for actual circles and false for every
fake (item 4) exists and is elementary, but it returns exactly the original problem. Open: the regime
`P_+(N) < M^(2 - delta)`, where no fakes are known and Theorem B does not apply (§6).

## 1. Setting and the two kinds of analytic input

**Clusters.** As in `round4/outside.md` §1 and `round8/stability.md` §1: after normalisation,

```text
z_x = eps_x prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj),   N' = prod_j p_j^(e_j),   W = log N (>= log N'),
arg z_x = <2a_x - e, theta> + (pi/2) k_x   (mod 2 pi),     theta_j = arg pi_j,
```

and the cluster condition is that all `arg z_x` lie in one interval `theta_0 + [0, Delta]`, `Delta = C N^(-1/4)`.
`r` is the number of distinct primes of `N'`. Under the near-extremal hypothesis `H(eps)` of `stability.md`,
`r = Theta(M)` and every `p_j <= M^(4/eps - 3 + o(1))` (Cor. 3.6).

**Angle models (fakes).** `outside.md` §1: a profile `(p_j, e, a_x, delta_x, k_x)` and real numbers `phi_j`
(the *fake angles*) with `<2a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x`, `delta_x in [0, Delta]` distinct,
satisfying the gap inequalities of `outside.md` Lemma 0.1 for every `v in Z^r \ {0}`. Its "points" are
`N^(1/2) exp(i(theta_0 + delta_x))`; no Gaussian integer stands behind them. Every actual cluster is an angle
model with `phi = theta`. The refereed fakes used below are:
* `round8/lw.md` Theorem 6.1 (class `A_0`: profile, every monomial gap, Haar-generic angles): skew double core,
  `M = q + 1`, `r = 2q`, one prime per column, all in `[P', P' e^(1/(8q))]` with `P <= P' <= 2P`, for **every**
  `P >= P_1 = max((40M/C_0)^4, (64M)^(2/eps))`, and `W <= 2(M - 1) log(2P) + 1/4`. With `P = M^A`
  (`A >= 4`, `M` large) this is `W <= (2A + o(1)) M log M`.
* `round5/round6_sharp/sharp.md` Theorem 5.1 (class P3-char, residue data at moduli `<= X = 3M^A`): primes in
  `[Z, 2Z]`, `Z >= 3 M^A (log M)^2`, `W <= (b max(A,2) + o(1)) M log M`, `b >= 33`.

Both families allow the prime size to be any fixed power of `M` at the cost of the constant only. This
freedom is what Theorem B exploits.

**Gaussian primes and their angle function.** Let `G` be the set of Gaussian primes `varpi` with
`0 < arg varpi < pi/2` and `Norm varpi = p = 1 mod 4`. Each split `p` has exactly two elements of `G` above it,
`varpi_p` with `arg varpi_p = theta_p in (0, pi/4)` and `varpi'_p = i conj(varpi_p)` with argument `pi/2 - theta_p`.
The Hecke Grossencharacters of `Q(i)` are `lambda_k(alpha) = (alpha/|alpha|)^(4k)`, `k in Z`; on `G`,
`lambda_k(varpi) = e^(4 i k arg varpi)` and `lambda_k(varpi'_p) = conj lambda_k(varpi_p)`. Inert and ramified
primes have `lambda_k = 1` and play no role below. The *genuine angle function* is
`theta : G -> R/2piZ`, `varpi -> arg varpi`.

**Two kinds of analytic input.** An analytic statement used in an argument about a cluster is of one of two
kinds.
1. *Fibre inputs*: statements about the norm fibre `{z : |z|^2 = N}`, the Gaussian divisors of `N`, and Hecke
   characters evaluated on them (sector counts in the fibre, the sums `S_k(N) = sum_(|z|^2 = N) lambda_k(z)`,
   Riesz products, ...). They involve only the primes of `N`.
2. *Ensemble inputs*: statements about the angle function on **all** Gaussian primes up to some height:
   large sieve inequalities for Hecke characters, bounds for `sum_varpi f(Norm varpi) lambda_k(varpi)` (from
   zero-free regions, zero-density estimates or GRH via the explicit formula), counts of Gaussian primes in
   polar boxes `{x < Norm varpi <= x + y, arg varpi in I}` (Kubilius, Ricci, Harman-Lewis; sieve upper bounds),
   joint sector and residue-class counts. Applied to a cluster, they constrain the angles `theta_j` because
   `theta_j = theta(pi_j)` is one value of the genuine angle function.

Fakes have no Gaussian integers, so for them only the fibre inputs make sense directly. To ask whether an
ensemble input separates actual clusters from fakes one must say what a fake's ensemble is. Section 3 uses
the most favourable choice for the fake: keep every genuine angle, except at the cluster primes (the *graft*).
Section 5 treats the inputs that are not of this robust kind.

## 2. Fibre inputs are functions of the angle vector (Proposition A)

**Proposition A (proved; a re-reading of `outside.md` direction (a)).** Let `N'` be odd with only split prime
factors, `N' = prod_j p_j^(e_j)`. For every `k in Z`,

```text
S_k(N') := sum_(|z|^2 = N') lambda_k(z) = 4 prod_j sum_(a=0)^(e_j) e^(4ik(2a - e_j) theta_j)
         = 4 prod_j sin(4k(e_j + 1) theta_j) / sin(4k theta_j),                                   (2.1)
```

and more generally, for every function `f` on the circle, `sum_(|z|^2 = N') f(arg z)` depends only on `e` and on
`theta = (theta_j)`, through the angle GAP `{<2a - e, theta> + u pi/2}`. In particular every sector count in the
fibre, every Hecke sum over the fibre or over the Gaussian divisors of `N'`, and every Riesz product or Bohr-set
statement about them, is a function `F(e, theta)` that is defined for **all** real `theta`.

*Proof.* Unique factorisation in `Z[i]`: the elements of norm `N'` are `i^u prod pi_j^(a_j) conj(pi_j)^(e_j - a_j)`,
`0 <= a_j <= e_j`, `u in Z/4`, each once. `lambda_k(i) = i^(4k) = 1` gives the factor 4, and the inner sum is a
Dirichlet kernel in `8k theta_j`. QED

*Consequence.* Evaluate `F` at a fake's angles `phi`. Any inequality `F(e, theta) <= B` that is derived for
actual clusters from the fibre alone (that is, from identities and inequalities valid for every real angle
vector) holds for every angle model. These are exactly the GAP/Fourier/Bohr/Riesz-product/inverse
Littlewood-Offord arguments of `outside.md` direction (a), which Theorem A.3 and Corollary A.4 there cap, and the
angle models of `lw.md` Theorem 6.1 and `exact.md` cap them at `W = O(M log M)`. So **Hecke characters
restricted to the norm fibre are the GAP Fourier analysis, and add nothing**. The input that separates `theta`
from `phi` must involve primes outside `N`: the ensemble. (The trivial identity `theta_j = theta(pi_j)` is such an
input; it is all of the problem.)

*Remark A.1 (why fibre Fourier sees nothing).* A Fourier proof of `#{x : arg z_x in theta_0 + [0, Delta]} <= M_0`
needs `S_k(N')` for `|k| <~ 1/Delta = e^(W/4)/C`. By Parseval, `(1/(2K+1)) sum_(|k| <= K) |S_k(N')|^2 -> 4 r(N')` as `K -> infinity` for
`theta` generic (the factor 4 counts associates; `r(N') = 4 prod (e_j + 1)`), so `|S_k|` is typically `~ r(N')^(1/2)`, and the
absolute-value sector bound is `~ r(N')^(1/2) = 2^(Theta(M))` (heuristic, consistent with the rigorous cap above).

## 3. The graft barrier for ensemble inputs (Theorem B)

**Definition 3.1 (graft).** Let `S` be a finite set of primes `p = 1 mod 4` and `phi = (phi_p)_(p in S)` real numbers.
The *graft* `theta^S = theta^(S, phi)` is the angle function on `G` with
`theta^S(varpi_p) = phi_p`, `theta^S(varpi'_p) = pi/2 - phi_p` for `p in S`, and `theta^S = theta` elsewhere.
Write `lambda^S_k(varpi) = e^(4ik theta^S(varpi))`, `s = |S|`, `P_0 = min S`, `P_1 = max S`, and for a norm
range `J`, `s_J = |S cap J|`. For an angle model with primes `S` and angles `phi`, its graft is the genuine
ensemble with the fake's angles put in at the fake's primes. It agrees with the fake on the fibre and with the
genuine primes everywhere else.

**Lemma 3.2 (proved; the graft moves at most `2s` points).**
1. *(Large sieve.)* For `X >= 2`, `K >= 1` and complex `a` supported on `{varpi in G : Norm varpi <= X}`,

   ```text
   ( sum_(|k|<=K) |sum a_varpi lambda^S_k(varpi)|^2 )^(1/2)  <=  ( sum_(|k|<=K) |sum a_varpi lambda_k(varpi)|^2 )^(1/2)
                                                              +  ( 8 s_[1,X] (2K+1) )^(1/2) ||a||_2 .
   ```
2. *(Prime sums.)* For bounded `f` supported in `[X_0, X_1]` and every `k`,
   `|sum_varpi f(Norm varpi) log(Norm varpi) (lambda^S_k(varpi) - lambda_k(varpi))| <= 4 s_[X_0,X_1] ||f||_inf log X_1`.
   The same bound holds with prime powers included (`lambda_k(varpi^m) = lambda_k(varpi)^m`), with
   `log X_1` counting `sum_(m : p^m <= X_1) log p`.
3. *(Boxes.)* For every polar box `B = {varpi in G : Norm varpi in J, theta(varpi) in I}` (any interval `J`,
   any measurable `I`), optionally intersected with a congruence condition on `varpi`,
   `|#{varpi in B under theta^S} - #{varpi in B under theta}| <= 2 s_J`.

*Proof.* Only the `2 s_J` elements of `G` above `S cap J` change. (1) The difference
`D_k = sum_(varpi above S) a_varpi (lambda^S_k - lambda_k)(varpi)` has `|D_k| <= 2 sum_(varpi above S, Norm <= X) |a_varpi|
<= 2 (2 s_[1,X])^(1/2) ||a||_2` by Cauchy-Schwarz; sum `|D_k|^2` over `2K + 1` values of `k` and apply Minkowski's
inequality in `l^2`. (2) Each changed term has modulus `<= 2 |f| log Norm varpi`, two elements per prime.
(3) Each changed element leaves or enters `B` at most once; congruence classes of `varpi` are not changed by a
graft. QED

**Definition 3.3 (robust ensemble inputs).** Fix `c > 0` and `gamma in (0, 1]`. A *`(c, gamma)`-robust ensemble
input* is a family, indexed by any set, of inequalities about an angle function `vartheta`, each of one of the
forms
* (R1) `sum_(|k|<=K) |sum_(Norm varpi <= X) a_varpi e^(4ik vartheta(varpi))|^2 <= Lambda(X, K) ||a||^2` with
  `Lambda(X, K) >= c (X + K X^(1/2))`;
* (R2) `|sum_varpi f(Norm varpi) log(Norm varpi) e^(4ik vartheta(varpi)) - Main(f, k)| <= E(f, k)` with
  `E(f, k) >= c ||f||_inf X_1(f)^gamma`, `X_1(f) = sup supp f` (needed only when `X_1(f)` exceeds the height of
  the grafted primes; below it a graft changes nothing);
* (R3) `L(B) <= #{varpi in B under vartheta} <= U(B)` for polar boxes (optionally with congruence conditions), where
  for each box either `L(B) = 0` or `L(B) >= c x_B^gamma`, and either `U(B) = infinity` or `U(B) >= c x_B^gamma`,
  `x_B = sup` of the norm range of `B`.

*What this covers (quoted, from memory; only the shape of the error terms matters).* The large sieve for Hecke
Grossencharacters in the form `<< (X + K^2) ||a||^2` (which is `>= (1/2)(X + K X^(1/2))`), or with the lattice
count `X/K + X^(1/2)` per arc; Hecke prime number theorems from the classical zero-free region
(`E = X_1 exp(-c' (log X_1)^(1/2))`), from zero-density estimates, and under GRH for Hecke L-functions of `Q(i)`
(`E << X_1^(1/2) log^2(X_1(|k| + 2))`, the main terms `Main` being anything); Kubilius' and Ricci's theorems on
Gaussian primes in sectors, Harman-Lewis (sectors of angle `X^(-0.381)`; lower bounds `L(B) >> |I| X/log X`),
Brun-Titchmarsh and Selberg-sieve upper bounds in sectors, short norm intervals and polar boxes, and their
versions in residue classes. All have error terms or bounds at least a fixed power of the height. In every case
`c` and `gamma` are absolute.

**Theorem B (graft barrier; proved).** Fix `c > 0`, `gamma in (0, 1]` and let `E` be any `(c, gamma)`-robust ensemble
input that holds for the genuine angle function `theta`. Fix `A > max(2, 1/gamma)`, `A >= 4`, and `C_0 in (0, 1]`.
For every Paley order `M = q + 1` (`q = 3 mod 4` prime) with `M >= M_3(A, C_0, c, gamma)` there is an angle model
`F_M` (the `lw.md` Theorem 6.1 model with `eps = 1`, `P = M^A`) such that:
1. `F_M` has `M` distinct points on an arc of length `<= C_0 sqrt R`, satisfies every monomial gap of
   `outside.md` Lemma 0.1, and `W <= (2A + o(1)) M log M`;
2. the graft `theta^S` of `F_M` satisfies every inequality of `E` with `Lambda`, `E` replaced by
   `(1 + eta_M) Lambda`, `(1 + eta_M) E`, and `L`, `U` by `(1 - eta_M) L`, `(1 + eta_M) U`, where
   `eta_M = O_c(M^((2 - A)/4) + A M^(1 - A gamma) log M) -> 0`.

The same holds with `F_M` replaced by the P3-char members of `sharp.md` Theorem 5.1 at residue exponent `A`
(primes `>= 3 M^A`, `r = b(M-1)` columns, `W <= (b A + o(1)) M log M`): residue data at all moduli `<= 3M^A` used through characters,
together with any robust ensemble input, still allow `W = O(M log M)`.

*Proof.* `F_M` exists with item 1 by `lw.md` Theorem 6.1 (`P = M^A >= P_1` for large `M`); its primes `S`
satisfy `s = 2q < 2M` and `P_0 >= M^A`. Apply Lemma 3.2.
* (R1) If `X < P_0` nothing changes. If `X >= P_0`: `8s(2K+1) <= 24 s K <= (24 s/(c X^(1/2))) Lambda(X, K)`, so
  the graft satisfies (R1) with `(1 + t^(1/2))^2 Lambda`, `t = 24 s/(c P_0^(1/2)) <= 48 M^(1 - A/2)/c`. For `t <= 1`,
  `(1 + t^(1/2))^2 - 1 = 2 t^(1/2) + t <= 3 t^(1/2) = O_c(M^((2 - A)/4))`, which tends to `0` because `A > 2`.
* (R2) The change is `0` unless `X_1 >= P_0`; then it is `<= 4 s ||f|| log X_1 <= (4 s log X_1/(c X_1^gamma)) E`,
  and `log x/x^gamma` decreases for `x >= e^(1/gamma)`, so the factor is `<= 4 s A log M/(c M^(A gamma))
  <= 8 A M^(1 - A gamma) log M / c -> 0` because `A gamma > 1`.
* (R3) Same: the change is `<= 2 s < 4 M` and the bounds are `>= c x_B^gamma >= c M^(A gamma)` whenever the norm
  range of `B` meets `S`.

For P3-char the primes lie in `[Z, 2Z]` with `Z >= 3 M^A (log M)^2` and `r = b(M - 1)`, and the same three
estimates apply with `s = r`. QED

**Corollary B.1 (what Theorem B rules out).** Let an argument use, about an actual near-extremal cluster:
the exponent profile; the monomial gaps (or, for the second family, residue data through characters at moduli
`<= 3 M^A`); any Haar-generic property of the realisation; and any robust ensemble input, applied to the angle
function whose values at the cluster primes are the cluster angles, with its implied constants (any fixed
relaxation `1 + eta`, `eta > 0`). Then the argument cannot prove `W/(M log M) -> infinity`, i.e. it cannot prove
`M = o(log R/loglog R)`. The fakes of Theorem B satisfy all of its inputs at `W <= (2A + o(1)) M log M`.

*Remarks 3.4.*
1. *Selection is the reason.* The graft moves `2s = O(M)` points out of the `~ P/log P` Gaussian primes of norm
   `<= P`. Every statement whose error term exceeds `O(M)` at heights `>= M^A` cannot see the move. Square-root
   error terms (the GRH quality, and the truth for typical boxes, §7.1) are `>= M^(A/2)`.
2. *Not covered, and why.* (R1) with the sharp constant `Lambda = c(X + K)` for `K >> X^(1/2)` is **not** robust:
   it is the dual large sieve for points with spacing `>> 1/X`, and the spacing of Gaussian prime angles is the
   lattice fact `|sin(arg varpi - arg varpi')| >= (Norm varpi Norm varpi')^(-1/2)` (as `Im(varpi conj varpi')` is a nonzero
   integer). Likewise "this particular sector, of predicted count `>> 1`, is empty" is not a robust statement.
   These *pointwise* statements are the subject of §5, where they turn out to collapse the fakes back to actual
   clusters.
3. *Residues.* A graft changes only angles. If the argument also uses residues of the cluster primes, the
   fake's designed residues can be grafted too; Lemma 3.2(3) is unchanged. So Theorem B also covers joint
   sector-and-residue distribution statements, but only for those residue classes in which the fake class is
   capped (P3-char). For P3-all and P3-G it says nothing, because their caps are open.

## 4. Adversarial N: why counting and averaging cannot reach a single circle (Proposition C)

The growth statement is about every `N`. An analytic theorem about Gaussian primes can bound a single,
adversarially chosen `N` in only three ways: (a) a statement with an empty exceptional set; (b) counting the
possible configurations and showing the count is `< 1` (or `<` the number of copies one cluster produces);
(c) a pointwise statement about the specific angles. This section disposes of (a) and (b) as far as robust
inputs go; (c) is §5.

**4.1 The selection.** By `stability.md` Cor. 3.6, a near-extremal cluster (`M >= eps log R/loglog R`) uses
`r = Theta(M)` primes, all below `M^(h_0)`, `h_0 = 4/eps - 3 + o(1)`, and all but a fraction
`O(M^(-delta)/log M)` of `W'` sits on primes `>= M^(1 - delta)`. There are `~ M^(h)/(h log M)` Gaussian primes in
`G` of norm `<= M^h`. So the cluster's angles are an adversarially chosen subset of density
`<= M^(1 - h + o(1))` of the angle ensemble at the heights that carry its weight, and an adversary may move them
anywhere by changing `O(M)` values (this is the graft of §3). The one regime in which the adversary has no
choice is the sharp regime of `stability.md` Thm 5.2: at efficiency `A <= 3 + delta`, all but a `delta`-fraction of
the split primes `<= M` divide `N'`. That regime concerns the constant (`2/3`), not the growth rate, and the
fibre Fourier analysis on a complete prime set is still the capped GAP analysis of Proposition A.

**4.2 Averages.** Statements of Kubilius type (the prime factorisation and the Hecke angles of a random
Gaussian integer of norm `<= x` are close in total variation to an independent model), equidistribution of
`{arg varpi}` in sectors, and density statements all have the form
`|#{N <= x : N in E} - (model count)| <= Err(x)`. Such a statement excludes a single `N` only if
`Err(x) < 1` and the model count is `< 1`. The error terms of the known theorems of this kind grow with `x`
(quoted, from memory; I know of none that is `< 1` in counts), while
near-extremal `N` (if they exist) are already known to be rare: `#{N <= x : M_N(C) >= 2} << C x^(3/4) log x`
(`outside.md` Prop. D.1). So averaging statements are consistent with the existence of near-extremal `N` at every
density they can resolve.

**4.3 Counting with absolute values.** A counting proof would take, for each column `j`, a candidate set `Pi_j`
of Gaussian primes (for example all of `G` up to some height), a fixed integer matrix `Cmat in Z^(m x r)` (the
cluster's pair-difference forms `2(a_x - a_1)`, `m = M - 1`, or any other forms the cluster makes small), and
count

```text
T = #{ tau in Pi_1 x ... x Pi_r : Cmat theta(tau) in B (mod 2 pi) },     B = prod_x [-Delta', Delta'],
```

where `theta(tau) = (theta(tau_1), ..., theta(tau_r))`. With a majorant `g >= 1_B` on the torus `T^m` with
absolutely summable Fourier coefficients,

```text
T <= sum_tau g(Cmat theta(tau)) = sum_(k in Z^m) ghat(k) prod_j S_j((Cmat^T k)_j),    S_j(mu) = sum_(varpi in Pi_j) e^(i mu theta(varpi)),
```

and the analytic input (large sieve, GRH, zero density) bounds the `|S_j(mu)|`, `mu != 0`. The resulting bound is
`Phi(g) = sum_k |ghat(k)| prod_j |S_j((Cmat^T k)_j)|`.

**Proposition C (proved).** For all finite `Pi_1, ..., Pi_r`, every `Cmat` and every majorant `g >= 1_B`
(`g >= 0`, `sum |ghat| < infinity`, `0 in B`):
1. `Phi(g)` is unchanged if the angles of each `Pi_j` are rotated by an arbitrary column-dependent constant `s_j`
   (the columns are treated as separate copies of their sets, so `Pi_j = Pi_l` is allowed);
2. `Phi(g) >= sup_(s in R^r) sum_tau g(Cmat (theta(tau) + s)) >= g(0) >= 1`.

So this method never proves `T = 0`, whatever the sets, the forms and the quality of the bounds for `S_j`.

*Proof.* (1) Rotating `Pi_j` by `s_j` multiplies `S_j(mu)` by `e^(i mu s_j)`. (2) For each `s`,
`sum_tau g(Cmat(theta(tau) + s)) = sum_k ghat(k) prod_j e^(i (Cmat^T k)_j s_j) S_j((Cmat^T k)_j) <= Phi(g)`. Take any
`tau^0` and `s = -theta(tau^0)`: the term `tau = tau^0` equals `g(0) >= 1_B(0) = 1`, and all terms are `>= 0`. QED

*Meaning.* A column-wise rotation is exactly what turns genuine angles into fake ones: after rotating column `j`
by `s_j = phi_j - theta(tau^0_j)`, the tuple `tau^0` carries the angles `phi`. Proposition C is the fake barrier in
Fourier form: the absolute-value method cannot distinguish the genuine ensemble from one in which each column
has been rotated, and the rotated ensemble contains a solution. Only the **phases** of the `S_j(mu)`, at the
frequencies `mu = (Cmat^T k)_j` with `|k| <~ 1/Delta'`, can see the difference. For a cluster those frequencies are
up to `e^(W/4)`, astronomically beyond any height at which phases of Hecke sums are controlled.

*Amplification (heuristic).* If one cluster produced `Z` counted tuples (column automorphisms of the profile,
conjugations, sub-clusters), the method would need `Phi(g) < Z`, not `< 1`. Proposition C does not exclude this.
But under square-root behaviour `|S_j(mu)| ~ n_j^(1/2)` for `mu != 0` (the GRH size and, by Fejer positivity,
`sum_(|mu| <= K) (1 - |mu|/(K+1)) |S_j(mu)|^2 >= (K + 1) n_j` unconditionally, the best possible on average),
`Phi(g) >= sum_k |ghat(k)| prod_j |S_j| ~ (prod_j n_j)^(1/2) >= (M^(A-1))^(r/2)`, while the natural amplification
factors are `e^(O(M log M))` and the main term `ghat(0) prod n_j ~ Delta^(M-1) prod n_j = e^(-Theta(M W))`. For
`A` larger than an absolute constant this fails by `e^(Theta(M log M))`; this is the simultaneous version of
`stability.md` §8.4 ("square-root cancellation reaches the Liouville scale and no further").

**4.4 Summary for adversarial `N`.** Averages give exceptional sets, never empty ones (4.2). Counting with
absolute values is blind to column rotations, hence to fakes (Prop. C). Worst-case statements of robust type are
satisfied by grafted fakes (Thm B). What remains is a statement about **where the specific cluster angles sit
relative to the specific genuine ensemble**: a pointwise statement. That is §5.

## 5. Pointwise sector inputs: cross-gaps collapse fakes to actual clusters (Theorem D)

The pointwise facts that robust inputs miss (Remark 3.4.2) are separation statements between a cluster prime
and the rest of the lattice. They are exactly what distinguishes `theta_j` from a fake `phi_j`.

**Definition 5.1 (cross-gaps).** For real `phi` and an integer `n >= 2`:
* `LCG_n(phi; Y)` (lattice cross-gap): `|sin(phi - arg alpha)| >= (n Norm alpha)^(-1/2)` for every
  `alpha in Z[i] \ {0}` with `Norm alpha <= Y` that is not parallel to `e^(i phi)` (i.e. `alpha e^(-i phi) notin R`);
* `PCG_n(phi; Y)` (prime cross-gap): `|sin(phi - arg varpi)| >= (n Norm varpi)^(-1/2)` for every Gaussian prime
  `varpi` with `Norm varpi <= Y` and `Norm varpi != n` (parallel ones included).

`Y = infinity` is allowed.

**Lemma 5.2 (actual clusters satisfy both; proved).** If `pi` is a Gaussian prime of norm `p`, then
`LCG_p(arg pi; infinity)` and `PCG_p(arg pi; infinity)` hold.

*Proof.* For `alpha in Z[i]`, `Im(pi conj alpha)` is an integer and equals `|pi||alpha| sin(arg pi - arg alpha)`. It
vanishes iff `alpha` is parallel to `pi`. Otherwise `|sin(arg pi - arg alpha)| >= 1/(|pi||alpha|) = (p Norm alpha)^(-1/2)`.
A Gaussian prime `varpi` parallel to `pi` is a real multiple of `pi`, hence `varpi = +-pi` (both are primitive), so
`Norm varpi = p`; such `varpi` are excluded from PCG. QED

*Remark.* Lemma 0.1 of `outside.md` at `v = e_j` is `LCG_(p_j)` against the four directions of norm `<= 2`:
`alpha in {1, i}` gives `|sin phi_j|, |cos phi_j| >= p_j^(-1/2)` and `alpha = 1 +- i` gives the diagonal gap
`>= (2 p_j)^(-1/2)`. The fakes satisfy these (and all monomial versions). Theorem D shows that going up to external
directions of norm `p_j + 1` kills every fake that is not the shadow of an actual cluster.

**Lemma 5.3 (Minkowski; proved).** Let `n >= 2` and suppose `LCG_n(phi; n + 1)` holds. Then
`phi = arg beta (mod 2 pi)` for a primitive `beta in Z[i]` with `Norm beta <= n`.

*Proof.* Put `Q = n + 1/2` and `K = {w in C : |Re(w e^(-i phi))| <= Q^(1/2), |Im(w e^(-i phi))| <= Q^(-1/2)}`, a compact
convex set, symmetric about `0`, of area `4`. By Minkowski's convex body theorem (compact form) it contains some
`alpha in Z[i] \ {0}`. Then `Norm alpha <= Q + 1/Q < n + 1` (as `n >= 2`), so `Norm alpha <= n`. If `alpha` were not
parallel to `e^(i phi)`, then `|sin(phi - arg alpha)| = |Im(alpha e^(-i phi))|/|alpha| <= Q^(-1/2) |alpha|^(-1) < (n Norm alpha)^(-1/2)`,
contradicting `LCG_n(phi; n + 1)`. So `alpha e^(-i phi)` is real, `phi = arg(+-alpha)`, and `beta` = the primitive part of
`+-alpha` has `Norm beta <= Norm alpha <= n`. QED

**Theorem D (collapse; proved).** Let `(P, phi)` be an angle model (`outside.md` §1) with primes `p_j`, exponents
`e_j`, `N = prod p_j^(e_j)`, and `M` distinct points on an arc of length `<= C N^(1/4)` of the circle of radius
`N^(1/2)`. Suppose `LCG_(p_j)(phi_j; p_j + 1)` holds for every `j`. Then there are `M` distinct Gaussian integers
`y_1, ..., y_M` of one norm `N'' <= N` on an arc of length `<= C N''^(1/4)` of the circle of radius `N''^(1/2)`.
Consequently, for every `M` and `C`,

```text
inf { W : angle models with M points, constant C, satisfying LCG_(p_j)(phi_j; p_j + 1) for all j }
   =  inf { log N : actual M-point clusters with constant C on x^2 + y^2 = N }.
```

*Proof.* Lemma 5.3 gives primitive `beta_j` with `Norm beta_j <= p_j` and `phi_j = arg beta_j (mod 2 pi)`. Put
`y_x = i^(-k_x) prod_j beta_j^(a_xj) conj(beta_j)^(e_j - a_xj)`. Then `|y_x|^2 = prod_j Norm(beta_j)^(e_j) =: N'' <= N`
and, by the realisation equations (1.1) of `outside.md`,
`arg y_x = <2a_x - e, phi> - (pi/2) k_x = theta_0 + delta_x (mod 2 pi)`. The `delta_x` are distinct, so the `y_x`
are distinct; they lie in an arc of angular width `Delta <= C N^(-1/4)`, of length
`<= N''^(1/2) Delta <= C N''^(1/2) N^(-1/4) <= C N''^(1/4)`. For the identity: `<=` because actual clusters are angle
models satisfying LCG (Lemma 5.2 and normalisation, which only lowers `N`); `>=` by the first part. QED

So the lattice cross-gap is an input that is **true for actual circles and false for every angle model that is
not the shadow of an actual cluster of no larger size**. It is not a barrier and not a reduction: the class it
defines has exactly the actual minimum. (Compare `round7` route G, whose collapse (Thm 3.1 of `p3G.md`) is
conditional on a pinning closure; Theorem D is unconditional, because it uses the archimedean lattice directly.)

### 5.1 The analytic version: prime cross-gaps and Gaussian primes in very narrow sectors

`LCG` is a lattice statement. The analytic, sector-distribution version replaces the external lattice points by
external Gaussian **primes**. Whether `PCG` still collapses the fakes is a question about Gaussian primes in very
narrow sectors.

**Conjecture GPA(h) (Gaussian prime approximation at exponent 1/2).** Let `h >= 2`. For every real `phi` that is
not `arg beta (mod 2 pi)` for any `beta in Z[i]` with `0 < Norm beta <= h`, there is a Gaussian prime `varpi` with
`Norm varpi != h` and `|sin(phi - arg varpi)| < (h Norm varpi)^(-1/2)`.

**Proposition D.1 (proved, conditional on GPA).** If GPA(`p_j`) holds for every `j`, then every angle model with
`PCG_(p_j)(phi_j; infinity)` for all `j` satisfies the conclusion of Theorem D. So under GPA the class
(angle models + prime cross-gaps) has exactly the actual minimum.

*Proof.* GPA(`p_j`) and `PCG_(p_j)(phi_j; infinity)` together force `phi_j = arg beta_j` with `Norm beta_j <= p_j`; the
rest is the proof of Theorem D. QED

**What GPA asks, against known results.** Split by the nature of `phi`.
1. *`phi` not a lattice direction.* GPA asks for a Gaussian prime, at some norm `X`, in the sector of half-angle
   `~ (hX)^(-1/2)` around `phi`; that sector contains `~ (X/h)^(1/2)` lattice points of norm in `[X, 2X]` and
   heuristically `~ (X/h)^(1/2)/log X` primes, so heuristically GPA holds with room to spare. Known (quoted, from
   memory): Kubilius and Ricci proved primes in sectors of angle `X^(-gamma)` for small `gamma > 0`; Harman and
   Lewis (Mathematika 2001) reach `gamma = 0.381`. GRH for the Hecke L-functions of `Q(i)` gives, by the explicit
   formula with a smoothed sector, primes in every sector of angle `X^(-1/2) (log X)^(3 + o(1))` (my derivation:
   main term `~ eta X/log X`, error `<< X^(1/2) log^2 X` summed against `sum |ghat| << 1`). GPA needs angle
   `(hX)^(-1/2)`, a factor `h^(1/2) (log X)^3` narrower, at a single adversarial `phi`. So GPA is beyond GRH, in the
   way that primes in every interval `[x, x + c x^(1/2)]` are beyond RH. Being free to choose `X` does not help
   GRH, since the ratio main/error `~ 1/(h^(1/2) log^3 X)` is `< 1` at every scale.
2. *`phi = arg beta`, `beta` primitive, `Norm beta > h`.* The lattice points near the direction `beta` lie on the
   lines `Im(alpha conj beta) = m`; `|sin(arg beta - arg alpha)| = |m|/(|alpha||beta|)`, which is
   `< (h Norm alpha)^(-1/2)` iff `|m| < (Norm beta/h)^(1/2)`. If `beta` is itself a Gaussian prime, `m = 0` works.
   Otherwise GPA asks for one Gaussian prime on one of the lines `alpha = alpha_m + t beta`, `0 < |m| < (Norm beta/h)^(1/2)`,
   i.e. for a prime value of one of the quadratics `f_m(t) = Norm(alpha_m + t beta)` (discriminant `-4m^2`). That
   is a Bunyakovsky-type statement, open for every single quadratic, though trivially checkable for each given
   `beta` (§7.3).
3. *`phi = arg beta`, `Norm beta <= h`.* Excluded, and necessarily so: by the lattice bound (the proof of Lemma 5.2
   with `beta` in place of `pi`) no prime that is not parallel to `beta` is close enough. A witness exists only when
   the direction itself carries a prime of norm `!= h` (`beta = 1` carries the inert primes, `beta = 1 + i` carries
   `1 + i`, a prime `beta` of norm `< h` carries itself); such `phi` violate PCG, which only helps the collapse.

So GPA is a sharp "primes in every very narrow sector" statement plus a Bunyakovsky-type statement; it is not
implied by GRH, it is consistent with every heuristic I know, and the computations of §7.2-7.3 find no
counterexample.

**Lemma D.2 (rigidity at the genuine angle; proved).** Let `pi` be a Gaussian prime of norm `p >= 5`, and suppose
there are Gaussian primes `varpi_+`, `varpi_-` with `Im(varpi_(+-) conj pi) = +-1` and `Re(varpi_(+-) conj pi) > 0`
(then `Norm varpi_(+-) != p`, since the primes of norm `p` give `Im` in `{0, +-p, +-(a^2 - b^2), +-2ab}`, never `+-1`). Then every `phi` with `PCG_p(phi; Y)`,
`Y >= max Norm varpi_(+-)`, and `|phi - arg pi| < min_(+-) 2 arcsin((p Norm varpi_(+-))^(-1/2))` equals `arg pi`.

*Proof.* `sin(arg varpi_(+-) - arg pi) = +-1/(|varpi_(+-)||pi|)` and `cos(arg varpi_(+-) - arg pi) > 0`, so
`arg varpi_+ = arg pi + rho_+` with `rho_+ in (0, pi/2)`, `sin rho_+ = (p Norm varpi_+)^(-1/2)`, and PCG forbids the open arc `(arg pi, arg pi + 2 rho_+)` (where
`|sin(phi - arg varpi_+)| < sin rho_+`). Symmetrically `varpi_-` forbids `(arg pi - 2 rho_-, arg pi)`. QED

So under PCG the genuine angle is an **isolated** admissible point: a fake cannot sit near it. (Existence of
`varpi_(+-)` is the Bunyakovsky-type input of item 2; §7.3 finds them for every `p <= 10^5`.)

### 5.2 What the collapse means for growth

Theorem D (unconditionally) and Proposition D.1 (under GPA) show that the pointwise sector statements kill every
fake, and that what they leave is **exactly** the set of actual clusters. A proof of growth from them is a proof
of growth for actual circles from scratch: nothing has been reduced. The same happens for truncated versions:
`LCG(·; p_j + 1)` is already exact. For PCG truncated at a polynomial height `Y = M^K` the status is open.
*Heuristic:* PCG`(·; Y)` confines each `phi_j` to within `~ (Norm beta_j Y)^(-1/2)` of a lattice direction `beta_j` of
norm `<= p_j` (and to the point `arg pi_j` itself when `beta_j = pi_j`, Lemma D.2), so a surviving fake comes from
`M` actual Gaussian integers `prod beta_j^(a) conj(beta_j)^(e - a)` within angular width `~ r Y^(-1/2)`, squeezed to
width `Delta` by moving each angle by `<~ Y^(-1/2)`. Counting `beta`-tuples against the fixed profile
(`e^(W)` tuples, probability `~ (r^2/Y)^((M-1)/2)`) predicts such fakes for `Y <~ P^(2r/(M-1))` and none beyond;
proving either side is again a Diophantine problem about actual Gaussian integers.

## 6. The dichotomy, and what it leaves open

Putting §§2-5 together, every analytic input about Gaussian primes and Hecke characters that has been proposed
falls into one of three bins (proved statements, with the scope stated):

| input | example | separates actual from fake? | growth? |
|---|---|---|---|
| fibre (Hecke characters on the norm fibre, divisors of `N`) | `S_k(N)` (2.1), sector counts in the fibre, Riesz products | no: functions of the angle vector (Prop. A) | capped (`outside.md` A.3, `lw.md` 6.1) |
| robust ensemble (error `>= c X^gamma`, or large sieve `>= c(X + K X^(1/2))`) | Hecke large sieve, Hecke PNT (zero-free region, zero density, GRH), Kubilius/Ricci/Harman-Lewis sectors, sieve bounds in sectors and short intervals, joint residue-sector counts | no: grafted fakes satisfy them up to `1 + o(1)` (Thm B) | capped at `W <= (2A + o(1)) M log M` |
| absolute-value counting of prime tuples | majorant + bounds for `|S_j(mu)|` | no: invariant under column rotations (Prop. C) | cannot certify `T = 0` |
| pointwise lattice separation | `LCG(phi_j; p_j + 1)` | yes, exactly (Thm D) | equivalent to the original problem |
| pointwise prime separation | `PCG(phi_j; infinity)`; the sharp `(X + K)` large sieve | yes under GPA (Prop. D.1), conjecture beyond GRH; locally rigid at `theta_j` (Lemma D.2) | equivalent to the original problem (under GPA) |

So **the distribution of Gaussian primes in narrow sectors, used either statistically or pointwise, does not give
a growth improvement**: statistically it cannot see the `O(M)` grafted angles; pointwise it is the lattice itself.
The task asked which input is true for actual circles and false for the fakes. The answer here is the
separation `|sin(theta_j - arg alpha)| >= (p_j Norm alpha)^(-1/2)` against external lattice points of norm
`<= p_j + 1` (Theorem D). It is false for every fake that is not the shadow of an actual cluster, but it
reconstitutes the actual problem exactly.

**What is not covered.**
1. *The regime `P_+(N) < M^(2 - delta)`.* Theorem B needs the grafted primes to be `>= M^A`, `A > 2`. No fake with all
   primes `<= M^(2 - delta)` is known (the equal-weight windows `[P, P e^(1/(4r))]` of all known fakes need
   `P >~ r^2 log r`). In this regime the cluster's angles are a positive-power fraction (`>= P^(-1/2 + delta')`) of all
   Gaussian prime angles up to `P`, and the large sieve at `K ~ P^(1/2)` is non-trivial for them. It is open whether
   near-extremal actual clusters can live there and, if they can, whether a large-sieve statement constrains them
   (it would have to force angular concentration, which nothing in the profile structure suggests).
2. *Inputs that couple angles and residues of the same primes beyond P3-char.* Theorem B covers joint statements
   only where the residue class is capped. For P3-all and P3-G the caps are open (round 9).
3. *Truncated prime cross-gaps* `PCG(·; M^K)`: heuristically fake-killing for `K` beyond a constant depending on
   the prime sizes (§5.2), unproved either way, and reducing again to a statement about actual Gaussian integers.

## 7. Computations

All scripts are in `analytic_checks/`, outputs next to them. None is used in a proof except as an exact
re-check of a proved lemma.

**7.1 Window fluctuations (`window_fluct.py`; evidence for Remark 3.4(1)).** Genuine Gaussian primes in the fake
windows `[P, P e^(1/(8q))]`:

| `P` | `q` (`M = q+1`) | Gaussian primes `n` | `2s = 4q` | 64 sectors: mean, max/median `|dev|` | median `|S(k)|`, `k <= 2000` | `sqrt n` |
|---|---|---|---|---|---|---|
| `10^7` | 11 | 7078 | 44 | 110.6; 22.6 / 4.0 | 49.6 | 84.1 |
| `10^8` | 11 | 61958 | 44 | 968.1; 44.1 / 16.0 | 134.1 | 248.9 |
| `10^8` | 23 | 29470 | 92 | 460.5; 37.5 / 6.5 | 95.5 | 171.7 |
| `10^9` | 47 | 128144 | 188 | 2002.2; 48.8 / 21.5 | 192.9 | 358.0 |

At these sizes (`P` far below the `(40M)^4` of `lw.md` Thm 6.1) the coarse sector deviations are of the same order
as `2s`, and the typical Hecke sum is a few times `2s`; individual `|S(k)|` fall below `2s` for 17-49% of `k`.
Theorem B does not rely on this: it uses only the shape `>= c X^gamma` of the published error terms. The table
shows that at small heights statements finer than those error terms would be true and could see a graft; the
ratio `sqrt(n)/s ~ (P/(M^3 log P))^(1/2)` grows with `P = M^A`, `A > 3`. (The median `|S(k)|` is about `0.54 n^(1/2)`,
below the Rayleigh value `0.83 n^(1/2)` of independent angles: lattice regularity is visible.)

**7.2 GPA, random directions (`gpa_check.py` part A; evidence).** All Gaussian primes of norm `<= 2*10^6`
(148,945 classes mod units). For 2000 uniform `phi` and each `h`, the smallest-norm witness of GPA(`h`):

| `h` | witnesses found | median witness norm | max witness norm |
|---|---|---|---|
| 5 | 2000/2000 | 9 | 29 |
| 13 | 1999/2000 | 5 | 199,657 |
| 101 | 2000/2000 | 41 | 730,889 |
| 1009 | 1983/2000 | 1493 | 1,786,753 |

The misses are directions close to a low-height lattice direction, where the witness norm grows like
`~ eps^(-2)`: at `phi = arg(3 + 4i) + eps` (`3 + 4i = (2+i)^2`, composite, `h = 101`) the witness norms are 557, 1129,
10733, 118061 for `eps = 10^(-2), 3*10^(-3), 10^(-3), 3*10^(-4)`. This is the expected behaviour (GPA has no uniform
height bound) and not a counterexample.

**7.3 GPA at lattice directions and Lemma D.2 witnesses (exact).**
* `gpa_check.py` part B: every primitive direction `beta` with `h < Norm beta <= 20h` (`h = 5, 13, 101`: 44, 114, 914
  directions) has a GPA witness (a prime on a line `0 <= |m| < (Norm beta/h)^(1/2)`); worst witness norm
  `16.7 Norm beta`.
* `lines_check.py`: for every split `p <= 10^5` (4783 primes) and `pi = a + bi`, both lines `Im(varpi conj pi) = +-1`
  carry a Gaussian prime with `Re(varpi conj pi) > 0`; the worst needs 113 steps along the line
  (`p = 91577`, norm ratio `1.28*10^4`). So Lemma D.2 applies to every Gaussian prime of norm `<= 10^5`.

**7.4 Lemmas 5.2 and 5.3 (`lcg_check.py`; exact).** Lemma 5.2: 1,384,674 pairs (`p <= 2000`, `Norm alpha <= 3000`),
no violation. Lemma 5.3: for `n in {5, 13, 17, 29, 101}` and all 19,600 primitive directions `(u, v)` with
`|u|, |v| <= 40`, a violator of `LCG_n(·; n + 1)` exists iff `u^2 + v^2 > n` (exact rational comparison).

## 8. Open points

1. **Regime `P_+(N) < M^(2-delta)`.** Either construct fakes there (then Theorem B extends, with a spreading argument
   instead of grafting, because the graft would move a positive-power fraction of the primes) or show that
   near-extremal actual clusters cannot live there. The latter would be a restricted theorem; no mechanism is known.
2. **GPA.** Prove the metric version (almost every `phi`, by a second-moment argument with the large sieve and a
   Chung-Erdos lemma across dyadic scales). This would make PCG provably false for the random fakes of `lw.md`
   Theorem 6.1, though not for all fakes.
3. **Truncated PCG at polynomial height** (§5.2): does `PCG(·; M^K)` cap or collapse for large fixed `K`? The
   heuristic count puts the threshold near `Y ~ P^(2r/(M-1))`.
4. **Joint angle-residue inputs in P3-all / P3-G.** Theorem B transfers any cap of a residue class to that class
   plus robust ensemble inputs. It cannot help where the residue class is uncapped.
5. *Citations.* Harman-Lewis (Mathematika 48 (2001)), Kubilius, Ricci, and the Hecke large sieve are quoted from
   memory; only the shape of their error terms is used (Theorem B needs `>= c X^gamma`).
