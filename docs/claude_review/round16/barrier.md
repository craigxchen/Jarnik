# Round 16 — the barrier / counterexample side of AC (route 5, T = 1)

Status: final for this round. Unrefereed.

## 0. Summary and verdict

**No growth improvement, no uniform bound, route 5 neither capped nor unblocked.** What is new is a sharper map of
the barrier side of AC / MIN^pair at `T = 1`, `w = sqrt(8(n-1) ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`, regime
`n - 1 >= 316 ln 4r`, barrier scale `W = O(r log n)`.

Proved here (modulo quoted refereed results; H_light is a stated hypothesis where used):
1. **AC is box-free** (Prop 2.1, Cor 2.2): AC(u) holds iff no subset of the ternary cube of density `>= 1/2` has `u` as an
   `eps_0`-near tangent-cone point. A proof of AC would decide the barrier for all fake boxes and all half-cubes at once;
   AC is strictly stronger than MIN^pair and controls `kappa_V`, not `kappa_H`.
2. **Only light characters need Diophantine control** (Thm 2.3): the barrier-scale transfer of fake Thm 2.1 holds with
   (2.2) replaced by `|lambda.u*| >= (13 s*/C) e^{-tau sigma_lambda/2}` for characters with `sigma_lambda < r/2 + O(n/tau)`;
   heavy relations (even exact ones) are free. So the barrier test is exactly a MIN^char test on light characters.
3. **Collision strata never give barrier points for fakes** (Prop 2.4): `sdc_q` has `sigma_xy = 1` for every pair
   (exact, `q <= 43`); the fakes' own pair gaps then beat every transversal-vanishing ceiling once `r > 2m'`.
4. **Collision-free rational strata would be fatal** (Prop 2.5): a rational `L ⊂ TC_H` with no pair relation, light
   relations of slack `>= c_0 r` and transversal order `o > 1/(2c_0)` (`o >= 3` for `sdc`) yields barrier-scale fakes.
5. **Products/sections give nothing** (Prop 3.1, with a self-contained proof of additivity of `nu` on products).
6. **Cyclic-symmetry dichotomy** (Thm 3.2, Cor 3.3): if a row permutation of order `N > 2m` preserves the patterns
   (translations of `sdc_q`; Singer cycles of Sylvester), then on each isotypic circle `sup kappa_H <= sqrt(m+1) inf kappa_H`,
   and by Galois the circles of a whole Galois class are in `TC_H` or none is. The circles satisfy the Diophantine
   condition automatically (Lemma 3.4, under H_light), so **route 5 is capped by `sdc` fakes if one Fourier circle is an
   `eps_0/sqrt(m'+1)`-near zero** (Thm 3.5) — a single explicit algebraic question, not decided.
7. Non-linear base curves of low degree, and through all coordinate points, are excluded (§4 row 7).

Exact / evidence data: the M8 tangent cone has no collision-free real zero, including on all 42 additive-quadruple
hyperplanes (high-precision search, evidence); in `sdc7`, `sdc11` (`L = 2`, `m = 8`) all tangent forms are
translation-invariant and nonzero on the Fourier circles, so the circles are **not** in `TC` (exact mod two primes);
`sdc19` (`L = 2`, `m = 5, 6`, where `2m < q` so Thm 3.2 applies): circle not in `TC`. AC is consistent with all exact data: its only exact violations
(`sdc7`, `sdc11` zeros) occur in underdetermined toys of density `< 1/2` at `m = rho`.

**Answer to the three targets.** (1) No nondegenerate near-zero was constructed; every known mechanism yields
collision-bearing points (excluded at the barrier scale) except two explicitly identified open classes (Fourier
circles; collision-free rational strata), each of which would cap route 5. (2) Not capped; the property that rules
barrier points out in every analysed case is "**cheapness comes only from few-value (integral) projections, and few
values force light relations**" — (P1) for rational strata, (P2) for directions with no light relation; this is the
content AC must supply (§5). (3) Consistent (§6).

## 1. Setting, notation, what is quoted

**Labels.** *Proved* = proved here with all quantifiers (or a one-line consequence of a quoted refereed result, named).
*Quoted* = refereed record result, cited by file and number. *Exact* = integer / exact-rational computation.
*Evidence* = floating point or finite search; nothing proved uses it. *[M]* = literature quoted from memory;
*[L]* = literature located by search in this round (abstract / snippet level). Scripts: `barrier_checks/`
(`x.py` writes `x_out.txt`). Record = `docs/claude_review/` (read-only).

**Objects (record notation, `round12/README.md` §1).** Cluster of `n` points, `d = n - 1`, `r` split primes with
weights `log p`, `W = sum log p`, sign columns `beta_p` with patterns `P_p` (rows with entry `-1`, base row
omitted), so `<k,beta_p> = +-2 k(P_p)` on `Z^n_0`.
* `T = {-1,0,1}^d` (relative coordinates, base coordinate free), `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`,
  `A = {k' in T : |k(P_p)| <= j for all p}`; regime `d >= 316 ln 4r`, where `|A| >= |T|/2` and
  `m' = ceil(4w) <= m_max` (t1rigour §2, quoted). Note `m'/(2j) in [4, 4 + O(1/w)]`.
* `V(A,m)` = measures on `A` annihilating `Poly_{<m}|_A`; `V_H` = span of the integer `c in V(A,m')` with
  `|c|_inf <= H`; `Phi_c(u) = sum_k c_k <k',u>^m`; `kappa_c(u) = |Phi_c(u)|/(|c|_1 d^m)`;
  `kappa_H = sup_{c in V_H} kappa_c`, `kappa_V = sup_{c in V(A,m)} kappa_c`; `TC_H = {kappa_H = 0}`.
* `J(A)` = leading-form ideal of `I(A)`; `u in TC(A,m)` iff `l_u^m in J(A)` (`l_u = <.,u>`); tangent degree
  `nu_A(u) = max{i : l_u^i notin J(A)}` (zetc Lemma 1.1, quoted).
* **Duality** (nearstrata Prop 2.1, quoted; it is finite-dimensional LP duality and holds for every finite
  `S subset R^d`): `d^m kappa_V^{(S)}(u) = dist_inf(l_u^m|_S, Poly_{<m}|_S)`.
* Step 0 level `eps_0 = C d e^{-W/4}/(m'+1)`; pair excesses `eps_xy`; gap region
  `G = {u : |u|_inf = 1, |u_x - u_y| >= e^{-eps_xy/2} for all pairs (u_base = 0)}`.
* **MIN^pair** (referee_nearstrata §4, quoted): `min_{u in G} kappa_H(u) > eps_0`. With Step 0, rigidity and
  separation it gives `n = O(loglog R)` on all norms (nearstrata Thm 4.1 + referee; quoted).
* **AC(u)** (steelman S3): for every polynomial `Q` of degree `m'` whose top homogeneous part is `l_u^{m'}`,
  `#{k in T : |Q(k)| <= delta} < |T|/2`, with `delta := d^{m'} eps_0`. (AC_A: the same with `T` replaced by `A`.)

**Fake pattern data** (fake.md §1, quoted). Pattern set `T_pat`, multiplicities `mu_S`, windings `k`;
`L_T = {lambda in Q^n : sum lambda = 0, lambda(S) in Z for all S}`; `lambda` *relevant* if `lambda.k in Z`;
normalised slack `sigma_lambda(mu) = sum_S mu_S (|lambda(S)| - 1/2)`; (Mar): (AI), `sigma_lambda >= rho_1 ||lambda||_1` on
relevant `lambda`, irrelevant `lambda` have `||lambda.k|| >= 1/den_T`. Barrier profiles: `sdc_q` (`r = 2n - 2`,
`rho_1 = 1/2`), `pal_q` (`r = 9n - 8`, `rho_1 = 1`); `syl_k` is actual-type (zero pair margin), not a fake.
Angle models: primes in `[P, P e^{1/(4r)}]`, `tau = log P`, `W in [r tau, r tau + 1/4]`, positions in
`[0, Delta]^n`, `Delta = C e^{-W/4}`.

**Fake Theorem 2.1 (corrected; fake.md §2.1 + referee_fake §7(a), quoted).** If `u*` is nondegenerate and
`theta > 0`, then for `P` large there are gap-admissible angle models with `|hat u - u*|_2 <= theta e^{-W/4}` and
`kappa_c(hat u) <= kappa_c(u*) + m theta e^{-W/4}/sqrt d`; at the barrier scale `W = O(r log n)` this needs the
Diophantine condition (2.2) `|lambda.u*| >= n^{-a||lambda||_1}` for relevant `lambda`, `||lambda||_1 <= r/rho_1 + 1`.

## 2. Reduction: what a counterexample to AC must look like

### 2.1 AC is a statement about `u` alone (box-free)

**Proposition 2.1 (proved).** Fix `d, m >= 1`, `delta > 0`, `u in R^d`. The following are equivalent.
1. AC(u) fails: some `Q` with `deg Q <= m` and top part `l_u^m` has `#{k in T : |Q(k)| <= delta} >= |T|/2`.
2. There is `S subset T` with `|S| >= |T|/2` and `dist_inf(l_u^m|_S, Poly_{<m}|_S) <= delta`.
3. There is `S subset T` with `|S| >= |T|/2` and `kappa_V^{(S)}(u) <= delta/d^m`, where `kappa_V^{(S)}` is the
   `H = infinity` kappa of the finite set `S` (sup over all measures on `S` annihilating `Poly_{<m}|_S`).

*Proof.* (1)⇒(2): put `S = {|Q| <= delta}` and `R = l_u^m - Q in Poly_{<m}`; then `|l_u^m - R| = |Q| <= delta` on `S`.
(2)⇒(1): `Poly_{<m}|_S` is finite-dimensional, so a best approximant `R` exists; `Q := l_u^m - R` has top part `l_u^m`
and `|Q| <= delta` on `S`. (2)⇔(3) is the duality of §1 applied to the finite set `S`. QED

**Corollary 2.2 (proved).** In the regime (`delta = d^{m'} eps_0`, `m = m'`):
1. AC(u) holds iff **no** subset of the ternary cube of density `>= 1/2` has `u` as an `eps_0`-near point of its
   `H = infinity` tangent cone. The profile, the primes and the box enter only through `delta`.
2. Hence AC on a set `U` of directions implies, simultaneously for every profile (actual, Sylvester-type, every
   (Mar) fake pattern set) whose slab box has density `>= 1/2`, that `kappa_V^{(A)} > eps_0` on `U`. In particular a proof
   of AC on the gap regions would decide the barrier question of fake Thm 2.1 for **all** fake pattern boxes at once,
   and also for half-cubes that are not slab boxes at all.
3. AC feeds MIN^pair only through `kappa_V >= kappa_H` in the wrong direction: AC gives `kappa_V > eps_0`, MIN^pair
   needs `kappa_H > eps_0` (the skeptic's point 3; `V_H = V(A,m')` is open). AC_A (density inside `A`) is box-dependent
   and its failure set contains `{kappa_V^{(A)} <= eps_0}`.

So the barrier side of AC has two levels. **Level 1 (box-free):** find any half-cube `S` and a gap-region `u` with
`kappa_V^{(S)}(u) <= eps_0`; this refutes AC as a route but says nothing about clusters. **Level 2 (fake boxes):**
find such `u` for a (Mar) fake pattern box `S = A_f` with the Diophantine property of §2.2; this caps route 5.
§3-§4 attack level 2 and, through Prop 2.1, record what is known at level 1.

### 2.2 Only light characters need Diophantine control

The record states the barrier-scale transfer with (2.2) for all relevant `lambda` with `||lambda||_1 <= r/rho_1 + 1`.
The proof of fake Thm 2.1 gives more: heavy characters are handled by the random perturbation whatever `lambda.u*` is.

**Definition.** For pattern data `(T_pat, mu, k)` and `tau > 0` put
`sigma*(tau) := r/2 + 1 + (2/tau) [ n ln(2 den_T (r+2)/rho_1 + 1) + ln(1800 n^{3/2}/(theta C ell_0)) ]`
(`ell_0 = min ||lambda||_1` over `L_T \ 0`, `theta` as in Thm 2.1). A relevant `lambda != 0` is **light** if
`sigma_lambda(mu) < sigma*(tau)`, **heavy** otherwise. Light `lambda` have `||lambda||_1 <= sigma*/rho_1`. For `tau`
large, `sigma* -> r/2 + 1`; at the barrier scale `tau = K ln n` it is `r/2 + O(n/K)`.

**Theorem 2.3 (refined barrier transfer; proved, by re-running the proof of fake Thm 2.1).** Let `(T_pat, mu, k)`
satisfy (Mar); let `u* in R^d`, `|u*|_inf = 1`, `C in (0,1]`, `theta in (0,1]`, and let `tilde u`, `s* in [1,2]` be as in
fake Thm 2.1. Assume `r >= 2n - 2`, `n >= 11`, `tau >= max(4, tau_2)` (`tau_2` = the class-2 and prime-window
threshold of Thm 2.1), `sigma*(tau) <= r + 1`, and `||lambda||_1 Delta <= 1/2` for all light `lambda`. If

```text
|lambda.u*|  >=  s* max( 13 e^{-tau sigma_lambda(mu)/2} / C ,  2 eta ||lambda||_1 )      for every light lambda,     (2.2')
```

with `eta = theta e^{-W/4}/(8 sqrt n)`, then for every assignment of primes in `[P, P e^{1/(4r)}]` there is a
gap-admissible angle model with `|hat u - u*|_2 <= theta e^{-W/4}` and `kappa_c(hat u) <= kappa_c(u*) + m theta e^{-W/4}/sqrt d`
for every `c`. **No condition is imposed on heavy `lambda`, in particular none on `lambda` with `lambda.u* = 0`
if they are heavy.**

*Proof.* Use the positions `delta = (t/(1+2eta))(tilde u + xi + eta 1)` of the proof of Thm 2.1 and split the nonzero
`lambda in L_T` into: (1') light relevant; (2) irrelevant with `||lambda||_1 <= 1/(n Delta)`; (3') heavy relevant, and
irrelevant with `||lambda||_1 > 1/(n Delta)`.
* *Class (1').* `|lambda.xi| <= eta ||lambda||_1 <= |lambda.tilde u|/2` by the second term of (2.2') (`lambda.tilde u =
  lambda.u*/s*`). Hence `|lambda.delta| >= Delta |lambda.tilde u|/6` for every `xi`, and `|lambda.delta| <= pi/4` by the
  choice of `P`. Since `lambda.k in Z`, (2.1) holds once `Delta |lambda.tilde u|/12 >= arcsin(e^{-w/2})`; with
  `e^{-w/2} <= e^{-W/4 - tau sigma_lambda/2 + 1/16}` (Thm 2.1, weights) and `arcsin x <= 1.01 x` for `x <= 1/8`, this is
  implied by `C |lambda.tilde u| >= 12 * 1.01 * e^{1/16} e^{-tau sigma_lambda/2}`, i.e. by the first term of (2.2')
  (`12.9 < 13`).
* *Class (2)*: verbatim Thm 2.1 (uses only `tau >= tau_2`).
* *Class (3').* The class-3 estimate of Thm 2.1 never uses `lambda.tilde u`: the bad set of (2.1) is a union of at most
  `||lambda||_1 Delta + 1` intervals of length `<= 8 e^{-w/2}` in the variable `lambda.delta`, and `lambda.xi` has density
  `<= n/(2 eta ||lambda||_1)`. So `Prob(bad_lambda) <= 24 n (1 + 1/(Delta ||lambda||_1)) e^{-w/2}/eta`. For heavy relevant
  `lambda`, insert `e^{-w/2} <= e^{-W/4 - tau sigma_lambda/2 + 1/16}`, `Delta = C e^{-W/4}`, `1/eta = 8 sqrt n e^{W/4}/theta`,
  `W <= r tau + 1/4`, `||lambda||_1 >= ell_0`, `C <= 1`:

  ```text
  Prob(bad_lambda) <= (450 n^{3/2}/(theta C ell_0)) exp( -tau (sigma_lambda - r/2)/2 ) .
  ```

  The number of relevant `lambda` with `sigma_lambda < s + 1` is `<= N((s+1)/rho_1) <= (2 den_T (s+1)/rho_1 + 1)^n`
  (Thm 2.1). Group the heavy `lambda` by `s = floor(sigma_lambda) >= s_0 := floor(sigma*)`. Consecutive group bounds
  have ratio `<= (1 + 1/(s+1))^n e^{-tau/2} <= e^{2n/r - tau/2} <= 1/2` (`s >= r/2`, `r >= 2n - 2`, `n >= 11`,
  `tau >= 4`). The first group has `s_0 + 1 <= r + 2` and `s_0 - r/2 >= sigma* - 1 - r/2`, so its bound is
  `<= (450/1800) e^{-tau/2}`. Hence the heavy relevant sum is `<= 1/2`. The irrelevant part of class (3') is doubly
  exponentially small as in Thm 2.1 (and `< 1/4` for `P` large).

So some `xi` satisfies (2.1) for all `lambda`, and `exact.md` Prop 1.1 gives the gap-admissible model; the direction
and `kappa` estimates are those of Thm 2.1 (with the referee's constant). QED

**Reading.**
1. For pairs, (2.2') is `|u*_x - u*_y| >= (13 s*/C) e^{-tau sigma_xy/2}`: the pair-gap region of the fake's own data
   (pair excess `eps_xy ~ tau sigma_xy`). For light non-pair characters it is the **character-gap** inequality. So at
   the barrier scale the barrier test is exactly a test of fatpoint's `MIN^char` restricted to light characters: a
   (Mar) fake box caps route 5 iff (sufficient direction) it has `u*` with `kappa_H(u*) <= eps_0/2` in the light-character
   gap region `G_light(tau)` (take `theta <= C d^{3/2}/(2 m'(m'+1))`).
2. For random-like patterns, `E|lambda(S)|` for a balanced `lambda in {-1,0,1}^n` with `||lambda||_1 = 2, 4, 6, 8` is
   `1/2, 3/4, 15/16, 35/32`, so `sigma_lambda ~ r (E|lambda(S)| - 1/2)` is `0, r/4, 7r/16, 19r/32`: at large `tau` the
   light characters are the balanced ones with `||lambda||_1 <= 6` (plus borderline unbalanced ones such as
   `2e_x - e_y - e_z`, `E|lambda(S)| = 1`). Exact sample for `sdc_q` (`sigma_pairs.py`, 200 random `lambda` each,
   `q = 7..43`): pairs `sigma = 1` exactly; 4-term `sigma/r ~ 0.18..0.24`; 6-term `0.45..0.52` (borderline);
   8-term `0.53..0.65`. *Evidence-level reading*, not used in proofs.
3. Every relation `lambda.u* = 0` with heavy `lambda` (e.g. `sum_{x in Z_q} u_x = 0` in §3.3, whose `sigma ~ r q/4`) is
   harmless. This is what makes the symmetric candidates of §3.3 admissible.

### 2.3 Rational strata: which ones can produce barrier points

**Lemma 2.0 (transversal ceiling; proved).** If every `Phi_c`, `c in V_H`, vanishes to order `>= o` along a linear
subspace `L` (i.e. `sum_k c_k <k',y>^{m'-i} <k',v>^i = 0` for `i < o`, all `y in L`, all `v`), then
`kappa_H(y + eta v) <= 2^{m'} eta^o` for `y in L`, `|y|_inf, |v|_inf <= 1`, `0 <= eta <= 1` (unnormalised `kappa`).
*Proof.* The coefficient of `eta^i` in `Phi_c(y + eta v)` is `binom(m',i) sum_k c_k <k',y>^{m'-i}<k',v>^i`, of modulus
`<= binom(m',i) |c|_1 d^{m'}`; it vanishes for `i < o`, and `sum_{i>=o} binom(m',i) eta^i <= 2^{m'} eta^o`. QED
(t1rigour §5 has the sharper `(2+eta)^{m'} eta^{o}`-type ceiling for projection-cheap spans.)

**Proposition 2.4 (proved).** Let `(T_pat, mu, k)` satisfy (Mar) and let `lambda_0` be a light character, e.g. a pair.
Let `u* in G_light(tau)` (i.e. satisfy (2.2')) and let `L` be any linear subspace with `lambda_0.L = 0` (e.g. any `L`
with a collision `x ~ y`, `lambda_0 = e_x - e_y`). Then `dist_2(u*, L) >= (13/C) e^{-tau sigma_{lambda_0}/2}/||lambda_0||_2`.
Consequently the transversal ceiling of Lemma 2.0 (any version `kappa <= c_1^{m'} eta^o` with `c_1 >= 1`; `o <= m'`),
evaluated at `eta = dist(u*, L)`, exceeds `eps_0` as soon as

```text
o sigma_{lambda_0} tau/2  <  W/4 + m' ln c_1 + o ln(13/(C ||lambda_0||_2)) - ln(C d/(m'+1)) ,
```

which, with `W >= r tau`, `C <= 1`, `||lambda_0||_2 <= 13`, holds as soon as `2 o sigma_{lambda_0} < r - 4 ln(d)/tau`;
for all `L` with a collision it suffices that `r > 2 m' sigma_pair + 4 ln(d)/tau`, `sigma_pair = max_{pairs} sigma_xy`.
For `sdc_q`, `sigma_xy = 1` for **every** pair (exact, `sigma_pairs.py`, `q = 7..43`) and `r = 2n - 2`, so the condition
is `n > m' + O(1)`: the whole regime (`m' = O(sqrt(n log n))`).

*Proof.* `|lambda_0.u*| = |lambda_0.(u* - y)| <= ||lambda_0||_2 |u* - y|_2` for `y in L`; insert (2.2'). The display is
the logarithm of `c_1^{m'} eta^{o} > eps_0` with that `eta` (and `c_1 >= 1`). QED

*Remark (scope).* This is a statement about the only **proved** source of near-zeros (transversal vanishing along
rational strata: cheap spans, projection-cheap families, XOR/Fano blocks, parity sections, J5's coordinate spans). It
shows that for fake data with bounded pair margins this source never produces barrier points: J5-type failures
(refJ §1) need *concentrated* pair excess `eps_xy >> tau sigma_xy` inside the classes, which the actual-type data
allow but the fixed, uniform margins of a fake profile do not. It does **not** exclude near-zeros that are far from all
rational strata. (For `sdc_q`, `sigma_xy = #{S separating x,y} - q` over the `2q` patterns, and the computation gives
exactly `q + 1` separating patterns for every pair.)

**Proposition 2.5 (the rational strata that are NOT excluded are fatal; proved).** Let `(T_pat, mu, k)` satisfy (Mar),
`r >= 2n - 2`, and let `L` be a rational subspace contained in `TC_H(A, m')` along which every `Phi_c`, `c in V_H`,
vanishes to order `>= o` (e.g. `o = m' - nu_A(L)`). Let `Lambda_L = {lambda relevant integral : lambda.L = 0}` and
`sigma_min(L) = min{sigma_lambda : lambda in Lambda_L light}` (`+infinity` if `Lambda_L` has no light element). Suppose
1. `L` is **collision-free**: no pair `e_x - e_y` (and no `e_x`) lies in `Lambda_L`;
2. `o sigma_min(L) > r/2 + 1 + (2/tau)(m' ln 6 + ln(4 d/C) + o ln(26/C) + o (d-1) ln(2r+3))`;

3. H_light(`c_0`) (§3.3): light non-pair `lambda` have `sigma_lambda >= c_0 r`.

Let `b = dim L` and `h_L` the height of an integral basis `w_1..w_b` of `L`. Then for every
`tau >= tau_L := max( 2 ln(60 b (N_2+n) h_L/C) , (2/(c_0 r)) (2 b ln(8 (r+1) h_L) + ln(100 (N_2+n)/C)) ) + 4`
(and the thresholds of Theorem 2.3) there are gap-admissible fakes with `kappa_H(hat u) <= eps_0`: **route 5 is capped
by this pattern set** at `W ~ r tau`. Condition 2 at `tau = K ln n` reads `o sigma_min > r/2 + O(o n/K)`; with
`sigma_min >= c_0 r` it holds for `o > 1/(2c_0)` (`o >= 3` for `c_0 = 0.18`) and `K >= K_0(c_0, o)`; `tau_L = O(log n)` when
`h_L = n^{O(1)}`, so this is the barrier scale. If `Lambda_L` has no light element, condition 2 is void and points of `L`
itself are barrier points.

*Proof.* *Choice of `y in L`.* Take `t in {0..T}^b` with `T = 2(N_2 + n)`, avoiding the zero sets on `L` of the `N_2 + n`
pair forms (all nonzero on `L` by condition 1; each kills at most a `1/(T+1)` fraction of `t`'s along a coordinate line),
and put `y_1 = sum t_i w_i`; then `|lambda.y_1| >= 1` for all pairs. Put `y_2 = sum_i B^{i-1} w_i` with `B = 4(r+1)h_L + 1`:
for every integral `lambda` with `||lambda||_1 <= 2r+2` not vanishing on `L`, `lambda.y_2 = sum_i B^{i-1}(lambda.w_i)` is a nonzero
integer (base-`B` uniqueness, `|lambda.w_i| < B/2`). Let `y = (y_1 + epsilon y_2)/|y_1 + epsilon y_2|_inf` with
`epsilon = 1/(4 |y_2|_inf ... )` chosen so that `epsilon |lambda.y_2| <= 1/2` for pairs, i.e. `epsilon = 1/(4(2r+2) b B^{b-1} h_L)`.
Then pairs have `|lambda.y| >= 1/(2 |y_1 + epsilon y_2|_inf) >= 1/(3 b T h_L)`, and every light non-pair `lambda notin Lambda_L`
has `lambda.y_1 + epsilon lambda.y_2 != 0` with modulus `>= epsilon` (if `lambda.y_1 = 0`) or `>= 1/2` (else), so
`|lambda.y| >= epsilon/(2 b T h_L)`.
*Choice of the transversal perturbation.* `v = (1, B', ..., B'^{d-1})/B'^{d-1}`, `B' = 2r+3`, has `|lambda.v| >= B'^{-(d-1)}`
for every nonzero integral `lambda` with `||lambda||_1 <= 2r+2`. Put `u* = y + eta v`, `eta = 26 e^{-tau sigma_min/2} B'^{d-1}/C`.
*(2.2').* Light `lambda in Lambda_L`: `|lambda.u*| = eta|lambda.v| >= 26 e^{-tau sigma_min/2}/C`. Pairs: `|lambda.u*| >= 1/(3bTh_L) - 2 eta
>= 26 e^{-tau/2}/C` for `tau >= tau_L`. Light non-pair `lambda notin Lambda_L`: `|lambda.u*| >= epsilon/(2bTh_L) - eta(2r+2) >=
26 e^{-tau c_0 r/2}/C` for `tau >= tau_L` (H_light). The `eta`-terms of (2.2') are `O(n e^{-W/4})`.
*kappa.* Lemma 2.0 gives `kappa_H(u*) <= 2^{m'} eta^{o}` unnormalised; `|u*|_inf >= 1/2`, so `kappa_H(hat u*) <= 4^{m'} eta^o
<= 6^{m'} eta^o`, and `6^{m'} eta^o <= eps_0/2` is the logarithm of condition 2 with `W <= r tau + 1/4`.
Theorem 2.3 (with `theta` small) gives the fakes. QED

**Consequence (the barrier, located).** By Prop 2.4 and Prop 2.5, at the barrier scale and for fake data with bounded
pair margins:
* rational tangent-cone strata **with a collision** never yield barrier points;
* rational tangent-cone strata **without collisions**, whose relations are heavy or are light with
  `o sigma_lambda > r/2`, always do. For `sdc_q` the lightest non-pair characters have `sigma >= 0.18 r` (exact for
  `||lambda||_1 = 4`, `q <= 31`, §3.3), so **any collision-free rational stratum of transversal order `o >= 3`** in an
  `sdc_q` box would cap route 5. (For Sylvester data the additive quadruples have `sigma = 0` like pairs, so there
  collision-free strata inside quadruple hyperplanes are excluded too; this is the toy content of J4, §6.)
* everything else is a near-zero near an irrational real component of `TC_H`, or a non-Lojasiewicz near-zero.

The record's README-12 §3.2 already notes that "a collision-free rational stratum is fatal" for MIN^pair; Prop 2.5 shows it
is fatal for the **barrier** as well, i.e. it would cap route 5 and not merely refute a hypothesis. §3 tests the irrational
candidates (symmetric circles) and §4 lists the rational ones.

## 3. Fake pattern boxes: tangent cones and near-zeros

### 3.1 Products and sections give nothing new

**Proposition 3.1 (proved).** Suppose the non-base rows split as `X_1 ⊔ ... ⊔ X_g` and every pattern lies in one
block (block-diagonal profile). Put `A_b = {k in {-1,0,1}^{X_b} : |k(P)| <= j for the patterns P subset X_b}`.
1. `A = A_1 x ... x A_g`, `J(A) = sum_b J(A_b) S` and `nu_A(u) = sum_b nu_{A_b}(u_b)` (`u_b = u|_{X_b}`).
2. `nu_{A_b}(u_b) >= 2 min(j, s_b)`, `s_b = #{x in X_b : u_x != 0}`.
3. Hence every `u in TC(A,m)` with no collision with the base point (`u_x != 0` for all `x`) satisfies
   `sum_b 2 min(j, |X_b|) < m`. At `m = m' <= 8j + 4` (regime): at most 4 blocks have size `>= j`, and then each block
   `b` carries a base-collision-free point with `nu_{A_b}(u_b) <= m' - 1 - 2j(g-1)`.

*Proof.* (1) `A` is the product because each constraint involves one block. `S/I(A) = ⊗_b S_b/I(A_b)` with the degree
filtration of the left side equal to the tensor product of the degree filtrations, so `gr` commutes with `⊗`:
`S/J(A) = ⊗_b S_b/J(A_b) = S/sum_b J(A_b)S`; since `J(A) ⊇ sum_b J(A_b) S`, the two ideals coincide. For `g = 2`,
`l_u^m = sum_i binom(m,i) l_{u_1}^i l_{u_2}^{m-i}` is bihomogeneous componentwise and `S/J(A)` is bigraded with
`(S/J)_{(i,m-i)} = (S_1/J_1)_i ⊗ (S_2/J_2)_{m-i}`, where the pure tensor `l_{u_1}^i ⊗ l_{u_2}^{m-i}` vanishes iff one factor
does. So `l_u^m in J(A)` iff for every `i` either `i > nu_1` or `m - i > nu_2`, iff `m > nu_1 + nu_2`; induct on `g`.
(2) Choose `X subset supp(u_b)`, `|X| = min(j, s_b)`. The cube `B = {k : supp k subset X}` lies in `A_b` (`|k(P)| <= |X| <= j`),
so `J(A_b) subset J(B) = (k_x^3 : x in X) + (k_y : y notin X)`, and `l_{u_b}^i in J(B)` iff `i > 2|X|` (the monomials
`k^alpha`, `alpha <= 2` on `X`, of degree `i` have coefficient `i!/alpha! u^alpha != 0` iff `i <= 2|X|`). (3) follows. QED

*Reading.* Products (the record's `M8 x M8`, sums of blocks) reduce the barrier question to their blocks, and with
`>= 5` blocks of size `>= j` they cannot carry collision-free tangent-cone points at all. Fake profiles are not
block-diagonal; their small row-blocks have full traces (fatpoint Thm 6.4 discussion), so projections onto row
blocks are uncut cubes and give `nu = 2|X_b|`. Sections (sub-boxes) go the other way: `B subset A` gives
`TC(A) subset TC(B)`, i.e. obstructions, never constructions. **So neither products nor sections can produce barrier
points.**

### 3.2 Cyclic symmetry: near-zeros on isotypic circles are all-or-nothing

**Theorem 3.2 (proved).** Let `g` be a permutation of the non-base rows of order `N` that maps the pattern set to
itself (so `gA = A`, and `g` preserves `V(A,m)`, `V_H`, `|c|_1`). Let `p, p' in R^d` and `a in Z` with `gcd(a, N) = 1`
such that `u_theta := cos(theta) p + sin(theta) p'` satisfies `g^{-1} u_theta = u_{theta + 2 pi a/N}` for all `theta`.
If `2m < N`, then for every `c in V(A,m)` and every `theta`, with `f_c(theta) := Phi_c(u_theta)`,

```text
max_{0 <= s < N} |Phi_{g^s c}(u_theta)|  >=  ||f_c||_{L^2(d theta/2pi)}  >=  (m+1)^{-1/2} ||f_c||_inf .
```

Consequently, for `kappa = kappa_H` and for `kappa = kappa_V`,

```text
inf_theta kappa(u_theta)  >=  (m+1)^{-1/2} sup_theta kappa(u_theta) ,
```

and the circle `C = {u_theta}` lies entirely in `TC_H(A,m)` or is disjoint from it.

*Proof.* `Phi_{g^s c}(u) = Phi_c(g^{-s} u)` (permutations are orthogonal and `(g^s c)_k = c_{g^{-s} k}`), so
`Phi_{g^s c}(u_theta) = f_c(theta + 2 pi a s/N)`. `<k', u_theta> = alpha_k cos theta + beta_k sin theta`, so `f_c` is a
trigonometric polynomial with frequencies in `{-m, -m+2, ..., m}` (`m+1` of them) and `|f_c|^2` has frequencies in
`[-2m, 2m]`. For `0 < |l| <= 2m < N` and `gcd(a,N) = 1`, `sum_{s<N} e^{i l (theta + 2 pi a s/N)} = 0`; hence the
average of `|f_c|^2` over the `N` equispaced points equals its mean `||f_c||_2^2`, and the maximum is at least the
mean. `||f_c||_inf <= sum |hat f_c(l)| <= (m+1)^{1/2} ||f_c||_2` by Cauchy-Schwarz. For the kappa statement use
`|g^s c|_1 = |c|_1`, `g^s c in V_H` (resp. `V(A,m)`), and take sups over `c`. QED

*Normalisation.* `kappa` is computed at `hat u = u/|u|_inf`; on the circles below `|u_theta|_inf in [cos(pi/N), 1]`, so
the ratio bound holds for `hat u_theta` up to a factor `cos(pi/N)^{-m} <= e^{pi^2/(4N)}` (`m < N/2`).

**Corollary 3.3 (Galois; proved).** In the setting of Theorem 3.2 with `p + i p' = psi_a := (zeta^{a t(x)})_x` for a
primitive `N`-th root of unity `zeta` and a labelling `t : rows -> Z/N` on which `g` acts as `t -> t + 1` (one `g`-orbit
of size `N` on the support, zero elsewhere): if `C_a ⊂ TC_H(A,m)` then `C_b ⊂ TC_H(A,m)` for every `b` with
`gcd(b, N) = gcd(a, N)`. Let `K_H subset Sym^m` be the apolar annihilator of `{Phi_c : c in V_H}` (so
`u in TC_H` iff `l_u^m in K_H`, and `K_H = J(A)_m` when `V_H = V(A,m)`). Then `C_a ⊂ TC_H` iff every monomial
`l_{psi_a}^i l_{psi_{-a}}^{m-i}` lies in `K_H ⊗ C`; in particular it implies `psi_a in TC_H(A,m)(C)`, i.e. (when
`V_H = V(A,m)`) the function `k -> (sum_x zeta^{a t(x)} k_x)^m` coincides on `A` with a polynomial of degree `< m`.

*Proof.* `TC_H` is a cone defined over `Q` (the forms `Phi_c`, `c` integral). `C_a ⊂ TC_H` gives the real plane
`W_a ⊂ TC_H`, hence its Zariski closure `W_a ⊗ C = span_C(psi_a, psi_{-a})`. `sigma_s in Gal(Q(zeta)/Q)` maps
`psi_a -> psi_{sa}` and preserves `TC_H`; `{sa}` runs over all `b` with `gcd(b,N) = gcd(a,N)`. For the second claim:
`g` acts on `l_{psi_a}^i l_{psi_{-a}}^{m-i}` by the scalar `zeta^{a(m - 2i)}`, distinct for distinct `i` since `2m < N`;
`J_H ⊗ C` is `g`-stable, so it contains `(c_1 l_{psi_a} + c_2 l_{psi_{-a}})^m` for all `c_1, c_2` iff it contains each
eigencomponent (`K_H` is `g`-stable because `V_H` is). The last sentence is the duality of §1. QED

**Where it applies (exact checks).**
* `sdc_q`: translations `x -> x + 1` of `F_q` (base `inf` fixed) preserve the pattern set (verified for `q = 19, 23`
  in this round's session; for all `q` by construction, `S'_a = S_a xor {inf, a}` is translation-equivariant).
  `N = q`, circles `u_theta(x) = cos(theta + 2 pi a x/q)`, all `a in F_q^*` in one Galois class. Regime: `2m' < q`
  holds since `m' = O(sqrt(n log n))`.
* `syl_k` (actual-type): a Singer cycle `x -> gamma x` of `F_{2^k}^*` (base `0` fixed) maps hyperplane patterns to
  hyperplane patterns; `N = d = 2^k - 1`; circles `u_theta(gamma^t) = cos(theta + 2 pi b t/d)`.
* `pal_q` is **not** translation-invariant (the twin of row `inf` breaks it; verified `q = 19, 23`), so Theorem 3.2
  applies to it only through subgroups that fix that twin (none of order `> 2m'` were found; not pursued).

### 3.3 The circles pass the Diophantine test automatically

**Lemma 3.4 (proved, under the stated hypothesis).** Let `q` be prime, data `sdc_q` (`mu = 1`, `r = 2q`), and assume
**H_light(c_0)**: every relevant integral `lambda` that is not a multiple of a pair `e_x - e_y` and is light in the sense
of §2.2 has `sigma_lambda >= c_0 r`. Then there is `tau_D = O_{c_0}(log n + log(1/C))` such that for every `tau >= tau_D`
and every `a in F_q^*`, a set of `theta` of measure `>= pi` (out of `2 pi`) makes `u_theta / |u_theta|_inf` satisfy (2.2').

*Proof.* Write `z = (0, u_theta)` (base first). Since `|u_theta|_inf <= 1` and `s* <= 2`, it suffices that
`|lambda.z| >= max(26 e^{-tau sigma_lambda/2}/C, 4 eta ||lambda||_1)` for every light `lambda`; we bound the `theta`-measure
of the failure set of each term.
* *Pairs `x, y in F_q`*: `z_x - z_y = -2 sin(theta + pi a(x+y)/q) sin(pi a(x-y)/q)`, and `|sin(pi a(x-y)/q)| >= 2/q`. The
  set of `theta` with `|sin(theta + phi)| < epsilon` has measure `<= 2 pi epsilon` (`epsilon <= 1/2`). Base pairs:
  `z_x = cos(theta + 2 pi a x/q)`, same bound. Multiples `l(e_x - e_y)` have `sigma >= sigma_xy` (`|l lambda(S)| >= |lambda(S)|`)
  and the same hyperplane, so they are implied. With `sigma_xy = 1` (exact, `sigma_pairs.py` for `q <= 43`; in general
  `sigma_xy >= rho_1 ||e_x - e_y||_1 = 1`), the total pair measure is `<= 2 pi (binom(n,2) + n) (13 q/C) e^{-tau/2}`.
* *Other light `lambda`*: `lambda.z = Re(e^{i theta} beta_lambda)` with `beta_lambda = sum_{x in F_q} lambda_x zeta^{ax}` in
  `Z[zeta_q]`. `beta_lambda = 0` iff `lambda|_{F_q}` is constant (`q` prime), and then `lambda = c(1_{F_q} - q e_inf)` with
  `sigma_lambda >= sum_S (|c| |S cap F_q| - 1/2) >= r(|c|(q-1)/2 - 1/2)`: heavy. Otherwise `|N(beta_lambda)| >= 1` and every
  conjugate has modulus `<= ||lambda||_1`, so `|beta_lambda| >= ||lambda||_1^{-(q-2)}`; since `|Re(e^{i theta} beta)| < epsilon`
  on a set of measure `<= 2 pi epsilon/|beta|`, the two bad sets have measure `<= 2 pi ||lambda||_1^{q-2} (26 e^{-tau sigma_lambda/2}/C
  + 4 eta ||lambda||_1)`. Light `lambda` have `||lambda||_1 <= sigma*/rho_1 = 2 sigma* <= 2r + 2`, and there are at most
  `(4r+5)^n` of them. By H_light, `e^{-tau sigma_lambda/2} <= e^{-tau c_0 r/2}`; and `eta ||lambda||_1 <= (2r+2) e^{-r tau/4}`.
  So this part is `<= 2 pi (4r+5)^n (2r+2)^{q-1} (26 e^{-tau c_0 r/2}/C + 4 e^{-r tau/4})`.
* For pairs the `eta`-term is implied by the first term (`e^{-tau/2} >> n e^{-W/4}`).

With `r = 2q`, `n = q + 1`, both parts are `<= pi/2` once

```text
tau >= tau_D := max( (2/c_0)(2 ln(4r+5) + ln(1000/C)) ,  4 ln(4r+5) + 8 ,  2 ln(1000 q n^2/C) ) .   QED
```

*Evidence for H_light* (`hlight.py`, exact): for every integral `lambda` with `||lambda||_1 = 4` not supported on two
points, `min sigma_lambda / r = 0.214, 0.227, 0.184, 0.196, 0.210` for `q = 7, 11, 19, 23, 31` (type `e+e-e-e`;
type `2e-e-e` has `>= 0.53 r`, i.e. it is heavy). Six-term characters: sample means `0.45..0.52 r` (`sigma_pairs.py`),
minimum not computed. For Paley-type patterns Weil's bound predicts `sigma_lambda = r(E|lambda(S)| - 1/2) + O(||lambda||_1 r/sqrt q)`
with `E|lambda(S)| >= 3/4` for every balanced non-pair `lambda` with `||lambda||_1 <= 6` [M]; H_light is not proved here.

**Theorem 3.5 (sharp barrier test on the circles; proved modulo H_light).** Let the regime parameters be those of
the `sdc_q` data (`r = 2n - 2`, `w = sqrt(8 d ln 4r)`, `m' = ceil(4w)`, `2m' < q`), let H_light(`c_0`) hold, and let
`tau = O(log n)` be at least `tau_D` and the thresholds of Theorem 2.3. Then:
* (i) if some point of a circle `C_a` has `kappa_H(hat u) <= epsilon`, every point of `C_a` has
  `kappa_H(hat u) <= sqrt(m'+1) e^{pi^2/(4q)} epsilon` (Theorem 3.2);
* (ii) exact version: some point of some `C_a` is in `TC_H` iff all of them are, for all `a in F_q^*` at once, iff
  `l_{psi_1}^i l_{psi_{-1}}^{m'-i} in K_H ⊗ C` for all `i` (Corollary 3.3).

**Consequence.** If some point of some circle `C_a` has `kappa_H <= eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)})` (a loss of
`O(log n)` in `W`), then by (i) every point of `C_a` has `kappa_H <= eps_0/2`, Lemma 3.4 puts half of `C_a` in the
light-character gap region, and Theorem 2.3 produces gap-admissible `sdc_q` fakes at the barrier scale violating
MIN^pair (and MIN^char): **route 5 would be capped.** Conversely, route 5 can only survive the circle test if
`kappa_H` is bounded below on `C_a` by `eps_0/(2 sqrt(m'+1))`, uniformly (Theorem 3.2 makes "somewhere" and
"everywhere" the same).

*Proof.* (i) is Theorem 3.2 with the normalisation remark, (ii) is Theorem 3.2 plus Corollary 3.3. Lemma 3.4
supplies (2.2') for half of each circle, Theorem 2.3 the fakes, and the fakes satisfy `kappa_H(hat u) <= eps_0`, while
their pair data define a gap region containing `hat u` (they are gap-admissible). QED

**What this decides and what it does not.** It turns the barrier question for the most natural *irrational*
candidates of the `sdc` boxes into **one explicit algebraic question**: is `l_{psi_1}^{m'}` (and its mixed companions)
in the leading-form ideal of the regime `sdc_q` box, equivalently is the `m'`-th power of the first Fourier coefficient
`hat k(1)^{m'}` interpolable on `A` by polynomials of degree `< m'`? The Diophantine side is automatic (given H_light).
I could not decide it: lower bounds for `nu_A(psi_1)` beyond `2w` need non-zonotopal certificates (record), and upper
bounds need an identity in `J(A)` along a direction with no integral structure. *Heuristic* (not used): `<k', psi_1>`
takes `|A|` distinct values on `A` (only heavy relations), the opposite of the few-value structure that makes every known
tangent-cone direction cheap; so one expects `nu_A(psi_1) >= m'`, i.e. the circles are **not** in `TC`, and by Theorem
3.2 then `kappa_H >= (m'+1)^{-1/2} max_{C_a} kappa_H > 0` uniformly on each circle — but no quantitative lower bound at
the `eps_0` scale follows.

**Same test for the actual-type Sylvester boxes.** With the Singer cycle, Theorem 3.2 and Corollary 3.3 hold verbatim
(`N = d`). If some Singer circle lay in `TC_H(A_syl, m')`, MIN^pair would fail for Sylvester-type data at every `W`
(the circle points have pair gaps `>= c/n^2` for suitable `theta`), without any fake being involved; for those profiles
walsh-extraction Theorem B already bounds `M` unconditionally (README-13 J9), so this would only remove a hypothesis.

## 4. Explicit near-zero constructions (product / section / symmetry)

All entries at `T = 1`, `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)` (`m'/(2j) ~ 4`), regime `d >= 316 ln 4r`,
barrier scale `tau = K ln n`, fake data `sdc_q` / `pal_q` (bounded pair margins), unless stated otherwise.

| # | candidate source of tangent-cone points / near-zeros | collision / character structure | barrier verdict | where |
| --- | --- | --- | --- | --- |
| 1 | width-cheap spans (`<= 4` long directions, `<= 16` classes, `<= (m'-1)/2` points) | collisions | **excluded** (Prop 2.4) | t1rigour Thm 6.3; §2.3 |
| 2 | projection-cheap blocks: Fano / XOR blocks `Y_t` (`t <= 4`), parity sections, congruence families | `<= 2^t` classes, collisions | **excluded** (Prop 2.4); they are Sylvester structures anyway (actual-type) | refJ, loctheory; §2.3 |
| 3 | J5-type perturbations of coordinate spans | need concentrated in-class pair excess | **excluded for fakes** (uniform `sigma_xy = 1` in `sdc_q`); finite-`n` failures of coll-type hypotheses for actual-type data only | refJ §1; Prop 2.4 |
| 4 | products of blocks, sums of coset blocks | additive `nu` | **nothing new**; `>= 5` blocks of size `>= j` admit no base-collision-free TC point | Prop 3.1 |
| 5 | sections / sub-boxes | give `TC(A) ⊂ TC(B)` | obstructions only, never constructions | §3.1 |
| 6 | linear joins and base cycles of the fat-point (width) system | inside spans of `<= O(m'/j)` integer directions | collisions as in 1-2 | refJ Prop 2.3, fatpoint §2 |
| 7 | non-linear base curves (rational normal curves) through column / coordinate fat points | low degree: inside the span of `<= 9` columns; through all coordinate points: impossible | **excluded** in cases (a), (b); high-degree curves through many column points excluded only by a dimension count (c) | remark below |
| 8 | Fourier circles `C_a` of `sdc_q` (translation symmetry), Singer circles of `syl_k` | no light relation; only the heavy `sum_{F_q} u = 0` | **decisive and open**: all-or-nothing (Thm 3.2), Galois-uniform (Cor 3.3), Diophantine automatic (Lemma 3.4, H_light); caps route 5 iff the circle is a `eps_0/sqrt(m'+1)`-near-zero (Thm 3.5) | §3.2-3.3 |
| 9 | collision-free rational strata (relations only non-pair characters, `sigma >= c_0 r`) | e.g. 4-term relations `u_x + u_y = u_z + u_w` | **fatal if they exist** with `o >= 3` (Prop 2.5); none known in any fake box; toy analogue (Sylvester `M8`, quadruple hyperplanes) tested in §6 | §2.3, §6 |
| 10 | Kronecker / sum-form directions `u(x,y) = a(x) + b(y)` on product-indexed profiles | 4-term relations `u(x,y) - u(x,y') - u(x',y) + u(x',y') = 0` with `sigma ~ 0` (light like pairs) for Kronecker columns | **excluded** by (2.2') (exact relation on a light character) | remark below |
| 11 | non-Lojasiewicz near-zeros far from all real zeros | — | **open**; no construction known; nearstrata Prop 2.7 shows general Lojasiewicz bounds cannot exclude them | — |

*Remark (row 7).* (a) *Proved.* A degree-`m'` form with multiplicity `>= m' - 2j` at `N` points of a rational normal curve
of degree `e` (in its span `P^e`) restricts to a binary form of degree `m' e` vanishing to total order `>= N(m' - 2j)`; so the
curve is in the base locus once `N(m' - 2j) > m' e`. If the curve is determined by `N = e + 3` column points in general
position, this needs `3m' > 2j(e+3)`, i.e. `e <= 8` for `m'/(2j) <= 4 + O(1/w)`; such a curve lies in the span of `<= 9`
columns, whose points take `<= 2^9` values: collisions. (b) *Proved.* A rational normal curve of degree `d - 1` through
all coordinate points `e_x` has the form `t -> (c_x/(t - a_x))_x` (distinct `a_x`); it meets two `+-1` vectors `beta, beta'` at
`t_1, t_2` only if `(t_2 - a_x)/(t_1 - a_x) = +-mu` for every `x`, which allows at most two values of `a_x`: impossible for
`d >= 3`. (c) *Not proved.* Curves of degree `e` through `N` column points and no coordinate points are base curves once
`e < N(m' - 2j)/m' ~ 3N/4`, while the naive parameter count for rational curves through `N` points of `P^{d-1}` needs
`e >~ N(1 - 2/d)`. So generic position forbids them; special (Paley, Sylvester) position is not excluded.

*Remark (row 10).* For Kronecker columns `beta_{(p,p')}(x,y) = beta_p(x) beta'_{p'}(y)`, the 4-term `lambda` above has
`lambda(S) in {0, +-2}` with `|lambda(S)| = 2` exactly when both factors separate, frequency `1/4`, so `sigma_lambda = r(2/4 - 1/2) = 0`:
as light as a pair. A sum-form `u` has `lambda.u = 0` exactly, violating (2.2'). Product-indexed profiles therefore
cannot use sum-form directions; product-form directions `a(x)b(y)` have no reason to be cheap (no factorisation of
`<k,u>` through a cut image).

**Answer to (1) as far as it goes.** At regime ratios and growing `n`, every tangent-cone point or near-zero that any
known mechanism produces in the fake pattern boxes is either (i) collision-bearing, hence excluded at the barrier
scale by the fakes' own uniform pair margins (Prop 2.4; strictly more than the record's separation radius, which is
for arbitrary pair data), or (ii) one of the two genuinely open, explicitly identified candidate classes: the
symmetric Fourier circles (row 8; a single algebraic question by Thm 3.2/Cor 3.3/Thm 3.5) and collision-free
rational strata with only non-pair relations (row 9; fatal by Prop 2.5). I could construct no nondegenerate
near-zero and could prove for neither class that it is empty.

## 5. Separation from collision strata (or not), and what AC must use

**What is proved about separation.**
1. For every cluster and every fake (any pair data): `dist(hat u, L) >= (1/2) e^{-N_2 abar/(2P(L))}` (record).
2. For fake data with bounded pair margins (`sdc_q`: `sigma_xy = 1` for all pairs), the fakes' own (2.2') gives the much
   larger radius `(13/C) e^{-tau/2}/sqrt 2` around **every** collision subspace, and Lemma 2.0 then cannot certify
   `kappa <= eps_0` near any collision stratum (Prop 2.4). So, for the barrier question, the collision strata are
   *never* the issue; what is at issue is whether near-zeros exist **away** from them.
3. Near-zeros away from all collision strata exist in fake boxes if and only if (sufficient direction via Thm 2.3)
   one of the following exists in `TC_H` or near it: a collision-free rational stratum with heavy / `c_0 r`-light
   relations and order `o > 1/(2c_0)` (Prop 2.5), a whole Fourier circle (Thm 3.5), or an irrational real component /
   non-Lojasiewicz near-zero of a kind not yet identified.

**Was "every near-zero is within the separation radius of a collision stratum" proved?** No. It is not even the
right statement for the barrier: by Prop 2.5 the correct target is "every near-zero of `kappa_H` in a fake box lies
near a subspace cut out by **light** characters with `o sigma < r/2`", and pairs are only the lightest such characters.

**The structural property that would rule barrier points out (the content AC must use).** Every mechanism in §4 that
produces tangent-cone points or near-zeros produces them along directions `u` for which the projection `k -> <k',u>`
has **few values on `A`** relative to `m'`: integral combinations of `<= O(m'/j)` cheap integer directions (columns, Fano /
XOR / parity / congruence families, coordinates). Few values force (pigeonhole) many coincidences among the `u_x`,
i.e. light relations. The two open classes are exactly the places where this "few values ⇒ light relations"
implication is not known to hold:
* (P1) **Value-count rigidity on rational strata:** if `L ⊂ TC_H(A,m')` is rational, then `Lambda_L` contains a pair
  (or a character with `o sigma < r/2`). Equivalently: a collision-free integral family `F` (rows of `F` pairwise
  distinct) has `nu_A(span F) >= m'`. For column families of `sdc/pal` boxes with `|F| <= (1/2 - eps) log_2 n` this is
  what full-Venn images (Weil) predict; beyond, nothing is known.
* (P2) **No cheapness without integrality:** a direction `u` whose projection `<k',u>` takes `|A|` distinct values
  on `A` (no light relation) has `nu_A(u) >= m'`, quantitatively `kappa_H(u) >= e^{-o(W)}` uniformly on the
  light-character gap region. For the Fourier circles (P2) is exactly the open question of Thm 3.5.

(P1) is a statement about the box's integral structure (a strengthening of ZE-TC^rat with the light-character
refinement). (P2) is the genuinely anti-concentration statement: it is AC restricted to directions with no light
relation, and by Prop 2.1 a proof of AC would give it for **all** half-cubes simultaneously. Neither (P1) nor (P2)
follows from anything in the literature I could locate (§7); both are consistent with all exact data (§6).

**A caution on the AC formulation itself.** AC as stated (Prop 2.1) quantifies over all half-cubes `S`, so its
failure set contains every `eps_0`-near tangent-cone point of every half-cube. Half-cubes are much more flexible than
slab boxes (any integer direction `v` can be cut at its median width `~0.55 ||v||_2`, so integer directions with
`||v||_2 <~ m'/1.1 ~ 3.6 w` are cheap in some half-cube). The resulting failure points are again collision-rich
(integer vectors with `||v||_2 = O(sqrt(d log r))` take `O(sqrt(log r))` typical values per coordinate), consistent with AC on
the gap region, but AC must survive this adversarial choice of `S`, which the slab-box statement MIN^pair does not
have to. This is the skeptic's "strict strengthening", made concrete; it is evidence-level, not a refutation of AC.

## 6. Consistency of AC with the exact data (M8, M8', pal11, sdc7/sdc11)

All toy models are `T = 1` with ad hoc `L`, LP-subcritical, far below the regime (`n <= 12` vs `n >= 2900`).
Densities `|A|/3^d` from `round12/fake_checks/sizes_out.txt` (exact counts):

| model | `(L, m)` | `|A|` | density | `dim V` | real tangent cone (record + this round) | AC / MIN^pair reading |
| --- | --- | --- | --- | --- | --- | --- |
| `M8'` = Syl(8) | (3, 11) | 1877 | **0.858** | 7 | every complex zero on a collision hyperplane (R11 Thm 5.2(a), exact) | the only toy with density `>= 1/2`: consistent with AC (exact-zero form) on collision-free `u` |
| `M8` = Syl(8) | (2, 8) | 1037 | 0.474 | 42 | zeros on collision hyperplanes or on the 4-term hyperplanes of additive quadruples (J4, exact); **this round**: 300 + 42x25 high-precision searches, incl. one restricted to each of the 42 quadruple hyperplanes, found **no** zero with min gap `>= 1e-3`: all slid to strata (evidence) | consistent; for Sylvester, quadruples have `sigma = 0` (light like pairs), so even collision-free zeros there would be excluded by (2.2') |
| `sdc7` | (2, 8) | 757 | 0.346 | 1 (`H_A(8) = 1`) | certified nondegenerate zeros, `sep >= 0.2496` (fake Prop 3.1) | AC_A false there (exact zero); AC_T not applicable (density `< 1/2`). Cause: **underdetermination** (one octic). this round: the single octic is translation-invariant (`delta_0 = 1`), hence constant on each Fourier circle, and nonzero there: `C_1` is **not** in `TC(A,8)` (`sdc_circles.py`); also `psi_1 notin TC(A,i)` for `i <= 6` (`sdc_fourier.py`) |
| `sdc11` | (2, 8) | 9373 | 0.053 | 5 | certified regular collision-free zeros, `sep >= 0.106` (referee_fake) | AC_A false there; cause: underdetermination (5 octics in `P^10`). **This round**: all 5 tangent forms are translation-invariant (character deficiencies `delta_0 = 5`, others 0), so each `Phi_c` is **constant** on every Fourier circle, and the constant (`|hat k|^8` component) is nonzero: circles `C_1`, `C_2` are **not** in `TC` (`sdc_circles.py`, two primes) |
| `pal11` | (2, 8) | 12093 | 0.068 | `<= 436` (mod `p`) | only unreliable floating-point evidence (referee_fake) | no data; `pal_q` has no translation symmetry, so §3.2 does not apply |
| `sdc19` | (2, `m = 5, 6`) | 264291 | 2.3e-4 | — | — | circle test with `2m < 19` (Thm 3.2 applies): **no** `F_i` (`i = 0..m`) lies in `Poly_{<m}`, so `C_1 ∩ TC(A,m) = ∅` for `m = 5, 6` (exact mod two primes; the `m = 7` run was stopped for time) |

**Explanation, item by item.**
1. **sdc7, sdc11 (zeros exist).** These are *exactly* the cases fake Thm 2.1 / 2.2 predict: an underdetermined system
   (`H_A(m) <= 5 < d - 1`) at `m = rho(A)` has a real variety of positive dimension through nondegenerate points. AC
   (in any form) must fail there, and it does; nothing transfers to the regime, where `H_A(m')` exceeds `d` by a factor
   `e^{Theta(m' log(d/m'))}`. The circle computations show that even in these underdetermined toys the symmetric
   Fourier circles avoid `TC`: the nondegenerate zeros there are not symmetric ones.
2. **M8, M8' (no nondegenerate zeros).** Overdetermined; every real zero lies on a rational hyperplane cut out by a
   light character (pairs; for `M8` also additive quadruples, which are light for Sylvester data). This is the toy form
   of (P1): integrality forces light relations. No collision-free zero was found even on the quadruple hyperplanes.
3. **pal11.** Undecided; the only faithful overdetermined *fake* model, not touched here (no symmetry to exploit).
4. **Densities.** Only `M8'` reaches density `1/2`. So no toy model tests AC in its stated form (Prop 2.1); in the
   regime the density is `>= 1/2` by t1rigour, and AC's premise holds.

**Is AC consistent with all exact data?** Yes, in the only sense available: every exact zero that violates the
box-level statement (sdc7, sdc11) sits in an underdetermined toy where the density premise of AC fails and the order
is `m = rho`; every overdetermined toy has only light-relation zeros. No datum distinguishes AC from its negation in the
regime.

## 7. Literature: what anti-concentration theory says here

Located this round [L] (abstract / snippet level):
* **Meka-Nguyen-Vu**, *Anti-concentration for polynomials of independent random variables*, Theory of Computing 12
  (2016), art. 11 (doi 10.4086/toc.2016.v012a011). Degree-`m` Littlewood-Offord bounds in terms of a rank parameter;
  as steelman S3 computed, vacuous at `m' ~ sqrt(d log r)` (the rank is `<= d/m'` and the loss is `(log)^{O(m log m)}`).
  They give **upper** bounds on concentration, so they speak to AC's positive side only; nothing on the barrier side.
* **Kane**, *A structure theorem for poorly anticoncentrated polynomials of Gaussians*, Ann. Probab. 45(3) (2017)
  1612-1679 (arXiv 1204.0543): a poorly anticoncentrated degree-`m` polynomial is close to a function of `O_m(1)`
  normalised polynomials. This is the right *shape* for the barrier side ("concentration ⇒ structure"), and its
  discrete analogue would be exactly (P1)/(P2) of §5: concentration at `0` on half the cube ⇒ `<k',u>` factors through a
  few cheap (integral) directions ⇒ light relations. But it is Gaussian, and its bounds depend on `m` in an
  uncontrolled (effective but huge) way, so it is unusable at `m' -> infinity`.
* **Small sums of roots of unity** (Myerson's problem; Barber arXiv 2104.15057, Kalogirou arXiv 2607.06098, and a 2016
  Diophantine-approximation paper arXiv 1608.04547): polynomial lower bounds `f(k,N) >> N^{-c_k}` are known only for
  `k <= 4`; for `k >= 5` only the Galois-norm bound `k^{-N}` (and Konyagin-Lev `k^{-N/4}` for prime `N`, distinct roots).
  Relevance: Lemma 3.4 needs only the Galois-norm bound because non-pair light characters have `sigma ~ c_0 r`
  (exponentially weak Diophantine requirement); polynomial bounds would be needed only if some non-pair character
  had `sigma = O(1)`, which H_light excludes.

Quoted from memory [M], not re-checked this round: Carbery-Wright (2001); Mossel-O'Donnell-Oleszkiewicz invariance
(2010); Costello-Tao-Vu quadratic Littlewood-Offord; Kwan-Sauermann (quadratic inverse Littlewood-Offord);
Fox-Kwan-Sauermann (combinatorial anti-concentration); Razborov-Viola; Tao-Vu / Nguyen-Vu inverse Littlewood-Offord
for linear forms; Alon-Furedi / combinatorial Nullstellensatz. What they give here:
* **Linear inverse LO (Nguyen-Vu)** explains the *cheap directions*: if `<k,v>` has an atom `>= 1/(2(m'+1))` on half
  the cube, most coordinates of `v` lie in a small GAP. This is the "few values ⇒ additive structure" step of (P1), but
  only for a **single** direction and only for half-cubes, not for `nu_A`.
* **Quadratic inverse LO (Kwan-Sauermann, [M])**: structure for degree 2 only. Route 5 needs degree `m' ~ 11 sqrt(d ln 4r)`.
* **Invariance principles / regularity lemmas** lose `m' tau^{1/(8m')} >= 1` (steelman), and regularity decompositions
  for low-degree polynomials on the hypercube have bounds tower-type in the degree: vacuous at growing degree.
* **Alon-Furedi (min support of a degree-`m` function on `{-1,0,1}^s` is `>= 3^{s - m/2}`)** gives the exact-zero
  version of AC only against sets missing `>= 3^{d - m'/2}` points of the cube; slab boxes miss far more
  (`~ 3^d (4r)^{-6}`), so it cannot see the slab cuts (checked in §4-type arithmetic: needs `s < 0.34 j` but `s >= 4j`).

**Verdict on the literature.** Nothing located gives either (P1) or (P2). The barrier side has no literature analogue at
all: the known "structure from concentration" theorems (Tao-Vu, Nguyen-Vu, Kane, Kwan-Sauermann) are all of the
form "concentrated ⇒ structured", i.e. they would *exclude* barrier points if they existed at growing degree on the
discrete cube; none constructs concentrated polynomials, and the only constructions of concentrated degree-`m'`
polynomials with prescribed leading form on half-cubes are the integral (few-value) ones of §4, rows 1-7.

## 8. What AC must use / open points

**Growth / uniformity: unchanged.** Nothing here improves `M <= (2/3 + eps) log R/loglog R` or touches the uniform
bound. Route 5 is **not capped** (no barrier point constructed) and **not unblocked** (neither (P1) nor (P2) proved).

**Open points, in order of decisiveness.**
1. **The circle question (sharpest single test).** For the regime `sdc_q` box (`r = 2n - 2`, `w = sqrt(8 d ln 4r)`,
   `m' = ceil(4w)`), is `inf_{theta} kappa_H(u_theta) <= eps_0/(2 sqrt(m'+1))` on the Fourier circle `C_1`? Exactly:
   is `hat k(1)^i conj(hat k(1))^{m'-i}` in `Poly_{<m'}|_A ⊗ C` for every `i` (Cor 3.3)? A yes caps route 5 (Thm 3.5,
   modulo H_light). Toy evidence (sdc7, sdc11 at `m = 8`; sdc19 at `m = 5, 6`): **no**, and
   in sdc7/sdc11 the tangent forms are translation-invariant, hence constant and nonzero on the circles.
2. **Collision-free rational strata (P1).** Does any regime fake box have a rational `L ⊂ TC_H` with no pair in
   `Lambda_L` and transversal order `o >= 3`? A yes caps route 5 (Prop 2.5). Toy evidence (`M8`): no.
3. **H_light.** Prove `sigma_lambda >= c_0 r` for non-pair light characters of `sdc_q` / `pal_q` (exact for `||lambda||_1 = 4`,
   `q <= 31`: `c_0 = 0.18`; Weil should give it for `||lambda||_1 = O(1)`).
4. **(P2) in any non-trivial instance**: a lower bound `kappa_H(u) >= e^{-o(W)}` for one explicit direction with no light
   relation (e.g. `C_1`) on one regime box. Even this is beyond current tools: all known certificates are zonotopes,
   capped at `zeta <= 2w < m'`.
5. **pal_q**: no symmetry; the circle reduction does not apply. Its overdetermined toy `pal11` remains undecided.

**Recommendation.** If route 5 is resumed, the two concrete questions are 1 and 2; both are finite algebraic questions
about one explicit box and one explicit family of subspaces, and a "yes" to either caps the route. A "no" to both would
leave only non-symmetric irrational components and non-Lojasiewicz near-zeros, for which nothing is known. The
anti-concentration literature (§7) points the right way ("concentration ⇒ integral structure") but is Gaussian or
fixed-degree; a discrete, growing-degree analogue of Kane's structure theorem on slab-cut ternary cubes is what (P1)+(P2)
amount to.

## 9. Files

* Write-up: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/barrier.md`
* Notes: `.../round16/barrier_notes.md`
* Scripts and outputs (`.../round16/barrier_checks/`):
  * `m8_tc.py`, `m8_refine.py` (+ `_out`): M8 tangent cone, exact `V(A,8)` basis (dim 42), floating search + 60-digit
    Gauss-Newton refinement (evidence).
  * `m8_refine_quad.py`, `run_quads.sh`, `quads.txt`, `m8_quads_L2_m8_out.txt`: the same restricted to each of the 42
    additive-quadruple hyperplanes (evidence).
  * `sigma_pairs.py`, `hlight.py` (+ `_out`): exact pair slacks (`sigma_xy = 1` for all pairs of `sdc_q`, `q <= 43`) and
    exact minimum 4-term slacks (H_light evidence).
  * `sdc_fourier.py`, `sdc_circles.py` (+ `_out`): Fourier-direction and circle tests in `sdc7`, `sdc11`, `sdc19`
    (exact mod two primes `p = 1 mod q`; char-0 conclusion is evidence unless `rank_p = rank_Q`).
  * `flib.py`: copy of `round12/fake_checks/flib.py`.
