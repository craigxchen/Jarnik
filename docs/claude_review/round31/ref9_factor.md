# Referee report: Part IX.3, the factored residue construction (statement through the paragraph on pairs and offsets)

Referee block "factor". This covers the text of `astra_free_IX/answer.md` from the statement of IX.3 (line 2034)
through the paragraph that ends "These estimates are uniform in the cutoff" (line 2069), with the same-level
deletion at parent primes checked first, as Astra asked. The non-pair paragraph ("For nonpairs ...", (IX.7)), the
e-all failure paragraph and IX.4 are outside this block. Where this block uses them I say so.

Scripts (my own code, exact integer and Gaussian-integer arithmetic; Astra's scripts were read, not imported):
`ref9_factor/linalg_check.py` (output `linalg_out.txt`) and `ref9_factor/residue_check.py` (output
`residue_out.txt`), both in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9_factor/`.

## 0. Verdict

| item | claim | verdict |
|---|---|---|
| 1 | What is verified is exactly P3-char (sharp.md Def. 1.1 (a)-(d), (e-char) with the same-level clause), at every cutoff `X`, on one fixed model | **Correct.** The model, its windings and its residue data do not depend on `X`. Every collision modulus at every cutoff divides one global Gaussian integer. For pairs, Astra checks more than the record requires (all four units at the parent primes). (a)-(c) and the non-pair characters rest on paragraphs outside this block. |
| 2 | Paley exponent matrix, coverage of both endpoints, odd exponent only at `j = i` | **Proved and computed.** The exact counts are `#{a = 2} = (q-3)/4` and `#{a = 0} = (q-1+2chi(j))/4` on the full class. |
| 3 | `D = H mod 2`; `D` has an integral right inverse (all Smith invariants are 1) | **Proved and computed** (Smith normal form for `q <= 71`, column Hermite form for `q <= 103`, random `S`) |
| 4 | `V = U(H-D)/2` is integral and `DF = H` | **Proved** (one line); **computed** |
| 5 | Witnesses (IX.6): norms, `r/conj r`, mark ratio `gamma_i/gamma_i'`, unit bookkeeping | **Proved and computed.** Even winding differences are needed and are allowed. The local character products do not depend on the choice of `U`. |
| 6 | Same-level deletion at parent primes | **Proved and computed.** The record asks nothing at `p_l` for pairs at different levels or for non-pair characters. One phrase is vacuous: no parent prime lies below `2m` here. |
| 7 | Global collision cost `D_{u0} <= sqrt2 L_0^{n/2}`; `D_1 <= L_0 |beta_i - beta_i'|` | **Proved and computed** |
| 8 | Primes in one narrow sector and radial bin; (IX.5); offsets; pair inequality; uniform in the cutoff | **Proved.** I recomputed the constants: the unit-one pair inequality holds with a factor of about 35 to spare. Minor correction: `m/omega(N') -> 1/2` holds only along the diagonal `rho -> 1`. |
| 9 | Small instances reproduced | **Computed:** 0 failures in 5,712 exterior pair identities, 31,200 exterior random characters and 4,448 same-level identities. The negative control (odd winding difference) fails as predicted. |

**Most serious issue.** There is no error in this block. The points that need attention are wording and scope:
* (i) The same-level clause of the record is the unit-one clause of `lower.md` Def. 1.1(e). Astra's "every Gaussian unit" is more than P3-char asks, and her proof covers it too.
* (ii) "parent primes below `2m`": none occur in this construction, because every `p_j >= q^d > 2m`.
* (iii) `m/omega(N') -> 1/2` needs `rho -> 1` along a diagonal sequence (the convention of README_VIII §1.3). At fixed `rho` the ratio tends to `rho/2`.
* (iv) The witnesses have enormous heights for a generic right inverse `U` (entries above `10^500` at `q = 103` from Hermite elimination). This is irrelevant to P3-char, which does not see heights. It is exactly where e-all fails.

## 1. The record's definition (sharp.md §1.1, Definition 1.1), written out

**Cutoff.** Fix `M`, `C > 0` and `X >= 3`. For an odd prime `q` put `A_q = A_q(X) = max{a >= 0 : q^a <= X}`.
"At cutoff `X`" means that residue data are supplied, and collisions are counted, only at odd primes `q <= X`,
to the power `q^{A_q}`.

**Data.** A member is a tuple `(p, e, (a_x), phi, theta_0, delta, k, r)`:
* `p_j` are distinct primes `= 1 mod 4` (the labels);
* `e_j >= 1`;
* the rows `a_x in prod_j [0, e_j]` (the marks `x = 1..M`) are pairwise distinct, with `min_x a_xj = 0` and
  `max_x a_xj = e_j` (normalisation);
* `phi in R^r` (real angles) and `theta_0`;
* `delta in [0, Delta]^M` with `Delta = C e^{-W/4}` and `W = log N`, `N = prod p_j^{e_j}`;
* windings `k in Z^M`;
* prime residues `r`.

The conditions are:
* **(b)** `<2a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x`.
* **(c) (G)** Both gap inequalities of Lemma 0.1 hold for every `v != 0`.
* **(d) Factored residue data.** For every odd prime `q <= X` and every `j` with `p_j != q` there is a unit
  `r_j(q)` of `Z[i]/q^{A_q}` with `Norm r_j(q) = p_j mod q^{A_q}`. Lower levels are reductions.
  * *Point residues*, for `q` not dividing `N`: `rho_x(q) = i^{-k_x} prod_j r_j(q)^{a_xj} conj(r_j(q))^{e_j - a_xj}`.
  * *Level residues*, for `q = p_l`: the same product over `j != l`, including the unit `i^{-k_x}`.
  * "Factored" means exactly this. The point and level residues are the monomials above in per-prime residues
    `r_j(q)` of the individual norms `p_j`. The `r_j(q)` may be chosen independently for each `q`.
* **(e-char)** For every `c in Z^M` with `sum c = 0` and `v(c) = sum c_x a_x != 0`, and every unit
  `eta in {1, i, -1, -i}`:
  `|sin((c.delta - arg eta)/2)| >= 2^{-1/2} m_{(c,eta)} e^{-w(v(c))/2}`. Here
  `m_{(c,eta)} = prod_{q <= X, q not dividing N} q^{a_{q,eta}(c)}` and
  `a_{q,eta}(c) = max{a <= A_q : prod_x rho_x(q)^{c_x} = eta mod q^a}`.
* **The same-level clause.** sharp.md says "Add the same-level pair clause at the primes `p_l <= X`, exactly as
  `lower.md` Def. 1.1(e)". That clause reads: for `q = p_j <= X`, `a_q(c) = 0` *unless* `c = e_x - e_y` with
  `a_xj = a_yj`. In that case `a_q(c) = max{a : rho_x(p_j) = rho_y(p_j) mod p_j^a}`, and `p_j^{a_q(c)}` multiplies
  the modulus on the right of (1.1). As written this is a **unit-one** clause, and it concerns **only pairs at the
  same level**. Nothing is required at `p_l` for pairs at different levels or for non-pair characters.
* **(e-all)** For every `v in Z^r \ {0}` and `s in Z/4`:
  `|sin(<v,phi> - s pi/4)| >= nu_s m_{(v,s)} e^{-w(v)/2}`, where `m_{(v,s)}` counts the collisions
  `prod (r_j/conj r_j)^{v_j} = i^s` at `q <= X`, `q` not dividing `N`.

P3-char is (a)-(d) with (e-char). P3-all adds (e-all).

**Astra's claim against the definition.** She claims: "satisfying P3-char for **every** finite cutoff `X` with the
same prime labels, marks, real angles, windings and compatible prime residues". Her residues at an exterior odd
prime power are the reductions of fixed Gaussian rationals `r_j`. At `p_l` they are the reductions of the fixed
`r_j^{(l)}` (with `gamma_l` omitted). Hence:
* The data at cutoff `X` are the reductions mod `q^{A_q(X)}` of objects that do not depend on `X`.
* The point and level residues at cutoff `X'` reduce to those at `X <= X'`.
* Every collision modulus `m_{(c,eta)}(X)`, together with the same-level factors, divides the cutoff-free product
  `D_{u0}` of §7.

So one fixed model lies in `F^char_X` for every `X`. That is exactly P3-char at every cutoff, and slightly more,
because the record quantifies the residue data existentially for each `X` separately. **Verdict: Correct.**
* For pairs she also verifies a same-level clause with all four units. This is stronger than the record, and her
  `D_{u0}` bound covers it (§7).
* (a), (b), (c) and (e-char) for non-pair characters rest on the moment construction of VIII.3/VIII.4 and on
  (IX.7). Those are outside this block: VIII.3/VIII.4 are refereed, and (IX.7) belongs to another block.

## 2. The matrix (Proved, Computed)

Setup: `K_ij = chi(j-i)` on `F_q`, `q = 3 mod 4`, labels `I = F_q^*` (column 0 deleted), `a_ij = 1 + K_ij`,
`e_j = 2`. The full class is `A = {i in F_q^* : K_i0 = 1} = {i : chi(i) = -1}`, of size `(q-1)/2`.

**Coverage.** For `j != 0` and `epsilon = +-1`, the number of `i in A` with `K_ij = epsilon` is
`sum_{i != 0, j} (1+chi(-i))(1+epsilon chi(j-i))/4 = (q - 2 - epsilon - chi(-j) - epsilon chi(j))/4`.
This uses three facts:
* `sum_{i != 0,j} chi(-i) = -chi(-j)`;
* `sum_{i != 0,j} chi(j-i) = -chi(j)`;
* `sum_{i != 0,j} chi(-i) chi(j-i) = -1`, the Paley inner product.

With `chi(-j) = -chi(j)` this gives:
* endpoint `a = 2` (`epsilon = +1`): `(q-3)/4` occurrences;
* endpoint `a = 0` (`epsilon = -1`): `(q-1+2chi(j))/4` occurrences.

Both are `q/4 + O(1)`. The `i = j` entry has `a = 1`. The formula was checked exactly for every column at every
`q` in the list below.
* Under independent sampling at density `rho`, a column misses an endpoint with probability at most
  `2 rho^{(q-3)/4}`. The union over `q - 1` columns is `o(1)`. **Proved.**
* At small `q` coverage can fail (`q = 11`, `rho = 0.7`, `m = 3`). This is consistent with a statement that is
  only asymptotic.

**Parity.** `a_ij` is odd iff `K_ij = 0` iff `j = i`. Every mark has `i in F_q^* = I`, so mark `i` has an odd
exponent at exactly the one label `j = i`. **Proved; computed.**

## 3. `D`, `H` and the integral right inverse (Proved, Computed)

Definitions:
* `D` has rows `a_i - a_{i0} = K_i - K_{i0}` on `I`, for `i in S \ {i0}`.
* `H` has rows `e_i - e_{i0}` on `I`.

**`D = H mod 2`.** `K_ij = 1 - [j = i] mod 2`, so `K_ij - K_{i0 j} = [j=i] - [j=i0] mod 2`.

**Right inverse.** The argument has four steps.
* A rational combination `sum_{i != i0} lambda_i D_i` equals `v = K_S^T c`, restricted to `I`, with `c_i = lambda_i`
  and `c_{i0} = -sum lambda`. So `sum c = 0`.
* If `v` is integral on `I`, it is integral on all columns: `v_0 = sum c_i K_{i0} = sum c_i = 0`.
* Extend `c` by zero. Then `K K^T = qI - J` gives `K v = q c`, so `y = qc` is integral. Also `K^T y = q v = 0 mod q`.
* `(K^T y)_j = sum_i y_i (j-i)^{(q-1)/2} mod q` is a polynomial in `j` of degree `(q-1)/2 < q` that vanishes on
  `F_q`, so it is the zero polynomial. Its coefficients are `binom((q-1)/2, t)(-1)^t sum_i y_i i^t`, and the
  binomials are units. Hence the power sums vanish for `t = 0, ..., (q-1)/2`. The Vandermonde matrix on the at most
  `m < q/2` support points is injective, so `y = 0 mod q`, i.e. `c` is integral.

The same argument with `v = 0` gives full row rank (`qc = Kv = 0`). So the row lattice of `D` is primitive in
`Z^{q-1}`. All Smith invariants are 1, and an integral `U` with `DU = I_{m-1}` exists. **Proved.**

**Computed.** `q in {11,19,23,31,43,47,59,67,71,79,83,103}`, three random `S` each, at `rho = 1, 0.95, 0.7`:
* the sympy Smith normal form over `Z` has all invariants equal to 1 (`q <= 71`);
* my own column Hermite elimination gives `DQ = [L|0]` with unit pivots, and `U = Q[:, :m-1] L^{-1}` satisfies
  `DU = I` (all `q`).

`U` is not unique. Hermite elimination gives entries up to about `10^{500}` at `q = 103`, while Astra's flint HNF
gives maximum entry 9 at `q = 11` (mine: 18). The choice is immaterial for P3-char (§5).

## 4. `V` and `F` (Proved, Computed)

`H - D` is even (§3), and `U` is integral, so `V = U(H-D)/2` is integral. Then
`DF = D + 2DV = D + DU(H-D) = H`. Verified exactly in every instance of §3.

## 5. The witnesses (IX.6) (Proved, Computed)

Definitions: `gamma_l = pi_l/conj(pi_l)` and `u = U((k_i - k_{i0})_{i != i0}/2)`.

**The witnesses.** `r_j = i^{u_j} pi_j prod_l gamma_l^{V_jl}`.
* `conj(gamma) = gamma^{-1}`, so `r_j conj(r_j) = pi_j conj(pi_j) = p_j` and
  `r_j/conj(r_j) = i^{2u_j} gamma_j prod gamma_l^{2V_jl} = (-1)^{u_j} prod_l gamma_l^{F_jl}`.
* At an odd prime `q` not dividing `N`, every `pi_l` and `conj(pi_l)` is a unit, so these are units of
  `Z[i]/q^{A_q}` with the norms required by (d).

**Marks.** With `e_j = 2`, `r^a conj(r)^{2-a} = p_j (r/conj r)^{a-1}`, so
`rho_i = i^{-k_i} prod_j p_j prod_j (r_j/conj r_j)^{K_ij}`. Hence
`rho_i/rho_{i'} = i^{-(k_i-k_i')} (-1)^{sum_j (K_ij - K_i'j) u_j} prod_l gamma_l^{((D_i - D_i')F)_l}`.
* `(D_i - D_i')F = H_i - H_i' = e_i - e_i'`.
* `sum_j (K_ij - K_i'j) u_j = (Du)_i - (Du)_{i'} = (k_i - k_i')/2`, with `(Du)_{i0} := 0` for the anchor.
* So the sign is `i^{k_i - k_i'}` and cancels the winding factor. Hence `rho_i/rho_{i'} = gamma_i/gamma_{i'}`.
* For `sum c = 0`: `prod rho_i^{c_i} = prod gamma_i^{c_i} = A_c/conj(A_c)` with
  `A_c = prod_{c_i>0} pi_i^{c_i} prod_{c_i<0} conj(pi_i)^{-c_i}`. **Proved.**

**Windings.** The hypothesis "all `k_i - k_{i0}` even" is used once: to make `(k_i - k_{i0})/2` an integer, so that
`(-1)^{(Du)_i}` can cancel `i^{-(k_i - k_{i0})}`.
* With an odd difference the mark ratio is `+-i gamma_i/gamma_{i'}`. At `q = +-3 mod 8` this is the parity flip of
  sharp.md Cor. 2.2(4) / Rem. 2.3(1).
* The record permits the restriction. Windings are part of the model's data, constrained only by (b), and the
  record's own construction uses `k = 0`.
* Astra's sources supply even differences: VIII.3 has zero windings. For the VIII.4 spread variant she supplies the
  hypothesis (row differences of `2Kb` even). That is outside this block; the README_VIII spread referee checked
  the "winding shift".

**Independence of `U`.** The ratio `gamma_i/gamma_{i'}` does not involve `U`, `V` or `u`. So the e-char content of
the construction does not depend on the right inverse chosen. Only the witnesses' heights depend on it, which
matters for e-all.

## 6. The same-level deletion at parent primes (Proved, Computed)

At `q = p_l`, (d) asks for `r_j(p_l)` only for `j != l`. The level residue is the product over `j != l`. Astra
takes `r_j^{(l)} = i^{u_j} pi_j prod_{l' != l} gamma_{l'}^{V_{jl'}}`.
* `gamma_l` is not a unit mod `p_l`, since `p_l | pi_l conj(pi_l)`. It must be dropped.
* All other factors are units mod `p_l`, and `Norm r_j^{(l)} = p_j`.
* `r_j^{(l)}/conj(r_j^{(l)}) = (-1)^{u_j} prod_{l' != l} gamma_{l'}^{F_{jl'}}`.

**The ratio.** For marks `i, i'`:
`rho_i/rho_{i'} = i^{-(k_i-k_i')} (-1)^{sum_{j != l}(K_ij-K_i'j)u_j} prod_{l' != l} gamma_{l'}^{sum_{j != l}(K_ij-K_i'j)F_{jl'}}`.
* This differs from the exterior computation only by the missing `j = l` term, with coefficient `K_il - K_i'l`.
* For a **same-level** pair (`a_il = a_i'l`) that coefficient is 0. The `u`-sum then still equals
  `(k_i - k_i')/2`, and the `gamma` exponents equal `(e_i - e_{i'})_{l'}` for `l' != l`.
* The `l`-coordinate is dropped. It vanishes provided `l` is not in `{i, i'}`.
* That holds: `a_{i'l} = 1 + chi(l - i') = 1` iff `i' = l`. So when `l in S`, the owner `l` is the **unique** mark at
  level 1, and it is never part of a same-level pair. When `l` is not in `S`, no mark is at level 1.
* Same-level pairs are therefore exactly the pairs both at level 0 or both at level 2, and for them
  `rho_i/rho_{i'} = gamma_i/gamma_{i'}`. Both factors are units mod `p_l`, since `p_i, p_{i'} != p_l`.
* This holds at every power `p_l^a` and is compatible under reduction. **Proved.**

**What the record requires at `p_l` beyond this.** Nothing (§1). The `lower.md` clause sets `a_{p_l}(c) = 0` for
every `c` that is not a same-level pair, and `m_{(c,eta)}` runs over `q` not dividing `N`. So pairs at different
levels and non-pair characters impose no condition at parent primes, and the construction meets the record.
* Computed (information only): at different levels the identity sometimes holds and sometimes fails, e.g.
  396/972 at `q = 19`. That is irrelevant.
* The phrase "this handles parent primes above `2m` as well as below it" is harmless but vacuous here. Every
  parent prime is `>= P = q^d > q > 2m`. The argument is size-free, so it would also handle small parent primes.

**The unit question at `p_l`.** As written the record's clause is unit-one. If one reads e-char's unit `eta` into
the same-level clause, the construction still satisfies it. A collision `gamma_i/gamma_{i'} = u0 mod p_l^a` is
equivalent to `p_l^a | pi_i conj(pi_{i'}) - u0 conj(pi_i) pi_{i'}` (the denominators are units mod `p_l`), so it
enters the same `D_{u0}` of §7.

## 7. Global collision cost (Proved, Computed)

Fix a character `c` and a unit `u0`.
* Each maximal colliding power, exterior or (for a pair) same-level internal, divides `A_c - u0 conj(A_c)`. The
  powers are pairwise coprime, so their product `D_{u0}` divides it.
* `A_c - u0 conj(A_c) != 0`. The `pi_i` are split primes over distinct `p_i`, so `A_c` and `conj(A_c)` have
  different prime factorisations unless `c = 0`.
* It is divisible by `1+i`: odd-norm elements and units are `1 mod (1+i)`.
* Hence `sqrt2 D_{u0} <= 2|A_c| <= 2 L_0^{n/2}`, i.e. `D_{u0} <= sqrt2 L_0^{n/2}`, with no cutoff.

For a pair and `u0 = 1`:
* `A - conj(A) = 2i Im(pi_i conj(pi_{i'}))`. `D_1` is odd and divides it, so `D_1 | Im(...)`, and
  `1 <= D_1 <= sqrt(p_i p_{i'}) |sin(beta_i - beta_{i'})| <= L_0 |beta_i - beta_{i'}|`.
* `Im != 0` is a rational integer, so `|beta_i - beta_{i'}| >= |sin| >= 1/sqrt(p_i p_{i'}) >= 1/L_0`. **Proved.**

**Computed** (`q = 19`): 12,384 checks of "collision of `gamma_i/gamma_{i'}` with `u0` mod `Q` iff
`Q | A - u0 conj A`", over all pairs, all four units, and all odd `Q = l, l^2` with `l < 200` prime to
`p_i p_{i'}`. All passed. The maximum of `odd|Im|/(L_0|dbeta|)` was 0.84 <= 1.

## 8. Primes, (IX.5), offsets and the pair inequality (Proved)

**Pigeonhole.** There are `~ P/(2 log P)` split primes in `[P, 2P]`. Choose the Gaussian factor with argument in
`[0, pi/4]`. Use `ceil(pi/(4 b_*)) = O(1/min(1,C))` sectors of width `b_*` and
`B = c P/(q log P)` equal-ratio radial bins of ratio `e^sigma`, `sigma = log 2/B`.
* The average cell holds `~ q min(1,C)/(1571 c)` primes, which is `>= q` for `c <= min(1,C)/1600` and large `q`.
  Some cell therefore holds `q - 1` primes.
* `sigma = (log 2/c) q log P/P = O_C(q log q/P)`. **Proved.** No equidistribution in sectors is needed.

**(IX.5).** For marks `i != i'`, the coefficients `|K_ij - K_{i'j}| - 1` on `I` lie in `{-1, 0, 1}`. Their sum is
`(2q - 2)/2 + 2 - (q-1) = 2`: use `||K_i - K_{i'}||^2 = 2q`, the two zero positions, and `K_{i0} = K_{i'0}`.
* So `w - W/2 = 2 log P_0 + sum_j coef_j log(p_{pi(j)}/P_0)`.
* Under a random permutation the mean is `2 x mean in [0, 2 sigma]`. A transposition moves the sum by `O(sigma)`.
  The permutation bounded-difference inequality gives deviation `O(sigma sqrt(q log q))`, simultaneously over the
  `O(q^2)` pairs.
* This is `o(1)` iff `d > 3/2`, which is exactly where the hypothesis `d > 3/2` enters. **Proved.**
* README_VIII §1.3 notes that this `o(1)` is about `-1` at computable `q`. That is harmless here too (see below).

**Offsets.** Quarter-turn units, `delta = (pi/2) epsilon`, `Delta = min(1,C) e^{-W/4}/10`. The physical width is
`(pi/2) Delta < C e^{-W/4}`.
* Centres `Delta(1/4 + (beta_i - beta_min)/(2b_*))` lie in `[Delta/4, 3Delta/4]`, and the radius is
  `Delta/(100 L_0)`. So `epsilon_i in [0, Delta]`.
* `|epsilon_i - epsilon_{i'}| >= Delta |dbeta|/(2b_*) - 2Delta/(100 L_0) >= Delta(1/(2b_*) - 0.02)|dbeta|`, using
  `|dbeta| >= 1/L_0`.

**Unit one.** The left side is `|sin((pi/4) depsilon)| >= |depsilon|/2 >= (25 - 10^{-3}) |dbeta| e^{-W/4}`.
* By (IX.5) and `L_0 <= P_0 e^sigma`, the right side is
  `2^{-1/2} D_1 e^{-w/2} <= 2^{-1/2} (L_0/P_0) e^{o(1)} |dbeta| e^{-W/4} = (0.707 + o(1)) |dbeta| e^{-W/4}`.
* The factor of about 35 in hand also absorbs an `O(1)` loss up to about `2 log 35` in (IX.5).
* This is the reason for the narrow sector. With VIII.3's centres `Delta(1/4 + beta_i/2)` the left side would be
  about `0.024 |dbeta| e^{-W/4}`. That is too small, because here `D_1` comes from the physical primes, so
  `L_0/P_0 ~ 1`. In VIII.3 the auxiliary primes had `L_0 = o(P)`.

**Other conditions.**
* *Other units:* the sine is `>= sin(pi/4) - o(1)`, while `2^{-1/2} D_{u0} e^{-w/2} <= L_0 e^{-W/4}/P_0 = O(e^{-W/4})`.
* *The (G) gaps for the pair monomial* hold with the same room.
* *Uniform in the cutoff:* every bound used is `D_{u0}`, which does not involve `X`.
* *Quantifiers:* `b_*` depends on `C`; the threshold on `q` depends on `C`, `rho` (and `d`); nothing depends on `X`.
  **Proved.**

**Minor correction (statement).** `m/omega(N') -> 1/2` in (IX.3) holds along a diagonal sequence `rho_k -> 1`.
At fixed `rho` the limit is `rho/2`, since `m = (rho/2 + o(1)) q` and `omega = q - 1`. Corrected statement:
*for fixed `.99 < rho < 1`, `m/omega(N') -> rho/2`; along a diagonal sequence `rho_k -> 1`, `m/omega(N') -> 1/2`.*
This is the convention already accepted in README_VIII §1.3.

## 9. Reproduction (Computed)

`residue_check.py` builds `U`, `V`, `F` and the witnesses for random `S`. It uses small split primes, randomly
permuted onto the columns: heights are irrelevant to the identities. Windings are nonconstant and odd, with even
differences. For each instance it verifies:
* `Norm r_j = p_j` and `r_j/conj r_j = (-1)^{u_j} prod gamma^{F_j}`, mod each `Q`;
* every exterior pair identity `rho_i gamma_{i'} = rho_{i'} gamma_i`;
* 200 random zero-sum characters per modulus: `prod rho^c = prod gamma^c`;
* at every parent prime power `p_l, p_l^2, p_l^3` (with `gamma_l` omitted), every same-level identity, and that the
  owner `l` is the unique mark at level 1.

Exterior moduli: `l^a`, `l < 80` an odd prime not among the labels, `a <= 3` (`a <= 2` for `q >= 23`).

| instance | exterior pair identities | random characters | same-level identities | failures |
|---|---|---|---|---|
| `q = 11`, `m = 5` | 360 | 7,200 | 90 | 0 |
| `q = 19`, `m = 9` | 1,296 | 7,200 | 756 | 0 |
| `q = 19`, `m = 8` | 1,008 | 7,200 | 588 | 0 |
| `q = 23`, `m = 9` | 864 | 4,800 | 648 | 0 |
| `q = 31`, `m = 14` | 2,184 | 4,800 | 2,366 | 0 |

* **Negative control:** an odd winding difference at `q = 19` gives 1,656 failures. This confirms that the parity
  hypothesis is used.
* Astra's `q = 11` counts (330 exterior with 11 moduli primes, 90 same-level) are consistent with mine
  (12 primes, giving 360).

## 10. Summary of labels

* **Proved:** §§2-8 as stated, with the corrections in §0 (iii) and §8 to the `m/omega` limit. The phrase in §6
  about parent primes "below `2m`" is vacuous here.
* **Computed:** §§2, 3, 4, 7, 9.
* **Gap / False:** none in this block.
* **Not checked here:** the non-pair paragraph (IX.7) and its offset count, and the e-all failure paragraph. Both
  are other blocks.
