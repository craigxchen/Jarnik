# Round 7: lower bounds from factored residue data (P3-all, P3-gen)

Conventions are those of `round4/outside.md`, `round5/lower.md` and `round6/sharp.md`:

* `W = log N` and `R = N^(1/2)`. The arc has angular width `Delta = C e^(-W/4)`.
* `w(v) = sum_j |v_j| log p_j`, and `d_xy = w(a_x - a_y)`.
* The *pair margin* is `mu_xy := d_xy - W/2`. This is the slack `s_c` of the pair character
  `c = e_x - e_y`.
* `m_xy` is the unit-1 pair collision modulus of the class in question. Precisely, it is
  `m_c(rho)` of `lower.md` Def. 1.1(e) for `c = e_x - e_y`, including the same-level clause at
  primes `p_j <= X`.
* The *pair demand* is `D_xy := 2 log m_xy >= 0`. The *point total* is `D_x := sum_(y != x) D_xy`.
  Put `c_C := 2 log C - log 2`.
* `lambda_max(B)` is the largest eigenvalue of a real symmetric matrix `B`. For entrywise
  nonnegative `B` it is the Perron root. `J` is the all-ones matrix.

*Theorem*, *Lemma* and *Proposition* mean proved here with all quantifiers, or cited exactly. The
citations are classical facts, or statements of `round5`/`round6` that are re-derived where they are
used. *Evidence* means a computation that no proof uses.

Scripts live in `round7/p3lower_checks/`. Each runs from that directory with `python3 <script>` (or
`./anneal`), and its output is stored next to it as `*_out.txt`. The repository and the round-5/6
folders were not modified.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The target was not reached.** No growth improvement `W >= f(M) M log M` with `f -> infinity` is
proved for actual circles, for P3-all or for P3-gen. This round proves instead what any such proof
must avoid, and one conditional restricted theorem. Every statement below names its class. The
classes are nested: `P3-gen ⊂ P3-all ⊂ P3-char ⊂ A_X`. Actual clusters lie in all of them
(Prop. 1.2).

| # | statement | class | status |
|---|---|---|---|
| A | **The pair channel is a cut LP.** Margins are nonnegative combinations of cuts of total weight `W` (Lemma 1.4). The least `W` compatible with pair demands `D` satisfies `2 lambda_max(D) <= W_cut(D) <= pi lambda_max(D)` (Perron below, Goemans-Williamson above), with an exact LP dual (Theorem A) | every member of `A_X`, hence actual clusters | **proved** |
| A.1 | valid lower bound: `W >= 2 lambda_max(D) - 2(M-1) c_C^+` | `A_X`, actual | proved (contains the cut/Plotkin budget and every PSD consequence) |
| A.2 | **barrier:** no inequality derived from the pair inequalities and the cut structure with real weights exceeds `pi (lambda_max(D) + (M-1) c_C^-)` | the real-weight relaxation | proved |
| A.3 | consequence: growth through pairs needs `lambda_max(D) >> M log M`, so needs a **point total** `max_x D_x >> M log M`. The maximal **pair** demand enters only through `lambda_max(D) >= D_max` | -- | proved |
| A.5 | **for P3-all no pair-cut certificate proves growth.** The members of `sharp.md` Thm 6.4 have `lambda_max(D) <= 36 M log M + 111 M` | P3-all at `X = 3M^A` | proved |
| B.1 | `K`-uniform margins: `W >= (2(M-1)/K)(D_max - c_C)` | `A_X` | proved (elementary) |
| B.2 | column-balanced profiles: `sum_y mu_xy = W/2` exactly, so `W >= 2 max_x D_x - 2(M-1) c_C` | `A_X` | proved |
| B.3 | **forced Legendre parity** refines the pigeonhole to `E(n_0, m/2) + E(n_1, m/2)`. This gains at most a factor 2 per modulus, and only through parity bias. Over many moduli at once (CRT) the gain in `D_max` is `<= 2|S| log 2` | P3-char and smaller, actual | proved, exact check on actual points |
| B.5 | monomial heights plus pigeonhole give only `W >= r log min(r, X) - O(r)`, which the count of distinct primes already gives | P3-all | proved |
| C | **conditional restricted growth.** Suppose the margins are `K`-uniform and the collision counts satisfy the Poisson-type factorial-moment Hypothesis `H(eta, delta)`. Then `W >= ((8(1-eta) - eps)/K) M log^2 M / loglog M`, so `M <= (K/(4(1-eta)) + o(1)) log R logloglog R/(loglog R)^2` | `A_X` (actual, P3-gen, P3-all) with these two hypotheses | **proved, conditional on H** |
| C.2/C.3 | the round-6 design `D_A` violates `H` for large `M`. The independent random residue model satisfies `H(eta, delta)` with probability `-> 1` for small `delta` | -- | proved |
| D.1 | with **equal weights**, `r >= lambda_max(k)` columns are needed. Clique-partition demand patterns attain this exactly; uniform margins cost `(M-1) k_max` | `+-1` matrices | proved, exact checks |
| D | discreteness (prime weights, few columns) is the one pair-level channel left open. The real-weight LP optima use `binom(M,2)` cuts | -- | open; evidence in §5 |

**What this means for the growth question.**

1. *Forced high-demand pairs cannot give growth through pair-cut certificates.* These are all
   arguments that use the pair inequalities together with the cut/Plotkin/PSD structure of
   margins, with real weights (Cor. A.2). This would remain true even if the Poisson heuristic
   `max_xy D_xy ~ 4 log^2 M/loglog M` were proved for actual circles. The reason is Theorem A.
   Non-uniform margins absorb heavy pairs at total cost at most `pi lambda_max(D)`. For random-like
   demands (the same model that produces the Poisson tail) that is `~ 2 pi M log X`. Uniform
   margins would need `2(M-1) D_max ~ 8 M log^2 M/loglog M`. For actual circles `lambda_max(D)` is
   not known; the bound `lambda_max(D) <= max_x D_x` is all that is proved.
2. *The mechanism works only when margins are forced to be near-uniform.* There it does give
   growth, conditional on `H` (Theorem C).
3. *Uniform margins are an assumption, not a theorem.* The round-6 heuristic obstruction for
   coherent designs (`sharp.md` §6.6: "max pair demand forces `(b-4) tau/2 >~ log^2 M/loglog M`")
   applies to uniform-margin (Paley) profiles. At the real-weight level it disappears for profiles
   adapted to the demand pattern (Cor. A.2). Whether prime weights allow the adaptation is open
   (Question Q).
4. *Parity, heights and PSD give at most constant factors.* For PSD this is exact: PSD is the Perron
   bound.
5. *Three channels remain for P3-all\* and for actual circles:*
   * (i) discreteness of the weights, i.e. few columns of prime weight (Question Q, §5);
   * (ii) completion and covering by the monomial slabs outside `L` (P3-all\*, unchanged);
   * (iii) global arithmetic of the actual pair integers: Hypothesis `H` together with a proof
     that margins must be near-uniform, or that point totals must be irregular.

**The task's mechanisms, one by one.**

| mechanism | result | where |
|---|---|---|
| forced high-demand pairs | useless against non-uniform margins: the pair channel is the cut LP, capped by `pi lambda_max(D)` | Thm A, Cor. A.2, A.3 |
| cut/Plotkin budget `sum mu <= W M/4` | the Perron test vector `u = 1/sqrt M` | Cor. A.1 |
| PSD-ness of the margin Gram matrix | carries exactly the Perron bound and nothing more; cut realisability adds at most a factor `pi/2` | Cor. A.1, Thm A(3) |
| height condition for every monomial | its pigeonhole consequences reproduce the trivial count of primes | Prop. B.5 |
| forced Legendre parities (products of prime residues) | at most a factor 2 per modulus, only through parity bias; `+2|S| log 2` over many moduli at once | Prop. B.3, Rem. B.4 |
| restricted class: Hadamard-type, near-uniform margins | growth, conditional on the named hypothesis `H` about joint smooth parts of the pair integers | Thm C |
| restricted class: primes in a short interval (equal weights) | reduces to a sign-matrix question (Q); structured patterns are realisable at the Perron bound | §5, Prop. D.1 |
| consistency with the round-5/6 fakes | no contradiction; `D_A` violates `H` | §6, Lemma C.2 |

**Consistency with rounds 5/6** (§6):
* every statement here is either a valid inequality for all members or concerns a restricted class;
* the round-5/6 fakes are `K`-uniform with `K <= (b+4)/(b-4) + o(1)`;
* their `D_A` residues violate `H` (Lemma C.2);
* none of the statements contradicts the constructions.

---------------------------------------------------------------------------------------------

## 1. Classes, inputs, and the cut representation of margins

### 1.1 Classes

* `A_X(M, C)` is the class of angle models with point-residue data (`lower.md` Def. 1.1).
* **P3-char** `= F^char_X(M, C)` and **P3-all** `= F^all_X(M, C)` are the classes of factored
  residue data (`sharp.md` Def. 1.1). In P3-all every monomial `A_v`, `v in Z^r`, is strengthened
  by its own collision modulus. Relations: `F^all_X ⊂ F^char_X ⊂ A_X` (`sharp.md` §1.1).

**Definition 1.1 (P3-gen: genuine prime residues).** Fix `X >= 3`. A member of `F^all_X(M, C)` is
*genuine-residue* if there are Gaussian primes `pi_j` with `pi_j conj(pi_j) = p_j` such that

```text
r_j(q) = pi_j  mod q^(A_q)    for every odd prime q <= X with q != p_j.
```

The angles `phi` are not required to be the arguments of the `pi_j`. We write
`P3-gen_X(M, C) ⊂ F^all_X(M, C)`.

In P3-gen:
* the point residues are those of the genuine Gaussian integers
  `zhat_x = i^(-k_x) prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj)` (Def. 1.1(d) of `sharp.md`);
* pair collisions are divisibilities of genuine integers. With `g_xy = gcd(zhat_x, zhat_y)` and
  `zhat_x/g_xy = u B_xy`, `zhat_y/g_xy = u' conj(B_xy)` (units `u, u'`, `|B_xy| = e^(d_xy/2)` by
  `two_thirds` Lemma 2.1), the unit-1 collision `zhat_x = zhat_y mod q^a` means
  `q^a | (zhat_x - zhat_y)/g_xy`. That reduced chord is a nonzero Gaussian integer of modulus
  `<= 2 e^(d_xy/2)`; it is nonzero because distinct rows give non-associate `zhat`;
* the **height conditions** of `sharp.md` Prop. 6.1 hold automatically. The modulus `m_(v,s)`
  divides the nonzero rational integer `D_s(prod_j pi_j^(v_j^+) conj(pi_j)^(v_j^-))`, whose absolute
  value is at most `nu_s^(-1) e^(w(v)/2)` (proof of `sharp.md` Prop. 1.2). The slab conditions (1.2)
  themselves involve the free angles `phi`, and are *not* automatic.

**Proposition 1.2 (actual clusters are members).** Let `z_1, ..., z_M` be a normalised actual
cluster on an arc of length `<= C sqrt R` (`two_thirds` Lemma 1.1). Then for every `X >= 3` it
defines a member of `P3-gen_X(M, C)`, hence of P3-all, P3-char and `A_X`.

*Proof.* `sharp.md` Prop. 1.2 constructs from the cluster a tuple in `F^all_X(M, C)` for every `X`.
* It takes `phi_j = arg pi_j`, with `delta` and `k` from (1.1) of `outside.md`.
* It takes the residues `r_j(q) = pi_j mod q^(A_q)`. These are genuine by Definition 1.1.

That proof was re-read here. Its points are:
* `rho_x(q) = conj(eps_0) z_x mod q^(A_q)` for a common unit `eps_0`, which cancels in
  `prod rho_x^(c_x)` since `sum c = 0`;
* a collision `U = eta V mod q^a` forces `q^a | U/G - eta V/G`, which is a nonzero multiple of
  `(1+i)`;
* the identity `|U - eta V| = 2 N^(n/4) |sin((c.delta - arg eta)/2)|` holds;
* for monomials, `q^a` divides the nonzero rational integer `D_s(A_v)`, and
  `|D_s(A_v)| = nu_s^(-1) |A_v| |sin(<v,phi> - s pi/4)|`. QED

### 1.2 The pair inequality

**Lemma 1.3 (pair inequality).** Every member of `A_X(M, C)` satisfies
`mu_xy >= D_xy - c_C` for all `x != y`.

*Proof.* Take `lower.md` (1.1) with `c = e_x - e_y`. Then `w(v(c)) = d_xy`, and
`|sin((delta_x - delta_y)/2)| >= 2^(-1/2) m_xy e^(-d_xy/2)`. The left side is
`<= |delta_x - delta_y|/2 <= Delta/2 = (C/2) e^(-W/4)`. Hence

```text
log m_xy <= log C - (1/2) log 2 + mu_xy/2,   i.e.   mu_xy >= 2 log m_xy + log 2 - 2 log C.
```

QED. In P3-char and P3-all this is (e-char) with `eta = 1`.

### 1.3 Margins are cut combinations

A *cut* is a nonconstant `s in {+-1}^M`, taken modulo `s ~ -s`. For a function `lambda >= 0` on cuts,
with finite support, put

```text
W(lambda) := sum_s lambda_s,   mu(lambda)_xy := sum_s lambda_s ([s_x != s_y] - 1/2) = -(1/2) sum_s lambda_s s_x s_y,
sigma(s) := sum_x s_x.
```

**Lemma 1.4 (cut representation).** Take any profile `(p, e, (a_x))` normalised as in
`lower.md` Def. 1.1(a), with `w_j = log p_j`. For each column `j` and each `1 <= t <= e_j`, let
`s^(j,t)_x = +1` if `a_xj >= t` and `-1` otherwise. Let `lambda^prof` put weight `w_j` on each
`s^(j,t)`. Then:
1. every `s^(j,t)` is a cut (nonconstant);
2. `W(lambda^prof) = W` and `mu(lambda^prof)_xy = mu_xy` for all `x != y`;
3. *(point totals)* `sum_(y != x) mu_xy = W/2 - (1/2) sum_s lambda_s s_x sigma(s)` for every `x`;
4. *(cut identity)* `sum_(x<y) mu_xy = M W/4 - (1/4) sum_s lambda_s sigma(s)^2 <= M W/4`.

*Proof.*
1. Normalisation gives rows with `a_xj = 0 < t` and rows with `a_xj = e_j >= t`.
2. For integers `0 <= a, b <= e`, `|a - b|` is the number of `t in [1, e]` with exactly one of
   `a >= t`, `b >= t` (`two_thirds` Lemma 2.2). Summing `w_j` over `j` and `t` gives `d_xy`,
   and the total weight is `sum_j e_j w_j = W`.
3. `sum_(y != x) (-(1/2)) s_x s_y = -(1/2) s_x (sigma(s) - s_x) = 1/2 - (1/2) s_x sigma(s)`. Multiply
   by `lambda_s` and sum.
4. Sum item 3 over `x` and halve. QED

For threshold cuts `sigma(s)^2/4 = (S_jt - M/2)^2` with `S_jt = #{x : a_xj < t}`. So item 4 is the
`two_thirds` identity `sum_(x<y) d_xy = (M^2/4) W - Q`, rewritten.

---------------------------------------------------------------------------------------------

## 2. Theorem A: the pair channel

For a real symmetric `B` with zero diagonal put

```text
W_cut(B) := inf { W(lambda) : lambda >= 0 on cuts,  mu(lambda)_xy >= B_xy  for all x < y }.
```

**Theorem A.** Let `M >= 2`.
1. *(LP duality)* The infimum is attained, and

   ```text
   W_cut(B) = max { sum_(x<y) nu_xy B_xy :  nu >= 0,  sum_(x<y) nu_xy ([s_x != s_y] - 1/2) <= 1  for every cut s }.
   ```

   For `nu != 0` let `surplus(nu) := max_s cut_nu(s) - |nu|/2 > 0` (`cut_nu(s)` = weight of the edges
   separated by `s`, `|nu| = sum nu_xy`). Then `W_cut(B) = max_(nu >= 0, nu != 0) <nu, B>/surplus(nu)`
   whenever `W_cut(B) > 0`.
2. *(Perron lower bound)* `W_cut(B) >= 2 sup { u^T B u : u >= 0, |u| = 1 }`. In particular
   `W_cut(B) >= 2 lambda_max(B)` if `B >= 0`. More generally, for every integer vector `b >= 0`
   with `sum b` odd and `|b|^2 >= 2`, `W_cut(B) >= 2 b^T B b / (|b|^2 - 1)` (hypermetric form).
3. *(Goemans-Williamson upper bound)* `W_cut(B) <= pi lambda_max(B^+)`, where `B^+` is the entrywise
   positive part.

So for `B >= 0`: `2 lambda_max(B) <= W_cut(B) <= pi lambda_max(B)`.

*Proof.*
1. *LP duality.* This is a finite LP with `2^(M-1) - 1` variables and `binom(M,2)` constraints.
   * It is primal feasible. Let `lambda` be uniform on all cuts. Then
     `P(s_x != s_y) = 2^(M-2)/(2^(M-1) - 1) > 1/2`, so every margin is positive, and scaling
     `lambda` up gives `mu >= B`.
   * The dual is feasible with `nu = 0`.
   * Strong LP duality gives equality of the optima, and both are attained.
   * *Surplus form.* `surplus(nu) > 0` for `nu != 0`, since a uniformly random cut has
     `E cut_nu(s) > |nu|/2`. Scale an optimal `nu` so that its constraint is tight.
2. *Perron lower bound.* Let `u >= 0`, `|u| = 1`, and put `nu_xy := 4 u_x u_y`. For every cut,

   ```text
   sum_(x<y) 4 u_x u_y (-(1/2) s_x s_y) = -(sum_(x != y) u_x u_y s_x s_y) = |u|^2 - (u.s)^2 <= 1,
   ```

   so `nu` is dual feasible, with value `sum_(x<y) 4 u_x u_y B_xy = 2 u^T B u`.
   * For the hypermetric form, take `nu_xy := 4 b_x b_y/(|b|^2 - 1)` and use `(b.s)^2 >= 1`, which
     holds because `b.s = sum b mod 2` is odd.
3. *Goemans-Williamson upper bound.* Since `mu >= B^+` implies `mu >= B`, we may assume `B >= 0`.
   * If `Lambda := lambda_max(B) = 0`, then `B = 0` and `lambda = 0` works.
   * Otherwise put `rho := I - B/Lambda`. All eigenvalues of `B` are `<= Lambda`, so `rho` is
     positive semidefinite with unit diagonal. Also
     `B_xy <= (e_x + e_y)^T B (e_x + e_y)/2 <= Lambda`, so `rho_xy = -B_xy/Lambda in [-1, 0]`.
   * Write `rho = Gamma^T Gamma` with unit columns `gamma_x`. Let `xi` be a standard Gaussian vector
     and `s_x := sgn <gamma_x, xi>`; ties have probability 0.
   * For unit vectors at angle `theta`, a uniformly random hyperplane through 0 separates them with
     probability `theta/pi` (Sheppard). So `P(s_x != s_y) = arccos(rho_xy)/pi`.
   * Let `pi_0 = P(s is constant)`, and put `lambda_(+-s) := pi Lambda P(s in {s, -s})` on
     nonconstant `s`. Constant outcomes are dropped. Then `W(lambda) = pi Lambda (1 - pi_0) <= pi Lambda`.
   * Constant `s` separate nothing, so

     ```text
     mu(lambda)_xy = pi Lambda (P(s_x != s_y) - (1 - pi_0)/2) >= Lambda (arccos rho_xy - pi/2)
                   = Lambda arcsin(-rho_xy) >= Lambda (-rho_xy) = B_xy,
     ```

     using `arcsin y >= y` on `[0, 1]`.
   * `lambda` is a finite measure on the finitely many cuts. QED

**Corollary A.1 (valid lower bound).** Every member of `A_X(M, C)` satisfies

```text
W >= W_cut(D - c_C (J - I)) >= 2 lambda_max(D) - 2 (M-1) c_C^+ .
```

This covers every member of P3-char, P3-all and P3-gen, and every normalised actual cluster
(Prop. 1.2). For `c_C <= 0`, i.e. `C <= sqrt 2`, the bound is `W >= 2 lambda_max(D)`.

*Proof.*
* Lemma 1.4 gives `lambda^prof` with `W(lambda^prof) = W` and `mu(lambda^prof) = mu`.
* Lemma 1.3 gives `mu >= B := D - c_C(J - I)`. Hence `W >= W_cut(B)`.
* Apply Theorem A(2) with the Perron unit vector `u >= 0` of `D`:
  `u^T B u = lambda_max(D) - c_C ((sum u)^2 - 1)`, and `0 <= (sum u)^2 - 1 <= M - 1`. QED

Special cases:
* `u = M^(-1/2) 1` gives the cut/Plotkin budget `sum_(x<y) mu_xy <= M W/4`. Together with the
  pigeonhole this is the `two_thirds` mechanism.
* Since the strongest consequence of "`W I - 2 mu` is PSD and `mu >= D - c_C`" is exactly the
  Perron bound, **PSD-ness of the margin Gram matrix carries no information beyond Corollary A.1**.
  * The best a PSD argument can do: if `B >= 0`, then `W' I - 2B` is PSD exactly when
    `W' >= 2 lambda_max(B)`.
  * Test vectors with mixed signs give nothing, because only lower bounds on `mu_xy` are known.

**Corollary A.2 (barrier for pair-cut certificates).** Call an inequality `W >= Phi(D, M, C)` a
*pair-cut certificate* if it holds for every pair `(W, mu)` such that:
* `mu >= D - c_C` entrywise, and
* `(W, mu) = (W(lambda), mu(lambda))` for some real `lambda >= 0` on cuts.

This includes:
* the cut identity, the Plotkin bound and point totals;
* PSD-ness, hypermetric, triangle and all other cut-cone inequalities;
* every maximal-demand argument;
* any combination of these with the pigeonhole lower bounds on `D`.

Then, for every nonnegative symmetric `D`,

```text
Phi(D, M, C) <= W_cut(D - c_C(J - I)) <= pi ( lambda_max(D) + (M-1) c_C^- ).
```

*Proof.* `(D - c_C(J-I))^+ <= D + c_C^-(J - I)` entrywise. The Perron root is monotone on
nonnegative matrices and subadditive, and `lambda_max(J - I) = M - 1`. Then use Theorem A(3). QED

**Corollary A.3 (what pair arguments would need).** `lambda_max(D) <= max_x D_x` (row-sum bound for
the Perron root). Hence a pair-cut certificate can prove `W >= f(M) M log M` on a class only if
every member of the class has

```text
max_x D_x  >=  (f(M)/pi) M log M - (M-1) c_C^- .
```

So pair arguments need **irregular point totals**. The maximal pair demand
`D_max = max_xy D_xy <= lambda_max(D)` enters only through this weak inequality.
* The pigeonhole of `two_thirds` forces only the *average* point total, `~ 2 M log M` at moduli
  `<= 2M`.
* Parity (Prop. B.3) raises that average by at most a factor 2.

**Lemma A.4 (point totals of the round-6 design).** In the design `D_A` (`sharp.md` Def. 3.2),
for every parity colouring and every choice of the random lifts and injections, every row `x`
satisfies

```text
D_x <= 4 sum_(d=1)^(M-1) (Lambda'(d) + log d) <= 36 M log M + 111 M .
```

*Proof.*
* Only small `q` (`q < 4M`) at levels `a <= a*(q)` produce unit-1 pair collisions. Level `a*(q)`
  is injective, lifts preserve injectivity, and the large-`q` maps are injections (`sharp.md`
  Lemma 3.3(1) proof).
* `sharp.md` Lemma 3.3(1) gives `log m_xy <= Lambda'(d) + log d` with `d = |x - y|`. Each `d` occurs
  at most twice as `|x - y|`, `y != x`, so `D_x <= 4 sum_(d<M) (Lambda'(d) + log d)`.
* Bound the first sum:
  `sum_(d<M) Lambda'(d) = sum_q log q #{d < M : p'(q) | d} <= sum_(q < 8M) log q (M/p'(q))`.
* By `sharp.md` Lemma 3.1(1), `q < 8 p'(q)`, so the sum is `<= 8M sum_(q < 8M) log q/q`.
* The elementary bound `sum_(q <= n) log q/q <= log n + log 4` follows from
  `n log n >= log n! >= sum_(q <= n) log q (n/q - 1)` and `theta(n) < n log 4` (`two_thirds`
  Lemma 7.1).
* So the first sum is `<= 8M (log M + log 32)`, and `sum_(d<M) log d <= M log M`.
* Altogether `D_x <= 4(9 M log M + 8 M log 32) <= 36 M log M + 111 M`. QED

**Corollary A.5 (P3-all: no pair-cut certificate proves growth).** Fix `A >= 1` and `C > 0`.
Suppose a pair-cut certificate satisfies `W >= Phi` on all of `F^all_(3M^A)(M, C)`. Then on the
members constructed in `sharp.md` Thm 6.4 (and Thm 5.1 for P3-char) it gives at most
`pi(36 M log M + 111 M + (M-1) c_C^-)`, by Corollary A.2 and Lemma A.4.

Those members use the design `D_A`. In Thm 6.4 the random fibre changes no point class
(`sharp.md` Lemma 6.3(1)), so their pair demands are those of `D_A`. In particular **no pair-cut
certificate proves `W_min^all >= f(M) M log M` with `f -> infinity`.** Any growth proof for P3-all
must use one of:
* the discreteness of the weights (§5);
* the monomial heights and completion (`sharp.md` §6);
* or, for actual circles, inputs outside P3-all.

*Evidence* (`lp_gw.py`, exact LP for `M <= 11`):
* on 48 demand matrices (uniform, heavy pair, matching, random collisions at `X = 3M` and
  `X = M^2`, iid Poisson tails), `W_cut/(2 lambda_max)` lies in `[1.000, 1.224]`;
* LP duality holds numerically, and so does the cut identity at the optimum;
* on the non-uniform matrices the uniform-margin cost (at least `2(M-1) D_max`, by Perron) is 1.4-7.1 times `W_cut`;
* every optimum uses many cuts, typically `binom(M, 2)` (an LP vertex). This matters in §5.

Statistics of the residue models (`stats_H.py`; evidence):

| model | `M` | `D_max/log M` | `lambda_max(D)/(M log M)` | `max_x D_x/(M log M)` |
|---|---|---|---|---|
| independent random residues, `X = 3M` | 64 ... 1024 | 10.2, 10.9, 12.0, 13.5, 14.5 | 1.905 ... 1.926 | 2.51 ... 2.08 |
| round-6 design `D_A` | 64 ... 1024 | 9.6, 11.6, 10.9, 11.5, 12.4 | 2.66 ... 3.19 | 3.34 ... 3.53 |
| genuine residues of random Gaussian primes `~10^6`, random profile `r = 4M` | 32, 64, 128 | 11.3, 11.3, 10.7 | 1.79, 1.72, 1.76 | 2.45, 2.07, 2.10 |

So `D_max/log M` grows with `M`, as the Poisson heuristic predicts (`~ log M/loglog M`). The
spectral quantity that controls the pair channel stays at about `2 M log M`. Even with the Poisson
tail, the real-weight pair channel is capped at about `pi * 1.9 M log M ~ 6 M log M`.

---------------------------------------------------------------------------------------------

## 3. Restricted lower bounds

Throughout, `E(n, Q)` is the least number of colliding pairs when `n` objects go into `Q` classes
(`two_thirds` Lemma 3.4). It satisfies:
* `n^2/(2Q) - n/2 <= E(n, Q) <= n^2/(2Q)`. The upper bound holds because, with `n = Qb + r`,
  `E = Q b(b-1)/2 + r b = n^2/(2Q) - Qb/2 - r^2/(2Q)`.
* `E(n_0, Q) + E(n_1, Q) <= E(n_0 + n_1, Q)`. Split an optimal distribution of `n_0 + n_1`
  objects class by class, and use `binom(a + b, 2) >= binom(a, 2) + binom(b, 2)`.

**Proposition B.1 (`K`-uniform margins).** Let a member of `A_X(M, C)` satisfy
`max_(x<y) mu_xy <= K mubar`, where `mubar := binom(M,2)^(-1) sum_(x<y) mu_xy` and `K >= 1`. Then

```text
W >= (2(M-1)/K) (D_max - c_C),    D_max := max_(x<y) D_xy.
```

*Proof.* By Lemma 1.4(4), `mubar <= (M W/4)/binom(M,2) = W/(2(M-1))`. By Lemma 1.3,
`D_max - c_C <= max mu <= K mubar`. QED

*Remarks.*
1. *What the hypothesis depends on.* `mu_xy = d_xy - W/2` depends only on the exponent profile and
   the weights. So `K`-uniformity is a property of the profile.
2. *The round-5/6 Paley profiles satisfy it.* They have `b` copies of each label and weights in
   `[tau, tau + 1/r]`. For a pair, `tau_j = c^T S_j in {0, +-2}`, and
   `2B - r = b +- (flip terms)` lies in `[b - 4, b + 4]` (`sharp.md` Prop. 4.3 and its proof). The
   slack formula `s_c = (1/2) sum_j w_j (|tau_j| - 1)` then gives
   `mu_xy in [(b-4) tau/2 - 1/2, (b+4) tau/2 + 1/2]`. So they are `K`-uniform with
   `K <= ((b+4) tau + 1)/((b-4) tau - 1)`, which is `37/29 + o(1)` for `b = 33`.
3. *Alone, B.1 is weak.* Unconditionally one only knows `D_max >= (2 - o(1)) log M` (Remark B.4),
   so B.1 gives only `W >= (4/K - o(1)) M log M`.

**Proposition B.2 (column-balanced profiles).** For every member of `A_X(M, C)` and every `x`,

```text
W/2 + (1/2) sum_s lambda^prof_s |sigma(s)|  >=  D_x - (M-1) c_C .
```

If every threshold cut is balanced (`sigma(s) = 0`, `M` even), then `sum_(y != x) mu_xy = W/2`
exactly, and `W >= 2 max_x D_x - 2(M-1) c_C`.

*Proof.* Lemma 1.4(3), the bound `|s_x sigma(s)| <= |sigma(s)|`, and Lemma 1.3. QED

Balanced profiles convert the maximal point total into `W`. Uniform-margin profiles convert the
maximal pair demand into `W`. For random-like demands, `max_x D_x ~ 2 M log X` and
`D_max ~ 4 log^2 M/loglog M` (§2 table; Lemma C.3). So B.2 gives `W >~ 4 M log X`, while B.1 gives
`W >~ 8 M log^2 M/(K loglog M)`. Only uniform margins can give growth through pairs.

**Proposition B.3 (parity-refined pigeonhole: forced Legendre parities).** Let a member of
`F^char_X(M, C)` be given; this covers P3-all, P3-gen and actual clusters by Prop. 1.2. Let:
* `q <= X` be an odd prime not dividing `N`, and `Q = q^a` with `1 <= a <= A_q`;
* `m = m_Q = (q - chi(q)) q^(a-1)` and `sigma_q(x) = sum_j a_xj nu_j(q) + lambda_q k_x mod 2`,
  where `nu_j(q) = [(p_j/q) = -1]` and `lambda_q = [q = +-3 mod 8]` (`sharp.md` (2.1));
* `n_i = #{x : sigma_q(x) = i}`, and `beta_q = (n_0 - n_1)/M`.

Then:
1. a unit-1 collision of `(x, y)` modulo `Q` implies `sigma_q(x) = sigma_q(y)`;
2. `#{x < y : unit-1 collision mod Q} >= E(n_0, m/2) + E(n_1, m/2)`;
3. `E(M, m) <= E(n_0, m/2) + E(n_1, m/2) <= E(M, m/2) <= M^2/m`. Also
   `E(n_0, m/2) + E(n_1, m/2) >= (M^2/(2m))(1 + beta_q^2) - M/2`;
4. summing over the pairs, every member satisfies

   ```text
   W >= (8/M) sum_(Q = q^a <= X, q odd, q not dividing N) log q [E(n_0(q), m_Q/2) + E(n_1(q), m_Q/2)] - 2(M-1) c_C .
   ```

5. The gain of item 4 over the parity-free count `E(M, m_Q)` is at most
   `4 M sum_(m_Q < 2M) log q/m_Q + 4 sum_(m_Q < 2M) log q = O(M log M)`. So it is at most a constant
   factor, and it is positive only through parity bias `beta_q != 0`.

*Proof.*
1. *Parity of the classes.* By `sharp.md` Cor. 2.2(2),(4), at the top power the point classes are
   `kappa_x = <a_x, ell> - iota k_x`, with `ell_j = nu_j(q) mod 2` and `iota = lambda_q mod 2`.
   So `kappa_x = sigma_q(x) mod 2`. The class at level `a` is `kappa_x mod m_Q`, and `m_Q` is even,
   so its parity is still `sigma_q(x)`. A unit-1 collision mod `Q` means equal classes mod `m_Q`
   (Cor. 2.2(3) at level `a`, via the reduction of Prop. 2.1(6)).
2. *Pigeonhole.* Points of parity `i` lie in the `m/2` classes of that parity. Apply `two_thirds`
   Lemma 3.4 to each parity.
3. *The comparison of counts.* The parity-constrained distribution is one distribution into `m`
   classes, which gives the first inequality. The second is superadditivity of `E` (above), and
   `E(n, Q) <= n^2/(2Q)` gives `E(M, m/2) <= M^2/m`. The last bound is `E(n, Q) >= n^2/(2Q) - n/2`
   together with `n_0^2 + n_1^2 = M^2 (1 + beta^2)/2`.
4. *Summation.* `D_xy = 2 sum_q a_q(xy) log q`, where `a_q(xy)` is the number of levels `a` with a
   unit-1 collision mod `q^a`, since collisions are nested. So
   `sum_(x<y) D_xy = 2 sum_Q log q #{pairs colliding mod Q}`. Dropping primes dividing `N` lowers
   the right side. Lemma 1.3 and Lemma 1.4(4) give
   `sum_(x<y) D_xy - binom(M,2) c_C <= sum mu <= M W/4`.
5. *The gain.* For each `Q`, `E(n_0, m/2) + E(n_1, m/2) - E(M, m) <= M^2/m - (M^2/(2m) - M/2)`, and
   this vanishes when `M <= m/2`. Multiply by `8 log q/M`. QED

*Exact check* (`parity_check.py`). It uses actual lattice points of squarefree `N` with 5, 6 and 7
split primes, random and arc-shaped clusters of 40-80 points, and all odd `q <= 60, 80, 90` with
their powers.
* It verifies `t^(m_q/2) = (-1)^(sigma(x) - sigma(y))` for `t = z_x conj(z_y)/N mod q` in 790,910
  pairs.
* Every one of the 48,138 observed collisions joins equal parities.
* Every collision count is `>= E(n_0, m/2) + E(n_1, m/2)`.
* The `E`-inequalities of item 3 hold exhaustively for `Q <= 12`, `n_0, n_1 < 40`.

ALL PARITY CHECKS PASSED.

**Remark B.4 (pigeonhole over many moduli at once).** Let `S` be a set of odd primes not dividing
`N`.
* *Without parity.* If `M > prod_(q in S) m_q`, two points have equal classes at every `q in S`
  (CRT), so `D_max >= 2 sum_(q in S) log q`. Taking `S` the odd primes up to `y`, with
  `prod m_q < M`, gives `D_max >= (2 - o(1)) log M`.
* *With parity.* Collisions at all of `S` occur only inside a class of the parity vector
  `(sigma_q(x))_(q in S)`. Each such class meets at most `prod_(q in S) m_q/2` residue classes.
  * If the largest parity class has `alpha M` points, a common collision is forced as soon as
    `alpha M > prod (m_q/2)`.
  * Since `alpha >= 2^(-|S|)`, the largest forced CRT modulus grows by a factor `alpha 2^(|S|) <= 2^(|S|)`.
  * So the bound on `D_max` rises by at most `2|S| log 2`. That is `O(log M/loglog M)` for sets
    with `prod (m_q/2) < M`.

So pigeonhole over many moduli, with the forced parities, gives `D_max = (2 + o(1)) log M` at best.
That is far below the `log^2 M/loglog M` of the Poisson heuristic, and the round-6 design `D_A`
indeed has `D_max = O(log M)` while respecting parity.

**Proposition B.5 (monomial heights plus pigeonhole).** Every member of `F^all_X(M, C)` satisfies

```text
(2r - 1) sum_j w_j  >=  sum_(q <= X, q odd, q not dividing N) log q E(2r, m_q),
```

and hence, if no `p_j` is `<= min(r, X)`, `W >= sum_j w_j >= r log min(r, X) - O(r)`. But distinct
primes `= 1 mod 4` already give
`sum_j w_j >= sum_(i <= r) log p^((i)) = (1 + o(1)) r log r`, where `p^((i))` is the `i`-th prime
`= 1 mod 4` (prime number theorem in progressions). So this pigeonhole adds nothing.

*Proof.*
* *The collisions.* For odd `q <= X` with `q` not dividing `N`, let `t_j = r_j/conj(r_j) mod q`
  in `T_q`, which is cyclic of order `m_q`. Consider the `2r` elements `t_j^eps`,
  `(j, eps) in [r] x {+-1}`. By pigeonhole at least `E(2r, m_q)` unordered pairs
  `P = {(j,eps), (l,eps')}` have `t_j^eps = t_l^(eps')`. For such a pair,
  `v_P = eps e_j - eps' e_l != 0` collides with unit 1 at `q`.
* *The heights.* The height condition (`sharp.md` Prop. 6.1, `s = 0`) for `v_P` gives
  `sum_q a_(q,0)(v_P) log q <= w(v_P)/2`. Sum over all unordered pairs `P`.
* *Evaluating the left side.* For `j != l` there are four sign patterns, each with
  `w(v_P) = w_j + w_l`. For `j = l` there is one pair, with `w(v_P) = 2 w_j`. So
  `sum_P w(v_P)/2 = 2(r-1) sum_j w_j + sum_j w_j`.
* *The final bound.* Use `E(2r, m) >= 2r^2/m - r`, keep only `q <= Y := min(r, X)`, and use
  `sum_(q <= Y) log q/(q+1) >= log Y - O(1)` (Mertens) and `theta(Y) <= Y log 4 <= r log 4`. QED

The analogous count with monomials of norm `<= k` forces only `max_j w_j >= log(2r/k) - O(1)`.
Pigeonhole consequences of the heights are subsumed by the trivial count of primes. A lower bound
in P3-all beyond `two_thirds` would have to come from the completion/covering of the slabs
(`sharp.md` §6.1), for which no necessary condition beyond the heights is known.

---------------------------------------------------------------------------------------------

## 4. Conditional growth for uniform-margin clusters

Fix `0 < eta < 1`. For a member of `A_X(M, C)` with `X >= M` let

```text
Q_eta := { q prime : M^(1-eta) < q <= M,  q odd,  q not dividing N },     lambda_Q := sum_(q in Q_eta) 1/m_q ,
c_xy  := #{ q in Q_eta : unit-1 collision of (x, y) mod q },
F_k   := sum_(x<y) binom(c_xy, k),         E_k := binom(M, 2) e_k( (1/m_q)_(q in Q_eta) ),
```

where `e_k` is the `k`-th elementary symmetric polynomial. Note `D_xy >= 2 (1 - eta) c_xy log M`.

**Hypothesis H(eta, delta)** (Poisson-type collision statistics, `0 < delta <= 1`). For every integer
`k >= 1` with `delta^k E_k >= 1`, one has `F_k >= delta^k E_k`.

`H` concerns only the residue data at primes in `(M^(1-eta), M]`.
* For actual clusters it is a statement about joint divisibility of the pair integers (the reduced
  chords `c_xy` of `two_thirds`) by primes of that window. Equivalently it concerns the smooth parts
  of `Im(z_x conj(z_y) u)/|g_xy|^2`-type integers: the hypothesis the task asked for.
* For independent uniform residues, `E[F_k] = E_k` exactly (Lemma C.3).

**Lemma C.1.** For `k >= 1` with `k p_max < lambda_Q`, where `p_max = max_(q in Q_eta) 1/m_q`,
we have `e_k >= (lambda_Q - k p_max)^k/k!`. Also `e_k <= lambda_Q^k/k!`.

*Proof.* `k! e_k` is the sum over ordered `k`-tuples of distinct primes of the products. Choosing
them in turn, each inner sum is `>= lambda_Q - i p_max`. The upper bound is
`lambda_Q^k = (sum p)^k >= k! e_k`. QED

If no prime factor of `N` lies in `(M^(1-eta), M]`, then `lambda_Q = log(1/(1-eta)) + o(1)`. This
follows from Mertens' theorem, since `1/m_q >= 1/q - 1/q^2`. It holds in all round-5/6 classes,
where `p_j > X >= M`.

**Theorem C (uniform margins plus H give growth).** Fix `eta in (0,1)`, `delta in (0,1]`,
`lambda_0 > 0`, `K >= 1`, `C > 0` and `eps in (0,1)`. There is `M_0` such that the following holds
for every `M >= M_0`. Let a member of `A_X(M, C)` with `X >= M` satisfy:
* `lambda_Q >= lambda_0`;
* Hypothesis `H(eta, delta)`;
* `K`-uniform margins (Prop. B.1).

Then

```text
W  >=  ((8(1 - eta) - eps)/K) M log^2 M / loglog M .
```

This applies in particular to actual clusters, members of P3-gen and members of P3-all with these
properties.

*Proof.*
1. *Choice of `k`.* Put `L = log M`, `l = log L` and `k := floor((2 - eps/4) L/l)`. For large `M`,
   `k p_max <= L/(M^(1-eta) - 1) <= lambda_0/2`.
2. *`delta^k E_k >= 1`.* By Lemma C.1, `binom(M,2) >= M^2/4` and `k! <= k^k`,

   ```text
   log(delta^k E_k) >= 2L - log 4 - k log k + k log(delta lambda_0/2).
   ```

   Here `k log k <= (2 - eps/4)(L/l)(l + log 2) = (2 - eps/4) L + O(L/l)`, and
   `k log(delta lambda_0/2) = -O(L/l)`. So `log(delta^k E_k) >= (eps/4) L - O(L/l) -> infinity`.
3. *A heavy pair.* `H` gives `F_k >= delta^k E_k >= 1 > 0`, so some pair has `c_xy >= k`. Hence
   `D_max >= 2(1-eta) k L >= (4 - eps/2)(1-eta) L^2/l - 2(1-eta) L`.
4. *Conclusion.* Prop. B.1 gives `W >= (2(M-1)/K)(D_max - c_C)`. Since
   `2(4 - eps/2)(1 - eta) > 8(1 - eta) - eps`, the losses `(M-1)/M`, `O(L)` and `c_C` are absorbed for
   large `M`. QED

**Corollary C.4.** Under the hypotheses of Theorem C, with `M -> infinity` and `C` fixed, an actual
cluster on an arc of length `C sqrt R` has

```text
M <= (K/(4(1 - eta)) + o(1)) log R logloglog R / (loglog R)^2 .
```

*Proof.* Put `c = (4(1-eta) - eps/2)/K`, so that `log R = W/2 >= c M log^2 M/loglog M`.
* Fix `theta < 1`. Either `M <= (log R)^theta`, or `log M >= theta loglog R`.
* In the second case, `loglog M <= logloglog R` (since `M <= log R` by `two_thirds`), so
  `M <= (c theta^2)^(-1) log R logloglog R/(loglog R)^2`.
* Let `theta -> 1` and `eps -> 0`. QED

**Lemma C.2 (the round-6 design violates H).** Fix `eta in (0,1)` and `delta > 0`. For all large
`M`, the point classes of `D_A` (`sharp.md` Def. 3.2) violate `H(eta, delta)`. This holds for every
parity colouring and every random choice, and also for the random-fibre version used in Thm 6.4.

*Proof.*
* *Which `q` are involved.* For `q in Q_eta` (`7 < q <= M < 4M`, a small prime of the design), a
  pair `(x, y)` with `d = |x - y|` collides at level 1 only if `p'(q) | d` (`sharp.md` Lemma 3.3(1)).
  Here `p'(q) > q/8 > M^(1-eta)/8` (Lemma 3.1(1)).
* *Few large prime factors.* The distinct primes `> M^(1-eta)/8` dividing `d < M` number fewer than
  `log M/((1-eta) log M - log 8) <= 2/(1 - eta)` for large `M`.
* *Bounded collision counts.* Each prime has at most 9 preimages `q` (Lemma 3.1(2)), so
  `c_xy <= c_0 := floor(18/(1 - eta))`.
* *Violation.* For `k = c_0 + 1`, `F_k = 0`. But `delta^k E_k >= delta^k (M^2/4)(lambda_Q/2)^k/k! -> infinity`,
  because `lambda_Q = log(1/(1-eta)) + o(1)` (all `p_j > X >= M`).
* *Random fibre.* The random fibre changes no point class (`sharp.md` Lemma 6.3(1)). QED

At small `M` the design has *more* window collisions than the random model, since it uses only
about `2p'(q) <= m_q` classes. In `stats_H.py`, `F_k/E_k` is about 2 at `k = 1` and reaches 31
(`M = 512`) and 64 (`M = 1024`) for `k` between 5 and 7. The violation is
asymptotic: it appears at `k = c_0 + 1` (37 for `eta = 1/2`), which needs `delta^k E_k >= 1`, i.e.
very large `M`.

**Lemma C.3 (independent residues satisfy H).** Let the classes `kappa^(q)_x` (`q in Q_eta`,
`x <= M`) be independent and uniform on `Z/m_q`. Suppose `lambda_Q in [lambda_0, Lambda_0]` with
`lambda_0 > 0`, and `delta <= 1/(4(1 + Lambda_0))`. Then `P(H(eta, delta)) -> 1` as `M -> infinity`.

*Proof.*
1. *Pairwise independence.* For distinct pairs `P != P'`, the indicator vectors
   `(1[P collides at q])_q` and `(1[P' collides at q])_q` are independent.
   * Different `q` are independent.
   * For one `q`: if `P, P'` are disjoint this is clear. If `P = {x, y}` and `P' = {x, z}`, then
     `P(kappa_x = kappa_y = kappa_z) = 1/m^2`, the product of the two probabilities. Two events with
     product probability are independent, and so are their complements.

   Hence the `c_P` are pairwise independent, and `Var(F_k) = sum_P Var binom(c_P, k)`.
2. *Moments.* `E binom(c, k) = e_k`, so `E F_k = E_k`. Write `binom(c,k)^2` as a sum over pairs
   `(S, S')` of `k`-sets and split by `T = S cap S'`. This gives
   `E binom(c,k)^2 <= sum_(i=0)^k e_i e_(k-i)^2 <= sum_i lambda^(2k-i)/(i! (k-i)!^2)`.
3. *Chebyshev.* `P(F_k < delta^k E_k) <= Var F_k/((1-delta)^2 E_k^2)`, and

   ```text
   Var F_k / E_k^2 <= (lambda/lambda')^(2k) binom(M,2)^(-1) sum_(i=0)^k binom(k,i)^2 i! lambda^(-i),     lambda' = lambda - k p_max .
   ```

4. *Large `k`.* For `k` with `delta^k E_k >= 1`, Lemma C.1 gives `binom(M,2)^(-1) <= delta^k lambda^k/k!`.
   Since `binom(k,i)^2 i!/k! = binom(k,i)/(k-i)! <= binom(k,i)`, the bound becomes
   `<= (lambda/lambda')^(2k) (delta (1 + lambda))^k <= 2 * 4^(-k)` once `k p_max` is small, which
   holds for `k <= log M`. Larger `k` have `E_k < 1` for large `M`.
5. *Small `k`.* For each fixed `k <= k_1` the bound is `O_k(M^(-2))`.
6. *Union bound.* `P(H fails) <= sum_(k <= k_1) O_k(M^(-2)) + 2 sum_(k > k_1) 4^(-k)/(1-delta)^2`.
   Let `M -> infinity`, then `k_1 -> infinity`. QED

*Evidence* (`stats_H.py`, window `(sqrt M, M]`):
* *Independent model.* `F_k/E_k` is within statistical error of 1. For example, at `M = 1024`:
  `k = 1..6`: `1.002, 1.000, 0.989, 0.958, 0.932, 0.931`.
* *Genuine residues of random Gaussian primes (P3-gen-type data).* The same, e.g. at `M = 64`:
  `0.98, 0.98, 1.02, 1.07`.
* *`D_A`.* `F_k/E_k` is about 2 at `k = 1` and reaches 31 (`M = 512`) and 64 (`M = 1024`) for
  `k` between 5 and 7.
  This is the small-`M` regime described after Lemma C.2.

*Remarks on Theorem C.*
1. *Both hypotheses are needed.* Without uniform margins, `H` gives nothing beyond Corollary A.1
   at the level of pair-cut certificates (Corollary A.2: in the real-weight relaxation the heavy
   pairs are absorbed). Without `H`, uniform margins give `W >= (4/K - o(1)) M log M` only.
2. *Is uniformity plausible for actual circles?* It is the profile type of every known
   `O(M log M)` fake (Hadamard cores with boundedly many flips per column). For *actual* Gaussian
   primes those profiles are already excluded at small flip mass by `round4/flipped.md`
   Theorem L, conditional on a Lang-Waldschmidt two-logarithm bound. Both conditional results act
   only on near-uniform margins. Neither touches profiles with non-uniform margins, and those are
   exactly the ones Theorem A leaves open.
3. *Why `H` is open for actual clusters.* Proving `H` would require joint equidistribution of the
   points modulo products of `k ~ (2 - eps) log M/loglog M` primes of size `~ M`. That is moduli of
   size `exp((2 - eps)(1 - eta) log^2 M/loglog M)`.
   * Suppose the primes of a P3-gen member were chosen at random, of size `P`. Equidistribution
     of Gaussian primes in residue classes modulo `Q` (Hecke characters) is known only for `Q` up
     to a power of `log P` (Siegel-Walfisz range), or for `Q <= P^(1/2 - eps)` under GRH.
   * Even under GRH, reaching these moduli needs `log P >~ 2 k log M ~ log^2 M/loglog M`, i.e. `W`
     already of the size Theorem C wants to prove. So, even for random primes, `H` cannot be
     bootstrapped from the size of the primes. (Heuristic remark; nothing depends on it.)

---------------------------------------------------------------------------------------------

## 5. Discreteness: prime weights and few columns

Theorem A is about the **real-weight relaxation**. In every member, and in actual clusters, the
cuts come with weights `w_j = log p_j` of *distinct* primes `p_j = 1 mod 4`. This has two
consequences.

**Lemma D.0 (column count).** A profile with `r` distinct primes has
`W >= sum_j log p_j >= sum_(i <= r) log p^((i)) = (1 + o(1)) r log r`, where `p^((i))` is the
`i`-th prime `= 1 mod 4`. If all `p_j > X` (the round-5/6 classes), then also `W >= r log X`.
So `W = O(M log M)` forces `r = O(M)` columns, of average weight `>= log r`.

*Proof.* Distinctness, then the prime number theorem in progressions (`theta(x; 4, 1) ~ x/2`,
`p^((r)) ~ 2 r log r`). QED

The real-weight optima of Theorem A are LP vertices. In every experiment (`lp_gw.py`) they use about
`binom(M, 2)` cuts, of average weight `~ W/binom(M,2) ~ log M/M`, far below `log 5`. So Theorem A
does not show that the *prime-weight* pair system is solvable at `O(M log M)`. The pair channel
with prime weights is governed by the **sparse cut problem**

```text
W >= min_r max( f_D(r), r log X ),     f_D(r) := min { W(lambda) : mu(lambda) >= D - c_C, |supp lambda| <= r }.
```

The extreme case is equal weights: primes in `[P, P(1 + 1/r)]`, weights `tau + O(1/r)`. There the
pair condition becomes an integer condition on a `+-1` matrix:
`<s_x, s_y> <= -k_xy` with `k_xy = ceil(2(D_xy - c_C)/tau)` (up to `O(1/tau)`).

**Proposition D.1 (equal weights: the Perron bound is attained on clique-partition patterns).**
1. If `S in {+-1}^(M x r)` satisfies `<s_x, s_y> <= -k_xy` for all `x != y`, with `k` symmetric and
   nonnegative, then `r >= lambda_max(k)`.
2. Let `M = g t`, with points `(a, i)`, `a <= g` (groups) and `i <= t` (members). Let `L >= 1` and
   `m := g - 1 + L`. Suppose Hadamard matrices `H_g`, `C_t` (normalised: first column all ones) and
   `G_m` exist; `m >= g` holds automatically. The rows

   ```text
   s_(a,i) := ( H_g[a, 2..g] ,  C_t[i, 2..t] (x) G_m[a, 1..m] )
   ```

   have `r = M - t + (t-1) L` nonconstant columns. Their correlations are `<s_(a,i), s_(a,j)> = -L`
   for `i != j` and `<s_(a,i), s_(b,j)> = -1` for `a != b`.
3. This `r` equals `lambda_max(k)` for `k = (J - I) + (L - 1)[same group]`. Uniform margins would
   need `(M-1) L` columns.

*Proof.*
1. Let `u >= 0` be the unit Perron vector of `k`. Then
   `0 <= |S^T u|^2 = sum_(x,y) u_x u_y <s_x, s_y> <= r - u^T k u = r - lambda_max(k)`.
2. Normalised Hadamard rows satisfy `<H[a, 2..], H[b, 2..]> = g [a = b] - 1`, likewise for `C_t`,
   and `<G[a], G[b]> = m [a = b]`. Hence:
   * same group, `i != j`: `(g - 1) + (-1) m = -L`;
   * different groups: `-1 + (<C_t rows>)(0) = -1`;
   * diagonal: `(g - 1) + (t - 1) m = r`.

   *Nonconstant columns.* All `g` rows of `H_g` occur, and each `C_t` column has both signs.
3. `J` and the block matrix commute. On `1` the eigenvalue is `M - 1 + (L-1)(t-1) = r`. On
   block-constant vectors orthogonal to `1` it is `-1 + (L-1)(t-1)`, and on vectors with zero block
   sums it is `-1 - (L - 1)`. QED

*Exact check* (`clique_blocks.py`, Sylvester orders). Twelve cases from `M = 4` to `M = 128`, `t = 2`
(matchings) to `t = 16`, and `L` up to 97. All Gram matrices are as claimed, and `r = lambda_max(k)`
in every case. Example: `M = 128`, `t = 4`, `L = 97` needs `r = 415` columns, against `12,319` for
uniform margins. ALL CLIQUE-BLOCK CHECKS PASSED.

So discreteness does **not** force uniform margins for structured patterns. A single modulus
produces exactly a clique-partition pattern: the residue classes are the groups.

What is not covered is the **sum** over all moduli `q <= X` of such patterns, each with a small
weight `2 log q/tau < 1` per pair. Concatenating one block per modulus costs `>= M - t` columns
per modulus, i.e. `~ pi(X) M`, which is far too many. The demand matrices of actual clusters are
such sums: `D = sum_(q^a <= X) 2 log q (C_(q^a) - I)`, where `C_Q` is the class-equality matrix
mod `Q`.

*Evidence on the sparse and equal-weight problems* (small `M`, exact where stated):

* `sparse_lp.py`, Part 2 (support pruning, upper bounds on `f_D(r)`, random residues at
  `X = 3M`):

  | `M` | `f_D(r)/W_cut` at `r = 2M` | at `r = 3M` | at `r = 4M` | `2(M-1) D_max / W_cut` (uniform) |
  |---|---|---|---|---|
  | 10 | 1.60-2.02 | 1.06-1.11 | 1.00-1.01 | 2.72-3.65 |
  | 12 | 3.59-6.06 (one trial: none found) | 1.16-1.23 | 1.01-1.02 | 3.00-3.30 |
  | 16 | none found | 2.12-2.52 | 1.24-1.31 | 4.05-4.11 |

  The heuristic is crude: greedy removal stalls at about `2.4 M` cuts. Still, `3M-4M` cuts already
  come within a factor 1.0-1.3 of the real-weight optimum.
* `sparse_lp.py`, Part 1 (`M = 8`, MILP with a 240 s limit; upper bounds on `f_D(r)`). Two demand
  matrices, with `W_cut = 59.6, 63.1` and uniform cost `149.6, 197.9`. Results as
  `f_D(r)/W_cut`:
  * `r = 7 = M - 1`: `<= 1.71, 1.82`;
  * `r = 8`: `<= 1.52, 1.74`;
  * `r = 10`: `<= 1.28, 1.41`;
  * `r = 14`: `<= 1.07, 1.25`.

  So even with the minimum number `M - 1` of cuts, adapted margins cost at most 0.58-0.68 of
  uniform margins.
* `equal_weight_ip.py` (**exact integer optimum** `r*` of the equal-weight problem):

  | `M` | `tau` | `r*/LP` (= `W_equal/W_cut`) | `r_unif/r*` (uniform margins vs adapted, equal weights) |
  |---|---|---|---|
  | 8 | `log M` | 1.35, 1.16, 1.19 | 2.58, 2.84, 3.50 |
  | 8 | `2 log M` | 1.58, 1.37, 1.74 | 3.50, 2.69, 1.75 |
  | 10 | `log M` | 1.27, 1.28 | 3.98, 2.70 |
  | 10 | `2 log M` | `<= 1.68` (time limit; `>= 1.57`), 1.52 | 2.03, 2.57 |
  | 12 | `2 log M` | not solved in 400 s: `29 <= r* <= 66`, i.e. `r*/LP in [1.22, 2.77]` (the uniform solution, 66 columns, is feasible) | `>= 1` (no information) |

  Here `r_unif` is the exact least number of equal-weight columns giving *every* pair margin
  `>= D_max`. All values are exact integer optima except the two marked (`M = 10`, second
  `tau`, first trial; `M = 12`). Output files: `equal_weight_ip_out.txt` (`M = 8, 10`) and
  `equal_weight_ip_M12_out.txt` (`python3 equal_weight_ip.py 400`).

**Question Q (prescribed-correlation sign matrices).** Let `D` be the demand matrix of independent
random residues at odd prime powers `<= 3M`, `tau = kappa log M` and `k = max(1, ceil(2D/tau))`.
Does there exist, with high probability, `S in {+-1}^(M x r)` with `r = O_kappa(M)` and
`<s_x, s_y> <= -k_xy` for all pairs? More generally: is `min_r max(f_D(r), r log X) = O(M log M)`
for such `D`?

Here is what each answer would mean.
* *Yes, for the demand matrices of actual clusters.* Then no pair-level argument proves growth,
  even with prime weights. The heavy-pair phenomenon would be absorbable discretely as well.
* *No, in the strong form "`r >= c (M-1) k_max` for every realisation" for `H`-type patterns.* Then,
  under `H`, every equal-weight actual cluster has `W >= c' M log^2 M/loglog M`. For general prime
  weights a corresponding sparse statement would be needed.

The equal-weight data above point to "yes": `r* <= 1.75 x LP` while uniform margins cost up to
`4 x r*`. They are small-`M` evidence only. At these sizes `D_max/tau` is only 2-10, so the
asymptotic question is untouched.

---------------------------------------------------------------------------------------------

## 6. Consistency with the round-5/6 constructions

| statement | applies to the constructions? | outcome |
|---|---|---|
| Cor. A.1, B.2, B.3, B.5 | yes: valid inequalities for all members of the stated classes | satisfied automatically, since the fakes are members. Quantitatively: Cor. A.1 says `2 lambda_max(D) <= W` for the `D_A` demands. The fakes have `W >= 33(M-1) tau >= (66 - o(1)) M log M` (`sharp.md` Thm 5.1). Lemma A.4 bounds `2 lambda_max(D) <= 72 M log M + 222 M`, which is only an upper estimate. Measured: `2 lambda_max(D) ~ 6.4 M log M` for `D_A` (§2 table) |
| Prop. B.1 | the Paley fakes are `K`-uniform with `K <= ((b+4) tau + 1)/((b-4) tau - 1)` | gives `W >= (2(M-1)/K)(D_max - c_C)` with `D_max <= (10 + o(1)) log M`, i.e. `<= (20 + o(1)) M log M/K`. The fakes have `W >= 33 (M-1) tau >= 66 (M-1) log M` |
| Theorem C | its hypothesis `H` fails for `D_A`, for all large `M` (Lemma C.2) | not applicable, consistent |
| Cor. A.2 / A.5 | the `D_A` members of P3-all (Thm 6.4) | pair-cut certificates give them at most `O(M log M)`. Their actual `W ~ 66 A M log^2 M/loglog M` comes from monomial heights (`sharp.md` Prop. 6.5), not from pairs |
| `lower.md` Thm 5.1 (free residues, auxiliary primes) | pair moduli `m_xy <= Y^4/4` (`lower.md` Lemma 4.1), so `D_xy <= 8 log Y = O(log M)` | `lambda_max(D) <= (M-1) max D_xy = O(M log M)`, consistent with `W = Theta(M log M)` |
| `sharp.md` §6.6 heuristic for coherent designs ("genuine-random pair collisions force `(b-4) tau/2 >~ log^2 M/loglog M`") | -- | **refined.** It is correct for uniform-margin (Paley) profiles, by Prop. B.1. By Cor. A.2 it is not an obstruction at the real-weight level for profiles adapted to the demand pattern, which need only `~ pi lambda_max(D) ~ 2 pi M log X` there. Whether discrete prime weights allow the adaptation is Question Q |
| `two_thirds` (`3 M log M - O(M) <= W`) | contained in Cor. A.1 plus the pigeonhole (`u = 1/sqrt M`) together with the chain lemma at `p | N` | consistent; the constant 3 is not improved here in general. In classes with all `p_j > X >= 2M` (no `p | N <= X`), B.3(4) gives `W >= (4 - o(1)) M log M`, and up to `(8 - o(1)) M log M` under full parity bias |

*Sanity check of B.3 against `D_A`.* `D_A` places classes inside parity classes, so B.3(2) holds
for it by construction. It uses only `2 p'(q) <= m_q` classes at small `q`, so its point totals
exceed the parity-free pigeonhole average: measured `~ 3 M log M` against `~ 2 M log M` for
independent residues (§2 table). That is the direction B.3 predicts, and it is far below the
fakes' `W >= 66 M log M`.

---------------------------------------------------------------------------------------------

## 7. What a growth proof must do; open problems

For actual circles, a proof of `M = o(log R/loglog R)` from residue data that factor through prime
residues must do at least one of the following. Each item names the relevant results.

1. **Force irregular point totals:** `max_x D_x >> M log M` (Cor. A.3).
   * Pigeonhole forces only averages.
   * Parity raises averages by at most a factor 2 (Prop. B.3).
   * Random, `D_A` and genuine-residue data all have `max_x D_x ~ (2-3.5) M log M` (evidence, §2).

   So this route needs new arithmetic.
2. **Force near-uniform margins** (Prop. B.1) **and prove `H`** (Theorem C).
   * Uniform margins are what all known fakes use. `round4/flipped.md` Theorem L (conditional on
     Lang-Waldschmidt) already penalises actual Hadamard-core clusters.
   * But no valid input forces uniformity. Theorem A shows that the pair inputs themselves allow
     adapted margins at the real-weight level.
3. **Use discreteness:** prime weights and `r = O(M)` columns (Question Q, Lemma D.0).
   * This is the only pair-level channel not closed here.
   * Structured patterns are realisable at the Perron bound (Prop. D.1).
   * Small-`M` integer optima are within 1.2-1.8 of the LP.
4. **Use monomial completion** (P3-all\*, `sharp.md` §6.7).
   * Pair-cut certificates cannot decide P3-all\* (Cor. A.5).
   * Height pigeonhole is trivial (Prop. B.5).
   * A covering statement for the slabs `{phi : |sin(<v,phi> - s pi/4)| < nu_s m_(v,s) e^(-w(v)/2)}`,
     `v` outside `L`, would be needed. No necessary condition beyond the heights is known.

**Open problems (precise).**
* **(Q)** Prescribed-correlation sign matrices, §5.
* **(H-actual)** Does every normalised actual cluster satisfy `H(eta, delta)` for some fixed
  `eta, delta`? (Joint divisibility statistics of the reduced chords by primes in `(M^(1-eta), M]`.)
* **(U)** Is there any valid input that forces `max_xy mu_xy <= K mubar` for actual clusters, or
  that forces `lambda_max(D) >> M log M`?
* **(P3-all\*)** Unchanged: is `W_min^all = O(M log M)`?
  * The bracket is still `3 M log M - O(M) <= W_min^all <= (66 + o(1)) M log^2 M/loglog M` at
    `X = 3M`.
  * Cor. A.5 shows that the upper bound's extra `log M/loglog M` comes from monomial heights alone,
    and that no pair argument can confirm it as a lower bound.

---------------------------------------------------------------------------------------------

## 8. Evidence and files (`round7/p3lower_checks/`)

| script | what it checks | type | result |
|---|---|---|---|
| `lp_gw.py` | exact cut LP `W_cut(D)` for `M <= 11` (48 demand matrices); sandwich `2 lambda <= W_cut <= pi lambda`; LP duality; cut identity at the optimum; uniform cost | floating-point LP (evidence) | ALL LP/GW SANDWICH AND DUALITY CHECKS PASSED; `W_cut/(2 lambda) in [1.000, 1.224]`; uniform cost/`W_cut` = 1.4-7.1 on non-uniform `D`; optimal supports of size `binom(M,2)` |
| `clique_blocks.py` | Prop. D.1 constructions, 12 cases up to `M = 128` | exact integer | ALL CLIQUE-BLOCK CHECKS PASSED; `r = lambda_max(k)` in all cases |
| `parity_check.py` | Prop. B.3 on actual lattice points (Euler criterion for 790,910 pairs; 48,138 collisions; collision counts vs parity-refined pigeonhole); `E`-inequalities | exact integer | ALL PARITY CHECKS PASSED |
| `stats_H.py` | `D_max`, `lambda_max(D)`, point totals, window factorial moments `F_k/E_k` for independent, `D_A` and genuine residues, `M` up to 1024 | evidence | see §2 table and §4 |
| `sparse_lp.py` | sparse cut LP: MILP at `M = 8` (time-limited; `python3 sparse_lp.py 1`), support pruning at `M = 10, 12, 16` (`python3 sparse_lp.py 2`) | evidence | §5; outputs `sparse_lp_part1_out.txt`, `sparse_lp_part2_out.txt` |
| `equal_weight_ip.py` | exact equal-weight integer optimum `r*` vs LP, Perron and uniform, `M = 8, 10, 12` | exact MILP (evidence) | §5 |
| `anneal.c`, `discrete_realise.py` | simulated annealing for prescribed-correlation sign matrices | evidence | **discarded**: on the uniform calibration pattern (`M = 32`, `k = 1`, true optimum 31) annealing needed `r = 141`, so its outputs carry no information; it is kept only for the record |

All exact checks use integer arithmetic. No proof in §§1-4 depends on any script. The repository
`/home/user/Jarnik` and the folders `round5/`, `round6/` were not modified.
