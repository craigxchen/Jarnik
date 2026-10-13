# P3-all* construction: structured fibres (round 7, draft, unrefereed)

## 0. Verdict

**The target was not reached.** No fake circle with `W = O(M log M)` in P3-all\* was found. The best
rate stays `sharp.md` Theorem 6.4: `W = O(M log^2 M/loglog M)`. No growth improvement for actual
circles is claimed.

What this round proves instead are exact statements about *construction methods* for
P3-all\*. Each names the obstacle it removes or isolates. Together with two labelled heuristics
they say that every structured fibre we could analyse stalls at the same rate, for three
independent reasons.

| # | statement | class / scope | status |
|---|---|---|---|
| 1.4 | twin coordinates: monomial label `= <u, ltilde> + <c', kappa>`; the twin labels are the solved coordinates `t_x = eps_x(kappa_x - P_x)` | profile of `sharp` §4, any data | proved (identity) |
| **Thm 1** | **Poisson barrier, robust form.** If on a set of primes in `[M^(1-eta), X]` of reciprocal mass `>= (log M)^(-C)` the labels of `>= eps M` columns are conditionally independent and uniform in their parity class, then every completion has `W >= (eps (1-eta)^2/(2(C+1)) - o(1)) M log^2 M/loglog M`, whatever the data elsewhere | random designs in `F^all_X`, any profile and weights, `X <= M^(A_0)` | **proved** |
| Cor. 1 | hybrids (structured or coherent except on such a set, e.g. the classes `q = +-1 mod 2^K` with `2^K <= (log M)^C`) stay at `Theta(M log^2 M/loglog M)`; only a random part of mass `M^(-Omega(1))` is harmless | as Thm 1 | proved |
| Lemma 3.2, Prop. 3.3 | **root-of-unity parity repair.** At every `q` the forced parity is repaired coherently by an element of order `2^(v_2(m_q))`; collisions are then certified: `q^a | Im(gamma_v^(2^k))`. For `q != +-1 mod 16` (`k <= 3`) this costs `12 log|gamma_v|`, so parity is *no* obstruction to coherence there | factored data | proved; exact check on all 89 odd prime powers `<= 400` |
| §3.2 | every odd-label element has order divisible by `2^k`; non-root repairs make the certificate quadratic in `||v||_1`; matching parities instead costs `2^K >~ X/(log X log p_j)` | factored data | proved (exact statements) |
| H1 | parity balance on the classes `q = +-1 mod 16`: certified cost `2^K` against random cost `log M/K` gives `Theta(log M/loglog M)` per unit weight | -- | heuristic |
| **Thm 2** | **coherent designs on the Paley profile cannot certify pairs and columns at once.** If column heights are `<= C` and the average pair height is `Hbar`, then `W >= (M-1)(theta(X;8,3) - 2 log 2)/(8 ceil(8(C + 2Hbar)))` (under a mild size condition). So at `W = O(M log M)`: `C + 2 Hbar >= M^(1/2 - o(1))` | G-coherent designs on `sharp` §4 profile, `X >= M` | **proved** (uncertainty for Hadamard matrices + pigeonhole on label sums) |
| H2 | pair integers of size `exp(M^(Omega(1)))` with Dickman-like smooth parts give pair demands up to `~2 log^2 M/loglog M`, i.e. the same rate under uniform margins | -- | heuristic |
| §5 | single-column demands of three designs (`M <= 2048`) and pair demands of the coherent design (`M <= 2040`) | -- | evidence (does not separate `log M` from `log M/loglog M` at these sizes) |

**Answer to the task's suggestion.** "Prime residues equal to residues of fixed auxiliary Gaussian
integers on most columns, with only the twin/repair columns designed":
* the coherent part is compatible with the forced Legendre parity on the classes `q != +-1 mod 16`
  after the root repair (Lemma 3.2), at constant cost;
* on the classes `q = +-1 mod 16` no cheap repair exists on present methods (§3.2, H1);
* the designed twin columns are the solved coordinates (1.4). If they are random, Theorem 1
  applies to them. If they are coherent too, Theorem 2 applies: the pairs then have huge height and
  are uncertified (H2).

So the obstruction is not only the parity: it is the pair/twin tension of `sharp` §6.6, now proved in
two forms (Theorems 1 and 2), plus the parity on one quarter of the primes.

**Validity for actual circles.** Nothing here is an inequality for actual circles; all statements
concern members of `F^all_X` produced by stated design methods. Each uses only conditions that
every member satisfies (the per-unit height condition (1.2), `sharp` Prop. 6.1 as corrected), and the
classes are those of `sharp` Def. 1.1.

## 1. Setting and notation

Everything refers to `round5/round6_sharp/sharp.md` (cited as `sharp`). The class is
`F^all_X(M, C)` (P3-all) of `sharp` Def. 1.1: genuine primes `p_j = 1 mod 4`, a genuine profile
`(a_x)`, real angles, factored residue data `r_j(q) mod q^(A_q)` with `Norm r_j(q) = p_j`, the
character conditions (1.1) and the per-unit monomial conditions (1.2) at every odd prime power
`q^a <= X`. Actual clusters are members for every `X` (`sharp` Prop. 1.2). The open question
**P3-all\*** is whether `W_min^all_(3M)(M, C) = O(M log M)`. The best construction is `sharp`
Thm 6.4, `W = O(M log^2 M / loglog M)`.

**Labels.** For an odd prime `q` not dividing `N`, `Q = q^(A_q)`, `T_Q` is cyclic of order
`m_Q = (q - chi(q)) q^(A_q - 1)` (`chi(q) = (-1/q)`), and we fix a generator `g`. The *label* of a
column is `ell_j = log_g(r_j/conj r_j) in Z/m_Q`; `i = g^iota`. `sharp` Cor. 2.2: the data are
exactly label vectors with the forced parity `ell_j = nu_j(q) mod 2`, `nu_j(q) = [(p_j/q) = -1]`.
A monomial `v in Z^r` has label `<v, ell>`; it *collides* at `q^a` with unit `i^s` iff
`<v, ell> = iota s mod m_(q^a)`. Write `a_q(v) = max_s a_(q,s)(v)` and
`m^gen_v = prod_q q^(a_q(v))`. Since `a_q(v) <= sum_s a_(q,s)(v)`,

```text
log m^gen_v  <=  sum_(s in Z/4) log m_(v,s).                                                  (1.1)
```

**The height condition** (`sharp` Prop. 6.1, per unit, as corrected by the referee) is
`m_(v,s) <= nu_s^(-1) e^(w(v)/2)` for every `v != 0` and `s`. With (1.1) it gives

```text
log m^gen_v  <=  2 w(v) + 2 log 2       for every member and every v != 0.                     (1.2)
```

(1.2) concerns residues and weights only, so it is a necessary condition on the *design* (the
residue data) before any angle is chosen.

**2-adic classes of primes.** Put `k(q) = v_2(m_q)`. Then `k(q) = 2` iff `q = +-3 mod 8`,
`k(q) = 3` iff `q = 7, 9 mod 16`, and `k(q) >= 4` iff `q = +-1 mod 16`. Reduction between levels
preserves `v_2(m_(q^a)) = k(q)`.

**Hit probabilities.** If a label is uniform on the `m/2` residues of a fixed parity `nu`, the
probability that the single monomial `e_j` collides at level 1 (with some unit) is the share of the
4-torsion `{0, m/4, m/2, 3m/4}` in that parity class:

```text
p(q, nu) = 4/m_q                if k(q) = 2   (two 4-torsion elements of each parity),
p(q, nu) = 8/m_q * [nu = 0]     if k(q) >= 3  (all four 4-torsion elements are even).         (1.3)
```

So at `q = +-1 mod 8` a column with `(p_j/q) = -1` can never collide as a single monomial.

**Twin coordinates (exact).** For the profile of `sharp` §4 (twins `f(x), s(x)` with
`a_(.,f(x)) - a_(.,s(x)) = eps_x e_x`), write `ell_f(x) = t_x`, `ell_s(x) = lambda'_x - t_x`, collect
the free labels `ltilde = (ell_j (j non-twin), lambda'_x)`, and let `atilde_x` be the row `x`
restricted accordingly (`a_(x,s(y))` in the `lambda'_y` slot). Then

```text
kappa_x = P_x + eps_x t_x,   P_x = <atilde_x, ltilde>;      <v, ell> = <u, ltilde> + <c', kappa>,   (1.4)
c'_x = eps_x (v_f(x) - v_s(x)),   u = vtilde - sum_x c'_x atilde_x,
```

where `kappa_x` is the point class (`sharp` Cor. 2.2(2)) and `vtilde = (v_j (j non-twin), v_s(x))`.
In particular `v in L` iff `u = 0` and `sum c' = 0`; the twin monomial `e_f(x)` has `c' = eps_x e_x` and
`u = -eps_x atilde_x` (a vector of weight about `r/2`); and the pair label is
`kappa_x - kappa_y = (twin-difference label) + <atilde_x - atilde_y, ltilde>`. The *solved*
coordinates `t_x = eps_x(kappa_x - P_x)` carry the tension between pair structure and monomial
structure that `sharp` §6.6 describes; §§2-4 make it quantitative.

*Proof of (1.4).* Substitute `a_(y,f(x)) = a_(y,s(x)) + eps_x [y = x]` into `kappa_y = <a_y, ell>` and
into `<v, ell>`. QED

## 2. Theorem 1: the Poisson barrier is robust (any prime set, any weights, any other structure)

`sharp` Prop. 6.5 shows that the uniform random fibre violates the height condition below
`(2 - eps) log^2 M/loglog M` per column, using single columns at the primes `q = +-3 mod 8` in
`[M^(1-eta), M]`. The mechanism needs much less: independence on *some* set of primes of
reciprocal mass `(log M)^(-C)`, and on a positive proportion of the columns. Everything else may be
structured, coherent or adversarially chosen.

**Hypothesis (Ind).** Fix the primes `p_j` (hence weights and parities), `X <= M^(A_0)`, a set `Q*`
of odd primes in `[M^(1-eta), X]` not dividing `N`, and a set `J` of columns. A random design
satisfies (Ind) if, conditionally on all residue data of the columns outside `J` and of the primes
outside `Q*`, the labels `(ell_j(q) mod m_q)_(j in J, q in Q*)` are independent, and `ell_j(q) mod m_q`
is uniform on the residues of parity `nu_j(q)`.

Put `p_j(q) = p(q, nu_j(q))` from (1.3) and `lambda_j = sum_(q in Q*) p_j(q)`.

**Theorem 1.** Fix `A_0 >= 1`, `eta in (0, 1/2)` and `C >= 0`. Put

```text
T_M = (1 - eta)^2 log^2 M / (2 (C + 1) loglog M).
```

There is `M_0 = M_0(A_0, eta, C)` such that for `M >= M_0` the following holds. Let a random design
satisfy (Ind), and let `J' ⊂ {j in J : lambda_j >= (log M)^(-C), w_j <= T_M}` (a set fixed by the
primes). Then

```text
Pr[ (1.2) holds for every e_j, j in J' ]  <=  exp(-|J'| M^(-1 + eta/2)).
```

**Corollary 1.** Under (Ind) with `lambda_j >= (log M)^(-C)` for all `j in J` and `|J| >= eps M`:
if a design completes to a member of `F^all_X` with probability `>= exp(-M^(eta/4))`, then

```text
W  >=  (eps M - M^(1 - eta/4)) T_M  =  (eps (1 - eta)^2 / (2(C+1)) - o(1)) M log^2 M / loglog M.
```

*Proof of Theorem 1.*
1. *Violation threshold.* Let `j in J'` and let `H_j` be the number of `q in Q*` at which `e_j`
   collides at level 1. Each such `q >= M^(1-eta)` divides `m^gen_(e_j)`, so
   `log m^gen_(e_j) >= H_j (1 - eta) log M`. By (1.2), the height condition fails as soon as
   `H_j (1-eta) log M > 2 w_j + 2 log 2`. Since `w_j <= T_M`, it suffices that `H_j >= k` with
   `k := floor((1-eta) log M/((C+1) loglog M)) + 2`.
2. *One column.* Conditionally on the data outside `(J, Q*)`, `H_j` is a sum of independent
   Bernoulli variables with means `p_j(q) <= 8/(q-1) <= 1/2` (for `q >= 17`). Let
   `p_max = max_q p_j(q) <= 8/(M^(1-eta) - 1)` and `lambda = lambda_j`.
   * `Pr[H_j = k] >= e_k(p) prod_q (1 - p_j(q)) >= e_k(p) e^(-2 Lambda)`, where `e_k` is the
     elementary symmetric function and `Lambda := sum_(M^(1-eta) <= q <= X) 8/(q-1) <= 8 log(A_0/(1-eta)) + 1`
     for large `M` (Mertens), using `1 - p >= e^(-2p)` for `p <= 1/2`.
   * Expanding `lambda^k` over ordered `k`-tuples, the tuples with a repeated prime contribute at
     most `binom(k,2) p_max lambda^(k-1)`. So `k! e_k(p) >= lambda^k (1 - k^2 p_max/(2 lambda))`, and
     the bracket is `>= 1/2` for large `M`, because `lambda >= (log M)^(-C)` and `k = O(log M)`.
   * Hence `log Pr[H_j >= k] >= -2 Lambda - log 2 - k C loglog M - k log k`. With `log k <= loglog M`
     and the choice of `k`, this is `>= -(1-eta) log M - O((C+1) loglog M) >= -(1 - eta/2) log M`
     for `M >= M_0`.
3. *All columns.* By (Ind) the events `{H_j >= k}`, `j in J'`, are conditionally independent. The
   probability that none occurs is at most `(1 - M^(-1+eta/2))^(|J'|) <= exp(-|J'| M^(-1+eta/2))`. If
   one occurs, (1.2) fails for that `e_j`. QED

*Proof of Corollary 1.* Let `J' = {j in J : w_j <= T_M}`. If `|J'| >= M^(1-eta/4)`, Theorem 1 bounds
the completion probability by `exp(-M^(eta/4))`. Otherwise at least `|J| - M^(1-eta/4)` columns have
`w_j > T_M`. QED

**Scope and validity.** Theorem 1 is a statement about *construction methods* (random designs), not
about inputs. It is valid for every profile, every weight vector, every `X <= M^(A_0)`, and
arbitrary (structured, coherent, adversarial) data outside `(J, Q*)`. It says nothing about actual
circles, whose residues are not random.

**Consequences.**
1. `sharp` Prop. 6.5 is the case `Q* = {q = +-3 mod 8, M^(1-eta) <= q <= M}`, `C = 0`, `J` the
   non-twin columns (Lemma 6.3(2) there gives (Ind)).
2. *Partial structure does not help.* A design that is structured or coherent at all primes except
   a set `Q*` of reciprocal mass `>= (log M)^(-C)`, and uses the random fibre on `Q*`, stays at
   `Theta(M log^2 M/loglog M)`. Examples: `Q*` = the primes `q = +-1 mod 2^K` in the window, for any
   `2^K <= (log M)^C` (mass `~ 2^(2-K) log(1/(1-eta))`); or `Q*` = the primes in `[M/L, M]` with
   `L = exp((log M)^(1-C))` (mass `~ 4 log L/log M = 4 (log M)^(-C)`).
3. *Only a set of mass `M^(-Omega(1))` may be left random.* If the random part has mass
   `lambda <= M^(-delta)`, the computation in step 2 gives `k ~ 1/delta`, i.e. no loss in rate. So a
   construction at `O(M log M)` built this way must be non-random (structured or certified) on all
   primes of `[M^(1-eta), X]` outside a set of reciprocal mass `M^(-Omega(1))`; for the 2-adic
   classes this means every class `k(q) = k` with `2^k <= M^(o(1))`.
4. *Twins.* By (1.4), `t_x = eps_x(kappa_x - P_x)`. If the non-twin labels are random and `kappa` is
   any fixed (structured) point design, then `P` is uniform on the image of `atilde`. Whenever
   `(atilde_x)_(x in R)` is surjective modulo `m_q/2` for `q in Q*` (on a row set `R`, `|R| >= eps M`),
   (Ind) holds for the twin columns `f(x)`, `x in R`, and Corollary 1 applies to them. So "structure
   the points, randomise the columns" fails at the twin columns. (The surjectivity is a property of
   the 0/1 profile modulo the primes dividing `m_q/2`; for the QR-code structure of Paley matrices it
   fails modulo 2, and then only the odd part of the label is uniform. We do not pursue this; the
   coherent version of the same tension is Theorem 2.)

*Evidence* (`check_tail.py`, §5): the single-column tails of the random fibre and of the hybrid
design of Consequence 2 at `M = 64 ... 4096`.

## 3. Structured fibres: coherence, parity repair by roots of unity, and its price

By Theorem 1, in any design whose random part satisfies (Ind), randomness must be confined to a
set of primes of reciprocal mass `M^(-Omega(1))`. The only
known non-random way to control *all* monomials is coherence (`sharp` §6.6): labels that are
discrete logarithms of fixed Gaussian numbers, so that a collision of `v` is a divisibility of one
fixed integer whose size is exponential in `||v||_1`. This section makes coherence compatible with
the forced parity wherever that is cheap, and computes the price elsewhere.

### 3.1 Coherent designs with root-of-unity repair

**Definition 3.1 (repaired coherent column).** Let `zeta_j in Q(i)^*` with `zeta_j/conj(zeta_j)` not
a root of unity, and let `q` be an odd prime with `q` coprime to the numerator and denominator of
`Norm zeta_j`. Let `f_q in T_Q` be an element of order exactly `2^(k(q))` (a generator of the
2-Sylow subgroup; reductions of `f_Q` are used at lower levels). Put

```text
eps_j(q) = [ (Norm zeta_j / q) != (p_j / q) ],       r_j(q)/conj r_j(q) := (zeta_j/conj zeta_j) * f_q^(eps_j(q)).   (3.1)
```

**Lemma 3.2 (repair is valid).**
1. `f_q` has odd label: `nu_Q(f_q) = -1`. For `k(q) = 2`, `f_q = +-i`.
2. (3.1) has `nu_Q = (p_j/q)` at every level, so by `sharp` Prop. 2.1(4) there are `r_j(q) mod Q`
   with `Norm r_j(q) = p_j mod Q` realising it, compatibly across levels.

*Proof.* 1. `T_Q` is cyclic of order `m = 2^k o`, `o` odd. An element of order `2^k` is
`g^(o u)` with `u` odd, so its label `o u` is odd. If `k = 2` the elements of order 4 are `+-i`.
2. `nu` is a homomorphism, `nu(zeta/conj zeta) = (Norm zeta/q)` and `nu(f_q) = -1`. Compatibility:
the 2-Sylow subgroups of `T_(q^(a+1))` and `T_(q^a)` have the same order `2^k`, and reduction is onto
(Prop. 2.1(6)), so it maps `f` to an element of order `2^k`. QED

**Proposition 3.3 (certified collision moduli).** Let `v in Z^r`, `Z_v = prod_j zeta_j^(v_j) = alpha/beta`
with `alpha, beta in Z[i]` coprime to `q`, `gamma_v = alpha conj(beta)`, and `E_v(q) = sum_j v_j eps_j(q)`.
1. If `v` collides at `q^a` with some unit, then `q^a | Im(gamma_v^(2^(k(q))))`.
2. Hence, for every `k >= 2`,
   `sum_(q <= X, k(q) = k) a_q(v) log q <= log |Im(gamma_v^(2^k))| <= 2^k log|gamma_v|`, provided
   `Im(gamma_v^(2^k)) != 0`, which holds when `gamma_v/conj(gamma_v)` is not a root of unity (e.g. if
   the `zeta_j` are pairwise non-associate, non-conjugate Gaussian primes of odd norm and `v != 0`).
3. If `k(q) = 2`, collision with *some* unit is equivalent to `q^a | D_s(gamma_v)` for some `s`; it
   does not depend on `E_v(q)`.
4. If `Z_v` is a unit times a rational (e.g. `Z_v = 1`), then `v` collides at every `q^a` with
   `k(q) = 2`, with unit `i^(+-E_v(q))`.

*Proof.* 1. The monomial is `(Z_v/conj Z_v) f^E`. If it equals `i^s`, raise to the power `2^k`:
`f^(2^k) = 1` and `i^(2^k) = 1`, so `(Z_v/conj Z_v)^(2^k) = 1`, i.e. `gamma^(2^k) = conj(gamma)^(2^k)
mod q^a`, i.e. `q^a | 2i Im(gamma^(2^k))`, and `q` is odd. 2. Distinct primes are coprime, and
`|Im w| <= |w|`. 3. With `f = +-i`, the monomial is `(gamma/conj gamma) i^(+-E)`, which lies in
`{1, i, -1, -i}` iff `gamma/conj gamma` does; that is `sharp` Prop. 1.2's computation. 4. Immediate. QED

So the classes `k = 2, 3` (`q != +-1 mod 16`, three quarters of all primes) can be made coherent
*despite* the parity, at a certified cost of `4 + 8 = 12` times `log|gamma_v|`. Class `k` costs
`2^k log|gamma_v|`.

*Exact check* (`check_root.py`): all 89 odd prime powers `<= 400`; per modulus 60 trials of five
columns with random `zeta_j`, random primes `p_j`, realisation of `r_j` by Hilbert 90 and a square
root, `Norm r_j = p_j` and parity verified; 17038 collisions tested against item 1; item 4 at every
`q = +-3 mod 8`. ALL ROOT-REPAIR CHECKS PASSED.

### 3.2 What the root repair buys, and why it is not enough

**Corollary 3.4 (coherent columns at the classes `k <= K`).** Let the `zeta_j` be distinct Gaussian
primes with `X < Norm zeta_j <= Y`, and use (3.1) at every prime with `k(q) <= K`. Then for every
`v != 0`, the part of `m^gen_v` supported on those primes satisfies

```text
log m^gen_v |_(k(q) <= K)  <=  2^(K+1) log|gamma_v|  <=  2^K ||v||_1 log Y.
```

This is linear in `||v||_1`, as Prop. 6.6 of `sharp` needs, but with constant `2^K`.

**What remains: the classes `k > K`.** The primes with `k(q) > K` have reciprocal mass
`~ 2^(1-K) log(A_0/(1-eta))` in `[M^(1-eta), X]`. Two options exist for them.
* *Certified (root repair).* Cost `2^k log|gamma_v|`, unbounded in `k`. For `2^k` near `X` the bound
  is useless.
* *Random.* By Theorem 1 (Consequence 2), mass `>= (log M)^(-C)` left random costs
  `(1/(2(C+1)) - o(1)) log^2 M/loglog M` per column (on a positive proportion of columns).

*Heuristic H1 (the parity balance).* Pseudo-random behaviour of the classes `k > K` (power-residue
symbols of the fixed `zeta_j` against the arbitrary choice of `f_q`) costs, by the computation of
Theorem 1 with `lambda ~ 2^(-K)`, about `beta_rand ~ A log M/(K log 2)` per unit of `||v||_1 log M`;
the certified classes `k <= K` cost `beta_cert ~ 2^K log Y/log M`. The sum
`2^K + log M/K` is minimised at `2^K ~ log M/loglog M`, giving `Theta(log M/loglog M)`, i.e. again
`W = Theta(M log^2 M/loglog M)`. This is heuristic: it treats the 2-power residue symbols of fixed
Gaussian numbers at primes `q = +-1 mod 2^k` as independent, and it ignores possible cleverer
deterministic choices of `f_q`.

**Why other repairs are worse (exact).**
* *A non-root repair* `f_q = xi/conj xi` with `(Norm xi/q) = -1` (a small Gaussian integer, possible
  for every `q` with the least such norm `O(log^2 q)` under GRH) makes the collision condition
  depend on the integer `E_v(q) = sum_j v_j eps_j(q)`: `v` collides at `q` iff
  `q^a | D_s(gamma_v xi^(E_v(q)))` for some `s`. The certified bound is
  `sum_(E) log|D(gamma_v xi^E)| ~ (2||v||_1 + 1) log|gamma_v|`, quadratic in `||v||_1`. The integers
  `eps_j(q)` are Legendre symbols of `p_j Norm zeta_j`, which the design cannot control (next item).
* *Matching the parities instead of repairing them* requires `(p_j Norm zeta_j / q) = 1` for all
  primes `q` of the classes to be matched. For a set `S` of such primes this has density
  `2^(-|S|)` among primes `p_j` (or among `zeta_j`). Matching the classes `k > K` up to `X` needs
  `log p_j >~ pi(X) 2^(2-K) log 2`, i.e. `2^K >~ X/(log X log p_j)`; with `log p_j = O(log M)` the
  certified cost `2^K` of the remaining classes is then `>= M^(A - o(1))`. (`sharp` §6.6: with all
  primes matched this is the pseudosquare obstruction.)
* *A root of unity of smaller order cannot repair:* every odd-label element of `T_Q` has order
  divisible by `2^(k(q))` (Lemma 3.2(1) and cyclicity). So `2^(k(q))` in Prop. 3.3 is forced for
  repairs by roots of unity.

So parity is a genuine obstruction only on the classes `q = +-1 mod 16`, and there it costs, on
present methods, the same `log M/loglog M` factor. The second obstruction, independent of parity,
is the pair/twin tension, which is §4.

## 4. Theorem 2: coherent designs on the Paley profile cannot certify the pairs

The natural structured fibre is coherent on all columns, twin columns included, with the root
repair of §3. Then every monomial, *and every pair character*, collides only at divisors of an
explicit Gaussian integer. The cheapest certificate for the pairs is the size of that integer. This
section shows that on the Paley profile of `sharp` §4 such certificates cannot be small for the
pairs and the columns simultaneously. The proof is an uncertainty principle for the Hadamard
matrix plus a pigeonhole on the label classes.

**Definition 4.1 (G-coherent design).** Let `gamma_1, ..., gamma_D` be Gaussian primes of odd norm,
pairwise non-associate and non-conjugate, and `G in Z^(D x r)`. Put `zeta_j = prod_t gamma_t^(G_tj)`.
A *G-coherent design* uses (3.1) with these `zeta_j` at every prime `q = 3 mod 8`, `q <= X` (such `q`
are inert, so they never divide `N` or `Norm gamma_t`); at all other primes it is arbitrary. Define

```text
column heights  h_j = sum_t |G_tj|,   C = max_j h_j;
pair heights    H_xy = sum_t |<g_t, a_x - a_y>|   (g_t = row t of G),   Hbar = binom(M,2)^(-1) sum_(x<y) H_xy.
```

For `v in Z^r`, `gamma_v = prod_t gamma_t^(|(Gv)_t|)` up to conjugation of factors, so
`log|gamma_v| = sum_t |(Gv)_t| log|gamma_t|`. By Prop. 3.3(3), at the primes `q = 3 mod 8` the
collision moduli of a single column are certified by `4 h_j log max|gamma_t|`, and those of the pair
`(x, y)` by `4 H_xy log max|gamma_t|`. "Certified" means: bounded through the size of the integer,
whatever its factorisation.

**Theorem 2.** Let the profile be that of `sharp` §4: `[1 | H']` an `M x M` Hadamard matrix (normalised
Paley-I, `M = q0 + 1`), `M - 1` labels with `b >= 1` physical columns each, canonical twins. Let a
G-coherent design be given with `X >= M`, and put `s = ceil(8 (C + 2 Hbar))`. If

```text
(M - 1)^2  >  4 s^2 (2 s b C + 1),                                                           (4.1)
```

then every set of data satisfying the height condition (1.2) for all monomials (in particular every
member of `F^all_X` using this design) has

```text
W  >=  (M - 1) (theta(X; 8, 3) - 2 log 2) / (8 s),        theta(X; 8, 3) = sum_(q <= X, q = 3 mod 8) log q.
```

**Corollary 2.** Fix `K`. Suppose `W <= K M log M` and `X >= M`. Then for large `M` either
`C + 2 Hbar >= M/(330 K log M)`, or (4.1) fails, which forces `(C + 2 Hbar)^3 C >= M^2/(2^(14) b)`
(when `C >= 1` and `C + 2 Hbar >= 1`).
In both cases the column heights or the *average* pair heights are `>= M^(1/2 - o(1))`. If moreover
single columns are certified (`C <= K' log M`), then `Hbar >= M^(1/2 - o(1))`: the typical pair
integer `gamma_(a_x - a_y)` has `log`-size `>= 0.8 Hbar`, and its collision modulus is not certified
by its size.

*Proof of Theorem 2.*
1. *Row images.* Write `col_a = (1 + H'_(.,a))/2` for the 0/1 column of label `a`, and let
   `psi_(t,a) = sum_i G_(t,(a,i))` (sum over the `b` physical columns of label `a`, flipped ones
   included). Since `a_(.,f(x)) = col_(alpha(x)) + eps_x e_x` and the other columns are equal to their
   label column,

   ```text
   (A g_t)_x = <g_t, a_x> = c_t + F_(t,x),     F_t = (1/2) H' psi_t + phi_t,    phi_t = (eps_x G_(t,f(x)))_x,
   ```

   with `c_t` constant in `x`. Hence `<g_t, a_x - a_y> = F_(t,x) - F_(t,y)`.
2. *Pair budget.* For any `F in R^M` with median `mu`: if `F_x >= mu`, at least `M/2` indices `y` have
   `F_y <= mu`, so `sum_y |F_x - F_y| >= (M/2)|F_x - mu|`; similarly below. Hence
   `sum_(x<y) |F_x - F_y| >= (M/4) ||F - mu 1||_1`. Summing over `t`:
   `sum_t ||F_t - mu_t 1||_1 <= (4/M) sum_(x<y) H_xy <= 2 Hbar M`.
3. *Uncertainty.* Let `psi in Z^(M-1) \ {0}` and `m in R`. Then `y = (1/2)H' psi - m 1 = [1 | H'] z` with
   `z = (-m, psi/2)`. Since `[1 | H']` is Hadamard, `||y||_2^2 = M ||z||_2^2` and `||y||_inf <= ||z||_1`, so
   `||y||_1 >= ||y||_2^2/||y||_inf >= M ||z||_2^2/||z||_1 >= M ||z||_1/||z||_0`. Here
   `||z||_1 >= ||psi||_1/2 >= ||psi||_0/2` (integer entries) and `||z||_0 <= ||psi||_0 + 1 <= 2||psi||_0`.
   So `||(1/2)H' psi - m 1||_1 >= M/4`.
4. *Few rows see the label classes.* Let `T* = {t : psi_t != 0}`. For `t in T*`, by 3,
   `||F_t - mu_t 1||_1 >= M/4 - ||phi_t||_1`. Summing and using 2 and
   `sum_t ||phi_t||_1 = sum_x h_(f(x)) <= M C`:
   `|T*| M/4 <= M C + 2 Hbar M`, i.e. `|T*| <= 4(C + 2 Hbar) <= s/2`.
5. *A degenerate monomial.* For `beta in Z^(M-1)` put `v_beta = sum_a beta_a s_a`, `s_a = sum_i e_(a,i)`.
   Then `(G v_beta)_t = <beta, psi_t>`, which vanishes for `t` outside `T*`. Let `Lab'` be the labels
   with `W_a := sum_i w_(a,i) <= 2W/(M-1)`; by Markov `n' := |Lab'| >= (M-1)/2`. Consider the
   `binom(n', s) >= (n'/s)^s` vectors `beta in {0,1}^(Lab')` of support exactly `s`. Their images
   `(<beta, psi_t>)_(t in T*)` lie in a box with at most `(2 s b C + 1)^(s/2)` points, since
   `|psi_(t,a)| <= sum_i h_(a,i) <= b C`. By (4.1), `(n'/s)^2 > 2 s b C + 1`, so two distinct `beta, beta'`
   have the same image. Then `v := v_(beta - beta')` is nonzero, `G v = 0`, and `v` is supported on at most
   `2s` light labels with coefficients `+-1`.
6. *`v` is not in `L`.* If `v = sum_x c_x a_x`, compare the twin columns `f(x)` and `s(x)`, which carry
   the same label: `v_(f(x)) = v_(s(x))`, while `(sum_y c_y a_y)_(f(x)) - (sum_y c_y a_y)_(s(x)) = eps_x c_x`.
   So `c = 0` and `v = 0`, a contradiction.
7. *Collisions everywhere.* `G v = 0` gives `Z_v = 1`. By Prop. 3.3(4), `v` collides at every
   `q^(A_q)` with `q = 3 mod 8`, `q <= X`. So `log m^gen_v >= theta(X; 8, 3)`.
8. *Weight.* `w(v) <= sum_(a in supp) W_a <= 2s * 2W/(M-1)`. By (1.2),
   `theta(X; 8, 3) <= 2 w(v) + 2 log 2 <= 8 s W/(M-1) + 2 log 2`. QED

*Proof of Corollary 2.* `theta(X; 8, 3) >= (1 - o(1)) X/4 >= (1-o(1)) M/4`. If (4.1) holds, Theorem 2
gives `K M log M >= (1 - o(1)) M^2/(32 s)`, so `s >= (1 - o(1)) M/(32 K log M)` and
`C + 2 Hbar >= s/8 - 1/8`. If (4.1) fails, `4 s^2 (2 s b C + 1) >= (M-1)^2`; with `2sbC + 1 <= 3sbC`
and `s <= 9(C + 2 Hbar)` this gives `12 * 729 (C + 2Hbar)^3 b C >= (M-1)^2`. QED

**Validity and scope.** Theorem 2 is exact and concerns a *construction method*: G-coherent designs on
the Paley profile of `sharp` §4. It uses only the height condition at the primes `q = 3 mod 8`, which
every member satisfies. It does not say that pair demands are large, only that they are not bounded
by the size of the pair integers; the actual demands are the `X`-smooth parts of those integers.

### 4.2 What Theorems 1 and 2 leave open

Combining: a fake at `W = O(M log M)` in P3-all built from labels needs
* non-random labels on all but a set of primes of reciprocal mass `M^(-Omega(1))` (Thm 1);
* a repair of the parity on the classes `q = +-1 mod 16` that is neither random nor certified by
  roots of unity at bounded cost (§3.2; heuristic H1 puts the natural trade-off at the same rate);
* on the Paley profile, if coherent: pair integers of typical size `exp(M^(1/2 - o(1)))` whose
  `X`-smooth parts must nevertheless be small for *every* pair (Thm 2).

*Heuristic H2 (pairs of large height).* If the `X`-smooth part `S` of a pair integer of size
`exp(M^(Omega(1)))` behaves like that of a random integer, then
`Pr[log S >= mu log M] ~ rho(mu)` (Dickman), and over `binom(M, 2)` pairs the maximum is
`~ 2 log^2 M/loglog M`. With the uniform Paley margins (pair slack `(b-4) tau/2`) this forces
`tau ~ log^2 M/loglog M`, the same rate once more. Counting gives the same picture for a
*deliberate* choice of the `zeta_j`: there are `exp(O(b M log M))` choices, against
`binom(M,2) rho(mu)` independent-looking constraints, so `rho(mu) <~ b log M/M` is needed, which is
`mu ~ log M/loglog M`.

*The one escape visible here* is the one `p3lower.md` isolates: non-uniform margins. By its
Theorem A the pair channel costs at most `pi lambda_max(D)` at real weights, and for random-like
demands `lambda_max(D) ~ M E[log S] = O(M log X)`. A coherent design on a cut-adapted profile would
then need (i) Question (Q) of `p3data.md` (prime weights, `O(M)` columns) for that demand matrix,
(ii) an analogue of Theorem 2 to fail for the adapted profile, and (iii) the parity repair on the
classes `q = +-1 mod 16`. None of the three is settled.

## 5. Checks and evidence (`round7/p3Dconstruct_checks/`)

No proof depends on this section except where marked "exact check". Each script runs from its
directory with `python3 <script>`; `gauss.py` holds the exact `Z[i]/Q` arithmetic.

**Exact checks.**
* `check_root.py` (Lemma 3.2, Prop. 3.3): all 89 odd prime powers `<= 400`. For each: the 2-Sylow
  generator `f` has odd label; `k(q) = 2` iff `q = +-3 mod 8` and then `f = +-i`; `k(q) = 3` iff
  `q = 7, 9 mod 16`. For 60 trials of 5 random columns: (3.1) has the parity `(p_j/q)`, and an `r_j`
  with `Norm r_j = p_j mod Q` and `r_j/conj r_j` = (3.1) is constructed (Hilbert 90, square root,
  Hensel) and verified. 17038 collisions of random monomials all satisfy `Q | Im(gamma_v^(2^k))`.
  Prop. 3.3(4) at every `q = +-3 mod 8`. ALL ROOT-REPAIR CHECKS PASSED.
* `check_T3.py` (Theorem 2): Paley-I matrices `q0 = 43, 47, 59` are normalised Hadamard; the
  uncertainty bound of step 3 (`>= M/4`) on 3000+ integer vectors `psi` (observed minimum `M/2`); the
  row-image identity of step 1 on 50 random rows of `G` for the profile `q0 = 43`, `b = 6`; two explicit
  designs showing the dichotomy: `G = I` has `C = 1` and every pair height `132 = r/2 + 3` (`r = 258`);
  the within-label path design has `C = 2`, average pair height `2.05`, and the label sum `s_5`
  (weight 6) is degenerate (`G s_5 = 0`) and not in `L`. ALL T3 CHECKS PASSED.

**Evidence 1: single columns** (`check_tail.py`, `X = 3M`, `r = 8M` columns, genuine parities of
random primes `p_j ~ M^4`, auxiliary Gaussian primes of norm `> X`). Entries are
`log m^gen(e_j)/log M`, maximum and mean over the `r` columns. RAND = random fibre; HYB = coherent with
root repair on `k(q) <= 3`, random on `k(q) >= 4`; COH = coherent with root repair everywhere (exact).

| `M` | RAND max / mean | HYB max / mean | COH max / mean | COH maximiser: `k<=3` / `k>=4` part | `log M/loglog M` |
|---|---|---|---|---|---|
| 64 | 8.04 / 3.77 | 6.53 / 3.33 | 6.63 / 3.22 | 5.70 / 0.93 | 2.92 |
| 256 | 11.52 / 3.82 | 8.03 / 3.20 | 7.26 / 3.11 | 4.19 / 3.07 | 3.24 |
| 1024 | 11.16 / 3.88 | 8.55 / 3.13 | 7.36 / 2.99 | 4.66 / 2.70 | 3.58 |
| 2048 | 10.92 / 3.87 | 8.92 / 3.11 | 7.28 / 2.97 | 4.63 / 2.65 | 3.75 |

(Full table, `M = 64 ... 2048`, in `tail_out.txt`.) Coherence lowers the maximum by about
`3.5 log M`; the classes `k(q) >= 4` contribute up to `3 log M` to the worst coherent column. At these
sizes `log M/loglog M` varies only from 2.9 to 3.8, so the data cannot separate `O(log M)` from
`Theta(log^2 M/loglog M)`; they are consistent with both.

**Evidence 2: pairs** (`check_pairs_fast.py`, cross-checked against the direct `check_pairs.py` on
`q0 = 43, 59`, identical output). Coherent design `G = I` with root repair on the profile of `sharp` §4
(`b = 6`), all odd prime powers `<= 3M`. `D1` = pair demand with unit 1 (the one constrained by the
pair slack); "some unit" = `log m^gen_xy`. Control: `check_pairs_random.py`, same parities, point
labels uniform in their parity class.

| `M` | coherent `D1` max / mean | control `D1` max / mean | coherent some-unit max | `lambda_max(D1)/(M log M)` (coh.) |
|---|---|---|---|---|
| 44 | 4.18 / 0.919 | 5.36 / 0.946 | 8.84 | 0.922 |
| 132 | 5.76 / 0.954 | 6.65 / 0.951 | 12.69 | 0.954 |
| 504 | 6.24 / 0.960 | 8.15 / 0.964 | 11.26 | 0.960 |
| 1020 | 7.43 / 0.963 | 8.83 / 0.964 | 14.22 | 0.963 |
| 2040 | 7.40 / 0.967 | 7.15 / 0.967 | 14.21 | 0.967 |

(in units of `log M`; full tables in `pairs_fast_out.txt`, `pairs_random_out.txt`.) The coherent pair
demands are statistically indistinguishable from random point labels, which supports the premise
of H2. They are close to the structured design `D_A` (`sharp` §3: `5.7-6.8 log M` at
`M = 64 ... 16384`). `lambda_max(D1) ~ 0.96 M log M` for both, consistent with `p3lower.md` Thm A: on
a cut-adapted profile the coherent pair channel would cost `O(M log M)` at real weights.

## 6. Open points

1. **(Coherent + cut-adapted profiles.)** Theorem 2 needs the Hadamard uncertainty of the Paley
   profile. On a profile adapted to the demand matrix (`p3lower.md` Thm A), uniform margins are not
   needed, and random-like pair demands cost `~ pi lambda_max(D) = O(M log X)` at real weights. The
   evidence in §5 (pairs) has `lambda_max(D) ~ 0.94 M log M`. Open: an analogue of Theorem 2 for
   cut-adapted profiles, together with Question (Q) of `p3data.md` (prime weights, `O(M)` columns).
2. **(Parity on `q = +-1 mod 16`.)** Is there a deterministic choice of the 2-Sylow repairs `f_q`
   (or of the auxiliary integers) whose collision moduli on the classes `k(q) >= 4` are linear in
   `||v||_1` with a constant independent of `k`? A negative answer would make H1 a theorem for
   coherent designs.
3. **(Certified versus actual.)** Theorems 1 and 2 bound what *certificates* can do (randomness,
   size). The actual demands of coherent designs are `X`-smooth parts of explicit Gaussian
   integers. Any construction at `O(M log M)` must control smooth parts of `binom(M, 2)` integers of
   size `exp(M^(1/2 - o(1)))` at once (Cor. 2). Conversely, a proof that this is impossible would
   be a statement about smooth values of products of Gaussian primes, close to Hypothesis H of
   `p3lower.md`.
4. **(Lower bound.)** None of this bounds `W_min^all` from below. A lower bound
   `W_min^all >= g(M) M log M`, `g -> infinity`, would be a growth improvement for actual circles
   (`sharp` §6.7); the barriers here concern constructions only.
