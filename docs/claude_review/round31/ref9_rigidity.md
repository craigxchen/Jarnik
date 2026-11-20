# Referee report on Part IX.7 (rigidity for rational prime witnesses)

Block: section IX.7 of `astra_free_IX/answer.md` (lines 2267-2300), its notes `astra_free_IX/ix_rigidity_work.md`
(identical text), the COMPUTED negative control (`ix_factor_check.py::real_power_test`, read, not imported).
My script: `ref9_rigidity/check_rigidity.py`, output `ref9_rigidity/check_rigidity_out.txt`.

## 0. Verdict

* **The theorem is correct (Proved).** I checked every step. Two hypotheses are worth stating more sharply:
  * the norm condition `r_j conj(r_j) = p_j` is implied by the other hypotheses;
  * (G) enters only through the irrationality of `phi_j/pi` and `(phi_i+phi_j)/pi`, and it cannot be dropped.
* **The negative control is reproduced exactly (Computed).** I add a second control. In it the angle is genuine and
  the lift is fake, and the test fails at exponent 2 for cutoffs `X >= 31`.
* **What it means.** IX.7 characterises genuineness. It is not a bound. Its message is the one the record already
  has in other languages: an input that is multiplicative and holds at all heights is the problem itself
  (composite.md Prop. 3.1, round 31 sieve rigidity, sharp.md Remark 2.3(4)).
  * New: the statement in residue language with rational lifts, proved by elementary means. It is also the first
    rigorous all-height complement to the coherent designs of sharp.md §6.6.
  * Its hypothesis carries no information at a finite cutoff. At every finite cutoff, every factored datum (d) has a
    rational lift of exact norm `p_j` (Proposition 3 below, also checked numerically). So IX.7 says nothing about
    P3-all* (cutoff `3M`), nothing about arbitrary adelic witnesses, and nothing about the growth rate.

## 1. The theorem, with all hypotheses (Proved)

**Data.** An angle model in the sense of sharp.md Definition 1.1(a)-(c). It has distinct primes `p_j = 1 mod 4`, a
profile, angles `phi` with (1.1) of outside.md, and the monomial gaps (G) of outside.md Lemma 0.1 for every
`v != 0`. Fix a label `j`. Residue data at every odd prime power: for every odd prime `q != p_j` and every `a >= 1`, a
unit `r_j(q^a)` of `Z[i]/q^a`, compatible under reduction.

**Hypothesis (L) (one rational lift).** There are `r_j in Q(i)^*` and a finite set `E_j` of primes such that
`r_j(q^a) = r_j mod q^a` for all `q not in E_j` and all `a`. `E_j` depends arbitrarily on the instance.

**Hypothesis (A) (all-height strengthened test, one unit, coordinate rays).** For every finite cutoff `X` and every
`n >= 1`, inequality (1.2) of sharp.md Definition 1.1(e-all) holds with `v = n e_j` and `s = 2` (unit `-1`). By
sharp.md (1.2) that inequality reads

```text
|sin(n phi_j - pi/2)| = |cos(n phi_j)|  >=  nu_2 * m_(n e_j, 2) * e^(-w(n e_j)/2) = m_(n e_j,2) * p_j^(-n/2),
m_(n e_j,2) = prod_(q <= X odd, q not | N) q^(a_q),  a_q = max{a <= A_q : (r_j/conj r_j)^n = -1 mod q^a}.
```

Here `nu_2 = 1` because `s = 2` is even, and `w(n e_j) = n log p_j`. The constant is one.

**Theorem IX.7 (corrected form).** Assume (G), (L) and (A) for label `j`. Assume also either (d), meaning
`Norm r_j(q^a) = p_j mod q^a` for every `q not in E_j`, or directly `r_j conj(r_j) = p_j`. Then:
* `r_j = u pi_j` for a Gaussian unit `u` and a Gaussian prime `pi_j` of norm `p_j`;
* `phi_j = eps_j arg r_j (mod pi)` for a sign `eps_j`.

Suppose (A) also holds for `v = n(e_i + e_j)` on the edges of a spanning tree. Then all `eps_j` are equal. In every
case the model, with `phi_j` replaced by `arg rho_j`, where `rho_j = r_j` or `conj r_j` according to `eps_j`, is an
actual normalised cluster. The new common phase is `theta_0 - pi sum_j e_j m_j`, where `phi_j = arg rho_j + pi m_j`.
The `delta_x` and `k_x` are unchanged.

**Quantifiers.** The cutoff needed to see a violation depends on the instance through `E_j`, `J` and the
Diophantine type of `phi_j`. The angular step is ineffective. Nothing is uniform in the radius.

## 2. Line-by-line check of the proof

### (a) Structure of `r_j` (Proved)

Write `r_j = eps (1+i)^t prod_(q inert) q^(c_q) prod_(p split) pi_p^(a_p) conj(pi_p)^(b_p)`. The integers
`t, c_q, a_p, b_p` may be negative.
* The norm is `2^t prod q^(2c_q) prod p^(a_p+b_p) = p_j`. Hence `t = 0` and `c_q = 0`. At the split primes,
  `a_p + b_p = 1` at `p_j` and `0` elsewhere.
* So no ramified factor, no inert factor and no rational-integer factor `pi conj(pi)` occurs. This is the
  "norm condition" step, and it is correct.
* `r_j/conj r_j = eps^2 prod (pi_p/conj pi_p)^(d_p)` with `d_p = a_p - b_p`.
  * `d_p = s_p mod 2` with `s_p = a_p + b_p`. This is Astra's "difference has the parity of the sum".
  * So `d_(p_j)` is odd and the other `d_p` are even.
* Put `alpha = prod (pi_p or conj pi_p)^(|d_p|)`, choosing the side by the sign of `d_p`. Then
  `r_j/conj r_j = u alpha/conj alpha` with `u = eps^2 in {+-1}` (a unit in any case). Also
  `J = |alpha|^2 = p_j^(|d_(p_j)|) prod_(p != p_j) p^(|d_p|) >= p_j`, and `J` is odd.
* `arg(alpha)/pi` is irrational. Otherwise `alpha^n` is real for some `n >= 1`, so `alpha^n = +- conj(alpha)^n`. That
  contradicts the disjoint prime supports of `alpha` and `conj alpha`, since `J > 1`.

**Correction (redundant hypothesis; Proved).** Assume (L) and (d) hold at all prime powers outside the finite
set `E_j`. Then the fixed rational number `Norm(r_j) - p_j` is divisible by arbitrarily high powers of every
`q not in E_j`, so it is zero. Thus `r_j conj(r_j) = p_j` follows from (L) and (d). Only `r_j/conj r_j` enters (A).
So a lift `r_j` of norm `p_j t^2`, `t in Q^*`, can be rescaled to norm `p_j`.

### (b) Clean real parts and (IX.9) (Proved)

`E` contains 2, the parent primes, `E_j`, the primes of the numerator and denominator of `r_j`, and the primes of
`J`. Take `L` to be a multiple of 4 and of the orders of `alpha` in `(Z[i]/l)^*` for `l in E` with `l not | J`. Let
`n = Lk`.
* For `l not | J`: `alpha^n = 1 mod l Z[i]`, hence `Re alpha^n = 1 mod l`. This includes `l = 2`, where the unit
  group `{1, i} mod 2` has order 2.
* For `l | J`: `l` splits, and `alpha` is divisible by exactly one of the two primes above `l`, say `pi`. Then
  `2 Re alpha^n = alpha^n + conj(alpha)^n = conj(alpha)^n != 0 mod pi`.
* So `Re alpha^n` is coprime to `E`. In particular it is odd, so every prime it has is odd and not a parent prime.
* For `q^k || Re alpha^n`, `q not in E`: `alpha^n = -conj(alpha)^n mod q^k`, and `conj alpha` is a unit mod `q`. So
  `(alpha/conj alpha)^n = -1`. Also `u^n = 1` because `4 | n`, and `r_j mod q^k` is the residue by (L). So
  `(r_j/conj r_j)^n = -1 mod q^k`, i.e. `a_q >= k` as soon as `X >= q^k`.
* Taking `X >= |Re alpha^n|` gives `m_(n e_j,2) >= |Re alpha^n|`. Then (A) is exactly (IX.9):
  `|cos(n phi_j)| >= p_j^(-n/2) |Re alpha^n|`. The trigonometric function (`cos`), the constant (one) and the weight
  (`n log p_j`) all match sharp.md (1.2).

### (c) `J = p_j` (Proved)

* `Re alpha^n = J^(n/2) cos(n theta)` with `theta = arg alpha`, and `L theta/pi` is irrational. Hence
  `|cos(kL theta)| >= 1/2` for infinitely many `k`.
* There (IX.9) gives `1 >= (J/p_j)^(n/2)/2` along an unbounded sequence of `n`, so `J <= p_j`.
* With (a): `|d_(p_j)| = 1` and `d_p = 0` elsewhere. Then `(a, b) = (1, 0)` or `(0, 1)` at `p_j` and `a = b = 0`
  elsewhere, so `r_j = eps pi_j` or `eps conj(pi_j)`.
* (G) is not used in this step.

### (d) The angular congruence (Proved)

Take `z = r_j` (or `z = r_i r_j` on an edge), `beta = arg z`, and `phi = phi_j` (or `phi_i + phi_j`). The weight of
`v = n e_j` (or `n(e_i+e_j)`) is `n log|z|^2`. Step (b) with `alpha` replaced by `z` (now `J = |z|^2`) gives
`|cos(kL phi)| >= |cos(kL beta)|` for all `k >= 1`, at all sufficiently large cutoffs.
* `beta/pi` is irrational because `z` and `conj z` have disjoint supports.
* `phi/pi` is irrational by (G): `dist(n phi, (pi/2)Z) > 0` for all `n >= 1`.
* Put `a_k = cos(2kL phi) - cos(2kL beta) = 2cos^2(kL phi) - 2cos^2(kL beta)`. Then `a_k` lies in `[0, 2]`.
* The Cesàro mean of `cos(k x)` is zero for every `x not in 2 pi Z`, rational or irrational, by the geometric sum.
  So `mean(a_k) = 0`. Then `0 <= a_k^2 <= 2 a_k` gives `mean(a_k^2) = 0`.
* Expand `a_k^2 = cos^2(2kL phi) + cos^2(2kL beta) - cos(2kL(phi+beta)) - cos(2kL(phi-beta))`. The first two terms
  have mean 1/2 each, because `4L phi`, `4L beta` are not in `2 pi Z`.
* **The point in question.** The mixed term has mean zero whenever `2L(phi +- beta) not in 2 pi Z`, including the case
  where `2L(phi +- beta)` is a rational multiple of `pi` other than a multiple of `2 pi`.
  * If both `L(phi + beta)` and `L(phi - beta)` are not in `pi Z`, the mean of `a_k^2` is one. That is a
    contradiction.
  * Hence `L phi = eps L beta (mod pi)`, i.e. `phi = eps beta + pi t/L (mod pi)`.
  * Astra's "both frequencies have zero mean" is justified exactly in the case she treats, `L phi != +- L beta mod pi`.
    Rational non-zero frequencies cause no problem.
* All Cesàro means used exist, because the terms are finite trigonometric combinations.

### (e) Removing `t` (Proved)

* `Re z != 0`. Otherwise `z = iy` with `y^2 = |z|^2 = p_j` (or `p_i p_j`), which is not a square.
* `B` = E-part of `|Re z|`, positive and finite. Let `L'` be the lcm of `L` and the periods of `z` in
  `(Z[i]/l^(v_l(Re z)+1))^*` for `l in E`, `l not | |z|^2`. For `n = 1 mod L'`, `z^n = z mod l^(v+1)`, so
  `v_l(Re z^n) = v_l(Re z)`. This includes `l = 2`.
* At `l | |z|^2` the real part is a unit mod `l`, as in (b).
* The exterior part of `Re z^n` is charged at unit `-1`, because `(z/conj z)^n = -1` there. No unit factor occurs,
  since `r_j/conj r_j = z/conj z` exactly. This gives `|cos(n phi)| >= |cos(n beta)|/B` for `n = 1 mod L'`.
* `n phi mod pi` is dense along the progression, since `L' phi/pi` is irrational. Along `n phi -> pi/2`, the
  congruence `n beta = eps(n phi - pi t/L) mod pi` holds because `n t/L = t/L mod 1`. So
  `|cos(n beta)| -> |sin(pi t/L)|`, while the left side tends to 0. Hence `L | t`.
* The only use of (G) is the irrationality of `phi/pi`. `Re z^n != 0` for all `n` because `beta/pi` is irrational.

### (f) Edges, trees, recovery (Proved)

* With `phi_i = eps_i beta_i`, `phi_j = eps_j beta_j` and `phi_i + phi_j = eps_(ij)(beta_i + beta_j) (mod pi)`,
  suppose the three signs are not all equal. Then exactly one differs from the other two. That gives
  `2 beta_i = 0`, `2 beta_j = 0` or `2(beta_i + beta_j) = 0 (mod pi)`.
* Each is impossible by irrationality, or equivalently by (G) applied to `2e_i`, `2e_j` or `2(e_i+e_j)`. A spanning
  tree propagates the equality of signs.
* Recovery: `phi_j = arg rho_j + pi m_j`. Then `<2a_x - e, phi> = <2a_x - e, arg rho> + pi sum_j e_j m_j (mod 2pi)`,
  because `2a_xj m_j pi = 0 (mod 2pi)`. The shift is the same for all marks and is absorbed into `theta_0`.
  `delta_x` and `k_x` are kept, and the drifting common phase is not set to zero.
* `z_x = prod rho_j^(a_xj) conj(rho_j)^(e_j - a_xj)` are distinct Gaussian integers of norm `N` with
  `arg z_x = theta_0' + delta_x + (pi/2)k_x`. So `i^(-k_x) z_x` is an actual cluster on an arc of length
  `<= C sqrt R`.
* The tree is needed only to make the *given* residue data those of the cluster, up to one global conjugation.
  Realisation by an actual cluster holds without it.

### (G) cannot be dropped (Proved; my addition)

Take `r_j = pi_j` genuine and `phi_j = 0`. Then `m_(n e_j,2) <= |Re pi_j^n| <= p_j^(n/2)`, so (A) holds for all `n` and
all cutoffs, yet `phi_j != +- arg pi_j`. Step (c) survives without (G); steps (d) and (e) do not.

## 3. Computations (Computed, `check_rigidity.py`)

1. **Astra's control reproduced.**
   * `(2+i)^6 = -117 + 44i` and `117 = 3^2 * 13`.
   * `((2+i)/(2-i))^6 = -1` mod 9 and mod 13, but not mod 27 or mod 169, so the collision powers are exactly 9 and 13.
   * `|cos(6 * beta/3)| = |cos 2 beta| = 3/5 < 117/125`, so (A) fails for every cutoff `X >= 13` (`A_3 >= 2` needs
     `X >= 9`).
2. **Second control: genuine angle, fake lift.**
   * `r = (2+i)(3+2i)/(3-2i)` has norm 5, `alpha = (2+i)(3+2i)^2 = -2 + 29i` and `J = 845`. Take `phi = arg(2+i)`;
     (G) holds for this angle because it is genuine.
   * At `n = 2`: `Re alpha^2 = -837 = -3^3 * 31`, with `E = {2, 5, 13}`. The required bound is `837/5 > 1 >= 3/5`.
     So (A) fails for `X >= 31`.
   * With `phi = arg r`, the test already fails at `n = 1` with `m = 1`. That is a failure of (G) itself, so it is not
     a meaningful control.
3. **Clean real parts.** For `alpha = -2 + 29i` with `E = {2,3,5,7,13}` (`L = 48`) and `E = {2,3,5,11,13,17,29}`
   (`L = 1680`), `Re alpha^(Lk)` is coprime to `E` for `k <= 25`.
4. **Finite-cutoff lifts (Proposition 3 below).** At cutoff 40 and label 5, five random factored data sets satisfying
   (d) at all odd prime powers `<= 40` were each lifted exactly by `(2+i) gamma/conj(gamma)`. Each lift has norm
   exactly 5, and `Norm(gamma)` has 25-26 digits.
5. **Sanity check.** The genuine `2+i` with `phi = beta` passes the ray tests for `n <= 119`.

## 4. Meaning

**Proposition 3 (finite-cutoff lifts; Proved).** Fix `X`. Any data satisfying sharp.md Definition 1.1(d) at the odd
prime powers `q^(A_q) <= X`, `q != p_j`, are reductions of `r_j = pi_j gamma/conj(gamma)` for some `gamma in Z[i]`.
This lift has exact norm `p_j`.
* *Proof.* `r_j(q)/pi_j` has norm 1 mod `Q = q^(A_q)`. By sharp.md Prop. 2.1(2) it equals `g_q/conj(g_q)` for a unit
  `g_q`. Choose `gamma = g_q mod Q` for all `q` by the Chinese remainder theorem.
* The height of `gamma` grows like `prod Q`, i.e. like `e^((1+o(1))X)`.

Consequences:
* **Relation to P3-all\*.** P3-all\* (sharp.md §6.7) is posed at cutoff `X = 3M`. There the rational-lift hypothesis is
  automatic (Proposition 3), and IX.7's mechanism is unavailable.
  * The mechanism needs collision moduli `|Re alpha^n| ~ J^(n/2)` for unbounded `n`, and the angular step needs `n`
    along which `n phi` approaches `pi/2`, with no rate.
  * So IX.7 neither caps P3-all\* nor gives a route to a lower bound in it.
  * A proof of growth from P3-all\* does not have to pass through IX.7, and it cannot.
* **What (L) really says.** It is a hypothesis on the projective limit over all cutoffs: one global element for all
  heights.
  * Compatible all-height witnesses form `pi_j * prod_q T(Z_q)` (uncountable). Rational lifts are countable. So
    generic adelic witnesses have no lift.
  * IX.7 says nothing about them, as Astra says.
  * The answer to her closing question ("can the local witnesses be forced to share a rational lift?"):
    * at a finite cutoff, always yes, trivially, with height `e^(O(X))`;
    * at all heights, not by local data. A lift is forced only by a global input such as "the witnesses are residues
      of actual Gaussian integers", which is the problem itself.
* **Relation to the record.** This is the same phenomenon as:
  * composite.md Prop. 3.1 (entire-L collapse: completely multiplicative, all-height smooth sums force
    `phi_p = +- theta_p mod pi/2`, "the problem itself");
  * the round 31 sieve rigidity (all-height multiplicative sharp sieve forces genuine angles mod `pi/4`;
    round31/README.md §1.6);
  * sharp.md Remark 2.3(4) (full residues determine `pi_j`).

  IX.7 is the residue-language version.
  * New: the lift need not be `pi_j`, and the strengthened test at all heights forces it to be. The proof is
    elementary (no L-functions, no Ratner), and it gives the angle conclusion mod `pi`.
  * It turns the coherent-design discussion of sharp.md §6.6 (finite cutoff, heuristic) into a rigorous all-height
    statement: coherent designs that pass all-height e-all at unit `-1` on the coordinate rays are actual clusters.
  * It belongs in the record's list "characterisations of genuineness, not bounds" (round31/README.md §3.2).
* **For the growth rate it gives nothing.** Its conclusion is "the data is an actual cluster". Any bound for such data
  is the original problem. No cap is removed and no input in README_VIII §3(e) changes status. The open items Astra
  lists (parent-only rate at `P_+ <= M (log M)^A`, unrestricted rate, constant 2/3) are unaffected.
