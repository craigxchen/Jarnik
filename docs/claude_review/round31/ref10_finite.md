# Referee report, Part X, block "finite": X.1 (all four PROVED paragraphs) and the X.5 least-cutoff table

Referee scripts (own code, exact integer arithmetic, Astra's scripts read but not imported):
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref10_finite/cut.py` (table),
`.../ref10_finite/wide.py`, `.../ref10_finite/wide2.py` (wider denominator class),
`.../ref10_finite/x1check.py` (order lemma and X.1 mechanism sanity checks). Outputs: `r500.txt`, `r1000.txt`, `r1000f.txt` in the same directory.

## 0. Verdict

1. **X.1 finite rigidity: Proved (CORRECT).** Every step was checked. The order step is right as written: if `n = 2^k` and `u^n` has order 4, then `ord(u) = 4n` exactly. "No `l` has `ord = 4n`" and "`4n` does not divide `ord`" are therefore the same condition.
2. **Direction pinning: Proved (CORRECT).** The equality case gives `phi = theta (mod pi)` or `phi = pi/2 - theta (mod pi)`, with `theta = arg r`. That is the statement of round 14 Lemma P0 and Proposition P. Astra's "`phi = +-arg r (mod pi/2)`" is a correct but weaker consequence.
3. **Dense families: Proved (CORRECT).** The cutoff `2 sqrt2 P^{1+2h}` and the density `1 - 2^{1-T_0}` were both checked.
4. **Scope of the height parameter: Proved (CORRECT).** The coordinator's literal conclusion "`gamma/conj(gamma) = +-1`" is false on genuine circles. This corrects the question, not Astra.
5. **X.5 table: Computed (CORRECT; all four rows reproduced exactly).** The rows are 439, 3881, 8681 and 5783. The case counts are 190, 2556, 8911 and 2556. The survivor at 438 is `985 - 1158i`, denominator 113, label 181, angle `arg(2+i)`. Two scope remarks apply (§3):
   - The table covers only lifts whose reductions are defined at every exterior prime up to `X` (denominator made of parent primes). That class is narrower than the hypothesis of X.1.
   - The killing exponent is `n = 1` in almost every case. So the table measures a low-exponent mechanism, not the dyadic one in X.1. The computed cutoffs are smaller than the proved bound (*) by factors up to about `10^570`.
6. **Meaning for P3-all\*: nothing about `W`.** X.1 gives a lower bound on the height of any rational lift of a nongenuine datum that passes the two ray tests. Lift height is not charged to `W`, and every finite datum has a lift of height `e^{O(X)}`. So X.1 excludes no fake of weight `O(M log M)`. It excludes only fakes with polynomial-height coherent witnesses, and only on parents with bounded `v_2(l - 1)` and labels below `X^{1/2^{T-1}}` (§4).

## 1. The test being applied

**Definition (sharp.md Def. 1.1 (e-all), (1.2)).** Fix the profile, the angles `phi in R^r` and a cutoff `X`. Write `A_q = max{a : q^a <= X}` for odd `q`. For `v in Z^r \ {0}` and `s in Z/4`, condition (1.2) reads

  `|sin(<v,phi> - s pi/4)| >= nu_s m_(v,s) e^{-w(v)/2}`,

with:
- `nu_s = 1` for `s` even and `nu_s = 2^{-1/2}` for `s` odd;
- `w(v) = sum |v_j| log p_j`;
- `m_(v,s) = prod_{q <= X odd, q not dividing N} q^{a_(q,s)(v)}`;
- `a_(q,s)(v) = max{a <= A_q : prod_j (r_j/conj r_j)^{v_j} = i^s mod q^a}`.

Only exterior odd primes are charged. Parent primes (`q | N`) and 2 are never charged.

**The case `v = n e_p`.** Here `e^{w(v)/2} = p^{n/2}`. For `s = +-1` the condition reads

  `|sin(n phi_p -+ pi/4)| >= m_(n e_p, +-1) / sqrt(2 p^n)`.

Unit `-i` is `s = 3`, and `|sin(t - 3pi/4)| = |sin(t + pi/4)|`.

**Which integer is charged at which unit.** Put `z = alpha^n = x_n + i y_n`, where `conj(alpha)` is a unit mod `q`.
- Unit `i`: `z = i conj(z)` mod `q^a` if and only if `q^a` divides `z - i conj(z) = (x_n - y_n)(1 - i)`, that is, `q^a | x_n - y_n` (`q` odd).
- Unit `-i`: the condition is `q^a | x_n + y_n`.
- Units `1` and `-1` charge `y_n` and `x_n` respectively.

So `x_n - y_n` belongs to unit `i` and the angle `n phi - pi/4`, and `x_n + y_n` belongs to unit `-i` and `n phi + pi/4`. This is exactly Astra's pairing.

**The charged modulus.** When every prime power exactly dividing `x_n -+ y_n` is exterior and at most `X`, the charged modulus is at least `|x_n -+ y_n|`; it equals it if the residues are the reductions of `alpha/b`. Her two displayed inequalities are then exactly (1.2) at units `i` and `-i`, with the full charged modulus and the record's normalisation `nu_s p^{-n/2} = 1/sqrt(2 p^n)`. **Proved.**

## 2. X.1, checked step by step

**Quantifiers.** First the parent `S` is fixed (all labels split, `p in S`). Then the height bound `H` on `Norm gamma`. Then the cutoff `X >= sqrt2 (p H^2)^{n/2}` with `n = 2^{T-1}`, `T = max_{l in S} v_2(l - 1)`. The conclusion holds uniformly over every lift `r = pi gamma/conj(gamma)` with `Norm gamma <= H` and over every real `phi_p`. It does not depend on `M`, `C`, the profile, (G), or the residues of the other labels. Only the residues of label `p` at exterior prime powers up to `X` enter, and only at the two tests `(n e_p, +-1)`.

**Set-up (Proved).**
- *Lowest terms and the bound on `b`.* `r = pi gamma^2 / Norm(gamma) = alpha/b` with `b > 0` minimal. Then `b | Norm gamma`, so `b <= H`, and `J = Norm alpha = p b^2`.
- *`gcd(x, y) = 1`.* Let a rational prime `q` divide `alpha`. If `q | b`, then `b` is not minimal. If `q` does not divide `b`, then `q^2 | p b^2` forces `q^2 | p`, which is impossible. This also covers `q = 2`.
- *`J` is odd, and `b` is never even.* If `b` were even, then `4 | J`, so `x` and `y` would both be even, a contradiction. An even `Norm gamma` cancels; for example `gamma = 1 + i` gives `gamma/conj(gamma) = i` and `b = 1`.
- *Disjoint supports.* An inert `q | J` divides both `x` and `y`, which is excluded. So every prime of `J` splits as `q = pi_q conj(pi_q)`, and exactly one of the two factors divides `alpha` (both would give `q | alpha`). Hence `alpha` and `conj(alpha)` have disjoint supports.
- *`gamma` sharing a factor with `pi` or `conj(pi)`.* `gamma = conj(pi)` gives `r = conj(pi)` (`b = 1`). `gamma = pi` gives `r = pi^3/p` (`b = p`, `J = p^3 = p b^2`). Here the parent prime `p` divides `b`, which is harmless because parent primes are never charged and the residues there are not used.

**Key step (Proved).**
- *Parity and coprimality.* `x_n^2 + y_n^2 = J^n` is odd, so `x_n` and `y_n` have opposite parity. A common prime factor would divide `J`, so it would split, and `conj(pi_q)` would divide `alpha^n`, which is impossible. So `gcd(x_n, y_n) = 1`.
- *Nonvanishing.* `x_n -+ y_n = 0` would mean `(alpha/conj alpha)^n = +-i`, that is, `alpha/conj(alpha)` is a root of unity. With disjoint supports and `J >= 5` this is impossible.
- *Size.* `|x_n -+ y_n| <= sqrt2 |alpha|^n = sqrt2 J^{n/2}`.
- *Coprime to `J`.* If `q | J` and `q | x_n - y_n`, then the factor of `q` dividing `alpha` divides `alpha^n - i conj(alpha)^n`, hence divides `conj(alpha)^n`, which is impossible. The same argument works for `+`. In particular `b` is invertible at every prime that is charged.
- *Parent primes `l` not dividing `J`.* `l | x_n -+ y_n` gives `u^n = +-i` in the norm-one group `T_l`, where `u = alpha/conj(alpha) mod l`. `T_l` is cyclic of order `l - 1`, and `i` has order 4 there.
- *Order lemma (Proved).* Let `n = 2^k`, `ord(u) = 2^j m'` with `m'` odd, and suppose `ord(u^n) = ord(u)/gcd(ord(u), n) = 2^{j - min(j,k)} m' = 4`. Then `m' = 1` and `j = k + 2`, so `ord(u) = 4n` exactly. The odd part cannot survive because it would appear in `ord(u^n)`. Since `u^{4n} = 1` always, `ord(u) | 4n`; "`ord(u) = 4n`" and "`4n | ord(u)`" are therefore equivalent. Astra's formulation is correct as written. A brute-force check found no counterexample in `Z/N` for `N < 400` and `k <= 5`.
- *Choice of `n`.* `4n = 2^{T+1}` does not divide `l - 1` for any `l in S`, so no parent prime divides `x_n -+ y_n`.
- *Full charging.* Every prime power exactly dividing `x_n -+ y_n` is odd, exterior, prime to `J` (so the residue is the reduction of `alpha/b`) and at most `sqrt2 J^{n/2} <= sqrt2 (p H^2)^{n/2} <= X`. Both integers are therefore fully charged, at units `i` and `-i` respectively.
- *Conclusion.* Square and add:
  - left side: `sin^2(n phi - pi/4) + sin^2(n phi + pi/4) = 1`, since `sin(t + pi/4) = cos(t - pi/4)`;
  - right side: `((x_n - y_n)^2 + (x_n + y_n)^2)/(2 p^n) = J^n/p^n = b^{2n}`.

  Hence `b = 1`, `Norm r = p`, and `r in {u pi, u conj(pi)}`. No choice of `phi` is involved.
- *Sanity check (Computed, `x1check.py`).* About 32,000 random lifts on three parents (`n = 2`, `32`, `128`): no parent-prime divisor appeared, and the parity, coprimality, size and sum-of-squares identities held.

**Arbitrary residues at the primes of `b`.** These primes divide `J`, so they never divide `x_n -+ y_n`. The remark that arbitrary residues are allowed there is therefore correct.

**The alternative `n = min(2^{T-1}, 2^{|S|})` (Proved).**
- The `|S| + 1` orders `2^{k+2}` (`k = 0, ..., |S|`) cannot all be attained by `|S|` primes, so some `n = 2^k <= 2^{|S|}` works for the given lift.
- By monotonicity of (*) in `n`, `X >= sqrt2 (p H^2)^{2^{|S|-1}}` is a uniform cutoff.
- As Astra says, the exponent `n` itself depends on the lift.

**Direction pinning (Proved, with the precise conclusion).**
- *Full charging at `n = 1`.* Take `b = 1`, `r = x + iy`, and the hypothesis `min S > sqrt(2 max S)`. Every prime factor of `|x +- y| <= sqrt(2p)` is smaller than every parent prime. At `X >= sqrt(2p)` both integers are fully charged.
- *Equality case.* The sum of squares is now `1 >= 1`, so both inequalities are equalities. Since `(x - y)/sqrt(2p) = -sin(theta - pi/4)`, this gives `sin^2(phi - pi/4) = sin^2(theta - pi/4)`, that is, `sin 2phi = sin 2theta`.
- *Precise conclusion.* `phi = theta (mod pi)` or `phi = pi/2 - theta = arg(i conj r) (mod pi)`. This is exactly the conclusion of Lemma P0 and Proposition P in round 14 (`round14/c1construct.md`). Proposition P assumes `X >= a_j + b_j` and no parent prime dividing the `D_s`, and works test by test. Astra's route by sum of squares uses only the units `+-i`, and her hypothesis `min S > sqrt(2 max S)` gives condition (i) of Proposition P for those `s`. Her "`phi = +-arg r (mod pi/2)`" is weaker than what is proved but true, and her description of it as the unpolluted case of Proposition P is accurate.
- *The `b > 1` remark (Proved).* For a candidate with `sqrt(2p) b < min S`, the same `n = 1` argument applies at `X >= sqrt(2p) b`.

**Dense families (Proved).**
- *Primes `5 mod 8`.* `l = 5 mod 8` gives `v_2(l - 1) = 2`, so `T = 2`, `n = 2`, and `(p H^2)^{n/2} = p H^2 <= 2P * P^{2h}`. The cutoff is `2 sqrt2 P^{1+2h}`.
- *Pinning for `P > 4`.* `min S >= P > 2 sqrt P >= sqrt(2 max S)` exactly when `P > 4`, and the cutoff exceeds `sqrt(4P)`.
- *Density.* By the prime number theorem in progressions, the share of primes `1 mod 4` with `v_2(l - 1) = k` is `2^{1-k}` (`k >= 2`). So `v_2 <= T_0` has relative density `1 - 2^{1-T_0}`, and "half" for `5 mod 8`. These densities are asymptotic in `P`.
- *Cutoff exponent.* It is `n/2 * (1 + 2h) = 2^{T_0-2}(1 + 2h)`, which depends only on `(T_0, h)`.

**Scope (Proved).**
- Genuine witnesses arise from twists that are not `+-1`: `gamma = conj(pi)` gives `r = conj(pi)`, and `gamma = 1 + i` gives `r = i pi`.
- `gamma/conj(gamma) = +-1` holds only for `gamma in Z u iZ`. The literal conclusion asked by the coordinator is false. The correct conclusion is `r in {u pi, u conj pi}`, which is what Astra proves.

## 3. The X.5 table (Computed, reproduced)

**What "least cutoff" means.**
- Fix the parent `S` (primes `1 mod 4`, or `5 mod 8`, in `[P, 2P]`).
- *Class of lifts.* Every nonintegral `r = pi gamma/conj(gamma)` with `Norm pi = p in S` and `Norm gamma <= 2P` whose reduction is defined at every exterior odd prime up to `X`.
- *Tests.* (1.2) for `v = n e_p`, all `n >= 1` and all four units, with the charged modulus `gcd(D_s(alpha^n), L(X))`. Here `L(X) = prod_{q <= X odd, q not in S} q^{A_q}`, and `D_s = y, x - y, x, x + y` for `s = 0, 1, 2, 3`.
- *Survival.* A lift survives at `X` if some real `phi` passes all `(n, s)`.
- *Least cutoff.* The least `X` at which no lift in the class survives.

**Why the experiment is exact.**
- *Killing.* At `X`, each case is killed by a single exponent `n`. For that `n` the four squared inequalities form a rectangle for `p^n (cos 2n phi, sin 2n phi)` that misses the circle of radius `p^n`. That kills the case for every real `phi`.
- *Survival.* At `X - 1`, one explicit lift with `phi = arg(2+i)` passes every `n`. Each comparison is the integer inequality `E_s^2 p^n >= m^2 5^n` with `E_s = D_s((2+i)^n) != 0`. Beyond the first `n` with `p^n >= 5^n L^2`, the comparisons hold automatically.
- *Minimality.* The charged moduli are monotone in `X`.

**Reduction of the class.**
- When `X - 1 >= 2P`, every prime of `b <= 2P` is at most `X`, so it must be a parent prime. Since `P^2 > 2P`, `b = l in S`.
- Up to units and conjugation, `alpha = pi_p pi_l^2` or `pi_p conj(pi_l)^2`. The single non-primitive case at `l = p` is dropped. This gives `s(2s - 1)` cases.
- *Independent enumeration.* I did not rely on this reduction. I enumerated all `gamma` with `Norm gamma <= 2P` and both primes of norm `p`, reduced, and kept the local-unit ones. This recovers exactly `s(2s - 1)` classes, with denominator set exactly `S`.

**Reproduced rows (own code `cut.py`).**

| parent | P | s | cases | least cutoff | survivor at X - 1 (label, b, alpha, phi; last n checked) | killing n at X | T | n of (*) | (*) |
|---|---:|---:|---:|---:|---|---|---:|---:|---|
| all split | 100 | 10 | 190 | **439** | 181, 113, 985 - 1158i, arg(2+i); 211 | n=1: 187, n=2: 3 | 6 | 32 | about 10^{111} |
| all split | 500 | 36 | 2556 | **3881** | 997, 557, -17482 - 1923i, arg(2+i); 1374 | 1: 2549, 2: 4, 5: 3 | 8 | 128 | about 10^{576} |
| all split | 1000 | 67 | 8911 | **8681** | 1601, 1381, -39869 - 38260i, arg(2+i); 2826 | 1: 8891, 2: 12, 3: 6, 4: 1, 5: 1 | 7 | 64 | about 10^{317} |
| 5 mod 8 | 1000 | 36 | 2556 | **5783** | 1621, 1493, -49085 - 34698i, arg(2+i); 1903 | 1: 2545, 2: 9, 3: 2 | 2 | 2 | about 1.1 x 10^{10} |

All four least cutoffs agree with Astra's, and the first-row survivor is exactly hers. The genuine primes pass at each tabulated `X`: the single-exponent rectangles are nonempty for `n <= 30` or `40`, and the genuine primes pass with equality structure at `phi = theta`.

**Scope remarks (Computed).**
- *(a) The proved bound is not what the table tests.* The least cutoffs are 4 to 9 times `P`. The bound (*) is astronomically larger, except on the `5 mod 8` parent, where it is still about `2 x 10^6` times larger. Almost every case is killed at `n = 1`, by partial charging of the coordinates of `alpha = pi_p pi_l^2` itself, not by the dyadic exponent of X.1. The table is consistent with X.1 but is not evidence for its mechanism or for its sharpness. Astra says (*) is "sufficient, not optimal", which is correct.
- *(b) The class is narrower than X.1's.* The table's class ("all exterior local reductions under test must be defined") excludes denominators with exterior prime factors. X.1 itself allows them, with arbitrary residues at those primes. I tested the wider class on the first parent (630 classes of nonintegral lifts with `Norm gamma <= 200`, `wide.py`). I charged nothing at the primes of `b`, which is the choice most favourable to survival and is not necessarily realisable. Two lifts are not killed by any single exponent `n <= 40` at `X = 439`: label 157 with `b = 185 = 5 * 37`, and label 173 with `b = 145 = 5 * 29`. Under that model they die at 769 and 443 (`wide2.py`). So 439 is the least cutoff for the local-unit class only. For X.1's full class it is not certified; it is at most 769 under the zero-charge model.

## 4. Meaning for P3-all\* (what is and is not implied)

**Proved (contrapositive of X.1).**
- On a fixed parent with `T = max v_2(l - 1)` and `n = 2^{T-1}`, consider a nongenuine label datum at cutoff `X` that passes the two tests `(n e_p, +-1)`. Every rational lift `pi gamma/conj(gamma)` of it has `Norm gamma >= b > ((X/sqrt2)^{2/n}/p)^{1/2}`.
- On `5 mod 8` parents this reads `b > (X/(sqrt2 p))^{1/2}`, which is item 8 of her summary on that class.
- At `X = 3M^A` with labels `p ~ M^{h_0}`, the bound is a positive power of `X` only when `A > 2^{T-2} h_0`, that is, `A > h_0` for `T = 2`. Otherwise it is vacuous.
- The `n = 1` variant needs `sqrt(2p) b < min S`. It applies only for `b < sqrt(P/2)` on comparable parents, at `X >= sqrt(2p) b` (regime `A > h_0/2`).
- For parents with large `T`, the uniform cutoff is exponential in `2^{T}` or `2^{|S|}` and says nothing at polynomial cutoffs.

**Not implied (Proved, by the record).**
- *Nothing about `W`.* The weight of a model is the physical parent weight. The height of a rational lift of its residue data is an auxiliary quantity that is never charged to `W`.
- *Every datum has some lift.* Every finite-cutoff factored datum has a lift `pi_j gamma/conj(gamma)` of height `e^{O(X)}` (README_IX §1.8, finite-cutoff lift proposition). So X.1 forbids no P3-all\* fake whose residues have only high-height lifts, and the generic independent residue choices of sharp.md §6 are of that kind.
- X.1 therefore shows: on comparable parents with bounded `v_2(l - 1)` and labels below `X^{1/2^{T-1}}`, a fake cannot be built from polynomial-height coherent witnesses. (A coherent witness is a global `r_j` with small twist; IX.3-type witnesses already have astronomical height.) It gives no lower bound on `W_min^all`.
- It also bears only on singletons and says nothing about the signed words that X.2 and X.4 concern.
- P3-all\* at `X = 3M` or `3M^A` remains **open**, as Astra states.

**Consequence when the hypothesis does hold (Proved).** Suppose a member of P3-all\* on a comparable `5 mod 8` parent has, for every label, residues that are reductions of a lift of twist norm at most `H`, at a cutoff `X >= 2 sqrt2 P^{1+2h}` (`H <= P^h`). Then:
- every `r_j` is genuine and every angle is pinned;
- with (G) and the profile, Proposition P of round 14 makes it an actual cluster, so `W >= 3M log M - O(M)` by two_thirds.

This is a conditional statement about witnesses of bounded height. It is not a bound on `W_min^all`.
