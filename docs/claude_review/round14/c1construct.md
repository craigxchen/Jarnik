# Round 14 — C1: a fake circle for P3-all* in regime (1.4)? (construction side)

Status: unrefereed draft (round 14). Scripts: `c1construct_checks/`. Notes: `c1construct_notes.md`.

## 0. Summary of results

**Verdict: no C1 construction and no obstruction to the class.** No fake circle in regime (1.4) was built, and no
lower bound on `W_min^all` beyond `3 M log M - O(M)` was proved; the growth status is unchanged
(`M <= (2/3 + eps) log R/loglog R`; P3-all\* open and uncapped in regime (1.4)). What was proved (unrefereed) is a
complete account of where every known construction breaks once `p_j <= X`, a set of barriers that exclude every
known mechanism of cross-modulus coordination as a C1 design, an `A`-uniform localisation of the monomial
obstruction to the deep 2-adic classes, and a norm-1 pinning theorem that sharpens the rigidity side.

1. **Break analysis (§2; proved, with exact finite instances).** In regime (1.4) a monomial of norm
   `k < 2A/h_0` dies from one collision at one large modulus.
   * `sharp.md` Thm 6.4's random fibre fails at **norm 1**: each non-twin column survives with probability
     `<= ((h_0/2 + o(1))/A)^2` (Thm B1); exact survival products 0.015-0.051 at `M = 30, 50`, `A = 3-4`, overshoot
     `~ (A - h_0/2) log M`.
   * Torsion-avoiding fibres fail at **norm 2** (Thm B2).
   * The bounded-norm systems that refute L4(K) (`L4K.md` Thm 2) fail at **norm `K + 1`, deterministically for
     every choice of permutations**: `>= r^3/36 - O(r^2)` norm-4 monomials vanish identically at every modulus above
     `q_T`, so some monomial has demand `>= (2/(9r)) theta(X)`, a factor `M^(A-1)` over budget (Thm B3; exact
     instances: 630/630 and 5460/5460 monomials over budget).
   * Coherent designs with root repair are certified on the shallow classes and not on the deep ones; the
     flexible tops of round 10 cannot seed a fake over a genuine core (Prop. P).
2. **Norm-1 pinning (Prop. P, Cor. P1-P2; proved).** If `X >= a_j + b_j` (`pi_j = a_j + i b_j`), labels that are
   genuine per unit at the prime powers dividing `a_j`, `a_j +- b_j` (all `<= sqrt(2 p_j) <= M^(h_0/2)`) pin
   `phi_j` exactly to the angle of a Gaussian integer of norm `p_j`; if this holds for all columns (and no `p_l`
   divides those numbers), the member **is an actual cluster**. So P3-G collapses in regime (1.4) at norm 1, a C1
   fake must be non-genuine below `M^(h_0/2)`, and a non-actual member is a union of actual clusters of a
   sub-circle glued by its unpinned columns. This is the fake-violated "label-angle" input that round 10's referee
   showed route 3 needs, at norm 1, given per-unit genuine collisions at the pinning moduli.
3. **Coherence (Lemma L, Prop. G, Thm C; proved given Linnik).** Parity matching on any deep class
   `q = -1 mod 2^(E+1)` forces `p_j Norm zeta` to be a square or twice a square (unconditional form of a GRH
   statement of round 7); exact coherence plus the norm-1 height forces `zeta_j = i^t pi_j` or `i^t conj(pi_j)`;
   hence **every exactly coherent member in regime (1.4) is an actual cluster** (generic cofactors).
4. **Pairs (Thm R; proved).** On the Paley profile, a coherent design meeting the height condition at
   `X >= M^4` has `C + 2 Hbar > (M - 1)/8` (round 7 had `M^(1/2 - o(1))`): pairs and columns cannot both be
   certified.
5. **Shallow classes (Thm S, Cor. S1; proved, exact checks).** With auxiliary Gaussian primes whose norms are
   `= 1 mod 2^(E_0+1)`, a root-repaired coherent design serves **every** shallow modulus `k(q) <= E_0` for **every**
   `X` at cost `(2^(E_0) - 2) sum_j |v_j| log l_j`, `l_j <= (1+o(1)) 2^(E_0) r log r`: for the monomials outside the
   character lattice, the parity obstruction lives only on the deep classes `q = +-1 mod 2^(E_0+1)`, uniformly in
   `A` (a construction-side localisation, not a cap: pairs and characters are not served).
6. **Residual problems (§3.5) and barrier logic (§5).** C1 reduces to DeepBal (balanced, non-coherent labels on
   the deep classes) and PairCompat; no mechanism is known (evidence: on the deep classes alone, root-repaired
   coherence and uniform labels already exceed norm-1 and norm-2 per-unit budgets at `X = 2*10^5`, genuine labels sit
   at ratio exactly 1). The implication "lower bound in regime (1.4) => growth" is correct with the quantifiers of
   (T'). A positive rigidity answer (Rig_pin suffices) would make P3-all\* **exact** in regime (1.4), turning the
   barrier question into L7 for genuine data; it would not by itself give growth.


## 1. Setting: the class P3-all*, regime (1.4), parameters

Conventions are those of `round7/README.md`, `round5/round6_sharp/sharp.md` (Def. 1.1, Prop. 2.1, Cor. 2.2) and
`round10/L4K.md` §1. Labels: *proved* (all quantifiers, here, unrefereed), *exact* (finite computation, script
named), *evidence* (computation no proof uses), *heuristic*. Nothing in this note is refereed.

* **Class.** P3-all\*`_X` `= F^all_X(M, C)` (`sharp.md` Def. 1.1): distinct primes `p_j = 1 mod 4` (`j <= r`),
  exponents `e_j`, distinct rows `a_x in prod [0, e_j]` (`x <= M`), real angles `phi` realising the profile
  (`<2a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x`, `delta_x in [0, Delta]`, `Delta = C e^(-W/4)`), the gap
  inequalities (G), factored residues `r_j(q) mod q^(A_q)` with `Norm r_j(q) = p_j` at every odd prime `q <= X`,
  `q != p_j` (lower levels are reductions), the character conditions (e-char) and the slab conditions (e-all)

  ```text
  |sin(<v, phi> - s pi/4)|  >=  nu_s m_(v,s) e^(-w(v)/2)      for every v in Z^r \ {0} and s in Z/4,        (1.2)
  ```

  `m_(v,s) = prod_(q <= X odd, q not dividing N) q^(a_(q,s)(v))`, `nu_0 = nu_2 = 1`, `nu_1 = nu_3 = 2^(-1/2)`.
  By `sharp.md` Cor. 2.2 the data at `q` are a **parity-respecting label vector** `ell(q) in (Z/m_Q)^r`,
  `ell_j(q) = nu_j(q) mod 2`, `nu_j(q) = [(p_j/q) = -1]`, `m_Q = (q - chi_(-4)(q)) q^(A_q - 1)`; `v` collides at
  `q^a` with unit `s` iff `<v, ell> = s iota mod m_(q^a)` (`iota = log i`). Taking `|sin| <= 1` in (1.2) gives the
  **per-unit budget** `log m_(v,s) <= B(v,s) := w(v)/2 + [s odd] (log 2)/2` (`sharp.md` Prop. 6.1, per-unit form
  of `referee_sharp.md`).
* **Classes of moduli.** `k(q) := v_2(m_q)` (`round7/p3Dconstruct.md` §3): `k(q) = 2` iff `q = +-3 mod 8`,
  `k(q) = 3` iff `q = +-7 mod 16`, and in general `k(q) >= E` iff `q = +-1 mod 2^E` (`q = chi(q) mod 2^E`).
  *Shallow* moduli: `k(q) <= E_0`; *deep* moduli: `k(q) > E_0`. Among odd primes the deep ones have relative
  density `2^(1 - E_0)`.
* **Regime (1.4)** (`round9/README.md` §0 item 1, §1): `W <= c_0 M log M`, all `p_j <= M^(h_0)` with
  `h_0 = 2 c_0 - 3 + o(1)`, `r <= rho M`, and `X = 3 M^A` with **`A` chosen after `c_0`**. In particular
  `X > sqrt(2 p_j)` for every `j` as soon as `A > h_0/2`, and `X` exceeds any fixed power of `p_max`.
* **C1** (the target; `round13/README.md` §7.2): for some fixed `c_0` and **every** `A`, for all large `M`, a
  member of P3-all\*`_(3M^A)` with `W <= c_0 M log M` and all `p_j <= M^(h_0)` that is **not** an actual cluster.
  (A weaker target: such members with `W = o(M log^2 M/loglog M)` and `p_j <= X`.)
* **Genuine data.** `pi_j = a_j + i b_j` with `a_j` odd `> 0`, `b_j` even `> 0`, `theta_j = arg pi_j in (0, pi/2)`.
  `D_0, D_1, D_2, D_3 (x + iy) = y, x - y, x, x + y`. For genuine labels, `v` collides at `q^a` with unit `s` iff
  `q^a | D_s(Z_v)` (`L4K.md` §1; `q` odd, `q` not dividing `N`).

## 2. Where the round-5/6 constructions break when p_j <= X

Every construction of rounds 5-7 has `p_j > X`; in regime (1.4) `X` is a large power of `p_max`. The single new
feature is then: **a monomial of norm `k` has budget `<= (k h_0/2) log M + 1`, which is below `log X` for
`k < 2A/h_0`, so one collision at one large modulus kills it.** This section records exactly which monomials
fail, and by how much, for each design family. (Moduli `q | N`, i.e. `q = p_l <= X`, are not in (1.2); every
statement below concerns only `q` not dividing `N`.)

### 2.1 The random fibre of `sharp.md` Thm 6.4 fails at norm 1

**Theorem B1 (proved).** Fix `h_0 > 0`, `A > h_0/2`, `b >= 33`. Take the profile and design of `sharp.md`
Thm 6.4 (Paley profile of §4, design `D_A` with the uniform random fibre, Def. 6.3), but with arbitrary primes
`p_j <= M^(h_0)`, and `X = 3M^A`. By `sharp.md` Lemma 6.3(2) the labels of the `r - 2M >= (b-2)M - b` non-twin
columns are, at each `q`, independent and uniform in their parity classes, independently over `q`. Then for
every non-twin column `j`,

```text
P( e_j meets its norm-1 per-unit budget )  <=  prod_(q = +-3 (8), sqrt(2 p_j) < q <= X) (1 - 4/m_q)
                                           <=  (1 + O(1/log M)) ( log sqrt(2 p_j) / log X )^2  <=  ((h_0/2 + o(1))/A)^2,
```

independently over `j`. Hence with probability `1 - exp(-Omega(M))` at least
`(1 - 2 (h_0/(2A))^2)(b-2)M` non-twin columns violate the height condition, so the design is not completable,
for every `A > h_0/2` (and the conclusion is uniform in `A` in the sense that the failing fraction tends to 1 as
`A -> infinity`).

*Proof.* At `q = +-3 mod 8` we have `m_q = 4 mod 8`, `iota = m/4` odd, and the 4-torsion `{0, m/4, m/2, 3m/4}`
has exactly two elements of each parity; a label uniform on its parity class (`m/2` elements) is torsion with
probability `4/m_q`, and then `e_j` collides at level 1 with some unit `s`. A collision at `q > sqrt(2 p_j)`
gives `log m_(e_j, s) >= log q > (1/2) log p_j + (1/2) log 2 >= B(e_j, s)`. Independence over `q` gives the
product. Mertens in progressions, `sum_(q <= x, q = +-3 (8)) 1/q = (1/2) loglog x + B + O(1/log x)`, and
`m_q <= q + 1` give `prod (1 - 4/m_q) <= exp(-4 sum 1/(q+1)) = (log Y/log X)^2 (1 + O(1/log Y))` with
`Y = sqrt(2 p_j)`; `log Y <= (h_0/2) log M + 1` and `log X >= A log M`. Independence over columns and Chernoff give
the count. QED

The excess is large: the last colliding modulus is typically of size `X^(1 - o(1))`, so the typical overshoot is
`(A - h_0/2) log M`. *Exact* (`break_random.py`, `break_random_out.txt`; exact survival products including the
`q = +-1 mod 8` factors computed from the actual Legendre symbols; the simulation is seeded and lists explicit
witnesses `(p_j, q, ell_j(q))` with `ell_j(q)` in the 4-torsion):

| `M` | `h_0` | `A` | `X` | exact mean `P(survive)` (all `q`) | `q = +-3 (8)` only | `(log sqrt(2p)/log X)^2` | violating columns (simulated) | median / max overshoot | `(A - h_0/2) log M` |
|---|---|---|---|---|---|---|---|---|---|
| 30 | 3 | 3 | 81000 | 0.051 | 0.223 | 0.232 | 193/200 | 4.48 / 6.17 | 5.10 |
| 30 | 3 | 4 | 2.43e6 | 0.018 | 0.132 | 0.137 | 198/200 | 6.89 / 9.46 | 8.50 |
| 50 | 2.5 | 3 | 375000 | 0.027 | 0.157 | 0.167 | 194/200 | 6.00 / 7.90 | 6.85 |
| 50 | 2.5 | 3.5 | 2.65e6 | 0.015 | 0.118 | 0.125 | 196/200 | 7.35 / 9.78 | 8.80 |

(The `q = +-1 mod 8` moduli roughly square the survival probability: there the torsion is even and is hit with
probability `8/m_q` when `nu_j(q) = 0`.)

### 2.2 Torsion-avoiding fibres fail at norm 2

The referee counterexample of `round7/p3Dconstruct.md` Thm 1 (labels uniform on the parity class **minus the
4-torsion**) removes all single-column collisions.

**Theorem B2 (proved).** Same setting, but at every `q` the non-twin labels are independent, each uniform on its
parity class minus the 4-torsion, independently over `q`. Let `A > h_0`. For disjoint pairs `{j, k}` of non-twin
columns, `P(e_j - e_k meets its budget) <= ((h_0 + o(1))/A)^2`, independently over disjoint pairs; so w.h.p. a
fraction `>= 1 - 2(h_0/A)^2` of `>= (b-2)M/2 - 1` disjoint pairs violate.

*Proof.* At `q = +-3 mod 8`, given `ell_k`, the monomial `e_j - e_k` collides iff `ell_j in ell_k + {0, m/4, m/2,
3m/4}`. Of these four values exactly two have the parity `nu_j(q)` (the shifts `0, m/2` are even, `m/4, 3m/4`
odd), and neither is torsion because `ell_k` is not and the torsion is a subgroup. So the conditional
collision probability is `2/(m/2 - 2) >= 4/m_q`. A collision at `q > sqrt(2 p_j p_k)` exceeds
`B(e_j - e_k, s) <= (1/2) log(2 p_j p_k)`. The rest is the proof of B1 with `Y = sqrt(2 p_j p_k) <= sqrt 2 M^(h_0)`. QED

The same argument kills any law that is, at each `q`, independent over columns and uniform on the parity class
minus a bounded set, at some bounded norm `K` with `A > K h_0/2`: removing `f` elements from the class leaves a
conditional collision probability `>= (2g - f)/(m/2)` for a family of `g` monomials `e_j - u e_k`.

### 2.3 The bounded-norm systems of `L4K.md` Thm 2 fail at norm `K + 1`, deterministically

`L4K.md` Thm 2 puts genuine labels at `q <= q_T = 4A_K + 5` and, at every `q > q_T`, the vector (2.1) there,
`a_j(q) = 2(Lambda + (K+1) beta_(sigma_q(j))) + nu_j(q)`, with any permutation `sigma_q` (`Lambda = K(K+1)B + K`).
Its Remark 2 says monomials of norm `> K` "collide about as often as for random labels". In fact they collide
far more often, for every choice of the permutations:

**Theorem B3 (proved; deterministic).** Let `K in {2, 3}` (`h = 1`, `beta_j = j`, `B = r`), and let
`F = {e_i + e_j - e_k - e_l}` (unordered disjoint pairs `{i,j}`, `{k,l}`; `|F| = r(r-1)(r-2)(r-3)/8`). For every
prime `q > q_T` not dividing `N`, every permutation `sigma_q` and every parity pattern `nu(q)`, at least

```text
Z_r := (1/6) sum_s N_s^2 - (1/2) binom(r, 2)  =  r^3/36 - O(r^2),     N_s = #{ {x, y} subset [r], x != y, x + y = s },
```

monomials `v in F` have `<v, a(q)> = 0` exactly, hence collide with unit `0` at **every** level `q^a <= X`.
Consequently `max_(v in F) log m_(v,0) >= (Z_r/|F|)(theta(X) - theta(q_T) - W)`, with `Z_r/|F| = (2/(9r))(1 + O(1/r))`,
while `B(v, 0) <= 2 log p_max` on `F`. In regime (1.4) (`r <= rho M`, `W <= c_0 M log M`, `X = 3M^A`, `A > 1`)
the maximal demand on `F` is `>= (2/(3 rho) - o(1)) M^(A-1)` against budgets `<= 2 h_0 log M`: the system misses
the height condition by a factor `M^(A - 1 - o(1))`.

*Proof.* `<v, a> = 2(K+1)(sigma_i + sigma_j - sigma_k - sigma_l) + (nu_i + nu_j - nu_k - nu_l)` and the second
term has absolute value `<= 2 < 2(K+1)`, so `<v, a> = 0` iff the `sigma`-sums and the `nu`-sums of the two pairs
agree. Group the unordered pairs of `[r]` by (`sigma`-sum `s`, `nu`-sum `tau in {0,1,2}`), with sizes `N_(s,tau)`;
since `sigma_q` is a bijection of `[r]`, `sum_tau N_(s,tau) = N_s` does not depend on `sigma_q`. Two distinct pairs
with the same `sigma`-sum are disjoint, so each group of size `N` yields `binom(N, 2)` monomials of `F` with
`<v, a> = 0`. By convexity `sum_tau binom(N_(s,tau), 2) >= 3 binom(N_s/3, 2) = N_s^2/6 - N_s/2`, and
`sum_s N_s = binom(r, 2)`. Since `ell = a mod m_Q`, `<v, ell> = 0` in every `Z/m_(q^a)`: unit `0`, every level.
Summing `A_q log q >= log q` over `q in (q_T, X]`, `q` not dividing `N` (whose `log`s sum to at most `W`), and
averaging over `F` gives the bound. `sum_s N_s^2 = r^3/6 + O(r^2)`. QED

*Exact* (`break_bh.py`, `break_bh_out.txt`): for `K = 2, 3`, `r = 10, 16`, random permutations at every `q` in
`(q_T, X]`, `X <= 60000`: **every** monomial of `F` violates its budget from the moduli above `q_T` alone
(630/630, 630/630, 5460/5460), e.g. `v = e_1 + e_9 - e_2 - e_8` (`K = 3`, `r = 10`, primes 3041, 5881, 3613, 5333)
has `log m_(v,0) >= 2010.9` against budget `16.74` and collides at 198 moduli. The minimum over `q` of the number
of vanishing `v in F` (13 for `r = 10`, 79 for `r = 16`) is above the proved bound `Z_r` (1.67 and 44.0 there).

So the systems that refute L4(K) are very far from fakes: they are built to be collision-free at norm `<= K`,
and the price is a rigid integer structure that makes norm-`(K+1)` relations exact integer identities, which
collide at every modulus. This is the Siegel trap (`sharp.md` Prop. 7.1) at norm `K + 1`.

### 2.4 Coherent designs with root repair: shallow classes certified, deep classes not

`round7/p3Dconstruct.md` §3 repairs the parity of a coherent design `zeta_j/conj zeta_j` by an element of order
`2^(k(q))`; collisions are then certified by `q^a | Im(gamma_v^(2^(k(q))))` (its Prop. 3.3). In regime (1.4) this
works uniformly in `X` on the shallow classes (Theorem S below, with a choice of the `zeta_j` that keeps their
norms off the shallow moduli). On class `k` the certified cost is `2^k log|gamma_v|`, which exceeds the budget per
unit of `||v||_1` as soon as `2^k >~ h_0` and certifies nothing at all once `2^k >~ X^(1/2)`; round 7's heuristic
H1 puts the actual cost on the deep classes at the random rate, and `deep_demand.py` (§4, evidence) finds norm-1
and norm-2 budgets exceeded on the deep classes alone. (The greedy systems of round 9's Prop. H, if completed
randomly above `Q_K`, fall under the argument of Theorems B1-B2 at norm `K + 1` once `A > (K+1) h_0/2`.) Two
further obstructions apply to *any* coherent design in regime (1.4): exact coherence with the right parity on a
deep class forces the genuine norm (Lemma L below; unconditional), and genuine-type exact coherence collapses to
actual clusters (Theorem C). On the Paley profile coherence also cannot certify pairs and columns at once
(Theorem R, linear in `M`).

### 2.5 The round-10 flexible tops cannot seed a fake as they stand

The flexible tops of `L4K.md` §4.4 are deformations of the genuine system **above** a genuine core. In regime
(1.4) every column is pinned by its norm-1 collisions at the prime factors of `a_j`, `a_j + b_j`, `a_j - b_j`, all
`<= sqrt(2 p_j) <= M^(h_0/2)` (Proposition P below). A system that is genuine (per unit) at those moduli has only
actual clusters as completions. So a C1 fake must deviate from the genuine labels **below** `sqrt(2 p_j)` on the
pinning moduli of its columns, i.e. inside the norm-`h_0` forced range, where round 10 found rigidity
(exhaustively, for `r = 4`, tiny primes, `K = 4`). A flexible top is useful to C1 only together with a
non-genuine core.

### 2.6 Summary of the break analysis

| design (source) | first failing monomials | norm | by how much (regime (1.4)) | status |
|---|---|---|---|---|
| random fibre (`sharp.md` Thm 6.4) | single non-twin columns `e_j` | 1 | one collision at `q ~ X^(1-o(1))`: overshoot `~ (A - h_0/2) log M`; failing fraction `>= 1 - (h_0/2A)^2` | Thm B1 |
| torsion-avoiding fibre (`p3Dconstruct.md`) | `e_j - e_k` | 2 | overshoot `~ (A - h_0) log M`; fraction `>= 1 - 2(h_0/A)^2` | Thm B2 |
| free `B_h` vectors above `q_T` (`L4K.md` Thm 2) | `e_i + e_j - e_k - e_l` | `K + 1 = 3, 4` | demand `>= (2/(9r)) theta(X)`, factor `M^(A-1)` | Thm B3, deterministic |
| coherent + root repair (`p3Dconstruct.md` §3) | none on shallow classes; uncertified on deep classes | -- | certified cost `2^k log|gamma_v|` on class `k` | Thm S; Lemma L; H1 heuristic |
| G-coherent on Paley profile | pairs or columns uncertified | -- | `C + 2 Hbar > (M-1)/8` | Thm R |
| genuine-type exact coherence | -- (no fake at all) | -- | collapses to actual clusters | Thm C |
| genuine core + flexible top (`L4K.md` §4.4) | -- (no fake at all) | -- | pinned at norm 1 | Prop. P |

## 3. Main results

### 3.1 Norm-1 pinning: genuine labels at the pinning moduli force an actual cluster

**Lemma P0 (one-sided pinning; proved).** Let `theta in (0, pi/2)`, `theta != pi/4`, and `S(theta) = {2, 3}` if
`theta < pi/4`, `S(theta) = {1, 2, 3}` if `theta > pi/4`. If `x in R` satisfies

```text
|sin(x - s pi/4)|  >=  |sin(theta - s pi/4)|        for all s in S(theta),
```

then `x = theta (mod pi)`, or `theta > pi/4` and `x = pi/2 - theta (mod pi)`.

*Proof.* All conditions are `pi`-periodic; work mod `pi`. `|sin y| >= |sin y_0|` iff
`dist(y, pi Z) >= dist(y_0, pi Z)`. `s = 2` reads `|cos x| >= cos theta`, i.e. `x in [-theta, theta] (mod pi)`.
*Case `theta < pi/4`.* `s = 3`: `dist(x + pi/4, pi Z) >= theta + pi/4` (note `theta - 3pi/4 = theta + pi/4 - pi`
and `theta + pi/4 < pi/2`). For `x in [-theta, theta]`, `x + pi/4 in [pi/4 - theta, pi/4 + theta] ⊂ (0, pi/2)`, where
`dist = x + pi/4`; so `x >= theta`, hence `x = theta`.
*Case `theta > pi/4`.* `s = 1`: `dist(x - pi/4, pi Z) >= theta - pi/4`. For `x in [-theta, theta]`,
`y = x - pi/4 in [-theta - pi/4, theta - pi/4] ⊂ (-pi, pi/4)`; on `[-pi/2, pi/4)` the condition is
`|y| >= theta - pi/4`, i.e. `x = theta` or `x <= pi/2 - theta`; on `[-theta - pi/4, -pi/2)` it reads
`pi - |y| >= theta - pi/4`, which holds, but there `x < -pi/4`. So `x = theta` or `x in [-theta, pi/2 - theta]`.
`s = 3`: `dist(x + pi/4, pi Z) >= 3pi/4 - theta` (as `theta + pi/4 > pi/2`). For `x in [-theta, pi/2 - theta]`,
`x + pi/4 in [pi/4 - theta, 3pi/4 - theta] ⊂ (-pi/4, pi/2)`, where `dist = |x + pi/4|`. Since
`theta - pi/4 < 3pi/4 - theta`, the condition forces `x + pi/4 = 3pi/4 - theta`, i.e. `x = pi/2 - theta`. QED

`|sin(theta_j - s pi/4)| = nu_s |D_s(pi_j)|/sqrt(p_j)` for every `s` (direct: `sin theta_j = b_j/sqrt p_j`,
`cos theta_j = a_j/sqrt p_j`). *Exact* (`pin_check.py`): this identity, as an identity of rationals for the squares,
for all 4783 split primes `<= 10^5`; and a 40000-point grid check of Lemma P0 for all 1125 split primes
`<= 2*10^4` (evidence; the proof does not use it).

**Proposition P (norm-1 pinning; proved).** Let a member of P3-all\*`_X` be given (any `M`, `C`, exponents,
rows). For each column `j` let `pi_j = a_j + i b_j` (`a_j` odd `> 0`, `b_j` even `> 0`), `theta_j = arg pi_j`,
`S_j = S(theta_j)`. Assume `X >= a_j + b_j` for all `j` and, for every `j`:

* (i) no `p_l` (`l <= r`) divides `prod_(s in S_j) D_s(pi_j)` (these are odd: `a_j`, `a_j + b_j`, `a_j - b_j`);
* (ii) for some `eps_j in {+1, -1}` and `t_j in Z/4`: for every `s in S_j` and every prime power `q^a`
  exactly dividing `D_s(pi_j)`, the monomial `e_j` collides at `q^a` with unit `eps_j s + 2 t_j`.

Then for each `j`, `phi_j = arg xi_j (mod pi)` for a Gaussian integer `xi_j` of norm `p_j`, and the member is an
**actual cluster**: its points `R i^(-k_x) exp(i <2a_x - e, phi>)` are `M` distinct Gaussian integers of norm
`N = prod p_j^(e_j)` in an arc of length `<= C sqrt R`.

(ii) holds with `(eps_j, t_j) = (1, t)` for labels coherent with `i^t pi_j` at those prime powers, and with
`(-1, t)` for `i^t conj(pi_j)`: `(i^t pi)/conj(i^t pi) = i^(2t) pi/conj(pi)`, and `pi/conj(pi) = i^s mod q^a` iff
`q^a | D_s(pi)`.

*Proof.* Fix `j`, `s in S_j`, and put `u = eps_j s + 2 t_j` (so `nu_u = nu_s`). Every prime power exactly
dividing `D_s(pi_j)` is odd, prime to `N` by (i), and `<= |D_s(pi_j)| <= a_j + b_j <= X`; by (ii) `e_j` collides
there with unit `u`, so `m_(e_j, u) >= |D_s(pi_j)|` (trivially if `D_s = +-1`). (1.2) for `(e_j, u)` gives
`|sin(phi_j - u pi/4)| >= nu_s |D_s(pi_j)|/sqrt(p_j) = |sin(theta_j - s pi/4)|`. With
`x = eps_j (phi_j - t_j pi/2)` we have `|sin(phi_j - u pi/4)| = |sin(x - s pi/4)|`, so Lemma P0 applies:
`x = theta_j` or `x = pi/2 - theta_j (mod pi)`. Unwinding, `phi_j = arg xi_j (mod pi)` with `xi_j` one of
`i^t pi_j`, `i^(t+1) conj(pi_j)` (`eps = 1`), `i^t conj(pi_j)`, `i^(t-1) pi_j` (`eps = -1`), all of norm `p_j`.
Put `z_x = prod_j xi_j^(a_xj) conj(xi_j)^(e_j - a_xj)`: `|z_x|^2 = N` and
`arg z_x = sum_j (2a_xj - e_j) arg xi_j = <2a_x - e, phi> (mod pi)` (integer coefficients). So the point
`P_x = R i^(-k_x) e^(i <2a_x - e, phi>) = +- i^(-k_x) z_x` is a Gaussian integer of norm `N`, and by the
realisation equation `P_x = R e^(i(theta_0 + delta_x))` with `delta_x in [0, Delta]`, `Delta = C e^(-W/4) = C/sqrt R`.
If `P_x = P_y` then `(z_x) = (z_y)` as ideals, and since `xi_j`, `conj xi_j` are non-associate primes this forces
`a_x = a_y`. QED

**Corollary P1 (P3-G in regime (1.4); proved).** If `X >= max_j (a_j + b_j)` (in regime (1.4): `A > h_0/2`),
every member of P3-G`_X` (genuine residues `r_j(q) = varpi_j mod q^(A_q)`, `varpi_j` any Gaussian prime of
norm `p_j`) whose prime set satisfies (i) is an actual cluster. (i) can fail only if some `p_l < sqrt(2 p_j)`
divides `a_j (a_j^2 - b_j^2)`. *Exact* (`pin_check.py`): for 20 random sets of `r = 50-200` split primes in each of
`[10^2, 10^3]`, ..., `[10^5, 10^6]`, no column fails (i); with primes spread over `[10^2, 10^6]`, 0.1 of 200
columns on average.

This is the norm-1 case of `round7/p3G.md` Thm 3.1 (collapse when the pinning closure contains `2 e_j`), made
automatic by `X > sqrt(2 p_j)`; at `X = 2M` (the p3G setting) the norm-1 relations are rarely tight, which is why
p3G needed long closures and why selected prime sets could defeat them there.

**Corollary P2 (structure of fakes; proved).** In the setting of Prop. P without (i)-(ii), let `U` be the set of
columns for which (i) or (ii) fails (for every choice of `eps_j, t_j`). Then for every pattern
`tau in prod_(j in U) [0, e_j]`, the rows with `a_x|_U = tau` give (after a sign per point) distinct Gaussian
integers `prod_(j not in U) xi_j^(a_xj) conj(xi_j)^(e_j - a_xj)` of norm `N_(U^c) = prod_(j not in U) p_j^(e_j)` whose
arguments lie in an interval of length `Delta = C e^(-W/4) <= C e^(-W_(U^c)/4)`: an actual cluster on the circle
`x^2 + y^2 = N_(U^c)` with constant `C N_U^(-1/4) <= C`. In particular a member that is not an actual cluster has
`U != {}`, and it is a union of at most `prod_(j in U)(e_j + 1)` actual clusters of a sub-circle, glued by the
unpinned columns.

*Proof.* As in Prop. P for the columns outside `U`; the `U`-columns contribute the same constant
`sum_(j in U)(2 tau_j - e_j) phi_j` to every row of the class. QED

**Consequence for C1.** A C1 fake must, for its unpinned columns `U`, carry labels that are **not** genuine
(per unit, up to the symmetries (ii)) at some prime power dividing `a_j (a_j^2 - b_j^2)`, i.e. at a modulus
`q <= sqrt(2 p_j) <= M^(h_0/2)`, and `|U|` must be large enough that the classes of Cor. P2 are not themselves
the whole cluster. Every deformation of the genuine system that keeps the labels genuine below `M^(h_0/2)` (in
particular every flexible top of `L4K.md` §4.4 over a genuine core of that size) is excluded. Conversely,
**rigidity of the norm-1 genuine collisions at the pinning moduli** (a per-unit, norm-1 statement, much narrower
than Rig(K)) would show that P3-all\* has no non-actual members in regime (1.4), with no label-angle step needed
(see §5).

### 3.2 Coherent designs: deep classes cannot be matched, and exact coherence collapses

**Lemma L (deep-class parity matching forces the genuine square class; proved, given Linnik).** Let `p = 1 mod 4`
be prime, `zeta` a Gaussian integer, `n = p Norm(zeta)`, `E >= 2`, and `X >= C_L (2^(E+1) n)^L` (Linnik's
constants as in `round9/L4.md` §2). If `(Norm zeta / q) = (p/q)` for every prime `q <= X` with
`q = -1 mod 2^(E+1)` and `q` not dividing `n`, then `n` is a square or twice a square.

Every such `q` has `k(q) = v_2(q + 1) >= E + 1`: it is a deep modulus, inert, never dividing `N`. So a coherent
column (any `zeta`, any odd units `u_q`) whose parity is right on the deep classes `k(q) > E` up to `X` has
`p_j Norm(zeta)` a square or twice a square, i.e. it is of genuine type.

*Proof.* This is `round9/L4.md` Lemma 2.1 with `8` replaced by `2^(E+1)` and the class `7 mod 8` by
`-1 mod 2^(E+1)`. Write `n = 2^e m s^2` with `m` odd squarefree, and suppose `m > 1`. For `q = -1 mod 2^(E+1)`
(so `q = 7 mod 8`, `(2/q) = 1`, `(q-1)/2` odd) reciprocity gives `(n/q) = (m/q) = sigma_m (q/m)`,
`sigma_m = prod_(l | m) (-1)^((l-1)/2)`. The Jacobi symbol `(./m)` is non-principal, so some `b` prime to `m` has
`(b/m) = -sigma_m`. With `k = lcm(2^(E+1), rad(n)) <= 2^(E+1) n`, choose `a mod k` with `a = -1 mod 2^(E+1)`,
`a = b mod m`, `a = 1 mod l` for the odd primes `l | s`, `l` not dividing `m`; then `gcd(a, k) = 1` and every prime
`q = a mod k` has `q` not dividing `n` and `(n/q) = -1`. Linnik gives one with `q <= C_L k^L <= X`, a
contradiction. So `m = 1`. QED

In regime (1.4) with `Norm zeta <= M^h` the hypothesis holds for every depth with
`2^(E+1) <= X^(1/L) M^(-(h_0 + h))/C_L^(1/L)`, i.e. for all classes up to depth `(A/L - h_0 - h) log_2 M - O(1)`. This
is the unconditional version, in regime (1.4), of the GRH statement of `round7/p3Dconstruct.md` §3.2 ("any other
matching of a whole 2-adic class is impossible").

**Definition.** A column `j` is *exactly coherent* with a primitive Gaussian integer `xi_j` (not divisible by a
rational prime) on a set `S` of odd primes if `r_j(q)/conj(r_j(q)) = xi_j/conj(xi_j) mod q^(A_q)` for every
`q in S` with `q` not dividing `N Norm(xi_j)`. (A rational factor of `zeta` does not change `zeta/conj zeta`; units
`u_q` and torsion shifts, allowed in `round9/L4.md` (2.1), are **not** allowed here.)

**Proposition G (exact coherence plus the norm-1 height force the genuine prime; proved, given Linnik).** Let
column `j` be exactly coherent with `xi_j` on all odd primes `q <= X`, with
`X >= max(C_L (8 p_j Norm xi_j)^L, sqrt 2 |xi_j|)`. Then:
(a) `xi_j = u varpi^(2w+1) sigma^2` with `u` a unit, `varpi in {pi_j, conj(pi_j)}`, `w >= 0`, and `sigma` primitive of
odd norm prime to `p_j`;
(b) if the per-unit height condition holds for `e_j` at `s = 1, 3`, then
`|xi_j|^2 <= 2 p_j g_j^2`, where `g_j` is the largest `N`-part (product of the prime powers `p_l^c`) of
`D_1(xi_j)`, `D_3(xi_j)`. In particular, if no `p_l` divides `D_1(xi_j) D_3(xi_j)`, then `xi_j = u varpi`:
the column is exactly coherent with a genuine prime.

*Proof.* (a) Exact coherence gives the parity of `log(xi_j/conj xi_j)`, i.e. `[(Norm xi_j/q) = -1]`, which must be
`nu_j(q)`; so `(p_j Norm xi_j / q) = 1` for every prime `q <= X` not dividing `N Norm xi_j`. Theorem A of
`round9/L4.md` (classes `7 mod 8`, which never divide `N`) gives `n = p_j Norm xi_j in {t^2, 2t^2}`. The class
`q = 3 mod 8` (never dividing `N`, and the least such prime not dividing `n` is `<= C_L (8n)^L <= X`) has
`(2/q) = -1`, which excludes `2t^2`. So `n` is a square. A primitive `xi` is `(1+i)^delta prod varpi_k^(f_k)` with
`delta in {0, 1}` and, for each split prime, only one of `varpi_k`, `conj varpi_k` occurring; `v_2(n) = delta` even
forces `delta = 0`, `v_(p_j)(n)` even forces an odd power of a prime above `p_j`, and every other exponent is even.
(b) `Norm xi_j` is odd, so `Re xi_j`, `Im xi_j` have opposite parities and `D_1, D_3 = Re -+ Im` are odd. A prime
`q | D_s(xi_j)` does not divide `Norm xi_j`: if `varpi_k | xi_j` (and `conj varpi_k` does not, by primitivity) and
`varpi_k | xi_j - i^s conj(xi_j)`, then `varpi_k | conj(xi_j)`, i.e. `conj(varpi_k) | xi_j`, a contradiction.
`|D_s(xi_j)| <= sqrt 2 |xi_j| <= X`. By exact coherence, `e_j` collides with unit `s` at `q^a` iff
`xi/conj xi = i^s mod q^a` iff `q^a | D_s(xi)` (`L4K.md` §1); so `m_(e_j, s) >= |D_s(xi_j)|/g_j` for `s = 1, 3`, and
the budget gives `|D_s(xi_j)| <= sqrt(2 p_j) g_j`. Since `D_1^2 + D_3^2 = 2 |xi_j|^2`, `|xi_j|^2 <= 2 p_j g_j^2`.
If `g_j = 1`: `|xi|^2 = p_j^(2w+1) Norm(sigma)^2 <= 2 p_j` forces `w = 0` and `Norm sigma = 1`. QED

**Theorem C (exactly coherent members in regime (1.4) are actual clusters; proved, given Linnik).** Fix `h`.
Let a member of P3-all\*`_X` with all `p_j <= M^(h_0)` have every column exactly coherent (on all odd `q <= X`)
with a primitive `xi_j`, `Norm xi_j <= M^h`. Suppose `X >= C_L (8 M^(h_0 + h))^L` (e.g. `X = 3M^A`,
`A >= L(h_0 + h) + 2`, `M` large), and the cofactor conditions: no `p_l` divides `D_1(xi_j) D_3(xi_j)`, and
condition (i) of Prop. P holds for every `j`. Then the member is an actual cluster.

*Proof.* Prop. G gives `xi_j = i^(t_j) varpi_j` with `varpi_j in {pi_j, conj pi_j}`; exact coherence with it at
the prime powers dividing `D_s(pi_j)` (prime to `N` by (i), and `<= X`) is condition (ii) of Prop. P. QED

So in regime (1.4) the only known mechanism of cross-modulus coordination, exact coherence with Gaussian
integers of polynomial height, produces no fakes at all (generic cofactors). Round 9 had reached "coherent
designs are of genuine type, then L7". Theorem C shows that on this branch the class is exact: exactly coherent
members are actual clusters, so there is no label-angle step left, and what remains for growth on this branch is
a lower bound for genuine (actual) data, i.e. L7/(I3) itself (§5 item 3). What it does **not**
cover: coherence up to units `u_q` or torsion shifts (system-level coherence, `L4K.md` §1), coherence on a
proper subset of moduli, and twisted or approximate coherence. Those are exactly where a C1 design would have to
live.

### 3.3 Coherence on the Paley profile cannot certify pairs and columns: a linear bound

Setting of `round7/p3Dconstruct.md` Thm 2: the profile of `sharp.md` §4 (`[1 | H']` normalised Paley-I of order
`M = q_0 + 1`, `M - 1` labels with `b` physical columns each, canonical twins `f(x)`, `s(x)`, `r = b(M-1)`), Gaussian
primes `gamma_1, ..., gamma_D` of odd norm (pairwise non-associate, non-conjugate), `G in Z^(D x r)`,
`zeta_j = prod_t gamma_t^(G_tj)`, and labels (3.1) of p3Dconstruct (root repair, `f = +-i`) at every prime
`q = 3 mod 8`, `q <= X`. Column heights `h_j = sum_t |G_tj|`, `C = max_j h_j`; pair heights
`H_xy = sum_t |<g_t, a_x - a_y>|`, `Hbar` their average.

**Theorem R (proved).** If the height condition (1.2) (per unit) holds for every monomial at the moduli
`q = 3 mod 8`, and

```text
theta(X; 8, 3)  >  2 (M - 1) b C W + log 2,                                                        (3.1)
```

then `C + 2 Hbar > (M - 1)/8`. In regime (1.4) (`W <= c_0 M log M`, `X = 3M^A`), (3.1) holds for every `C <= M`
once `A >= 4` and `M >= M_0(b, c_0)`; so every G-coherent design meeting the height condition has
`C + 2 Hbar > (M-1)/8`.

This improves p3Dconstruct's Cor. 2 (`C + 2 Hbar >= M^(1/2 - o(1))`) to a linear bound, using that `X` is large.

*Proof.* Let `psi_t in Z^(M-1)`, `psi_(t,a) = sum_i G_(t,(a,i))`, and `T* = {t : psi_t != 0}`. Steps 1-4 of
p3Dconstruct Thm 2 give, for `t in T*`, `sum_(x<y) |<g_t, a_x - a_y>| >= (M/4)(M/4 - ||phi_t||_1)` (row image
`A g_t = c_t 1 + (1/2) H' psi_t + phi_t`, Hadamard `l^1` uncertainty `||(1/2) H' psi - m 1||_1 >= M/4`, median
inequality), with `sum_t ||phi_t||_1 <= M C`. Summing over `T*`:
`Hbar M(M-1)/2 >= (M/4)(|T*| M/4 - M C)`, i.e. `C + 2 Hbar >= |T*|/4` (as `M - 1 < M`).
Suppose `|T*| <= (M-1)/2`. The `|T*|` linear forms `beta -> <beta, psi_t>` on `Z^(M-1)` have coefficients
`|psi_(t,a)| <= b C`; Siegel's lemma gives `beta != 0` with all of them zero and
`||beta||_inf <= ((M-1) b C)^(|T*|/(M - 1 - |T*|)) <= (M-1) b C` (any `beta != 0` if `T*` is empty). Put
`v = sum_a beta_a s_a` (`s_a` the sum of the `b` columns of label `a`); `v != 0`, `(G v)_t = <beta, psi_t> = 0` for
every `t`, so `Z_v = 1` and by p3Dconstruct Prop. 3.3(4) `v` collides at every `q^(A_q)`, `q = 3 mod 8`, with unit
`i^(+-E_v(q))`. Summing over the four units, `sum_s log m_(v,s) >= theta(X; 8, 3)`, while the per-unit budgets sum
to `2 w(v) + log 2 <= 2 ||beta||_inf W + log 2 <= 2 (M-1) b C W + log 2`, contradicting (3.1). Hence
`|T*| > (M-1)/2` and `C + 2 Hbar > (M-1)/8`. Regime (1.4): `theta(X; 8, 3) >= X/5 >= (3/5) M^A` for large `X`, and
`2(M-1) b C W <= 2 b c_0 M^3 log M` when `C <= M`. QED

Reading: columns certified at norm 1 need `C = O(h_0)` (their certificates `4 h_j log max|gamma_t|` must fit
the budget `(h_0/2) log M`), so then `Hbar > M/16 - O(1)`: the pair integers have `log`-size `Theta(M)` and the
pair demands (which must be `O(log M)` on average, `round8/README.md` §3) are not certified by size. The pair
side would have to be uncertified-random, which in regime (1.4) costs `~ 4 log X` per pair on average for
random-like collisions (heuristic H2 of round 7 at `X = M^A`) against a mean margin `<= (c_0/2) log M`.

### 3.4 The shallow classes are harmless for monomials, uniformly in `A`

**Theorem S (proved).** Let `E_0 >= 2`, let `p_1, ..., p_r` be distinct primes `= 1 mod 4`, and let
`l_1 < ... < l_r` be the `r` smallest primes `= 1 mod 2^(E_0 + 1)` not among the `p_j`, `zeta_j` a Gaussian prime
of norm `l_j`. For **every** `X` there are factored residues `r_j(q) mod q^(A_q)` (`Norm r_j(q) = p_j`, levels
compatible) at every odd prime `q <= X` with `q` not dividing `N` and `k(q) <= E_0`, such that for every
`v in Z^r \ {0}` and every unit `s`

```text
sum_(q <= X, q not | N, k(q) <= E_0)  a_(q,s)(v) log q   <=   sum_(k=2)^(E_0) log |Im(gamma_v^(2^k))|   <=   (2^(E_0) - 2) sum_j |v_j| log l_j,
```

where `gamma_v = prod_(v_j > 0) zeta_j^(v_j) prod_(v_j < 0) conj(zeta_j)^(|v_j|)`. By the prime number theorem in
progressions `l_r <= (1 + o(1)) 2^(E_0) r log r` (`E_0` fixed, `r -> infinity`).

*Proof.* Every `l_j` has `k(l_j) = v_2(l_j - 1) >= E_0 + 1`, so no shallow `q` equals an `l_j`, and every `zeta_j`
is a unit modulo every shallow `q`. Define `r_j(q)` by p3Dconstruct (3.1) at the top level `Q = q^(A_q)`, with
`f_Q` of order `2^(k(q))` (`f_Q = +-i` if `k(q) = 2`), and lower levels by reduction; Lemma 3.2 there shows these are
valid factored residues (parity `nu_j(q)` at every level, realisable with norm `p_j`, compatible), and reduction
preserves the order of `f`. By Prop. 3.3(1) there, a collision of `v` at `q^a` (any unit) implies
`q^a | Im(gamma_v^(2^(k(q))))`. `gamma_v/conj(gamma_v)` is not a root of unity (distinct non-associate,
non-conjugate primes of odd norm; `v != 0`), so `Im(gamma_v^(2^k)) != 0`, and distinct primes of the class `k`
contribute coprime factors: `sum_(k(q) = k) a_(q,s)(v) log q <= log |Im(gamma_v^(2^k))| <= 2^k log|gamma_v|`, with
`log |gamma_v| = (1/2) sum_j |v_j| log l_j`. Sum over `k = 2, ..., E_0`. QED

**Corollary S1 (proved).** In the setting of Theorem S, if `sqrt(p_j) >= 40 r l_r^(2^(E_0) - 2)` for all `j`, then the
shallow part of the union-bound completion sum is small:
`sum_(v != 0) sum_s m^(shallow)_(v,s) e^(-w(v)/2) <= 4 (e^(1/10) - 1) < 1/2`. In regime (1.4) with `r = Theta(M)`
this holds when every `p_j >= M^(2^(E_0 + 1) - 2 + eps)` (compatible with `p_j <= M^(h_0)` iff `h_0 > 2^(E_0 + 1) - 2`), for every `A`.

*Proof.* The sum is at most `4 (prod_j (1 + 2 sum_(t >= 1) x_j^t) - 1)` with `x_j = l_j^(2^(E_0) - 2)/sqrt(p_j) <= 1/(40 r)`,
which is `<= 4 (exp(4 sum_j x_j) - 1) <= 4 (e^(1/10) - 1)`. QED

*Exact* (`shallow_cert.py`, `shallow_cert_out.txt`): `r = 6`, `E_0 = 3`, `X = 3000`, `zeta`-norms
`17, 97, 113, 193, 241, 257`, primes `p_j ~ 10^4`; and `r = 8`, `E_0 = 4`, `X = 2000`. At all 324 (resp. 266) shallow
moduli and all 347 (286) levels, constructed at the top level and reduced: the parity is the forced one at every
level, an explicit `r_j(q)` with `Norm r_j(q) = p_j` realises the level-1 label for every `q <= 150` (`100`), every
one of the 1211 (487) collisions of a monomial with `||v||_1 <= 3` (`<= 2`) and a unit is certified by
`q^a | Im(gamma_v^(2^k))`, and the total shallow demand is at most `0.28` (`0.12`) of the bound of Theorem S.

**What Theorem S means.** For the monomials **outside the character lattice** (whose budgets grow linearly in
`||v||_1` and whose angles are Haar-uniform on the realisation torus, so that the union-bound completion applies),
the parity obstruction lives only on the deep classes: shallow moduli `k(q) <= E_0` can be served, uniformly in
`X` (hence in `A`), at cost `2^(E_0) log(2^(E_0) r log r)` per unit of `||v||_1`, which fits the budgets when
`h_0 > 2^(E_0+1) - 2`. The bound is useless for pairs and characters (`v in L`, `||v||_1 ~ r/2`, margins
`O(log M)`), and by Theorem R no coherent design on the Paley profile can serve them. So this is a localisation
for construction purposes, not a cap and not a statement about proofs: a C1 design may take Theorem S on the
shallow classes for the monomials outside `L`; it must supply the deep classes `k(q) > E_0` (DeepBal below) and the
pair/character side by other means.

### 3.5 Anatomy of a C1 fake: necessary conditions and the residual problems

Collecting §2-§3, a member of P3-all\*`_(3M^A)` in regime (1.4) that is not an actual cluster must have labels
that are, simultaneously:

| | necessary property of a C1 fake | why | status |
|---|---|---|---|
| N1 | non-genuine (per unit, up to the symmetries of Prop. P(ii)) at a pinning modulus `q <= sqrt(2 p_j)` of every column in a set `U` with `U != {}` and every class of Cor. P2 a proper part of the cluster | Prop. P, Cor. P2 | proved |
| N2 | not exactly coherent on all moduli with Gaussian integers of norm `<= M^h` (generic cofactors) | Thm C | proved (Linnik) |
| N3 | on the deep classes `k(q) > E`, `2^(E+1) <= X^(1/L) M^(-(h_0+h))`, not exactly coherent with the forced parity for any `zeta` of norm `<= M^h` unless `p_j Norm zeta` is a square or twice a square | Lemma L | proved (Linnik) |
| N4 | not independent over moduli on any set of moduli in `(M^(h_0 k/2), X]` carrying collision probability `>= c/q` for norm-`k` monomials (single-hit form); not independent on any set of reciprocal mass `>= (log M)^(-C)` in `[M^(1-eta), X]` | Thms B1, B2; `p3Dconstruct.md` Thm 1 | proved (probabilistic) |
| N5 | not a fixed integer vector, or a family of permuted small integer vectors, on a positive-mass set of large moduli | Thm B3; `sharp.md` Prop. 7.1 | proved |
| N6 | on the Paley profile: if coherent at `q = 3 mod 8`, then `C + 2 Hbar > (M-1)/8` (pairs uncertified) | Thm R | proved |

Every mechanism of cross-modulus coordination that rounds 5-14 know (independence with repairs, integer/Sidon
vectors, Gaussian coherence with any repair) is excluded by N2-N6, and the one structure that does satisfy all
budgets everywhere (genuine labels) is excluded by N1. What remains open is precisely:

> **DeepBal(E_0).** In regime (1.4), are there parity-respecting labels on the deep moduli `k(q) > E_0`,
> `q <= X = 3M^A`, with `log m^(deep)_(v,s) <= beta ||v||_1 log M` for all `v` and `s` (`beta` independent of `A`),
> that are not exactly coherent? (With Theorem S on the shallow classes this would give the monomial part of C1 at
> `h_0 ~ 2^(E_0+1) + 2 beta`.)
>
> **PairCompat.** A profile with `W = O(M log M)` and residue data whose pair demands are `O(log M)` on average
> uniformly in `A` (necessarily non-coherent point classes at the moduli `q >~ M^(c_0)`, by the Plotkin bound on
> margins) whose column labels still satisfy DeepBal and Theorem S's bound. Theorem R excludes coherent column
> labels on the Paley profile; with twins, the twin labels inherit the point classes.

*Heuristic assessment (labelled).* DeepBal asks for a family of lattices `Lambda_q` (index `~ q`, with the fixed
parity sublattice `{<v, nu(q)> even}` at `q = +-1 mod 8`) whose intersections over **every** set `S` of deep moduli
have no nonzero vector of `l^1`-norm below `log(prod_S q)/(beta log M)`. Random lattices fail this by the Poisson
tail (the union over `S` is too large); Gaussian coherence satisfies it with room (one certificate controls all
`S`) but is excluded by parity (N2, N3); evidence below (`deep_demand.py`) shows root-repaired coherence and
uniform labels already violate norm-1 and norm-2 per-unit budgets on the deep classes alone at `X = 2*10^5`, while
genuine labels sit exactly at ratio 1.00. We know of no third mechanism. This is evidence (no more) that C1 is
impossible and that the regime (1.4) is where the barrier stops; it is not a lower bound.

## 4. Finite verification

Scripts in `c1construct_checks/` (run from that directory with `python3 <script>`; outputs `*_out.txt` next to
them). `lib14.py` holds exact Gaussian arithmetic modulo `q^a`; `lib10.py` is a copy of `round10/L4K_checks/lib.py`
(used for monomial enumeration and the Legendre symbol). Random choices are seeded.

| script | output | what it checks | label |
|---|---|---|---|
| `pin_check.py` | `pin_check_out.txt` | (a) `sin^2(theta - s pi/4) = nu_s^2 D_s(pi)^2/p` as rationals, all 4783 split `p <= 10^5`; (b) Lemma P0 on a 40000-point grid, 1125 primes; (c) frequency of failures of Prop. P(i) on random prime sets | (a) exact; (b), (c) evidence |
| `break_random.py` | `break_random_out.txt` | Thm B1: exact survival products `prod (1 - 4/m_q) prod (1 - 8/m_q)` from actual Legendre symbols, 200 columns, four `(M, h_0, A)`; seeded simulation with explicit witnesses `(p_j, q, ell_j(q))`, `ell_j(q)` in the 4-torsion | exact products; witnesses exact; proportions evidence |
| `break_bh.py` | `break_bh_out.txt` | Thm B3: the `L4K.md` (2.1) vectors with random `sigma_q` at every `q in (q_T, X]`; exact count of norm-4 monomials with `<v, a> = 0` per modulus against `Z_r`; every monomial of `F` over budget; explicit witness `v` with its colliding moduli | exact (per instance) |
| `shallow_cert.py` | `shallow_cert_out.txt` | Thm S: parity at every level, explicit level-1 realisations, certificate `q^a | Im(gamma_v^(2^k))` for every collision (1211 and 487), total shallow demand `<= 0.28` / `0.12` of the bound | exact |
| `coh_check.py` | `coh_check_out.txt` | Lemma L: least prime `q = -1 mod 2^(E+1)` with `(n/q) = -1` for random `n` (evidence on the size of Linnik's prime); Prop. G(b): 12600 twisted coherent columns `xi = pi sigma^2` all violate the norm-1 per-unit height at `s = 1` or `3` (min ratio `3.545 = 5/sqrt 2`), `D_1^2 + D_3^2 = 2|xi|^2` | Prop. G(b) instances exact; Linnik sizes evidence |
| `deep_demand.py` | `deep_demand_out.txt` | per-unit deep-class (`k(q) > E_0`) level-1 demand of norm-1, 2 monomials, `X = 2*10^5`, `r = 10`, `p_j ~ 2*10^4`, for coherent+root repair, uniform labels and genuine labels, as a ratio to the per-unit budget | evidence |

`deep_demand_out.txt` (evidence): with `E_0 = 2` the deep classes alone push 5 of 40 norm-1 (and 29 of 400 norm-2)
monomial-units of the coherent+repair design over their per-unit budget (max ratio 1.70, 2.14), uniform labels 17
and 79 (max 3.25, 2.75), genuine labels none (max ratio exactly 1.00, as Prop. 6.1 requires). With `E_0 = 3`:
3, 9 (coherent+repair), 6, 19 (uniform), 0 (genuine).


## 5. Implications for the barrier logic

1. **The implication "lower bound in regime (1.4) => growth" is correct, with the quantifiers of (T').** By
   `round8/stability.md` Cor. 3.6 (refereed), an actual cluster with `M >= eps log R/loglog R` has
   `W = 2 log R <= (2/eps + o(1)) M log M`, all primes `<= M^(4/eps - 3 + o(1))` and `r <= (2/eps + o(1)) M`; with
   its genuine residues it is a member of P3-all\*`_X` for **every** `X` (`sharp.md` Prop. 1.2). So it lies in
   regime (1.4) with `c_0 = 2/eps + o(1)`, `h_0 = 2c_0 - 3`. Hence any theorem of the form
   *"for every `h_0` there are `A` and `F > c_0 = (h_0 + 3)/2` such that every member of P3-all\*`_(3M^A)` with all
   `p_j <= M^(h_0 + o(1))` has `W >= F M log M` for `M >= M_0`"* gives (T') and so `M = o(log R/loglog R)`.
   A bound `W >= g(M) M log M` with `g -> infinity` at one fixed `A` (e.g. `g ~ loglog M`) gives growth with an
   explicit rate. Two cautions: the bound must hold for **all** members of the restricted class (not just for
   members built by some method), and `A` may depend on `h_0` but not on `M`.
2. **Nothing here proves such a bound or a cap.** Theorems B1-B3, R and Lemma L are barriers to construction
   methods; Prop. P, Prop. G and Theorem C say that members with special labels (genuine at the pinning moduli,
   or exactly coherent) are actual clusters, which bounds nothing; Theorem S is a partial design (monomials, shallow
   moduli), not a member. The bracket `3 M log M - O(M) <= W_min^all <= O(M log^2 M/loglog M)` and the bound
   `M <= (2/3 + eps) log R/loglog R` are unchanged. P3-all\* stays open and uncapped in regime (1.4); no fake of any
   kind is known there.
3. **What the round adds to the map.** In regime (1.4):
   * every coordination mechanism known from rounds 5-13 is excluded as a C1 design (N1-N6 of §3.5);
   * for constructions, the parity obstruction for the monomials outside `L` is confined to the deep classes
     `k(q) > E_0`, uniformly in `A` (Theorem S), and the pair obstruction for coherent designs on the Paley
     profile is linear in `M` (Theorem R);
   * **the dichotomy is sharper than Rig(K) suggests.** Either (a) a C1 fake exists — it needs DeepBal and
     PairCompat by a mechanism that is neither coherent nor independent — and then P3-all\* is capped where growth
     needs it; or (b) P3-all\* has no non-actual members in regime (1.4). By Prop. P and Cor. P2, (b) follows (up to
     the columns failing condition (i), which are rare, Cor. P1) from a narrow, norm-1, per-unit rigidity statement
     (**Rig_pin**: in every good system the collisions of `e_j` at the prime powers dividing `a_j`, `a_j +- b_j` are
     the genuine ones, up to the symmetries of Prop. P(ii), for every column satisfying (i)), with no label-angle
     step: the angles are then pinned exactly, and a member is an actual cluster, or (if the set `U` of columns
     failing (i) is non-empty) a union of at most `prod_(j in U)(e_j + 1)` actual clusters of a sub-circle.
   * **The label-angle step is settled at norm 1 in regime (1.4).** `round10/referee_L4K.md` §7 showed that every
     fake, given the genuine labels of its own primes, satisfies all hypotheses and conclusions of L4(K), L4^forced
     and L4\*, so route 3's fake-separating power had to come from the label-angle step (L7 / slabs). Prop. P is
     such a fake-violated input: a fake whose labels are genuine (per unit) at the pinning moduli of all its
     columns violates (1.2) for some single column `e_j`. So route 3 in regime (1.4) is now "Rig_pin + Prop. P".
     Caveat: Prop. P needs **per-unit** genuine collisions; the conclusions of L4\*-type statements are
     system-level coherence (one odd unit `u_q` per modulus), which scrambles the units and is not enough as it
     stands (open point 3).
   * **What case (b) does and does not give.** In case (b) the class is exact in regime (1.4): its members are the
     actual near-extremal clusters, so W_min^all there equals the actual minimum and (T') for P3-all\* is literally
     the growth statement for actual circles. Rigidity therefore answers the **barrier** question (P3-all\* is not
     a relaxation there, so no fake barrier exists for it) but proves no growth by itself; a growth proof through
     the class still needs a lower bound for genuine data on near-extremal profiles, i.e. L7/(I3) of
     `round9/README.md` §4 (Poisson tails of genuine collision data), which is a statement about actual circles.
     This should be kept in mind when weighing Rig(K) against C1.


## 6. Open points

1. **DeepBal(E_0)** and **PairCompat** (§3.5): the two residual design problems of C1. A solution of DeepBal by a
   non-coherent mechanism would be new even without pairs; a proof that DeepBal fails (every family of lattices with
   the forced parity on the deep classes has a short vector in many intersections) would be the first rigidity
   mechanism, but by §5 item 3 it would bound only non-actual members.
2. **Rig_pin** (§5 item 3): decide whether good systems (all norms, all moduli `<= X`) must reproduce the genuine
   norm-1 collisions at the pinning moduli. A finite exact test is the one-modulus deformation of `L4K_checks/deform.py`
   restricted to the pinning moduli, with **all** norms imposed (feasible for `r = 4, 5`: the vectors with genuine
   slack `< log q` are finitely many once `||v||_1 > 2(psi(X) + log q)/log p_min`). Not run here. Note that a
   single-modulus deformation does not unpin a column in regime (1.4): the norm-2 relations `e_j +- e_k` (whose
   `D_s(pi_j pi_k^(+-))` are `<= sqrt(2 p_j p_k) <= X`) pin `phi_j +- phi_k` in the same way, so unpinning needs
   non-genuine labels at the prime factors of `D_s(Z_v)` for many short `v`.
3. **Theorem C beyond exact coherence**: system-level coherence (one odd unit `u_q` and torsion shifts per modulus,
   `L4K.md` §1), coherence on a proper subset of moduli, twisted coherence. Lemma P0 needs per-unit information;
   with scrambled units the one-sided arcs move and the solution set of Lemma P0 can be an interval.
4. **Theorem R for other profiles.** The proof uses the Paley `l^1` uncertainty. Near-extremal profiles are
   Plotkin-balanced (columns of weight `~ M/2`, pair distances `~ W/2`); an uncertainty principle for all such
   profiles would make "coherent designs cannot certify pairs" profile-independent.
5. **Coherence over other number fields** (heuristic): labels from elements of a larger field `K ⊃ Q(i)` meet the
   same two walls in regime (1.4) — parity matching on deep classes forces `|Norm_(K/Q)| >= p_j` up to squares
   (Linnik), and the norm-1 height then needs a large `X`-rough part of `Norm_(K/Q)(alpha_j - i^s conj alpha_j)`,
   impossible when `X` exceeds that norm. Not written out.
6. **Referee items.** Linnik's theorem (constants `C_L, L`) and Siegel's lemma are quoted; Mertens' theorem in
   progressions mod 8 is used with an unspecified `O(1/log x)`; the prime number theorem in progressions bounds
   `l_r` in Theorem S. Theorems B1, B2 are probability statements about the stated laws (independence over `q` and
   over columns, as in `sharp.md` Lemma 6.3(2)); B3 is deterministic. Prop. P needs the normalisation of
   `pi_j` (`a_j` odd, `b_j` even), which fixes which `D_s` are odd.
