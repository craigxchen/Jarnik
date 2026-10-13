# Referee report: `round7/p3Dconstruct.md` (P3-all\* structured fibres)

Adversarial referee, round 7. Scripts are in `round7/referee_p3Dconstruct_checks/`. They are my own
code and import nothing from the author's `p3Dconstruct_checks/`. Each runs with `python3 <script>`
from that folder.

## 0. Verdict

**The result survives as a partial barrier for construction methods, but its scope has to be
narrowed in five places.**

The formal statements check out on re-derivation:
* Theorem 1, once its hypothesis (Ind) is repaired (see M3);
* Lemma 3.2 and Proposition 3.3;
* Theorem 2 with Corollary 2;
* identity (1.4).

My own exact checks agree. The author's negative headline is also right: no `O(M log M)` member of
P3-all\* was found, and the rate stays `O(M log^2 M/loglog M)` (`sharp` Thm 6.4).

The barrier *claims* built on these statements go further than the proofs. In order of importance:

1. **(M1, major: a gap in the main answer.)** The answer to the task's suggestion is not proved.
   The suggestion is fixed auxiliary Gaussian integers on the non-twin columns, with the twin columns
   designed. The write-up says that if the twins are random Theorem 1 applies, and if they are
   coherent Theorem 2 applies. That dichotomy does not cover the suggested design. There the twin
   labels `t_x = eps_x(kappa_x - P_x)` are deterministic, so they are not random in the sense of (Ind).
   They are either not coherent (structured point classes such as `D_A`), or coherent only with
   column height `>= M^(1/2 - o(1))` (by Corollary 2 itself, when the pairs are certified), where
   Theorem 2 adds nothing. Neither theorem applies. The case is open, and it is supported only heuristically.
2. **(M2, major: overstated summary of Theorem 1.)**
   * "Every completion has `W >= ...`" is false as written. Corollary 1 is a probability statement:
     a random design with smaller `W` completes with probability `< exp(-M^(eta/4))`.
   * "Independent randomness on any prime set ... is fatal" is false as written. Take labels that are
     independent and uniform on the parity class *minus the 4-torsion*. They have no single-column
     collisions at all, so Theorem 1's mechanism disappears.
   * The summary replaces the per-column quantity `lambda_j` by the reciprocal mass of `Q*`. For the
     classes `k(q) >= 3` (for example the advertised `Q* = {q = +-1 mod 2^K}`) these differ through
     Legendre symbols that nobody controls.

   Two-column monomials repair the second and third points (sketch in §3).
3. **(M3, major in statement, trivial to repair.)** The hypothesis (Ind) asks for independence
   *conditionally on all data of the columns outside `J`*. That fails in every application the
   write-up lists:
   * the random fibre of `sharp` (Consequence 1, "Prop. 6.5 is the case ...");
   * the twin columns (Consequence 4).

   In the first, the twin data determine `M - 1` linear combinations of the non-twin labels. In the
   second, the non-twin data determine the twin labels outright. The proof only needs
   *unconditional* independence, so the fix is to restate (Ind) that way.
4. **(M4, major: a wrong "exact" claim.)** §3.2, "matching the parities instead costs
   `2^K >~ X/(log X log p_j)`", is labelled exact. It is a pseudosquare-type density heuristic.
   It also omits the choice `Norm zeta_j in p_j Q*^2` (`zeta_j = pi_j` up to rational and square
   factors), which matches every parity at no cost. `sharp` §6.6 lists that choice explicitly. What
   is true, and exact, is weaker and different:
   * matching through `Norm zeta_j in {1, 2} p_j Q*^2` forces `|zeta_j|^2 >= p_j`, and then
     certification by size is lost (exact);
   * any other matching of a whole 2-adic class on `[M^(1-eta), X]` is impossible under GRH, but
     unconditionally it is excluded only polynomially.

   So "parity is a genuine obstruction on `q = +-1 mod 16`" is a statement about one family of
   designs, not a theorem.
5. **(M5, scope.)** Open point 3 and §4.2 say that *any* construction at `O(M log M)` must control
   smooth parts of `binom(M,2)` integers of size `exp(M^(1/2-o(1)))`. Corollary 2 proves this only
   for G-coherent designs on the Paley profile, and only for the *average* pair height.

The minor points are in §6. One of them goes in the author's favour: the surjectivity hypothesis of
Consequence 4, which the author leaves open ("fails mod 2 for Paley; not resolved"), holds on a row
set of size `M/2` for every Paley order I tested (exact, `q0 <= 107`). The consequence itself is
redundant, though (§6).

## 1. What I re-derived and what I checked

### 1.1 Re-derivations (line by line)

* **(1.1), (1.2).** `a_q(v) = max_s a_(q,s)(v) <= sum_s a_(q,s)(v)` gives (1.1). The per-unit
  condition `m_(v,s) <= nu_s^(-1) e^(w(v)/2)` (`sharp` Prop. 6.1 / Def. 1.1 (e-all) with `|sin| <= 1`)
  then gives `log m^gen_v <= 2 w(v) + log 2`. The stated `2 log 2` is weaker, so (1.2) is a valid
  necessary condition for every member. It is aggregated *after* the per-unit bound, so it is not
  the false aggregated form that the `sharp` referee rejected. **Valid for actual circles:** yes,
  through `sharp` Prop. 1.2.
* **(1.3).** `T_Q` is cyclic and the 4-torsion is `{0, m/4, m/2, 3m/4}`. If `k(q) = 2`, `m/4` is odd,
  so the 4-torsion has two elements of each parity. If `k(q) >= 3`, all four are even. Correct.
* **(1.4).** Substituting `a_(y,f(x)) = a_(y,s(x)) + eps_x [y = x]` into `kappa_y = <a_y, ell>` and into
  `<v, ell>` gives exactly the stated formulas. I also checked the converse of
  "`v in L` iff `u = 0` and `sum c' = 0`" (with `L` the zero-sum combinations, as in `sharp` Lemma 6.2):
  if `u = 0`, then `v = sum_x c'_x a_x` coordinate by coordinate (non-twin, `s(x)`, `f(x)`). Correct.
  It is a change of variables, nothing more.
* **Theorem 1.**
  * Step 1: `k = floor((1-eta) log M/((C+1) loglog M)) + 2` gives
    `k(1-eta) log M > 2 T_M + 2 log 2`.
  * Step 2: `Pr[H = k] >= e_k(p) prod (1 - p_q)`. The repeated-index bound gives
    `k! e_k >= lambda^k (1 - binom(k,2) p_max/lambda)`. Then
    `k (C loglog M + log k) <= (1-eta) log M + 2(C+1) loglog M`.
  * Step 3: product over independent columns.
  * Corollary 1 follows. One strictness nit: the corollary assumes probability `>= exp(-M^(eta/4))`
    and the theorem gives `<= exp(-M^(eta/4))`. The bound `(1-x)^n < e^(-nx)` is strict, so the
    contradiction stands.
  * **Correct, under (Ind) as the proof actually uses it** (unconditional independence; see M3).
* **Lemma 3.2.**
  * Item 1: the elements of order `2^k` are `g^(o u)` with `u` odd, so their labels are odd. When
    `k = 2` they are `+-i`.
  * Item 2: `nu(z/conj z) = (Norm z/q)`. I re-derived this for split `q` (`T ≅ {(u, u^(-1))}`) and for
    inert `q` (`z^(1-q)` is a square in `T` iff `z` is a square in `F_(q^2)^*` iff `Norm z` is a square
    in `F_q^*`).
  * Realisation follows from `sharp` Prop. 2.1(4). Correct.
* **Proposition 3.3.**
  * Item 1: raise to the power `2^k`. Since `f^(2^k) = 1` and `i^(2^k) = 1`, a collision gives
    `q^a | 2i Im(gamma^(2^k))`. A rational integer divisible by `q^a` in `Z[i]` is divisible by `q^a`
    in `Z`.
  * Items 2-4: correct.
  * Nit: "Gaussian primes of odd norm" must read "*split* Gaussian primes". An inert prime `q'` has
    `q'/conj q' = 1`, a root of unity, and it violates Definition 3.1.
* **Corollary 3.4.** `sum_(k=2)^K 2^k = 2^(K+1) - 4` and `log|gamma_v| <= ||v||_1 log Y/2`. Correct.
* **Theorem 2**, all eight steps:
  * Step 1, the row-image identity, matches the profile of `sharp` §4: the flip column is
    `col_alpha(x) + eps_x e_x`, with `eps_x = -H(x, alpha(x))`.
  * Step 2, the median bound: `sum_(x<y)|F_x - F_y| >= (M/4)||F - mu 1||_1`.
  * Step 3, the Hadamard uncertainty: `||Hz||_1 >= M||z||_2^2/||z||_1 >= M||z||_1/||z||_0 >= M/4`.
  * Step 4: `|T*| <= 4(C + 2 Hbar) <= s/2`.
  * Step 5, the pigeonhole. It needs `n' > s`, which (4.1) implies, and `|psi_(t,a)| <= bC`.
  * Step 7: `G v = 0` makes `Z_v = 1`, so `v` collides at every `q = 3 mod 8` with the unit
    `i^(+-E_v(q))`.
  * Step 8: `w(v) <= 4sW/(M-1)`.
  * Step 6 (`v` not in `L`) is correct but not needed, since (1.2) holds for every `v != 0`.
  * **Correct.**
* **Corollary 2.** Correct. It can be sharpened: when `C <= K' log M`, both cases of the dichotomy
  give `Hbar >= M^(2/3 - o(1))`. (If (4.1) holds, `Hbar >= M^(1-o(1))`. If it fails,
  `(C+2Hbar)^3 C >= M^2/(2^14 b)`.) The stated `M^(1/2 - o(1))` is true but weaker.

### 1.2 My exact checks

| script | what | result |
|---|---|---|
| `root_ref.py` | Lemma 3.2 and Prop. 3.3 on all 89 odd prime powers `<= 400`, using my own `Z[i]/Q` arithmetic. Checked: `|T_Q| = m`; `f` has odd label and order `2^k`; `f = +-i` iff `k = 2`; `k = 2` iff `q = +-3 mod 8`; `k = 3` iff `q = 7, 9 mod 16`; every odd label has order divisible by `2^k`; `nu(z/conj z) = (Norm z/q)` on random `z`; (3.1) has parity `(p_j/q)`; `r_j` realised by explicit Hilbert 90 (`u = c + t conj c`) and Hensel, with `Norm r_j = p_j mod Q` verified; Prop. 3.3(3) (equivalence at `k = 2`); 60 x 400 random monomials per modulus. | **168562 collisions, all with `Q | Im(gamma_v^(2^k))`; 0 failures.** 13 monomials had `Im = 0`. Each came from the same Gaussian prime being drawn twice (`Z_v = 1`), which is consistent with item 4. (`root_ref_out.txt`) |
| `paley_snf.py`, `paley_snf2.py`, `snf_cross.py` | the surjectivity hypothesis of Consequence 4 for the Paley 0/1 label matrix `P = (1 + H')/2` (`M x q0`) | `rank_F2 P = M - 1` for `q0 = 3 mod 8` and `M/2` for `q0 = 7 mod 8`. For **every** `q0` in {43, 47, 59, 67, 71, 79, 83, 103, 107} there is a row set `R` with `|R| = M/2` and `P_R` **surjective over `Z`** (Smith form trivial; confirmed with sympy for `q0 = 43, 47`). See §6. |
| `torsion_avoid.py` | M2(b): labels uniform on (parity class) minus the 4-torsion | exact enumeration for `m` in {12, ..., 400}: single-column hit probability **0**; same-parity pair hit probability exactly `4/(m-4)` (`k = 2`), and `8/(m-8)` or `8/m` (`k >= 3`). Always `>= p(q,0)`. |

### 1.3 Reruns of the evidence (my own code, other seeds; evidence only)

| script | quantity | mine | author |
|---|---|---|---|
| `tail_ref.py` | single columns, `r = 8M`, `X = 3M`: RAND max/mean vs COH max/mean, in units of `log M` | `M = 64`: 8.96/3.73 vs 6.28/3.22; `M = 256`: 10.72/3.82 vs 6.83/3.12; `M = 1024` (2 seeds): 11.14/3.86 vs 7.33/3.00 and 10.83/3.87 vs 7.24/3.00 | 8.04/3.77, 6.63/3.22; 11.52/3.82, 7.26/3.11; 11.16/3.88, 7.36/2.99 |
| `pairs_ref.py` | pair demand `D1` of the coherent `G = I` design, `b = 6`: max/mean, in units of `log M` | `M = 44`, 5 seeds: max 4.22 to 6.13, mean 0.911 to 0.977. `M = 132`, 2 seeds: max 6.01 and 5.56, mean 0.952 and 0.961 | `M = 44`: 4.18/0.919; `M = 132`: 5.76/0.954 |

The means reproduce. Across seeds the maxima vary by about `2 log M` at `M = 44`. That is as large as
the coherent-versus-control differences in the author's table. So "statistically indistinguishable"
is consistent with the data, but it rests on one seed per size (minor).

## 2. M1: the answer to the task's suggestion is not covered by Theorems 1-2

**The claim** (§0, "Answer to the task's suggestion", and the result summary): "the twin labels are
the solved coordinates `t_x = eps_x(kappa_x - P_x)` of (1.4). If they are random, Theorem 1 applies;
if they are coherent, Theorem 2 applies."

**The suggested design.** The non-twin columns (and the `lambda'_x` slots, if desired) are coherent:
`ell_j = log_g((zeta_j/conj zeta_j) f_q^(eps_j(q)))` with fixed auxiliary `zeta_j` of small norm. The
point classes `kappa` are *designed* to make pair demands `O(log M)`, for example the class design
`D_A` of `sharp` §3. The twin labels are then solved from (1.4).

**Why neither theorem applies.**
1. *Not (Ind).* Every label is a deterministic function of the fixed data, and there is no
   probability space. Theorem 1 is a statement about the *distribution* of a random design. It says
   nothing about a fixed design, however pseudo-random its labels look.
2. *Not G-coherent with small height.* Definition 4.1 requires *every* column, the twin columns
   included, to be coherent at every `q = 3 mod 8`, `q <= X`. In step 1 of the proof the twin
   exponents `G_(t,f(x))` enter through `phi_t`, and in step 5 the vector `v_beta` uses all physical
   columns of a label. Two cases arise.
   * (a) `kappa` is not coherent (e.g. `D_A`, whose classes are residues `x mod p'`). Then `t_x` is not
     the label of a fixed Gaussian number, and Definition 4.1 fails.
   * (b) `kappa_x` is the label of a fixed Gaussian number `K_x` at every `q = 3 mod 8`. Then
     `zeta_(f(x)) = (K_x / Z_(atilde_x))^(eps_x)` up to units. If the pairs are certified (`Hbar` small),
     Corollary 2 itself forces `C >= M^(1/2 - o(1))`. So the twin (or other designed) columns have large
     height. Corollary 2's conclusion is then met by `C` alone and says nothing more.

   Without the twin columns the argument cannot be rerun. If only the non-twin columns are coherent,
   step 4 has no pair heights with which to bound `|T*|`. For example `G = I` on the non-twin columns
   has no degenerate monomial at all.
3. *Consequence 4 does not apply either.* It assumes *random* non-twin labels. In the suggestion they
   are coherent.

**What is true.** For the suggested design:
* in case (b), Theorem 2 shows that the columns cannot all be certified by size;
* in case (a), nothing is proved.

In both cases, whether the single twin columns and the two-column monomials `e_f(x) +- e_j` exceed
the height condition depends on the `X`-smooth parts of explicit Gaussian integers of height
`>= M^(1/2-o(1))` (case b) or on a pseudo-random comparison of `kappa_x` with `P_x` (case a). That is
the H2-type heuristic, not a theorem.

A further caution concerns weights. Single uncertified columns need not force a uniform weight
increase: with column-dependent weights, a column with `H` hits needs only about `H log X` of
weight, and `sum_x H_x = O(M)` heuristically. Any barrier for this design has to come from
*multi-column* monomials, of which there are about `M^2` (see §3, item (c)).

**Severity: major.** The answer to the task's own suggestion is presented as settled by Theorems 1-2.
It is not. Correct to: "the suggested design is not covered by Theorem 1 or 2. Theorem 2 shows that,
if the point classes are coherent with certified pairs, some columns have height
`>= M^(1/2 - o(1))`. Whether such a design reaches `O(M log M)` is open (heuristically no, by H2)."

## 3. M2: the summary of Theorem 1 claims more than Theorem 1

The formal Theorem 1 and Corollary 1 are correct (with the M3 repair). Three statements in the
summary and in Consequences 2-3 are not.

**(a) "Every completion then has `W >= ...`" (summary (1); table row Thm 1).** Corollary 1 says: *if*
the random design completes with probability `>= exp(-M^(eta/4))`, *then* `W >= ...`. On each parity
class the uniform law has full support. So every deterministic label vector with the forced parities
(coherent, structured, anything) is a possible outcome of the random design. Read literally, "every
completion has `W >= c M log^2 M/loglog M`" would therefore exclude *every* member of `F^all_X` below
that rate on the given primes and profile. That would be a negative answer to P3-all\*, and hence a
growth improvement for actual circles. It is not what is proved. Correct to: "for fixed primes with
`W < (eps(1-eta)^2/(2(C+1)) - o(1)) M log^2 M/loglog M`, such a random design completes with
probability `< exp(-M^(eta/4))`."

**(b) "Independent randomness on any prime set ... is fatal" ("Proved" list; Consequence 3).** This
fails for independent labels with a non-uniform law.
* *Counterexample to the mechanism.* For each `q in Q*` and each `j`, draw `ell_j(q)` independently and
  uniformly from the parity class `nu_j(q)` *minus* the 4-torsion `T4 = {0, m/4, m/2, 3m/4}`. The
  labels are independent, random and parity-respecting (that set is nonempty for `m >= 8`), and
  `H_j = 0` for every `j` and `q`. Theorem 1 and Consequence 3 ("only a set of mass `M^(-Omega(1))`
  may be left random") then say nothing.
* *What survives: two-column monomials.* Suppose that for each `q in Q*` and each parity `nu`, the
  columns of parity `nu` have i.i.d. labels with a common law `mu_(q,nu)`. Then for `j != k` with
  `nu_j(q) = nu_k(q)`,

  ```text
  Pr[ell_j - ell_k in T4] = sum_(C a T4-coset) mu(C)^2  >=  p(q, 0)   (= 4/m if k(q)=2, 8/m if k(q)>=3),
  ```

  by Cauchy-Schwarz. A coset meets a parity class in 2 elements when `k(q) = 2` (`m/4` cosets), and
  lies inside one class when `k(q) >= 3` (`m/8` cosets per class). For the law above (`k = 2`) the
  value is exactly `4/(m - 4)`.

  `Pr[H >= k]` is monotone in each Bernoulli mean. So Theorem 1's step-2 computation applies verbatim
  to `H_(jk)`, the number of `q in Q*` at which `e_j - e_k` collides, with threshold `w_j + w_k <= T_M`.
  Pick disjoint pairs to keep independence; see (c) for the matching. The result is again
  `W >= (c eps - o(1)) M log^2 M/loglog M`, now with `c = (1-eta)^2/(24(C+1))`. *This is a sketch,
  checked by me line by line against the proof of Theorem 1 but not written out in full.*
* If the laws differ between columns, even this fails. Columns given disjoint `T4`-cosets never
  collide in pairs, but such a design is structured rather than random. "Independent randomness" has
  to be quantified, for example by exchangeability within parity classes.

**(c) Reciprocal mass of `Q*` versus `lambda_j`.**
* The theorem needs `lambda_j = sum_(q in Q*) p(q, nu_j(q)) >= (log M)^(-C)` for each column used.
* When `k(q) >= 3`, `p(q, 1) = 0`: a column with `(p_j/q) = -1` never collides alone.
* Consequence 2's example `Q* = {q = +-1 mod 2^K}` lies entirely in `k(q) >= K`. There `lambda_j` is
  `8 sum_(q in Q*, (q/p_j) = 1) 1/m_q`. A lower bound on it for *every* light column
  (`log p_j <= T_M`, `q <= M^(A_0) = p_j^(o(1))`) requires equidistribution of the character
  `(./p_j)` on primes far below its conductor. That follows from GRH, but it is not unconditional.
* So "on a set of primes with reciprocal mass `>= (log M)^(-C)`" (summary (1), Consequence 2) is
  unproved as stated.
* *Repair (parity-agnostic, same pair argument).* For each `q`, at least
  `(1/2 - o(1)) binom(|J|, 2)` pairs have equal parity. So `lambda_(jk) := sum_q p(q,0)[nu_j = nu_k]`
  has average `>= (1/2 - o(1)) Lambda*`, where `Lambda* = sum_(q in Q*) p(q,0) >= 4 mass(Q*)`. Since
  `lambda_(jk) <= Lambda*`, at least a fraction `1/3 - o(1)` of pairs have
  `lambda_(jk) >= Lambda*/4`. Any graph with that edge density has a matching of size
  `>= (1/12 - o(1))|J|`, because a maximal matching of size `mu` meets every edge, so
  `#edges <= 2 mu |J|`. This gives the corrected Corollary 1 with mass in place of `lambda_j` and
  constant `eps/12` instead of `eps`.

**(d) Weights must not depend on the labels (scope; minor).** Theorem 1 fixes the primes *before*
the labels are drawn, and (Ind) requires uniformity *given* the primes. A design that picks the
primes *after* the labels is outside it.
* At `k(q) = 2` such a design is easy. Draw `u` uniform on `Z/m` and set
  `ell = u + (m/4)[parity(u) != nu]`. This is uniform on the class `nu`, and it collides iff
  `u in T4`, whatever `nu` is. One can therefore choose `p_j` after seeing `H_j`, and single columns
  then cost only `sum_j O(H_j log X) = O(M log M)` in total.
* The single-column mechanism is thus not robust to adaptive weights. The pair argument of (b)
  plausibly is: each column has a partner with about `log M/loglog M` pair hits, and a vertex cover
  of these `r` edges has size `>= r/3`. That is a sketch.
* The phrase "for any profile and weights" should read "for any profile and any weights fixed
  independently of the random labels".

**Severity: major** for (a)-(c) taken together. These sentences are the headline "proved" claims of
the round, and as written each is false or unproved. (d) is minor.

## 4. M3: hypothesis (Ind) fails in the write-up's own applications

(Ind) asks that, *conditionally on all residue data of the columns outside `J`* (and of the primes
outside `Q*`), the labels in `J x Q*` be independent and uniform in their parity classes.

* *Consequence 1 ("`sharp` Prop. 6.5 is the case ..., `J` the non-twin columns; Lemma 6.3(2) there
  gives (Ind)").* In `sharp` Lemma 6.3(2) the twin data are
  `k_(f(x)) = eps_x (t - (A k)^(non-twin)_x - sum_y a_(x,s(y)) u_y)` and `k_(s(x)) = u_x - k_(f(x))`.
  Conditioning on them reveals `u_x = k_(f(x)) + k_(s(x))`, and hence
  `(A k)^(non-twin)_x - t` for every `x`. These are `M - 1` independent linear constraints (the
  differences over `x`) on the non-twin labels. Conditionally, the non-twin labels are therefore
  uniform on an affine subspace of codimension `M - 1`, not independent.
* *Consequence 4 (twin columns `J = {f(x)}`).* Conditionally on the non-twin labels (columns outside
  `J`) and on `kappa`, `t_x = eps_x(kappa_x - P_x)` is deterministic.
* *Consequence 2* inherits the same problem when the random part is `sharp`'s fibre.

**Repair (trivial, and it strengthens the theorem).** The proof uses only two facts:
* labels of one column at distinct `q in Q*` are independent and uniform in parity;
* labels of distinct columns in `J` are independent.

Both hold *unconditionally* in all three applications. Restate:

> **(Ind').** The family `(ell_j(q))_(j in J, q in Q*)` is independent, and `ell_j(q)` is uniform on
> the residues of parity `nu_j(q)`. The data outside `J x Q*` may be arbitrary random variables,
> possibly dependent on this family.

Theorem 1 then holds verbatim with (Ind'). "Whatever the data elsewhere" becomes true in a stronger
sense, since (1.2) is checked only for the monomials `e_j`, `j in J'`.

**Severity: major in statement** (the flagship example does not satisfy the stated hypothesis);
**trivial to repair**.

## 5. M4: "matching the parities" is labelled exact, but it is a heuristic, and it omits the obvious match

**The claim** (§3.2, "Why other repairs are worse (exact)", second item; summary (2)): matching
requires `(p_j Norm zeta_j/q) = 1` on the matched primes. That has density `2^(-|S|)` among `p_j`
(or among `zeta_j`). So matching the classes `k > K` needs `2^K >~ X/(log X log p_j)`.

**Problems.**
1. *Density is not cost.* A density `2^(-|S|)` among random `p_j` or `zeta_j` does not bound the least
   admissible `p_j` or `zeta_j` from below. The designer chooses `zeta_j`.
   * With `Norm zeta_j in p_j Q*^2` (`zeta_j = pi_j` up to rational and square factors), every parity
     is matched at every `q`, with no repair.
   * With `Norm zeta_j in 2 p_j Q*^2` (e.g. `zeta_j = (1+i) pi_j`), every class `k(q) >= 3` is matched,
     because `(2/q) = 1` for `q = +-1 mod 8`. Only `k(q) = 2` then needs the `+-i` repair, at cost 4.

   `sharp` §6.6 lists the first option explicitly. §3.2 drops it.
2. *What is true is different.* Put `n = p_j Norm(zeta_j) > 0`.
   * **(Exact.)** If `n in {1, 2} Q*^2`, then `Norm(alpha) Norm(beta) in {1,2} p_j Z^2` for
     `zeta_j = alpha/beta`. So `|gamma_(e_j)|^2 >= p_j`, and the size certificate
     `m_(e_j, s) <= |D_s(gamma)| <= nu_s^(-1) |gamma|` is no better than the trivial necessary bound
     `nu_s^(-1) e^(w_j/2)`. Matching through these `n` is free, but it loses certification by size,
     which is the only reason coherence was introduced (`sharp` §6.6, "Advantage"). The coherent
     design with genuine `pi_j` is the "global" class of `sharp` Remark 2.3(4).
   * **(Under GRH.)** The only real quadratic characters of 2-power conductor are those of `-1`, `2`
     and `-2`, and `n > 0`. So if `n` is not in `{1, 2} Q*^2`, then `q -> (n/q)` twisted by any
     character mod `2^(K+1)` is non-principal. Under GRH it takes the value `-1` on about half of the
     primes of every class `q = c mod 2^(K+1)` in `[M^(1-eta), X]`, as soon as `X` exceeds a fixed
     power of `2^K log n`. (The GRH error `O(sqrt X log^2(n X))` must be small against
     `X/(2^K log X)`.) So under GRH, matching any *whole* class `k(q) = k` on that window
     forces `n in {1,2} Q*^2`, and hence the loss of certification above.
   * *Unconditionally* only Burgess/Linnik-type polynomial statements are available. `sharp` §6.6
     already says so.
3. *Consequence for §3.2's conclusion.* "Parity is a genuine obstruction only on the classes
   `q = +-1 mod 16`" holds for coherent designs with `|zeta_j|^2 < p_j`, and under GRH. It is a
   statement about one family of certified designs, not about P3-all\*.

**Severity: major** (a statement labelled exact is a heuristic, and it is false without the
restriction `|zeta_j|^2 < p_j`). Corrected statement:

> For coherent designs certified by size (`|zeta_j|^2 < p_j`), parity matching on a whole 2-adic
> class over `[M^(1-eta), X]` is impossible under GRH. So a repair is needed at a positive proportion
> of the primes of every class. The root repair costs `2^k log|gamma_v|`; whether cheaper repairs
> exist is open (H1). Without certification (`Norm zeta_j in {1,2} p_j Q*^2`), parity costs nothing
> and the obstruction moves entirely to the pairs (`sharp` §6.6, heuristic).

## 6. M5 and the minor points

**M5 (scope, minor-to-major).** Open point 3 says "*Any* construction at `O(M log M)` must control
smooth parts of `binom(M, 2)` integers of size `exp(M^(1/2 - o(1)))` at once (Cor. 2)". §4.2 adds
"... small for *every* pair". Corollary 2 covers only G-coherent designs on the `sharp` §4 Paley
profile, and it bounds only the *average* pair height (`log|gamma_(a_x - a_y)| >= 0.8 H_xy`, with
`Hbar >= M^(1/2 - o(1))` on average). Correct to: "a G-coherent design on the Paley profile with
certified columns has pair integers of *average* log-height `>= M^(1/2 - o(1))` (indeed
`>= M^(2/3 - o(1))`, see item 6)."

**Minor points.**
1. *Consequence 4: the open caveat holds, and the consequence is redundant.* The author writes
   "surjectivity ... fails modulo 2 for Paley; not resolved". By exact computation, for every Paley
   order `q0` in {43, 47, 59, 67, 71, 79, 83, 103, 107} there is a row set `R`, `|R| = M/2`, on which
   the label matrix `(1 + H')/2` is surjective over `Z`, hence modulo every `m_q/2`. The full row set
   is never surjective, since there are `M` rows and only `M - 1` distinct columns. So with (Ind')
   (M3), Consequence 4 holds with `eps = 1/2` for these orders. For general `q0` it is not proved.
   It is redundant anyway: if the non-twin labels are random, Theorem 1 already applies to the
   `r - 2M = (b-2)M - b` non-twin columns.
2. *Constant in "Prop. 6.5 is the case `C = 0`".* Theorem 1 goes through `m^gen` and the aggregated
   (1.2), and loses a factor 4. At `C = 0` it gives the threshold
   `(1-eta)^2 log^2 M/(2 loglog M)`. The repaired `sharp` Prop. 6.5, which counts hits with one
   fixed unit and uses `m_(v,s) <= nu_s^(-1) e^(w/2)`, gives `2(1-eta)^2 log^2 M/loglog M`. The two
   agree in rate only. Harmless for growth, but the sentence should say "in rate".
3. *Consequence 2, "stays at `Theta(...)`".* Only the lower bound `Omega` is proved for hybrid
   designs. The upper bound `O(M log^2 M/loglog M)` is `sharp` Thm 6.4 for the *fully* random fibre.
   A hybrid with the same rate is not constructed.
4. *Prop. 3.3(2) and Def. 4.1, "Gaussian primes of odd norm".* Read "split Gaussian primes". An inert
   prime has `gamma/conj gamma = 1`, and `gamma_t = q` would divide the modulus at `q = 3 mod 8`.
5. *"Why other repairs are worse (exact)", first item.* The existence of `xi` with
   `Norm xi = O(log^2 q)` is GRH (it is stated, but the item is labelled exact). The "quadratic in
   `||v||_1`" bound is an *upper* bound on one certificate. It does not show that non-root repairs
   are worse. Label the item "remark".
6. *Corollary 2 is weaker than its proof.* With `C <= K' log M`: if (4.1) holds,
   `Hbar >= M^(1 - o(1))`; if it fails, `Hbar^3 >~ M^2/(b C)`. Either way `Hbar >= M^(2/3 - o(1))`.
7. *Theorem 2, step 6* is not needed, since (1.2) is imposed for every `v != 0`.
8. *Attribution.*
   * Step 3 of Theorem 2 is the standard `l^1` support-uncertainty argument
     (`||Hz||_1 >= ||Hz||_2^2/||Hz||_inf`, then Cauchy-Schwarz), as in Donoho-Stark (1989) and
     Elad-Bruckstein (2002). Cite it.
   * Identity (1.4) is `sharp` Lemma 6.3(2) and Cor. 2.2(5) in other coordinates
     (`k_f(x) = eps_x(t - (Ak)^(non-twin)_x - sum_y a_(x,s(y)) u_y)`). It should be attributed, not
     listed as a result.
   * The root-of-unity repair (Lemma 3.2, Prop. 3.3) is elementary group theory. I found no earlier
     occurrence in `docs/claude_review` or in the round-7/8 drafts (grep for Sylow / root of unity).
     As far as I can tell it is new to this project.
9. *Heuristic H1* treats the classes `k > K` through the *single-column* Poisson tail at uniform
   weights. By M2(d), single columns can be absorbed by column-dependent weights at total cost
   `O(M log M)`. So the binding constraint in H1 should be the multi-column monomials. H1 is
   labelled heuristic, so this is a remark on its premise.
10. *Evidence.* The maxima in §5 come from one seed per size. My reruns show a seed-to-seed spread of
    about `2 log M` in the `D1` maxima at `M = 44` (4.22 to 6.13), which is as large as the
    coherent-versus-control differences quoted. The means do reproduce. Also,
    `lambda_max(D1)/(M log M)` is quoted as `~0.96` in §5 and the summary, but as `~0.94` in §6
    item 1.
11. *Corollary 1, `>=` against `<=`.* The boundary case is excluded only because
    `(1-x)^n < e^(-nx)` is strict. Say so, or assume `> exp(-M^(eta/4))`.

## 7. Corrected statement

> **Status: partial barrier for construction methods; no construction; no statement about actual
> circles.** No member of P3-all\* with `W = O(M log M)` was found. The best rate stays `sharp`
> Thm 6.4, `W = O(M log^2 M/loglog M)`.
>
> **Theorem 1 (proved, with (Ind') = unconditional independence).** Fix the primes first. Suppose
> that on a set `Q*` of primes in `[M^(1-eta), X]`, `X <= M^(A_0)`, the labels of a set `J` of
> columns are independent and uniform in their parity classes, and that `>= eps M` columns of `J` have
> `lambda_j >= (log M)^(-C)`. Then for `M >= M_0` the design completes with probability
> `< exp(-M^(eta/4))` unless `W >= (eps(1-eta)^2/(2(C+1)) - o(1)) M log^2 M/loglog M`. It covers the
> random fibre of `sharp` (in rate) and hybrids that are random in this sense on such a `Q*`.
>
> **Extension (sketched, not proved here).** Two-column monomials replace `lambda_j` by the
> reciprocal mass of `Q*` (constant `eps/12`), and replace uniform labels by labels that are i.i.d.
> within parity classes. The theorem does *not* cover adaptively chosen weights, labels with
> column-dependent laws, or deterministic designs.
>
> **Lemma 3.2 / Prop. 3.3 (proved; checked exactly on all 89 odd prime powers `<= 400`).** Repairing
> parity by a 2-Sylow generator `f_q` makes collisions divide `Im(gamma_v^(2^(k(q))))`. So coherent
> designs with `|zeta_j|^2 = M^(O(1))` are certified with constant `2^(K+1)` on the classes
> `k(q) <= K`. On the classes `k(q) <= 3` the parity costs a constant factor. Matching parities
> instead of repairing them either loses certification by size (`Norm zeta_j in {1,2} p_j Q*^2`,
> exact) or, under GRH, is impossible on any whole 2-adic class. Whether the classes `k(q) >= 4`
> admit a repair of bounded cost is open (H1, heuristic).
>
> **Theorem 2 (proved).** In a G-coherent design on the `sharp` §4 Paley profile, under (4.1),
> `W >= (M-1)(theta(X;8,3) - 2 log 2)/(8 ceil(8(C + 2Hbar)))`. Hence at `W = O(M log M)`,
> `C + 2Hbar >= M^(1/2 - o(1))`, and `Hbar >= M^(2/3 - o(1))` if `C = O(log M)`. Pairs and columns
> cannot both be certified by size.
>
> **Not covered by Theorems 1-2:** the task's suggested design (coherent auxiliary columns with
> designed twin/repair columns and structured point classes). There the twin labels are
> deterministic. They are coherent only with height `>= M^(1/2 - o(1))` (when the point classes are
> coherent with certified pairs), or not coherent at all. Its fate is open; heuristically (H2, and
> pair-of-columns Poisson tails) it stalls at the same rate.

## 8. Classification of the statements

| statement | class | valid inputs for actual circles? | status after review |
|---|---|---|---|
| (1.2) | necessary condition on every member of `F^all_X` | yes (`sharp` Prop. 1.2, 6.1) | correct |
| (1.4) | identity for the `sharp` §4 profile | n/a (algebra) | correct; it restates `sharp` Lemma 6.3(2) |
| Thm 1, Cor 1 | method barrier: random designs, primes fixed first | n/a (no actual input) | correct with (Ind'); the summary overstates it (M2) |
| Consequences 1-4 | applications of Thm 1 | n/a | need (Ind') (M3); Cons. 2 example needs `lambda_j` or the pair repair; Cons. 3 overstated; Cons. 4 holds for tested `q0` (exact), but is redundant |
| Lemma 3.2, Prop 3.3, Cor 3.4 | algebra of factored data | yes, the per-unit condition | correct; exact check passed |
| §3.2 "matching" | labelled exact | -- | heuristic; false without `|zeta_j|^2 < p_j` (M4) |
| §3.2 "smaller roots cannot repair" | exact | -- | correct |
| H1, H2 | heuristics (labelled) | -- | labelled correctly; H1's premise is weak (minor 9) |
| Thm 2, Cor 2 | method barrier: G-coherent on the Paley profile | uses (1.2) only | correct; Cor 2 can be sharpened to `2/3` |
| Answer to the suggestion | claimed to follow from Thms 1-2 | -- | **not proved** (M1) |
| §5 | evidence | -- | reproduces (means); the maxima are noisy |
