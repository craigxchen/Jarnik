# Referee report on `tau/L8.md` ((L_6) is false; a Liouville barrier on Mbar_{0,k}; (L_5) via integer tuples)

## Verdict

**The result survives.** I found no fatal or major error. Every load-bearing step of Theorem B was
re-derived, and every exact computation it needs was redone with independent code. That code shares
nothing with the author's scripts. It uses a different extraction of the 25 contact integers, a
symbolic Jacobian, and a non-torsion certificate that does not use Mazur's theorem.

1. **Theorem B ((L_6) is false) is correct.** For every proper closed `Z` it produces anchored
   six-point configurations off `Z` with `eta -> 0`, `w -> infinity`, `Y_j = +-1` and bounded
   residues. The negation of `(L_6)` follows with the quantifiers in the right order. As far as the
   notes show, the refutation is new: `referee_ptolemy.md` §3.4 left `(L_6)` open.
2. **The family itself is not new.** The 3-parameter rational family of "disjoint mixed pencils"
   (four rational common arguments; `lambda_b`, `lambda_c` as the second roots of the branch
   discriminants) is already in `research/docs/cm_disjoint_mixed_pencil_search.md`, eqs. (1)–(6),
   in the normalization `f_0 = x^2/(x-1)`. L8 presents §3.1 and Remark 2 ("why the family is
   rational") without citing it. Three things are new:
   * the generic positive rank (one non-torsion specialization plus Mazur);
   * dominance in `Mbar_{0,6}` (rank-2 Jacobian of a contact point);
   * the use of both to refute `(L_6)`.
3. **Theorem A(a) is correct and correctly cited.** `beta.E >= c_k ||E||` for every psef `E` follows
   from Keel–McKernan Thm 1.3(1): the effective cone of `Mbar_{0,n}/S_n` is simplicial, spanned by
   the `B_i`, for every `n`. Hassett–Tschinkel and Moon et al. quote it the same way. It is routine
   given KM, but it is not in the notes.
4. **Theorem A(b), A(c) and Corollary A are correct only as a formal barrier.** The class of
   "Liouville arguments" they exclude is never defined; "every Liouville estimate obtainable from F
   is equivalent to ..." is not a mathematical statement. A precise version is given in §2.3 below,
   and it holds. Two scope claims in §2.3 and §7 overreach (§2.4).
5. **Proposition C is correct.** The invisible top and singleton blocks need the same inflation as
   in Theorem B. `referee_ptolemy.md` §3.3 had already proved `(L_5)` false; C is a simpler route.
6. **The Section 6.1 data are confirmed exactly.** The interpretive claims built on them are
   overstated (§4).

Checks are in `tau/referee_L8_checks/`. All are exact (Python integers and `Fraction`). Floats
appear only in `log` for the balance columns and in the Euler-product estimate of R7.

---------------------------------------------------------------------------------------------

## 1. Theorem B, step by step

### 1.1 Lemma B.1 (pencil identities)

**R1 (`r1_identities.py`)** works symbolically in `Q[x_1..x_4, z, lam, mu, T]` and verifies:

* `P_lam - z^2 Q_lam = lam prod(z - x_j)`;
* `P_lam Q_mu - P_mu Q_lam = (lam - mu) prod(z - x_j)`;
* `disc_z P_lam = lam(lam(e_3^2 - 4e_2e_4) - 4e_4)`;
* `disc_z Q_lam = lam(lam e_1^2 + 4)`;
* for `Delta_lam(T) = disc_z(P_lam - T Q_lam)`: the constant term is `disc P_lam` and the
  `T^2`-coefficient is `disc Q_lam`.

Hence `Delta_b = kappa_b T(T - s)`, `Delta_c` is linear, and `Delta_a = 4T`. All correct.

B.1(2) on `P^1` needs one more remark, which is easy: the homogenized form of degree 4 has its four
zeros at the `x_j`, so it does not vanish at `infinity`. Consequently any coincidence
`u = v = z` of two moving coordinates forces `f_u(z) = f_v(z)`, hence `z = x_j` and `T = y_j`,
including at branch fibres and at `infinity`.

### 1.2 The open set `U` and Lemma B.2

All the conditions are non-vanishing conditions on rational functions of `t`, so `U` is open.

**R2 (`r2_family.py`)** re-implements the family from the formulas in L8 alone. At
`t_0 = (2,3,-5,1)` and at 60 random rational `t` (50 of which lie in `U`) it checks:

* every listed condition;
* a **direct enumeration of the 32 rational fibre points** over `y_1, ..., y_4`, with the
  collision partition of each 6-tuple. Exactly 25 points carry a collision, each carries exactly
  one cluster, and the clusters are exactly the 25 sides with at most one anchor, each occurring
  once.

The author's own `verify_L6.py`, re-run, gives 36 of 40 random points in `U`, as claimed.

*Genus (B.2(1)).* The three square classes ramify over `{0,inf}`, `{0,s}`, `{t',inf}`. Every
nontrivial product ramifies somewhere: `T Delta_b Delta_c` over `{s,t'}`, and `Delta_b Delta_c`
over all four points. So the cover is a connected `(Z/2)^3` cover with four branch points, each with
inertia of order 2, and Riemann–Hurwitz gives `g = 1`. Correct.

*Transversality (B.2(2)).* The author's derivative argument is correct. There is also a cleaner
argument, which I checked and which shows that the derivative condition in `U` is redundant:

* For moving `u != v`, the poles of `u` and `v` lie over different `T`. R2 checks this, and it also
  follows from B.1(2).
* So `u - v` has degree `4 + 4 = 8`, and it vanishes at the 8 rational fibre points where
  `u = v = x_j`.
* Hence every zero is simple. Similarly, `u - x_k` has degree 4 and 4 zeros.

Pairwise simple zeros mean contact order 1 with `D_C`.

### 1.3 Lemma B.3 and Corollary B.3 (positive rank off a closed set)

R2 reproduces exactly the model quoted in L8,
`E'_{t_0}: Y^2 = X^3 + (912465360/14641) X^2 - (513022816051200/1771561) X`, and `Q_1(t_0)` up to
sign. It also checks `W_j^2 = Delta_b(y_j) Delta_c(y_j)`. It verifies `phi(O) = +-Q_1` and
`phi(P) = -phi(O)` for `O = (x_1,x_1,x_1)` and `P = (x_1, b_1', x_1)`, so `psi(P) = -+2Q_1`.

**Non-torsion, two independent ways.**

1. `nQ_1 != O` for `n` in `{1..10, 12}` (Mazur), at `t_0` and at all 50 random points of `U`.
2. **A certificate that does not use Mazur.** On an integral model, the reductions of `Q_1(t_0)` at
   the good odd primes 23, 31, 41, 43, 47, 53 have orders 14, 16, 22, 22, 28, 26. For odd primes of
   good reduction, reduction is injective on rational torsion, because the formal group has no
   torsion when `p > 2`. A torsion point would therefore have the same order at every such prime,
   so `Q_1(t_0)` has infinite order. The same certificate passes at all 50 random points.

**Corollary B.3.** `E'_U -> U` is an elliptic scheme: `kappa_b kappa_c != 0` and `0, s, t'` are
distinct on `U`. `Q_1` is a section over `U`, `{[n]Q_1 = O}` is closed, and `t_0` lies outside
every such set. Mazur then gives infinite order for every `t` in `U(Q) \ Z_tor`. A morphism of
genus-one curves is an isogeny composed with a translation, so `P` has infinite order on `E_t`.
Correct.

### 1.4 Lemma B.4 (dominance)

**R3 (`r3_dominance.py`)** uses a different method from the author's dual numbers. It restricts to
lines `t_0 + s*delta` and computes `u(s)`, `v(s)` as **reduced univariate rational functions**,
using exact polynomial gcd. It then differentiates symbolically, checks that the directional
derivative is linear in `delta` (consistency), and takes ranks over `Q` and modulo two large
primes. Rank 2 is found for `C = {x_j, a}`, `j = 1, 2, 3`, at three parameter points: `t_0`, the
item-109 member `(1/5, 9/5, 6/5, 1)`, and `(7/2, -4, 9/5, 1)`. Correct.

The geometric argument is correct: the closure `W` of `Psi(E_U)` is irreducible, contains the
2-dimensional `D_C`, and meets `M_{0,6}`, so `dim W = 3`. Two justifications need tightening
(minor):

* *Why the bad set of parameters is closed.* The set `{t : E_t subset Psi^{-1}(Z)}` is closed
  because its complement `pi(E_U \ Psi^{-1}Z)` is open, and that holds because `pi` is smooth,
  hence flat and open. The text instead cites properness, which gives closedness of images, not
  this.
* *`Psi` is a morphism on all of `E_U`.* This deserves one sentence. The collision loci are
  sections over `U` with first-order separations that are uniform in `t`. So one blow-up gives a
  stable family, and the contact points map into `D_C` as used.

### 1.5 Lemmas B.5, B.6 and the actual configurations

*B.5.* I re-read `disjoint_mixed_elliptic_integer_contact_profile.md`. Its proof uses only:

* a smooth model with disjoint sections `B_C` off `Sigma`;
* compact recurrence of `NP`;
* `div(s_ij) = sum_{C containing i,j} [B_C]`, which follows from B.2.

So it applies to every `E_t`, `t` in `U`.

*B.6.* I re-derived it. For `C` containing the anchor the exponent of `d_C` is
`1 + [C contains i,j] - [i in C] - [j in C]`. This equals 1 exactly when `i, j` are both outside
`C`, and 0 otherwise. So `det(P_i,P_j) = (Sigma-part) * prod_{T contains i,j} n_T`, with
`n_T = d_C` for `T = C` or `T = C^c`. The `Sigma`-part is an integer, since its denominator is
`Sigma`-supported and divides an integer that is prime to `Sigma`. `Y_j = +-1`. Correct.

**R4 (`r4_configs.py`): an independent extraction of the 25 contacts on actual points.** For a
rational point `(a,b,c)` with `T = a^2`:

* the good primes dividing `num(T - y_j)` are the primes where the point meets the fibre over
  `y_j`;
* which of the 8 fibre points it meets is read off from `p | a -+ x_j`, `p | b - x_j` and
  `p | c - x_j`;
* the multiplicity is `v_p(T - y_j)`, because the fibre is étale.

Everything is done by gcd-splitting, without factoring and without using any determinant.
`Sigma` is computed from the fixed data of the curve: the discriminant of `E'`, differences of
fibre points, derivatives, leading coefficients and branch values. The determinant identity is
then tested on the Lemma B.6 rows.

The orbit is `Q_1 + N(Q_2 - Q_1)`, for `N <= 14`, at `t_0` and at a second point
`t = (7/2, -4, 9/5, 1)`:

| `t` | `N` | `w` | min/w | max/w | `max|Y_j|` | every `t_ij` `Sigma`-supported | max gcd of blocks |
|---|---|---|---|---|---|---|---|
| `t_0` | 6 | 80.6 | 0.775 | 1.221 | 1 | yes | 1 |
| `t_0` | 10 | 242.3 | 0.862 | 1.123 | 1 | yes | 1 |
| `t_0` | 14 | 492.4 | 0.920 | 1.081 | 1 | yes | 1 |
| second | 10 | 350.8 | 0.881 | 1.119 | 1 | yes | 1 |
| second | 14 | 711.8 | 0.922 | 1.081 | 1 | yes | 1 |

The full table is in `r4_configs.txt`. It agrees with L8 §3.3 (for example `w = 492.4` against
`492.5` at `N = 14`); the small differences come from the different `Sigma`. The `Sigma`-parts
fluctuate with `N` (`log max|t_ij|` up to about 45), as expected. B.5 bounds them only on a
recurrent subsequence. This is data, not proof, and the proof does not need it.

The author's `verify_L6.py`, re-run (`author_verify_L6_rerun.txt`), passes. It runs only to
`N = 12`; the `N = 14` row of L8's table is reproduced by R4.

### 1.6 Theorem B: quantifiers

For a given proper closed `Z`:

1. Choose `t` in `U(Q)`, avoiding `Z_tor` and the proper closed set of B.4. `U(Q)` is Zariski
   dense, so this is possible.
2. Apply B.5 and B.6.
3. `Psi_t(E_t)` is not contained in `Z`, so it meets `Z` in finitely many points.
4. `Psi_t` is birational onto its image: `deg Psi_t^* D_C = 1` forces degree 1. Hence infinitely
   many distinct `Phi_N` avoid `Z`.
5. Inflate the top block by `diag(n, 1)`, which keeps `Phi` and keeps rows primitive because
   `Y = +-1`. Choose the singleton blocks freely, since they enter no determinant.

This yields, for every `(c, eta_0, Z, C_0)`, a configuration with `H_Sigma < c w - C_0`. That is
exactly the negation of `(L_6)` in both `ptolemy.md` §7 and L8 §1. The constant `B_Z` depends on
`Z` through `t`, which is allowed.

### 1.7 Remark 1 (item 109 is a member)

**R5 (`r5_item109.py`)** confirms that `t = (1/3, 3, 2, 5/3)` lies in `U`, with
`s = -53062201/14082120` and `t' = 52156/5643`. More strongly, it identifies the three maps
exactly: `T/(T-4)` composed with `f^{109}_k(2z/(z-1))` equals `f^{fam}_k(z)`, for `k = a, b, c`, at
146 rational test values. Correct.

---------------------------------------------------------------------------------------------

## 2. Theorem A

### 2.1 Part (a)

The proof is correct.

* The psef cone is pointed, so the symmetrization `E_bar` of a nonzero psef `E` is nonzero.
* Symmetrizing effective approximants puts `E_bar` in the closure of the `S_k`-invariant effective
  cone. By KM Thm 1.3(1) that cone is `cone(B_2, ..., B_{floor(k/2)})`. It is simplicial, hence
  closed.
* `beta.E = beta.E_bar / k! = sum a_j |B_j| / k! > 0`, and compactness gives `c_k`.

Citation: a web search confirms that Hassett–Tschinkel ("for each n the cones of S_n-invariant
effective divisors are generated by boundary divisors [KeMc]") and Moon et al. ([KM96, Thm 1.3])
quote KM this way, for all `n`. I could not download the paper itself: arXiv is blocked by the
egress policy here.

*Independent check at `k = 6` (R6, `r6_classes.py`).* I rebuilt the class
`KV = 2H - sum_{i<=5} E_i - E_13 - E_14 - E_23 - E_24` in the 25-dimensional boundary space modulo
the Keel relations, where `N^1` has dimension 16. Checks:

* the `H`-formula is independent of the chosen triple, modulo Keel;
* the class is invariant modulo Keel under the **whole order-48 stabilizer** of `{12|34|56}`,
  including `(56)`, which moves the psi-point. A deliberately wrong quadric class fails this test.

The 15 images give `beta.KV = 5` for all of them. Together with Hassett–Tschinkel, this proves
(a) at `k = 6`. I also confirmed `K.beta = 2^{k-3}(k-8)+k+2` for `k = 4..9`.

### 2.2 Parts (b), (c)

*(c), the reduction to U-invariants.* Each step is correct:

* the ideal of `2 x 2` minors of `[P_i]_{i in T}` is `SL_2`-stable;
* its powers are its symbolic powers (maximal minors of a generic `2 x n` matrix; De Concini,
  Eisenbud, Procesi);
* every coefficient `F_r` of `F(u_c P)` lies in `I^{nu}`, by a Vandermonde argument in `c`;
* the top coefficient is `U`-invariant;
* `U`-invariants are bracket polynomials in `P_0 = e_1, P_1, ..., P_m` (first fundamental theorem
  plus transfer);
* coefficient growth is at most `2^{deg}`.

*(b), the height computation.* The computation is right:

* `h_L = sum nu_T h_{D_T} + h_{E'} + O(1)`;
* `E' == sum a_T D_T` in `Pic = N^1`;
* Lemma 6.2 of `ptolemy.md` gives `h_{D_T} = log n_T + O(eta w + H_Sigma) + O(1)`.

Together these give `h_{E'} >= (beta.E') w - C||E'|| (eta w + H_Sigma) - O_F(1)`.

### 2.3 What (b) actually proves (minor; rewording requested)

"Every Liouville estimate obtainable from F is equivalent to `h_{E'} >= -c`" is not a defined
statement. Read literally, `h_{E'}(Phi) >= -c` holds at every point outside `supp E'`, so "holds
automatically at balanced points" says nothing. Also, "equivalent" should be "implied by", because
the boundary archimedean terms `lambda_{D_T,inf} >= -O(1)` enter as well.

The correct, rigorous content is this. For every nonzero section `F` and every anchored
configuration `Phi` outside `supp div F`,

```text
sum_T nu_T(F) log n_T  <=  h_L(Phi) - (beta.E') w + (2 sum_T nu_T + C_k||E'||) H_Sigma + C_k||E'|| eta w + O_F(1),
```

with `beta.E' >= c_k ||E'||`. That is, the forced divisibility is already bounded by the
archimedean size, with slack `(beta.E')w`. The barrier is the meta-statement that an argument whose
only inputs are (i) `F(Phi) != 0`, (ii) forced divisibility by `prod n_T^{nu_T}`, and (iii) a size
bound through `h_L` and the dictionary bounds of Lemma 6.2 cannot reach a contradiction when
`H_Sigma <= c w`. The scheme should be stated as a definition. With that definition, Theorem A and
Corollary A are correct.

The growing-degree clause ("`O_{E'}(1)` is `o(w)` per unit degree when coefficient heights are
`o(w deg F)`") is asserted, not proved. It is plausible, since all comparison constants are linear
in the degree for a fixed boundary basis.

### 2.4 Scope overreach (minor)

* §7 says "Liouville is excluded (Theorem A). That includes several shells, N-dependent
  coefficients ...". Theorem A is form-blind and contains no `N`. Several shells with
  `N`-polynomial coefficients are the Gaussian setting of `construct.md`, Corollary C3, and the
  attribution should go there. Coefficients of height comparable to `W deg F` lie outside the
  hypothesis of A(b).
* Corollary A's "several tuples at once" needs the union of the tuples to be balanced. That is true
  for a random union, by item 532 at the union size, but not for an arbitrary union chosen by the
  argument.
* "Equivalent to `beta.E' >= 0`" (for the non-strict theorem of `plucker_compatibility_route.md`
  §1) is asserted, not shown. The notes average over all `S` in `[m]`, including the full diagonal.
  Nothing depends on it.

---------------------------------------------------------------------------------------------

## 3. Proposition C

The proof is correct.

* With `infinity` unmarked and the points integral, a prime dividing exactly one of the 10 numbers
  gives a single pair edge of `Mbar_{0,5}`.
* The sieve (two of the 10 forms divisible by `p` is a codimension-2 condition, `p <= 4X`) gives a
  positive proportion. Check: the empirical density 0.139 at `X = 10^6` matches the Euler product
  `prod_{p>=5} (p-1)(p-2)(p-3)(p+6)/p^4 = 0.137` (`r7_propC.py`).
* The tuples are dense in a box, hence Zariski dense.

R7 also verifies the anchored identity `det(P_i,P_j) = +-(x_i - x_j) x_k x_l` with `Y = +-1`
exactly, so `t_ij = +-1` and no top block appears. If `(L_5)` is read with balance on all nonempty
`T`, the top block is `1` and must be inflated by `diag(n, 1)` as in Theorem B. L8 omits this, but
the fix is identical (minor).

`(L_5)` was already shown false in `referee_ptolemy.md` §3.3, which L8 acknowledges. Proposition C
is a simplification, and C1 is its linear subfamily: `X_i = C(s)/u_i` with `u_i = (s - tau_i)/a_i`.

---------------------------------------------------------------------------------------------

## 4. Sections 5–7

* **Section 6.1, numbers.** R6 independently confirms all of the following:
  * `gamma_hat` satisfies every Keel relation for `k = 2..7`;
  * `pi_* gamma_hat = 2^{1-k} beta` for `k = 3..9`;
  * the `K.gamma_hat` table: `-3/8, -3/16, 1/32, 17/64, 65/128, 193/256`;
  * `gamma_hat.KV` takes only the values `1/8` and `1/4`, over all 720 labellings and all 15 KV
    divisors;
  * `gamma_hat.D_S >= 0`.

  Hence `gamma_hat` lies in the movable cone of `Mbar_{0,6}` (Hassett–Tschinkel plus BDPP duality).
* **Section 6.1, interpretation (minor).**
  * The claim "Liouville needs `gamma_hat.(D H - m D_PQ) < 0`, i.e. `m > 2D`" uses a class `H` that
    is never defined in L8, so it cannot be checked as written.
  * "Liouville with forced Gaussian divisibility works iff `gamma_hat` is not movable" is a
    heuristic dictionary; only the direction "movable implies barrier" is meaningful.
  * So "independently confirms `construct.md`'s Corollary C3 at `k = 4`" should read "is consistent
    with".
  * That the normalized height vector of circle points converges to `gamma_hat` is asserted, not
    proved. It is plausible: the short arc gives `W/4` and each Gaussian block `W/2^k`.
* **Section 5.** These are correctly labelled as open or heuristic.
  * `k7/cubicnet.py` computes its rank over `Q` with `Fraction`, so no second prime is needed.
    Its docstring says the degree is `d = 2^{k-4}`, but it uses `d = 6` at `k = 7` (not 8) and
    `d = 3` at `k = 6`. The write-up's "sextic" agrees with the code, so only the docstring is
    wrong.
  * The "overdetermined by 10" and "e^{8w}" counts are heuristics and are not used in any proof.
* **Section 4, integral chart.** The "2-(6,3,2) design is forced" paragraph is presented with
  "must", but no proof of the balance-implies-design step is given. It is explanatory only.

---------------------------------------------------------------------------------------------

## 5. Novelty

* **New:** the refutation of `(L_6)`, meaning dominance plus generic positive rank of the
  disjoint-mixed family; Theorem A(a) in this project, with its strict margin; the U-invariant
  reduction; and the `gamma_hat` computations.
* **Not new:** the rational 3-parameter family and the mechanism that makes it rational. Both are
  in `cm_disjoint_mixed_pencil_search.md`: there the split quartic `R(x)`, eq. (1); `lambda_b`,
  eq. (3); `lambda_c`, eq. (4); and the branch values, eq. (6).
* **Not new:** `(L_5)` false (`referee_ptolemy.md` §3.3) and the non-strict collision-order barrier
  (`plucker_compatibility_route.md` §1). L8 cites both.
* Grep finds no dominance or Zariski-density statement for six-curves in the notes. Item 112
  covers only seven rows.

## 6. Corrected statement

**Theorem B (correct as stated).** `(L_6)` is false. For every proper Zariski-closed `Z` in
`Mbar_{0,6}`, there are anchored configurations off `Z` with the following properties:

* pairwise coprime blocks, including inflated top and singleton blocks;
* `|log n_T - w| <= eta w` with `eta = O(w^{-1/2})` and `w -> infinity`;
* `Y_j = +-1` and `|t_ij| <= B_Z`.

They come from the 3-parameter rational family of disjoint mixed pencils of
`cm_disjoint_mixed_pencil_search.md`, here normalized by `f_a = z^2`. That family should be cited;
it contains item 109 at `t = (1/3, 3, 2, 5/3)`.

**Theorem A.**

* *(a) (correct).* `beta.E >= c_k ||E||` for every psef `E` on `Mbar_{0,k}`, `k >= 4`.
* *(b, c) (correct as reformulated in §2.3).* For every nonzero section, and every nonzero
  polynomial in the matrix entries via the U-invariant reduction, the forced boundary divisibility
  at an anchored configuration is at most `h_L(Phi) - (beta.E') w + O(deg (H_Sigma + eta w)) + O_F(1)`.
  Hence a Liouville argument, defined as in §2.3, cannot prove `(L_k)` or any bound
  `H_Sigma >= f(w)` with `f -> infinity`. This covers fixed `F`, and growing degree with coefficient
  height `o(w deg)`, the latter asserted rather than proved.
* The form-blind barrier does **not** by itself cover `N`-dependent multi-shell forms. Those are
  `construct.md`, Corollary C3.

**Proposition C (correct, with top-block inflation).** `(L_5)` is false via integer 5-tuples
`(0, x_1, ..., x_4)` whose 10 differences are pairwise coprime away from 2 and 3.

## 7. Files (`tau/referee_L8_checks/`)

| File | Content |
|---|---|
| `rl.py` | sparse multivariate polynomials, elliptic-curve arithmetic over `Q` and `F_p`, ranks over `Q` and `F_p` |
| `fam_ref.py` | the family rebuilt from L8's formulas; `U`; fibre enumeration; `E'`; `phi` |
| `r1_identities.py` / `.txt` | B.1, symbolic |
| `r2_family.py` / `.txt` | `U`; 25 contacts by enumeration; disjoint poles; `E'`, `Q_j`; non-torsion by Mazur and by a Mazur-free reduction certificate; `phi` relations; 50 random points |
| `r3_dominance.py` / `.txt` | rank-2 Jacobian by reduced univariate symbolic differentiation, checked over `Q` and mod two primes |
| `r4_configs.py` / `.txt` | actual orbits for `N <= 14` at two parameter points; contacts extracted from fibres (not determinants); determinant identity with `Sigma`-supported `t_ij`; `Y = +-1`; coprimality; balance |
| `r5_item109.py` / `.txt` | item 109 is the member `(1/3, 3, 2, 5/3)`: maps identified exactly |
| `r6_classes.py` / `.txt` | Keel relations for `beta` and `gamma_hat`; pushforward; `K.gamma_hat`; KV class invariant under its order-48 stabilizer; `beta.KV = 5`; `gamma_hat.KV` in `{1/8, 1/4}` |
| `r7_propC.py` / `.txt` | Proposition C: anchored identity; sieve density |
| `author_verify_L6_rerun.txt` | re-run of the author's `verify_L6.py`: passes |

Sources consulted for the KM citation:

* [Keel–McKernan, arXiv alg-geom/9607009](https://arxiv.org/pdf/alg-geom/9607009), via search
  summary; the PDF itself is blocked here.
* [Hassett–Tschinkel, math/0110231](https://arxiv.org/pdf/math/0110231).
* [Moon et al., 1406.2196](https://arxiv.org/pdf/1406.2196).
