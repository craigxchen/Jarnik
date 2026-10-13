# Round 10: additive and lattice-geometric mechanisms

*Unrefereed.* Labels:
* **proved**: all quantifiers, proof given here;
* **exact**: finite computation in exact integer arithmetic (scripts in `additive_checks/`);
* **evidence**: finite computation that no proof uses;
* **heuristic**: no proof.

Every statement says which inputs it uses, and whether those inputs hold for actual circles and for
fake circles (angle models, `round4/outside.md` §1).

## 0. Summary and verdict

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` stands. What this
round proves is a map of the additive inputs. This is a reading of the results below, not a theorem.
* Each additive identity examined is either local (already capped) or uses the existence of the
  relation lattice. Existence of the relation lattice is lattice-ness (Thm 2.1).
* Lattice-ness, used in full, is the problem itself (Thms 2.2, 2.3, 2.6).
* One height-graded input, `AG(M^A)`, separates actual circles from almost every fake. Turning it into
  a bound needs a missing lemma.

**Proved results.**
1. **Collapse theorems (§2).**
   * Thm 2.1: the relation lattice `K` has rank `M - 3` iff the points lie in a lattice coset.
   * Thm 2.2 (rigidity): `M >= 5` points of a circle that satisfy all relations of an actual cluster are
     that cluster up to similarity. So a fake with the actual relation lattice is the actual cluster.
   * Thm 2.3: a lattice-coset configuration of `>= 5` points on a circle is a norm-form problem for some
     imaginary quadratic order, with no bound on its discriminant from the hypotheses.
   * Thm 2.6: residues modulo `N` together with "small lifts" of the differences reconstruct the cluster
     exactly. Every `c` in `K_N = {c : sum c z = 0 mod N}` with `||c||_1 < N/L` is an exact relation.
   * Prop 2.7: the invariants of `K` are local data. The index `I` of the difference lattice is even, its
     odd primes are those where the cluster has `<= 2` residues, and projections of `K` have index
     `A_123/I`.
2. **The additive gap (§2.2, §2.5).** `AG(H)`: every zero-sum `c` with `||c||_1 <= H` has `sum c z = 0` or
   `|sum c z| >= 1`.
   * Actual circles satisfy it at every height.
   * Every configuration either has an exact relation or violates `AG` at height `M^(c_0/2 + o(1))`
     (Thm 2.4(1)). Almost every fake of a positive-measure family has no relation (Cor 2.5).
   * `AG` at height `e^(W/2 + O(M))` forces a lattice coset (Thm 2.4(2)).
   * `AG(M^A)` forces difference-module rank `<= 2 + (c_0/(2A)) M + o(M)` (Prop 2.8).

   So **`AG(M^A)` with `A > c_0/2` (and `M` large) is an input that is true for actual circles and
   provably false for almost every fake.** A dimension count (heuristic) puts fakes with independently imposed relations
   below `A = c_0`.
3. **Pair gcds at a power of `N` (Prop 3.1).** Every pair outside one matching has
   `Norm g_xy >= N^(1/4)/(sqrt2 C)`, i.e. `mu_xy <= W/4 + log C + (1/2) log 2`. This answers
   `round9/L1L2.md` open point 1 with `delta = 1/4`. **It uses only monomial gaps, so it holds for the
   fakes and is capped.**
4. **Effective height of the arc's position (Thm 3.2; actual only).** Let
   `h* = min_v |v|(1 + N^(1/4) eta(v)/C)` over primitive lattice directions `v`, where `eta(v)` is the
   angular distance of the arc from the line `R v`. Then:
   * `M <= C^2 h* + 2`;
   * every pair whose chord is not parallel to `i v*` has `mu_xy <= 2 log(sqrt2 C h*)` and
     `|B_xy|/t_xy <= C h*` (nearly sharp on the data);
   * `h* <= (2 pi/C)^(1/2) N^(1/8) + 1`.

   Cor 3.3: in **low position** (`h* <= M^b`) round 9's (R1-arch)/L1 holds off one matching. Inputs:
   lattice lines and lattice directions (false for fakes). But fakes in generic position meet every
   conclusion, so this is a restricted (R1-arch), not a growth step. It contains Lemma 1.3 of round 9 and
   Prop C.2 of round 4 as special cases.
5. **p-adic structure of relations (§4).** The Newton polygon and leading-residue cancellation of any
   relation are given (Lemmas 4.1-4.2). An isolation bound `p_j <= lambda_(M-3)(K)` holds for columns with
   a unique extreme row (Prop 4.3). Applied to four-point relations it gives `p_j <= C^3 N^(1/4)/2` for
   `M >= 5`, which is the triangle identity (profile data).
6. **Other additive tools (§5)** reduce to known classes:
   * sums and Sidon properties are pair data;
   * the `sigma`-cocycle is triangle areas;
   * the `F_q` polynomial method sees exactly `K ⊗ F_q` plus residue partitions (proved);
   * determinant and Bombieri-Pila variants are dominated by the pair data or capped at `tau* < 2`.

**Exact checks** (`additive_checks/`):
* Prop 3.1 and Thm 3.2: 1,192 runs, including all runs of the 7-, 8- and 9-point record circles, with no
  violation;
* Prop 2.7 and Thm 2.6(1): 468 clusters with `M >= 4`, with no failure.

**Barrier status of the mechanism** (§6). The only input found that is true for actual circles and
provably false for the fakes, other than full lattice-ness, is `AG(M^A)` with `A > c_0/2`. The missing
step is Lemma (AG) of §7 (`forall c_0 exists A`, the quantifier order of round 9's (T')), which would
imply growth. Its extreme case `rho = 2` is the norm-form problem for an arbitrary imaginary quadratic
order. Lattice lines and directions (Thm 3.2) separate actual circles from fakes only in low position,
and there they give (R1-arch), which the fakes satisfy anyway.

## 1. Setting

Conventions are those of `round9/L1L2.md` §1 and `round7/README.md`.

* `z_1, ..., z_M` are distinct lattice points on `|z|^2 = N = R^2`, on an arc of angular width
  `D = L/R` and length `L <= C sqrt R`. So `D <= C N^(-1/4)`. When a statement uses the cluster's own
  constant, `C := D N^(1/4)`.
* *Normalised*: the Gaussian gcd of the points is 1. Then `N` is odd with all prime factors
  `= 1 mod 4`, and
  `z_x = i^(k_x) prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj)`, `min_x a_xj = 0`, `max_x a_xj = e_j`.
* `W = log N`, `g_xy = gcd(z_x, z_y)`, `d_xy = sum_j |a_xj - a_yj| log p_j`, `Norm g_xy = N e^(-d_xy)`.
  The pair margin is `mu_xy = d_xy - W/2`, so `Norm g_xy = N^(1/2) e^(-mu_xy)`.
* `beta_xy = (z_x - z_y)/g_xy = i^(k_x) gamma_s B_xy` with `B_xy` a nonzero rational integer and
  `|gamma_s| in {2, sqrt 2}` (`round9/L1L2.md` Thm 1.1).
* `w(v) = sum_j |v_j| log p_j`, `A_v = prod_j pi_j^(v_j^+) conj(pi_j)^(v_j^-)`, `Norm A_v = e^(w(v))`.
* `Lambda_0 = {c in Z^M : sum_x c_x = 0}`. The **relation lattice** is
  `K = {c in Lambda_0 : sum_x c_x z_x = 0}`. The **difference lattice** is
  `Lambda = Z-span{z_x - z_y}`, a subgroup of `Z[i]`, and `I = [Z[i] : Lambda]`.
* A *fake circle* (angle model) has genuine primes `p_j`, a genuine exponent profile `(a_x)`, windings
  `k_x`, positions `delta_x in [0, D]` and real angles `phi` with
  `<2 a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x`, satisfying the monomial gap inequalities
  (`round4/outside.md` Lemma 0.1) and any further residue data of its class. Its "points" are
  `w_x = R e^(i(theta_0 + delta_x))`. They are not lattice points. For fixed combinatorial data the
  positions `delta` range over a subset of `[0, D]^M`, and the barrier constructions choose them.

**The question.** Which exact additive facts about actual clusters (integer combinations of the points,
chords, sums, lattice lines) are (i) true for actual clusters, (ii) false for the fakes, and (iii)
strong enough to beat `log R/loglog R`? The barrier map (`docs/claude_review/README.md`) says that
anything a fake satisfies cannot beat `log R/loglog R`.

## 2. The additive content of a cluster is lattice-ness

This section proves that every exact additive fact about a cluster is a consequence of one property:
the points lie in a coset of a lattice. It then measures how much of that property each weaker
additive input retains. `round9/L1L2.md` Prop 5.2 already says that the content of exact relations is
"lattice-ness". Here that sentence is made into theorems, with the quantifiers and heights needed for
the barrier check.

### 2.1 Relations and lattice cosets

**Theorem 2.1 (relation lattice = lattice coset; proved).** Let `w_1, ..., w_M in C`, `M >= 3`, not all on
one real line. Put `Gamma = Z-span{w_x - w_1}` and `K(w) = {c in Lambda_0 : sum c_x w_x = 0}`. Then
`rank K(w) = M - 1 - rank_Z Gamma <= M - 3`, and the following are equivalent:
1. `rank K(w) = M - 3`;
2. `Gamma` is discrete, i.e. a lattice in `C = R^2`;
3. the points lie in one coset `w_1 + Gamma` of a lattice.

*Proof.* The map `Lambda_0 -> Gamma`, `c -> sum c_x w_x = sum c_x (w_x - w_1)`, is onto (take
`c = e_x - e_1`), with kernel `K(w)`. So `rank K = (M - 1) - rank Gamma`. `Gamma` spans `R^2`, since
the points are not collinear, so `rank Gamma >= 2`. If `rank Gamma = 2`, `Gamma` is generated by two
elements, and these span `R^2`, so they are `R`-independent and `Gamma` is a lattice. Conversely a
lattice has rank 2. (3) is a restatement of (2). QED

For an actual cluster `Gamma = Lambda` is a sublattice of `Z[i]`, so `rank K = M - 3`. For a fake,
`K(w)` is in general `0` (Corollary 2.5).

**Theorem 2.2 (rigidity; proved).** Let `z_1, ..., z_M` (`M >= 5`) be distinct lattice points on
`|z| = R`, and let `w_1, ..., w_M` be distinct points on some circle `|w - w_0| = rho`. Suppose every
relation of the `z` is a relation of the `w`: `K(z) ⊂ K(w)`. Then there is `lambda in C^*` with
`w_x = w_0 + lambda z_x` for all `x`, or `w_x = w_0 + lambda conj(z_x)` for all `x`.

*Proof.* Since `K(z) ⊂ K(w)`, the map `Lambda -> C`, `z_x - z_1 -> w_x - w_1`, is a well-defined group
homomorphism. Indeed, two integer combinations of the `z_x - z_1` with the same value differ by an element
of `K(z)`. It extends `R`-linearly to `T : R^2 -> R^2`, and `w_x = T z_x + b` with `b = w_1 - T z_1`.
`T` is invertible, since otherwise the `w_x` would be collinear, and three distinct points of a circle
are not. So the `w_x` lie on the ellipse `E = T(circle) + b` and on the circle `S` around `w_0`. Two
distinct conics without a common component meet in at most 4 points, and `E` and `S` are irreducible.
Hence `E = S`. Comparing centres, `b = w_0`. `T` maps a circle centred at `0` onto a circle centred at
`0`, so `T` is a similarity: `T z = lambda z` or `T z = lambda conj z`. QED

So an angle model whose points satisfy all relations of an actual cluster is that cluster, rotated,
scaled and possibly reflected. Its positions `delta_x` are then the actual ones. By `round4/outside.md`
Thm A.1(i), its angles `phi` agree with the genuine ones on the span `L` of the exponent differences
(up to the reflection).

**Theorem 2.3 (lattice cosets on a circle are norm-form problems; proved).** Let `M >= 5` points of a
coset `w_1 + Gamma` (`Gamma` a lattice with basis `gamma_1, gamma_2`) lie on a circle `|w - w_0| = rho`.
Then:
1. the Gram matrix `G = (<gamma_a, gamma_b>)` is a positive real multiple of a rational matrix;
2. `w_0 - w_1 in Q gamma_1 + Q gamma_2`.

Hence there is a lattice `Gamma' ⊃ Gamma` of finite index with `w_0 in w_1 + Gamma'`, and after translating
`w_0` to `0` and scaling, the points are vectors `alpha` of a lattice with an integral positive definite
binary quadratic form `F` and `F(alpha) = n` for one integer `n`. Equivalently they are elements of fixed
norm in a fractional ideal of an order of an imaginary quadratic field.

*Proof.* Write `w_x = w_1 + m_x gamma_1 + n_x gamma_2`. The circle condition reads
`G_11 m^2 + 2 G_12 m n + G_22 n^2 + 2 l_1 m + 2 l_2 n + kappa = 0` with `l_a = <w_1 - w_0, gamma_a>` and
`kappa = |w_1 - w_0|^2 - rho^2`. So the vector `X = (G_11, G_12, G_22, l_1, l_2, kappa)` lies in the
kernel of an integer `M x 6` matrix with rows `(m^2, 2mn, n^2, 2m, 2n, 1)`. Its real kernel is spanned by
rational vectors. If it had dimension `>= 2`, the points `(m_x, n_x)` would lie on two independent conics.
They lie on the irreducible conic defined by `X`, so the other conic would share 5 points with it and
coincide with it, a contradiction. So the kernel is the line `R X` and contains a nonzero rational vector,
whence `X` is a real multiple of a rational vector. This gives (1). The centre `w_0 - w_1 = s gamma_1 +
t gamma_2` solves `G (s, t)^T = -(l_1, l_2)^T`, a rational system up to the common factor, which gives
(2). The rest is the standard dictionary between lattices with rational Gram matrix and ideals of
imaginary quadratic orders. QED

**Reading.** "The points lie in some lattice coset" is the generalised circle problem for an arbitrary
positive definite integral binary quadratic form, whose discriminant is not bounded by the hypothesis.
The Gaussian problem is the special case where `Gamma` is similar to `Z[i]` and the centre is a lattice
point. An argument that uses additive structure only through the existence of the relation lattice
treats every imaginary quadratic order alike.

### 2.2 How much of lattice-ness the weaker additive inputs carry

**Definition (additive gap).** A configuration `w` satisfies `AG(H)` if every `c in Lambda_0` with
`||c||_1 <= H` has `sum c_x w_x = 0` or `|sum c_x w_x| >= 1`. Actual clusters satisfy `AG(H)` for all `H`,
since `sum c_x z_x in Z[i]`.

**Theorem 2.4 (the additive gap; proved).** Let `w_1, ..., w_M` (`M >= 4`) satisfy `|w_x - w_1| <= L`
with `L >= 1`, and not all on one line.
1. *(Decision at polynomial height.)* Let `H >= 1` be an integer with
   `(H + 1)^(M-1) > (3(M-1) H L + 1)^2`. Then some nonzero `c in Lambda_0` with `||c||_1 <= 2(M-1)H`
   has `|sum c_x w_x| < 1`. So either the `w` have an exact relation of that height, or they violate
   `AG(2(M-1)H)`.
2. *(All heights force a lattice.)* If `H + 1 > 2^(M-2) M^2 L^2` and `w` satisfies `AG(2(M-1)H)`, then
   `rank K(w) = M - 3`, i.e. `w` is a lattice-coset configuration (Theorem 2.1).

For an arc of length `L = C N^(1/4)` and `W <= c_0 M log M`, the least admissible `H` in (1) is
`H_1 <= M^(c_0/2 + o(1))` (as `M -> infinity`, `C` fixed). The height in (2) is `e^(W/2 + O(M))`.

*Proof.* (1) Take the `(H+1)^(M-1)` vectors `c' in {0, ..., H}^(M-1)` and the values
`V(c') = sum_(x >= 2) c'_x (w_x - w_1)`. They lie in a disc of radius `(M-1) H L`, which is covered by
`(3(M-1)HL + 1)^2` squares of side `2/3` (diameter `< 1`). Two values share a square. Their difference
`c''` has entries in `[-H, H]`. Put `c_1 = -sum c''_x` and `c_x = c''_x` for `x >= 2`. Then
`||c||_1 <= 2(M-1)H`, `c != 0`, and `|sum c_x w_x|` is the distance between the two values, `< 1`.
(2) Let `kappa = rank K(w)` and let `K'` be the projection of `K(w)` to the coordinates `2..M`, which is
injective. Two vectors `c', c'' in {0..H}^(M-1)` have the same value iff `c' - c'' in K'`. A rank-`kappa`
lattice has at most `(2H+1)^kappa` points in `[-H, H]^(M-1)`, since some `kappa` coordinates determine its
points. So there are at least `(H+1)^(M-1)/(2H+1)^kappa` distinct values. By `AG(2(M-1)H)` distinct values
are at distance `>= 1`, and they lie in a disc of radius `(M-1)HL`, so there are at most
`(2(M-1)HL + 1)^2 <= 4 M^2 L^2 (H+1)^2` of them. If `kappa <= M - 4`, then, using `2H + 1 <= 2(H+1)`,
`(H+1)^(M-1) <= 2^(M-4) (H+1)^(M-4) 4 M^2 L^2 (H+1)^2`, i.e. `H + 1 <= 2^(M-2) M^2 L^2`. Hence
`kappa >= M - 3`, and `kappa <= M - 3` by Theorem 2.1. For the heights: in (1), `H^(M-3) >~ 9 M^2 L^2`
gives `log H_1 = (W/2 + O(log M + log C))/(M - 3)`, which is `(c_0/2 + o(1)) log M`. QED

**Corollary 2.5 (fakes violate the additive gap almost surely; proved).** Fix `M >= 4`, an angle `theta_0`
and integers `k_x`. Let `P_rel ⊂ R^M` be the set of position vectors `delta` for which the points
`e^(i(theta_0 + delta_x))` have an exact relation, i.e. `K != 0`. Then `P_rel` has Lebesgue measure 0.
Hence, for any family of angle models whose admissible positions form a set of positive measure, almost
every member violates `AG(2(M-1)H_1)` with `H_1` as in Theorem 2.4(1).

*Proof.* For fixed `c != 0` in `Lambda_0` let `F_c(delta) = sum c_x e^(i delta_x)`. Choose `x` with
`c_x != 0`. Then `dF_c/d delta_x = i c_x e^(i delta_x) != 0`, so at every zero either `Re F_c` or
`Im F_c` has nonzero gradient. By the implicit function theorem the zero set is locally contained in a
`C^1` hypersurface and is therefore null. `P_rel` is the countable union of these zero sets. The
consequence follows from Theorem 2.4(1). QED

So `AG` at polynomial height `M^(c_0/2 + o(1))` is an input that is **true for actual circles and false
for almost every fake** of any positive-measure family. Whether the barrier constructions can be modified
to satisfy it is discussed in §6.

### 2.3 The local-global form: residues modulo `N` plus small lifts

**Theorem 2.6 (proved).** Let `z_1, ..., z_M` be an actual cluster with chords `<= L`, and suppose
`2 R L + L^2 < N`, which holds when `L <= C N^(1/4)` and `N^(1/4) > 2C + C^2`.
1. *(Local-global relation principle.)* Let `K_N = {c in Lambda_0 : sum c_x z_x = 0 mod N Z[i]}`, a
   lattice of rank `M - 1` and index at most `N^2` in `Lambda_0`. Every `c in K_N` with
   `||c||_1 < N/L` lies in `K`. If moreover `N^(1/2) > 2 C^4`, these vectors span `K ⊗ Q`.
   Consequently `K_N` has at most `M - 3` successive minima (in `l_1`) below `N/L`, and exactly `M - 3`
   under that condition. The last two are `>= N/L = N^(3/4)/C`.
2. *(Small-lift reconstruction.)* Conversely, let `z_1` be a lattice point with `|z_1|^2 = N`, let
   `r_x in Z[i]/N` satisfy `r_x conj(r_x) = 0` in `Z/N`, and suppose `r_x - z_1` has a representative
   `w_x in Z[i]` with `|w_x| <= L`. Then `|z_1 + w_x|^2 = N` exactly.

*Proof.* (1) For `c in K_N`, `sum c_x z_x = sum c_x (z_x - z_1)` lies in `N Z[i]` and has absolute value
at most `||c||_1 L`. A nonzero element of `N Z[i]` has absolute value `>= N`. The index bound holds
because `Lambda_0/K_N` embeds in `Z[i]/N`. For the spanning statement, the `M - 3` four-point relations
`(A_23x, -A_13x, A_12x, -A_123)` on the points `1, 2, 3, x` (`x >= 4`, `A = 2 Area`, Cramer's rule) are
independent, lie in `K`, and have `||c||_1 <= 4 max 2 Area <= 2 C^3 N^(1/4) < N/L`, using
`2 Area = |z_a - z_b||z_b - z_c||z_a - z_c|/(2R) <= L^3/(2R)`. The successive-minima statement follows
because `K` has rank `M - 3`. (2) `|z_1 + w_x|^2 = N + 2 Re(conj(z_1) w_x) + |w_x|^2` lies in `(0, 2N)` because
`|2 Re(conj z_1 w_x) + |w_x|^2| <= 2RL + L^2 < N`. It is `= r_x conj(r_x) = 0 mod N`. So it equals `N`. QED

*Reading.* Residues modulo `N`, i.e. the level data at the primes of `N` that the P3 classes already
contain, together with the archimedean statement "the differences have small lifts", pin the cluster
down: one actual lattice point and the residues reconstruct all the others. For random level data the
small-lift property has probability about `(pi L^2/N^2)^(M-1)` (heuristic). So "level data + small lifts"
is not an intermediate input. It is the actual problem.

### 2.4 Exact invariants of the difference lattice

**Proposition 2.7 (proved; exact check).** For a normalised actual cluster with `M >= 3`:
1. `I = [Z[i] : Lambda]` is the gcd of all `2 x 2` minors of the differences, i.e. the gcd over triples
   of `2 Area`. It is even, because every point is `= 1 mod (1 + i)` (`N` odd).
2. For an odd prime `q` not dividing `N`, `q | I` iff the cluster has at most 2 residues modulo `q`.
3. For any three points `1, 2, 3`, the projection of `K` to the other `M - 3` coordinates has index
   `A_123/I` in `Z^(M-3)`, where `A_123 = 2 Area(z_1, z_2, z_3)`.

*Proof.* (1) The index of the lattice spanned by integer vectors is the gcd of the maximal minors, and
`det(z_y - z_x, z_w - z_x) = +-2 Area`. The minors of differences from one base point generate the
same ideal as all triples, since `2 Area` is a cocycle in the base point. `N` odd forces `Re z + Im z` odd
for every point, so all differences lie in `(1 + i) Z[i]`, of index 2. (2) `q | I` iff the differences
reduce to an `F_q`-line in `F_q^2`, i.e. the residues lie on an affine line. That line meets the
nondegenerate conic `x^2 + y^2 = N` over `F_q` in at most 2 points. Conversely 2 residues are collinear.
(3) `c -> (c_x)_(x >= 4)` maps `K` onto `{c'' : sum_(x>=4) c''_x [z_x - z_1] = 0 in Lambda/Lambda_123}`,
where `Lambda_123` is spanned by `z_2 - z_1` and `z_3 - z_1`. The map is injective because a relation
supported on three points would make them collinear. Its image has index `|Lambda/Lambda_123| = A_123/I`.
QED

*Exact check* (`collapse.py`, output `collapse_out.txt`). 468 actual normalised clusters with
`4 <= M <= 9` (random split `N < 10^14` at `C <= 4`, plus the record circles at `C <= 8`) pass:
* `rank K = M - 3`, from an exact kernel basis;
* items 1-3 of Proposition 2.7; the odd part of `I` exceeds 1 in 48 clusters;
* Theorem 2.6(1): in every cluster an LLL-reduced basis of `K_N` has exactly `M - 3` vectors with
  `||c||_1 < N/L`, and each of them is an exact relation.

Item 2 shows that the index `I` is residue data (class `A_X`). Item 3 shows that the covolume of `K` is
controlled by triangle areas, which are profile data (`round4/outside.md` Prop C.1). So the *invariants*
of `K` are known local data. What is not local is the existence of `K` itself (Theorem 2.1) and its
compatibility with the circle (Theorems 2.2, 2.3).

### 2.5 The additive-gap ladder `AG(M^A)`

Theorem 2.4 places the additive gap between two heights. The following quantifies what `AG` at a
polynomial height `M^A` forces. Write `rho = rank_Z Gamma = M - 1 - rank K(w)` for the rank of the
difference module of a configuration.

**Proposition 2.8 (rank bound; proved).** If `w` (as in Theorem 2.4) satisfies `AG(2(M-1)H)`, then

```text
(H + 1)^(rho - 2) <= 2^M M^2 L^2,     i.e.   rho <= 2 + (2 log L + M log 2 + 2 log M)/log(H + 1).
```

For `L = C N^(1/4)`, `W = c_0 M log M` and `H = M^A` this reads `rho <= 2 + (c_0/(2A)) M + O_C(M/log M)`.

*Proof.* This is the counting of Theorem 2.4(2) with `kappa = M - 1 - rho`:
`(H+1)^(M-1) <= (2H+1)^kappa (2(M-1)HL + 1)^2 <= 2^kappa (H+1)^kappa 4 M^2 L^2 (H+1)^2`. QED

**Heuristic 2.9 (where the fakes sit on the ladder).** Consider configurations with positions in a
positive-measure family (the fakes).
* *No relations (`rho = M - 1`).* `AG` fails at the height `H_1 = M^(c_0/2 + o(1))` of Theorem 2.4(1).
  This is proved, by Corollary 2.5, for almost every member. A Jacobian count suggests that it
  holds below `M^(c_0/4)`. For fixed `c`, the value `R sum c_x e^(i delta_x)` has first-order part of size
  `R|c|D` and curvature part of size `|c| C^2`, so `P(|sum c w| < 1) ~ 1/(|c|^2 C^3 N^(1/4))`. Summing
  over `||c||_1 <= H` gives `H^(M-3)/(C^3 N^(1/4))`, which is `< 1` iff `H < M^(c_0/4 + o(1))`.
* *Imposed relations.* `k` relations are `2k` real equations on `M - 1` free positions. For independent
  equations, solutions exist in a family only if `k <= (M-1)/2`, i.e. `rho >= (M-1)/2`. Then Prop 2.8 caps
  the attainable height at `log H <= (W/2 + O(M))/(rho - 2) = (c_0 + o(1)) log M`.
* So `AG(M^A)` with `A > c_0` cannot be met by any configuration whose relations are imposed
  independently on free positions. By Prop 2.8 it forces `rho <= (c_0/(2A)) M + o(M) < (M-1)/2`, an
  overdetermined system. Actual circles meet it because `rho = 2`.

The ladder therefore has this shape:

| height `H` | actual circles | relation-free fakes | relation-imposed fakes |
|---|---|---|---|
| `H < M^(c_0/4)` | yes | yes (heuristic) | yes |
| `H >= M^(c_0/2 + o(1))` | yes | **no** (proved, a.e.) | possibly |
| `H >= M^((1 + eps) c_0)` | yes | no | **no** if relations are independent (heuristic dimension count) |
| `H >= e^(W/2 + O(M))` | yes | no | no: forces a lattice coset (Thm 2.4(2)) |

The quantifier order matches round 9's (T') (`round9/README.md` §1 item 1): an argument may choose `A`
after `c_0`. `AG(M^A)` with `A = A(c_0)` is the cleanest input found here that is true for actual circles
and false for every generic fake. It is also strictly weaker than full lattice-ness. §7 states the
missing lemma.

## 3. Lattice geometry of the tangent range

All chords of a cluster point in directions inside the **tangent range**
`T = {theta + pi/2 : theta in [theta_lo, theta_hi]} mod pi`, an interval of length `D <= C N^(-1/4)`. The
chord of `(x, y)` is tangent to the circle at the mid-angle `(theta_x + theta_y)/2`. This section
extracts what the lattice can say about `T`. It first records what the monomial gaps (class `A_0`)
already give, and then what lattice lines and lattice directions add.

**Notation.** For a pair, `b_xy = a_x + a_y - e in Z^r` is the *bisector exponent*. The primitive
chord direction `u_xy` is `(z_x - z_y)/content(z_x - z_y)`, with `content(w) = gcd(Re w, Im w)`.

### 3.1 What the monomial gaps already give

**Proposition 3.1 (pair gcds are at least `N^(1/4)` off one matching; proved; class `A_0`).** Let an angle
model (in particular any normalised actual cluster) satisfy the monomial gaps of `round4/outside.md`
Lemma 0.1, with `D <= C N^(-1/4) < pi/4`. Then there is `b* in Z^r` such that every pair with
`b_xy != b*` has

```text
mu_xy <= W/4 + log C + (1/2) log 2,     equivalently  Norm g_xy >= N^(1/4)/(sqrt2 C).
```

The pairs with `b_xy = b*` form a matching.

*Proof.*
1. *Per column*, with `lo = min(a_xj, a_yj)` and `hi = max(a_xj, a_yj)` in `[0, e_j]`:
   `|lo + hi - e_j| <= e_j - hi + lo`. The two halves of this are `2 lo >= 0` and `2 hi <= 2 e_j`. Summing
   with weights `log p_j` gives `w(b_xy) <= W - d_xy = W/2 - mu_xy`.
2. Adding the realisation equations of `x` and `y`,
   `<b_xy, phi> = theta_0 + (delta_x + delta_y)/2 + (pi/4)(k_x + k_y)`. So for two pairs `P, P'`,
   `v = b_P - b_P'` has `dist(<v, phi>, (pi/4) Z) <= |delta_x + delta_y - delta_x' - delta_y'|/2 <= D`.
3. If `v != 0`, the gap inequality gives `dist >= arcsin(e^(-w(v)/2)/sqrt2) >= e^(-w(v)/2)/sqrt2`. Hence
   `w(v) >= W/2 - 2 log C - log 2`. With `w(v) <= w(b_P) + w(b_P') <= W - mu_P - mu_P'` this gives
   `mu_P + mu_P' <= W/2 + 2 log C + log 2`. So two pairs with different `b` cannot both exceed
   `W/4 + log C + (1/2) log 2`. Take `b*` to be the common value of `b` over the pairs that exceed it,
   or any value if there are none.
4. If `(x, y)` and `(x, y')` both have `b = b*`, then `a_y = a_y'`. Then `theta_y - theta_y' in (pi/2) Z`,
   and `|delta_y - delta_y'| <= D < pi/2` forces `theta_y = theta_y'`, i.e. `y = y'`. QED

*Remarks.*
1. This answers `round9/L1L2.md` open point 1 ("improve Lemma 1.3 to `Norm g_xy >= N^delta`") with
   `delta = 1/4`, for all pairs but one matching. That matching cannot be removed. Every complementary
   pair (`z_y = unit * conj z_x`, `mu = W/2`) has `b_xy = 0`, so all complementary pairs lie in one class,
   and the record clusters do have complementary pairs (§3.3).
2. Since it uses only monomial gaps, the fake circles satisfy it. It cannot be the source of growth.

### 3.2 Lattice lines plus lattice directions

**Definition.** For a primitive `v in Z[i]` let `eta(v)` be the angular distance from the arc to the line
`arg v + pi Z`. It is `0` if the arc meets that line, i.e. if `i v` lies in `T`. Put

```text
h(v) = |v| (1 + N^(1/4) eta(v)/C),      h* = min_v h(v),      v* a minimiser,
```

with `C = D N^(1/4)` the cluster's own constant. The minimum exists, since `h(v) >= |v|`. `h*` measures
how close the tangent range is to a lattice direction of small height. Call it the *effective height of
the arc's position*.

**Theorem 3.2 (effective height; proved; actual clusters only).** For every cluster of `M >= 3` lattice
points on an arc with `D < pi/2`, and every normalised cluster for (b):
* **(a)** `M <= C^2 h* + 2`.
* **(b)** Every pair whose chord is not parallel to `i v*` has `|u_xy| >= N^(1/4)/(C h*)`, hence
  `mu_xy <= 2 log(sqrt2 C h*)` and `|B_xy|/t_xy <= C h*`. If `eta(v*) > 0` there are no exceptional
  pairs. Otherwise the exceptional pairs form a matching.
* **(c)** `h* <= (2 pi/C)^(1/2) N^(1/8) + 1`.

*Proof.*
* (a) is `round4/outside.md` Prop C.2 applied to `v = v*`, with `L = RD`, `R D^2 = C^2` and
  `R D = C N^(1/4)`. If `eta = 0` the bound is `|v| C^2 + 2`. If `eta > 0` it is
  `|v| C N^(1/4) eta + |v| C^2/2 + 1 <= C^2 h(v) + 1`.
* (b) The chord of `(x, y)` is perpendicular to the radius at the mid-angle `theta_m`, which lies in the
  arc. So the angle between the lines `R u_xy` and `R i v*` equals `dist(theta_m, arg v* + pi Z)`, which is
  at most `eta + D`. If the two lines differ, `Im(u_xy conj(i v*))` is a nonzero rational integer, so
  `1 <= |u_xy| |v*| sin(angle) <= |u_xy| |v*| (eta + D)`. Since `eta + D = D (1 + N^(1/4) eta/C)` and
  `D = C N^(-1/4)`, this gives `|u_xy| >= N^(1/4)/(C h*)`.
  Next, `z_x - z_y = i^(k_x) g_xy gamma_s B_xy` with `B_xy` a rational integer (`L1L2.md` Thm 1.1), so
  `u_xy` is the primitive part of `g_xy gamma_s`. For `gamma_s = 2i` or `2` that is the primitive part of
  `g_xy`, and for `gamma_s = 1 +- i` its absolute value is at most `sqrt2 |g_xy|`. So `|u_xy| <= sqrt2 |g_xy|`
  and `Norm g_xy >= N^(1/2)/(2 C^2 h*^2)`, i.e. `mu_xy = W/2 - log Norm g_xy <= 2 log(sqrt2 C h*)`. The
  bound on `|B|/t` follows from the exact identity `mu = 2 log|B| + log nu_s - 2 log t` (`L1L2.md` (1.1))
  with `nu_s >= 2`.
  A chord parallel to `i v*` lies in `T`, so it forces `eta(v*) = 0`. Two chords from one point are never
  parallel, so the exceptional pairs form a matching.
* (c) Let `psi` be the mid-angle of the arc. Minkowski's theorem gives a nonzero lattice point in the
  symmetric box `|Re(w e^(-i psi))| <= Q`, `|Im(w e^(-i psi))| <= 1/Q` (area 4), and we may take it
  primitive. Its line makes an angle `alpha` with `psi`, where `sin alpha <= 1/(Q|v|)`, so
  `eta(v) <= alpha <= pi/(2 Q |v|)`. Hence `h(v) <= |v| + pi N^(1/4)/(2 C Q) <= Q + 1 + pi N^(1/4)/(2CQ)`.
  Take `Q = (pi N^(1/4)/(2C))^(1/2)`. QED

**Corollary 3.3 (polynomial reduced chords in low position; proved).** Let `b > 0`. Every normalised actual
cluster with `h* <= M^b` satisfies, for all pairs off one matching,

```text
mu_xy <= 2 b log M + log(2 C^2),    |B_xy| <= C^2 M^b,    t_xy >= 1/(C M^b).
```

So round 9's `(P_b)` (equivalently (R1-arch), i.e. L1) holds off one matching. The exceptional pairs are
chords parallel to the low-height direction `i v*`, and these are the only pairs that can carry a
super-polynomial `|B|/t`.

*Proof.* Theorem 3.2(b) with `t_xy <= C` and `|B_xy| >= 1`. QED

*Special cases.* Taking `v = i u_xy` (`eta = 0`) gives `h* <= |u_xy| <= sqrt2 |g_xy|`, so (a) contains
`L1L2.md` Lemma 1.3 (`M <= sqrt2 C^2 |g_xy| + 2`). With `v` arbitrary, (a) is `outside.md` Prop C.2
optimised over `v`.

**What Theorem 3.2 adds to Proposition 3.1, and what it does not.**
1. By (c), in generic position (`h* ~ N^(1/8)`) part (b) says `mu_xy <= W/4 + O(log C)`. That is
   Prop 3.1 again, except that the exceptional class is empty when `eta(v*) > 0`. Part (a) says
   `M = O(C^(3/2) N^(1/8))`, which is useless.
2. The gain is in **low position**: `h*` polynomial in `M`. There the margins off one matching are
   `O(log M)`, which is (R1-arch). By (a), `h* >= (M - 2)/C^2` always. So low position means that the
   tangent range sits within `~ N^(-1/4) M^b/|v|` of a lattice direction `v` of height `|v| <= M^b`.
3. Part (b) is an upper bound on margins. A growth proof needs a lower bound on `W`, i.e. large margins
   (`L1L2.md` §5.1: route 1 = (R1-arch) + (R1-arith)). For example the cut bound
   `sum_(x<y) mu_xy <= M W/4` is compatible with every conclusion of Theorem 3.2. Theorem 3.2 is a
   restricted (R1-arch). It is not a growth theorem.
4. *Inputs and fakes.* (a) uses lattice lines, i.e. integrality of `Re(z conj v)` for every `v`.
   (b) uses integrality of `Im(u conj v)` for one lattice direction `v` and the lattice chord directions.
   Both are false for fakes in general: a fake has no lattice, and its "chord directions" are real angles.
   But a fake whose tangent range is in generic position (for instance after a rotation, when its
   realisation coset allows one) satisfies every **conclusion** of Theorem 3.2 as soon as its margins are
   `<= W/4 - O(log C)`, which holds for the Hadamard-type fakes (margins `O(log M)`). So Theorem 3.2
   separates actual circles from fakes only in low position. There it gives (R1-arch), which those fakes
   satisfy anyway.

### 3.3 Exact check

`tangent.py` and `tangent_records.py` (outputs `tangent_out.txt`, `tangent_records_out.txt`) test
Prop 3.1 and Theorem 3.2(a)-(c) on every maximal run with `M >= 3`:
* 600 runs at `C <= 4` on 151 random split `N < 10^16` and the three record circles;
* all 592 runs of the three record circles at `C <= 8`, 5.3 and 5.9 (maximal `M = 9`, 7, 8 after
  normalisation).

The arithmetic is exact. Gcds, contents and parallelism use Python integers. `h*` is an exhaustive
minimum over primitive `|v| <= (2 pi/C)^(1/2) N^(1/8) + 3`, with angles to 60 digits.

There are no violations. Over the 600 random runs, `max(mu_xy - 2 log(sqrt2 C h*))` over
non-exceptional pairs is `-0.061`, so (b) is nearly sharp. On the 9-point record (`N = 1176852625`),
`h* = sqrt 2` with `v* = 1 - i` (the diagonal direction). Its 4 exceptional pairs are the chords
perpendicular to the diagonal (`record_detail.py`). Two of them are complementary (`mu = W/2`) and two have `Norm g = 25`. All other pairs have
`mu <= 2 log(sqrt2 C h*) - 2.08` with `h* = sqrt2`. This is the symmetric low-position mechanism that
`L1L2.md` Lemma 1.3 confines to `M <= sqrt2 C^2 + 2`. Theorem 3.2(a) is its exact generalisation to
every direction and every position.

## 4. p-adic structure of exact relations

Section 2 shows that the existence of `K` is the non-local input. This section records exactly what one
relation says at the primes of `N`. It extends `round9/L1L2.md` Lemma 5.3.

**Lemma 4.1 (Newton polygon of a relation; proved).** Let `c in K`, `S = supp(c)`, and let `j` be a column.
Write `v_x = v_(p_j)(c_x)`. Then `min_(x in S) (v_x + a_xj)` is attained at least twice, and
`max_(x in S) (a_xj - v_x)` is attained at least twice.

*Proof.* `v_(pi_j)(c_x z_x) = v_x + a_xj`, since `v_(pi_j)(c_x) = v_(p_j)(c_x)` for a rational integer. In
a vanishing sum the minimal valuation occurs at least twice. For `conj(pi_j)` the valuation is
`v_x + e_j - a_xj`. QED

**Lemma 4.2 (leading residues; proved).** In the setting of Lemma 4.1 let `m` be the minimum and
`S_m ⊂ S` the set where it is attained. Put `rho_x = z_x pi_j^(-a_xj) mod pi_j`, a nonzero element of
`Z[i]/pi_j = F_(p_j)`, and `c'_x = c_x p_j^(-v_x)`. Then

```text
sum_(x in S_m) c'_x conj(pi_j)^(v_x) rho_x = 0   in F_(p_j),
```

and symmetrically at the maximum with `conj(pi_j)`.

*Proof.* Divide `sum_(x in S) c_x z_x = 0` by `pi_j^m` and reduce modulo `pi_j`. The terms outside
`S_m` vanish. For `x in S_m`, `c_x z_x/pi_j^m = c'_x conj(pi_j)^(v_x) (z_x pi_j^(-a_xj))`. QED

The `rho_x` are the level residues at `p_j` (P3 level data at `q = p_j`). When `|S_m| = 2`, Lemma 4.2
says that a ratio of two level residues is a ratio of small integers in `F_(p_j)`.

**Proposition 4.3 (isolation; proved).** Let `lambda_k^inf(K)` be the successive minima of `K` in the sup
norm, and let column `j` have a unique row `x*` with `a_(x* j) = e_j` (or a unique row with `a = 0`). Then
`p_j <= lambda_(M-3)^inf(K)`.

*Proof.* Suppose `p_j > lambda_(M-3)^inf`. Take `M - 3` independent relations with entries of absolute
value `< p_j`, so every `v_x = 0` on their supports. If `x*` lay in a support, the maximum of `a_xj` over
that support would be `e_j`, attained only at `x*`, contradicting Lemma 4.1. So all `M - 3` relations
vanish at `x*`. But `K ∩ {c_(x*) = 0}` has rank `M - 4` (relations among the other `M - 1` points). QED

**Remark 4.4 (the four-point case is the triangle identity).** `K` contains the `M - 3` independent
four-point relations `(A_yzw, -A_xzw, A_xyw, -A_xyz)` (`A = 2 Area`), with entries
`<= L^3/(2R) <= C^3 N^(1/4)/2`. So Prop 4.3 applies whenever `p_j > C^3 N^(1/4)/2`. In fact more is
true. Lemma 4.1 applied to all four-point relations forces, for such a column, every four rows to show the
pattern `{a, a, b, b}` or `{a, a, a, a}`. For `M >= 5` this makes the column constant, which is
impossible after normalisation. So **every prime of a normalised cluster with `M >= 5` satisfies
`p_j <= C^3 N^(1/4)/2`.** This is the triangle identity again (`round4/outside.md` Prop C.1): `p_j`
divides some `2 Area` once some triple has range `< e_j` in column `j`. It is profile data, and the fakes
satisfy it. This matches `L1L2.md` §5 (four-point relations reduce to the profile).

**What short relations of large support would add.** By Theorem 2.4(1) and `L1L2.md` Prop 5.2, a
near-extremal cluster has a relation of height `M^(O(c_0))`. By Minkowski's second theorem it has `M - 3`
independent relations whose heights have geometric mean `M^(c_0/2 + 1/2 + o(1))`. Their supports are not
controlled. For a column with `p_j` above that height, Lemmas 4.1-4.2 force on
every such relation:
* doubled extremes on its support;
* one `F_(p_j)`-linear relation with small coefficients among the level residues at each extreme level.

For **random** level data each such relation holds with probability about `1/p_j`. Over a 0/1 profile
the product is `prod_j p_j^(-2) = N^(-2)`, against about `H^M` candidate vectors. That is
`exp(M log H - 2W) << 1` at `W = c_0 M log M` (heuristic). So random level data almost never admit even
one relation with doubled extremes and the residue cancellations everywhere.

This is not a mechanism. Actual level data are not random: by Theorem 2.6 they are exactly the data
with small lifts, and the relations are then automatic. A fake has no relations at all, so it satisfies
every statement of the form "each relation of height `<= H` obeys ...". A mechanism of this type would
have to **derive a contradiction from the existence of the relations**. That existence is Theorem
2.4(1) / Prop 5.2, i.e. lattice-ness at height `M^(O(c_0))`. §6 discusses why no such derivation is known.

## 5. Other additive tools: sums, Sidon sets, polynomial method mod q, determinants

The task lists further additive tools. Each one reduces, by the statements below, to a class already on
the barrier map.

1. **Sums.** `z_x + z_y = i^(k_x) g_xy (A + i^s conj A)` is again `g_xy` times a small Gaussian factor times
   a rational integer `B'_xy` (the complementary `D_s` of `L1L2.md` Thm 1.1), and
   `|z_x + z_y|^2 = 4N - |z_x - z_y|^2`. So sums carry the same data as chords: pair data, class `A_0`.
   Midpoint coincidences `z_x + z_y = z_u + z_w` force `{x, y} = {u, w}`, because two chords of a circle
   with the same midpoint `m != 0` coincide. So the cluster is additively Sidon (`round4/outside.md`
   §3), a 4-point statement.
2. **Sidon and sum-free properties of the `beta`'s.** The cocycle `g_xw beta_xw = g_xy beta_xy +
   g_yw beta_yw` is the 3-point identity. At odd `q` not dividing `N` its only content is transitivity of
   `q^a | B` (`L1L2.md` Prop 5.1(a)), i.e. point-residue partitions (class `A_X`). Chord lengths are
   pairwise distinct, so unit-distance and Sidon-type counts are trivial (`outside.md` Remark C.3).
3. **Integrality of `Im(z_x conj z_y)`.** `sigma_xy = Im(z_x conj z_y) = N sin(theta_x - theta_y)` is an
   integer, and `sigma_xy + sigma_yw + sigma_wx = 2 Area` (signed). So the additive defect of the integer
   cocycle `sigma` is the triangle area, which is profile data (Prop C.1). Longer cycles give sums of
   triangle areas.
4. **Polynomial method over `F_q` (reduction modulo `q`).** Let `q` be odd, not dividing `N`, and suppose the
   cluster has at least 3 residues modulo `q`. Then the `F_q`-space of affine relations among the reduced
   points equals `K ⊗ F_q`. *Proof.* `K` is saturated in `Lambda_0` (it is the kernel of a map to a
   torsion-free group), so `K/qK -> Lambda_0/q Lambda_0` is injective with image of dimension `M - 3`.
   Three distinct points of the nondegenerate conic `x^2 + y^2 = N` over `F_q` are not collinear, so the
   reduced points span the affine plane and their relation space also has dimension `M - 3`. It contains
   the image. QED
   So reduction modulo `q` sees exactly `K` and the residue partition. Any statement about the reductions
   alone, Stepanov-type counting included, is a statement about residue data (class `A_X`), because the
   arc condition does not survive reduction. Coupling the reduction to the archimedean place is the
   `q`-adic determinant method. By `growth/auxiliary.md` Prop A and §2.2 it is dominated by the pair data
   (Vandermonde times a symmetric factor, `q`-adic gain = collision counting).
5. **Bombieri-Pila with arithmetic input.** Interpolation determinants on the circle factor as Vandermonde
   times a symmetric Gaussian integer (`auxiliary.md` Prop A). Single-place auxiliary forms are capped at
   `tau*(M) < 2` (`threshold/upper.md`). Forced divisibility at the primes of `N` is exactly critical
   (`threshold/construct.md`). Low-height lattice conics through many points are auxiliary polynomials.
   Lattice **lines** are the degree-1 case, which §3 uses in full: Theorem 3.2(a) is the best line
   statement, and it is uniform only in low position.
6. **The sublattice covolume.** "The points span a sublattice whose covolume is forced small by the short
   arc" holds only in this form: `I | 2 Area` for every triangle, and `2 Area <= C^3 N^(1/4)/2`. The
   covolume of `K` is pinned by the triangle areas and `I` (Prop 2.7(3)). `I` is residue data and the areas
   are profile data. The short arc does **not** force `I` to be large. Prop 2.7(2) says that the odd part of
   `I` is supported on primes where the cluster has at most 2 residues.

## 6. Barrier map for additive inputs

Columns:
* "actual": the input holds for every actual cluster;
* "fakes": whether the fake circles of rounds 4-7 (or generic angle models) satisfy it;
* "strength": the best consequence known.

"Capped" means that the class contains fakes with `W = O(M log M)`, so the input alone cannot give growth.

| additive input | actual | fakes | strength | where |
|---|---|---|---|---|
| sums, midpoints, Sidon, chord lengths, `sigma_xy` cocycle, triangle areas, 4-point relations | yes | yes (profile/pair identities) | capped (`A_0`) | §5.1-5.3, `L1L2.md` §5, `outside.md` C.1 |
| index `I` of the difference lattice; reductions mod `q`; `F_q` polynomial method | yes | as residue data | class `A_X` (capped at each fixed `A`) | Prop 2.7, §5.4 |
| `q`-adic and global determinant methods, Bombieri-Pila | yes | n/a | dominated by pair data / `tau* < 2` | §5.5, `auxiliary.md`, `threshold/` |
| pair gcds `>= N^(1/4)` off one matching | yes | **yes** (monomial gaps) | capped (`A_0`) | Prop 3.1 |
| lattice lines (Prop C.2), Thm 3.2(a) | yes | not defined for fakes (no lattice); the conclusion fails only for fakes in low position | uniform bound iff `h* = O(1)`; `M <= C^2 h* + 2` | Thm 3.2(a) |
| lattice directions (`Im(u conj v) in Z`), Thm 3.2(b) | yes | not defined for fakes; the conclusions hold for small-margin fakes in generic position | (R1-arch) off one matching in low position `h* <= M^b` | Thm 3.2(b), Cor 3.3 |
| `AG(M^A)`, `A < c_0/4` | yes | yes for typical fakes (heuristic) | none | Heur. 2.9 |
| `AG(M^A)`, `A >= c_0/2 + o(1)` | yes | **no** for almost every relation-free fake (proved) | rank bound Prop 2.8; no bound on `M` known | Thm 2.4(1), Cor 2.5 |
| `AG(M^A)`, `A > c_0` | yes | **no** for any fake with independently imposed relations (heuristic) | open lemma (AG) of §7 | Prop 2.8, Heur. 2.9 |
| full relation lattice (`AG` at height `e^(W/2+O(M))`) | yes | no, unless the fake is a lattice-coset configuration | = norm-form problem for an arbitrary imaginary quadratic order (Thm 2.3) | Thms 2.1-2.4 |
| relations of an actual cluster | yes | no, unless the fake is that cluster rotated | = the cluster | Thm 2.2 |
| residues mod `N` + small lifts | yes | no, unless actual | = the actual problem | Thm 2.6 |

**Barrier check of the proposed mechanisms** (the requirement of the task: name an input that is true for
actual circles and provably false for the fakes).
* *Prop 3.1* uses monomial gaps only. It holds for fakes and cannot give growth. This is recorded so that
  "`Norm g_xy >= N^(1/4)` off a matching" is not mistaken for a new handle.
* *Theorem 3.2* uses lattice lines and lattice directions, which are false for fakes. But its conclusions
  are upper bounds on margins and on `M` in terms of `h*`. Fakes in generic position meet the conclusions,
  and near-extremal fakes have small margins anyway. So it separates actual circles from fakes only
  through (a) in low position, which is the known Lemma 1.3 mechanism.
* *`AG(M^A)` for `A >= c_0/2`* is provably false for almost every fake of a positive-measure family
  (Cor 2.5). This is the one input of this round that separates actual circles from fakes in the required
  sense. It is **not** known to bound `M`. Converting it into a bound needs the open lemma (AG) below.
* *The full relation lattice* separates, but using all of it is the norm-form problem for an arbitrary
  imaginary quadratic order (Thm 2.3). An argument that cannot tell `Z[i]` from another order must prove
  growth for all orders at once. The hypotheses bound the covolume of the lattice
  (`covol <= 2 Area <= C^3 N^(1/4)/2`) but not the discriminant of the primitive integral form.

**What no additive argument can avoid.** Every additive input above is either local (capped), or a
consequence of the existence of the relation lattice at some height. The existence of relations is
lattice-ness, and lattice-ness is the problem itself (Thms 2.3, 2.6). No intermediate additive class
was found whose inputs are strictly between "local" and "lattice coset", except the height-graded family
`AG(M^A)` of §2.5.

## 7. Open points

1. **Lemma (AG) (open; the missing step of the only separating mechanism found).**

   > For every `c_0` and `C` there is `A = A(c_0)` such that, for `M >= M_0(c_0, C)`, no angle model
   > (genuine primes, profile and windings, gap-admissible angles, positions in `[0, D]^M` with
   > `D <= C N^(-1/4)`) with `W <= c_0 M log M` satisfies `AG(M^A)`.

   * *(AG) implies growth (proved).* A normalised actual cluster is an angle model satisfying `AG(H)` for
     every `H`, and normalisation does not increase `W` or `C` (`two_thirds` Lemma 1.1). So (AG) gives
     `W > c_0 M log M` for every actual cluster with `M >= M_0(c_0, C)`. As `c_0` is arbitrary,
     `M = o(log R/loglog R)`.
   * *What is known.* `AG(M^A)` forces module rank `rho <= 2 + (c_0/(2A)) M + o(M)` (Prop 2.8). It fails
     for almost every relation-free angle model once `A >= c_0/2 + o(1)` (Cor 2.5). Imposing relations
     independently cannot reach `A > c_0` (Heuristic 2.9).
   * *What is missing.* A theorem that a module of rank `rho = o(M)` (or `rho <= eps M`) carrying `M`
     points of one circle on an arc of length `C N^(1/4)`, with the Gaussian profile structure, forces
     `W >= f(M) M log M` with `f -> infinity`. The extreme case `rho = 2` is the full problem (Thm 2.3).
   * *Caveat (it may be false for a reason unrelated to `Z[i]`).* Lattice-coset configurations for other
     imaginary quadratic orders satisfy `AG(infinity)`. If such a configuration also carried a Gaussian
     profile with gap-admissible angles at `W = O(M log M)`, (AG) would fail. So (AG) contains a
     statement about norm-form problems in all orders. No such hybrid is known, and none was looked for.
2. **A finite test of the ladder (proposed).** Take the round-5 skew double core fakes (`exact.md`) at the
   largest exactly solved `M`. Decide:
   * whether their admissible position set has interior;
   * the largest `A` with `AG(M^A)` attainable inside it, by LLL on the lattice
     `{(c, lambda Re sum c w, lambda Im sum c w)}`.

   If imposed relations reach `A > c_0` there, Heuristic 2.9 is wrong and the ladder must be redrawn.
3. **The modular relaxation `MOD_N` (open).** Drop the circle and keep only "Gaussian integers in a box of
   side `L = C N^(1/4)` whose norms are divisible by `N`". This class contains every actual cluster.
   * A Gaussian integer with `N | Norm z` lies in one of the `d(N)` ideals of norm `N`, all principal.
     Each is a square lattice of side `R`, and `R > sqrt2 L`, so the box meets each ideal at most once.
     Hence `M <= d(N)`, the divisor bound.
   * Whether `MOD_N` has members with `M >> W/log W` at `W = O(M log M)` is open. Equivalently: is
     curvature needed at all?
   * For a randomly placed box, the expected number of ideals met is about `d(N) L^2/N`, far below 1
     (heuristic).
4. **Theorem 3.2 in generic position.** In generic position `h* ~ N^(1/8)`, Theorem 3.2 is no stronger
   than the monomial gaps. Since fakes in generic position meet all its conclusions, no argument built
   from lattice lines and finitely many lattice directions can give (R1-arch) there. A lattice-geometric
   improvement must use the exactness of the circle (`|z|^2 = N`), not just lines.
5. **Relations of large support (Lemmas 4.1-4.2).** On near-extremal clusters, do the relations of
   height `M^(O(c_0))` have supports of size `Theta(M)`, and is their level-residue cancellation at the
   large primes compatible with genuine labels (P3-G) beyond what Theorem 2.6 makes automatic? Any
   contradiction must come from the existence of the relations (§4), i.e. through (AG).

## 8. Files

All scripts are in `additive_checks/` and run from that directory with `python3` (`mpmath` needed).
`glib.py` is a copy of `round9/L1L2_checks/glib.py`.

| script | checks | output |
|---|---|---|
| `tangent.py` | Prop 3.1 and Thm 3.2(a)-(c) on 600 runs (151 random split `N < 10^16` plus records, `C <= 4`) | `tangent_out.txt` |
| `tangent_records.py` | the same on all 592 runs of the three record circles (`C <= 8, 5.3, 5.9`) | `tangent_records_out.txt` |
| `record_detail.py` | exceptional pairs of the 9-point record | `record_detail_out.txt` |
| `collapse.py` | `rank K = M - 3` (exact kernel), Prop 2.7(1)-(3), Thm 2.6(1) via LLL on `K_N`, small-lift identity, on 468 clusters with `M >= 4` | `collapse_out.txt` |

`additive_notes.md` holds the working notes.
