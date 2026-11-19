# Referee report ref3_four: astra_lemma_III, sections III.1 to III.4 (and the III.9 tail they cite)

Scope: `astra_lemma_III/answer.md` lines 585-698: III.1 (corrections to Part I), III.2 (the signed
fixed-cofactor fourfold), III.3 (the arithmetic test and the empty-closure example), III.4 (moving planes, the
four-to-five count, fifteen cuts, specified singletons, the area formula and the bridge to 3-spaces). Summary items
10 and 11 also cite two results that are **not in III.1-III.4**: the interpolation degree bound is in III.9
(line 905-909) and the genus-four curve is in III.9 (lines 935-941). I refereed both here because the block names
them.

Scripts are in `ref3_four/`. They are my own; I did not rerun Astra's.

| script | what it checks | result |
|---|---|---|
| `common.py` | my real-sample construction (angles, then `log D = M_adj log n`), from the star incidence only | used by the others |
| `dim_check2.py` | analytic Jacobian of the incidence `|P_i|^2=S_i(n)`, `[P_i,P_j]=t_ij n_ij` (20 variables, 15 equations), 12 random real samples with mixed signs | rank 15 every time, so the incidence has dimension 5; the tangent map to `n` has rank exactly 4 (the fifth singular value is 0: the common rotation); the tangent map to the four angle differences has rank 4 (`dim_check.py`) |
| `dim_check.py` | at 6 samples: `Mc=0`, rank of the 5 linear forms, rank of `M`, all 30 triangle sextics | `|Mc|<1e-14`; 4 independent linear forms; rank `M` = 2; sextic residuals at most 1e-10 relative |
| `slice_exact.py` | exact over Q: a random affine 2-plane inside the 6-dim linear space `L0`; five Pfaffians, one triangle sextic, resultant | the 5 Pfaffians are constant multiples of one quadric on `L0`; the resultant has degree **12**, is squarefree, with no common root at infinity; all **12** points satisfy all 30 triangle identities (residual at most 1e-34) and have no zero coordinate. Three cofactor vectors, including the mixed-sign `t=(-4,5,3,1,-4,-3,5,3,-2,-1)` (`slice_out.txt`) |
| `monodromy.py` | degree from the star incidence alone: 9 equations in `(x_1..x_5,y_1..y_4)`, `y_5=0`, with 4 random complex affine equations on `n`; monodromy loops | see §2.3 |
| `param_check.py` | formula (III.2) `B^3=∏q/∏K`, `n_ij=d_ij q_i q_j/(K_iK_jB)` | relative error at most 7e-15 |
| `genuine_checks.py` | exact Gaussian arithmetic on 40 genuine pure ten-pair tuples (private odd split-prime blocks, random units and orientations) | parity lemma of III.3: 40/40; `∏_{i<j}(q_i-q_j)/∏q_i^2 = ±1024 τ/H^3`: 40/40; the true fifth mark is a root of the nonzero quartic of the four-to-five count: 40/40 |

## 1. III.1 Corrections to Part I: CORRECT (bookkeeping)

I compared each sentence with the round-31 Part I judge report and with ref_cover, ref_baker, ref_square and
ref_local, as listed in `reply3.txt`. Astra adopts all of them: the twisted toric cover with the same fan, real
torsion fibres for unit translates, the I-sheeted covering of `(S^1)^4`, uniformity for `I>=7` through (2) of §2.2,
Baker–Wüstholz checked only through a secondary source, the squareclass density fix, the 16 distinct projected
values, and withdrawal of the sentence on "different local obstructions". One new sentence, "neither additional
projected value can coincide with one of those fourteen ... distinct by their nonzero paired imaginary parameter",
supplies the distinctness that ref_local asked for. I accept it as a sketch; I did not recheck §5 itself, which
belongs to a different block. Nothing here is new mathematics.

## 2. III.2 The signed fourfold

### 2.1 Four linear equations, the Plücker quadric, one sextic suffices: CORRECT

Line by line:
* `Mc=0` with `c_i=1/∏_{j≠i}t_ij` (antisymmetric extension). From the quintic interpolation identity
  `Σ F(P_i)/∏_{j≠i}[P_i,P_j]=0` for binary cubics `F` and five points, with `F=|P|^2 P`, one gets
  `Σ_i c_i P_i = 0`. Then row `i` of `Mc` is `[P_i, Σ_j c_j P_j]=0`. Checked numerically. The kernel argument
  (bivectors annihilated by `c` form `Λ^2(c^⊥)`, which has dimension 6, and rank 2 there is one Plücker quadric) is
  standard linear algebra. Exact check: on `L0` all five 4x4 Pfaffians are constant multiples of one quadric, as
  `M·pf=0` forces.
* The converse ("one sextic suffices"). Factor `m_ij=[V_i,V_j]`. `Σc_iV_i=0` follows from `Mc=0` because the `V_i`
  span. The values `S_i` satisfy the two interpolation relations, so a unique form `Q` with `Q(V_i)=S_i` exists.
  Squaring the three-term Cramer identity gives `Q(V_i,V_j)=n_ij s_ij` with `s_ij=A_ijk/(2t_ik t_jk)`. I rederived
  `A_ijk` and got the same cubic. The Gram identity gives `det Q · t_ij^2 = D_ij - s_ij^2`. So the one triangle
  sextic is exactly `det Q=1`. Over R, `S_i>0` makes `Q` positive definite, and an `SL_2` change gives Euclidean
  stars; the incidence matrix of `K_5` has rank 5, so the edge phases exist. All steps hold. The complex version on
  `∏n≠0` is the same argument with `SO(2,C)`.
* Exact confirmation: on a generic affine slice, the 12 points cut out by the linear equations, the quadric and one
  sextic satisfy **all** 30 triangle identities (`slice_exact.py`).

### 2.2 Dimension four: CORRECT (independently reproduced)

The incidence has 20 variables and 15 equations. The analytic Jacobian has rank 15 at 12 random real samples with
mixed sign patterns, so its dimension is 5. The projection to `n` has rank exactly 4, and the remaining kernel
direction is the common rotation. So `V_t` has dimension 4. This confirms the coordinator's "about 4" count. It is
also the expected count: 10 star coordinates minus rotation, minus 5 star-norm equations, minus 0 for `n`, since
`n` is determined by the stars.

### 2.3 Degree exactly twelve and irreducibility: CORRECT

* **Upper bound and exact value on a slice.** A generic affine 6-plane meets `V_t` inside `L0∩(4 hyperplanes)`, a
  2-plane. There the variety is cut by a conic and a sextic curve, so Bézout gives at most 12. `slice_exact.py` finds
  exactly 12 distinct finite points over Q-random slices for three cofactor vectors: the resultant is squarefree of
  degree 12, with no common point at infinity and no zero coordinate. All 12 lie on `V_t`. This also rules out an
  extra 4-dimensional component inside a coordinate hyperplane or at infinity, which was the delicate boundary step
  in Astra's Bézout paragraph.
* **Irreducibility.** It follows from Astra's parametrization (III.2). I verified the formulas numerically
  (`param_check.py`). The cubic `B^3=∏q/∏K` has order −2 along each diagonal, so it is not a cube and the cover is
  irreducible of degree 3. The map from the cover to `n` is injective because `n` recovers `q` (through
  `(s+it)/(s-it)`) and then `B`.
* **Monodromy from the star incidence alone** (`monodromy.py`, `mono_audit.py`, `mono1.log`). This check does not
  use any of Astra's equations.
  * Setup: seed 1, `t=(5,-3,2,-5,2,-3,4,-3,1,1)`, 4 random complex affine equations on `n`, 80 monodromy loops
    started from one real sample.
  * Result: 110 star solutions. Of the 55 distinct `n` among them, exactly **12** satisfy `Mc=0` and rank `M`=2
    (residual at most 2e-16) and have no zero coordinate.
  * The other 43 lie on spurious components of the square star system inside `{n_e=0}`, where isotropic star
    vectors make the equations degenerate. The tracker reached them by path-jumping near those components; they are
    not points of `V_t`.
  * So the main component meets the generic slice in 12 points, all reached by monodromy from a single start. This
    supports irreducibility independently of (III.2).
  * Caveat: this is floating point and not certified. A run with a second seed hit its time limit and gave no
    result.

So the closure of the continuous norm vectors is an irreducible fourfold of degree 12. Real samples form an open
4-dimensional real set and are therefore Zariski dense.

### 2.4 Real dominance onto ratio space and "elimination gives zero": CORRECT

The angle construction is the one I used to generate samples, and every sample is a genuine continuous model. Rank 4
of the map to the four angle differences was checked. The consequence is the important one: **for a single fixed
`t` the ratio image is Zariski dense**, so no nonzero polynomial in the ratios vanishes on `V_t`. The fixed-`t`
closure is therefore not a proper subvariety in the coordinates where Schwartz–Zippel lives. This refutes the
coordinator's proposed route in `reply3.txt` item 2. It is not a mere failure of a proof: any polynomial in the
ratios that vanishes on all real-feasible points of one `t` is zero. The remark that an additional norm equation `G`
would transfer through a field norm down the cubic cover, with degree `O(deg G)`, is correct and elementary.

## 3. III.3 Arithmetic on the fixed-cofactor fourfold

### 3.1 Exact realization test (III.3): CORRECT

* Necessity holds: genuine reduced words are primitive, and with each prime in only one orientation a common rational
  prime divisor of `t` and `D` would divide `Z`.
* Sufficiency, step by step:
  * `s_ij` is rational by (II.6) and `s^2=D-t^2` is an integer, so `s` is an integer.
  * `gcd(t,D)=1` and odd `D` make `Z` primitive.
  * Inert primes are excluded.
  * The π-adic valuations of `q_j/q_i` form a coboundary that vanishes on non-crossing pairs and equals `±k` on the
    six crossing pairs.
  * The unit parity count holds: the unsigned incidence lattice is the even-sum lattice, and the common factor `i`
    flips parity because there are five stars.
* Uniqueness of the ratios is immediate.

Caveat: the hypothesis "private, each block coprime to its conjugate, pairwise coprime odd norms" is essential, and
the test says nothing once prime layers overlap. Astra says this in III.9.

### 3.2 Empty arithmetic closure for `t≡1`: CORRECT

Odd `|P_i|^2` forces `P_i ≡ 1` or `i (mod 2)`, so `t_ij ≡ [P_i,P_j] (mod 2)` is odd exactly across a bipartition,
and a triangle cannot have three crossing edges. `genuine_checks.py` confirms the parity rule on 40/40 genuine
tuples. In fact the conclusion is stronger: no genuine odd-norm tuple has odd `t` on any triangle. The real locus for
`t≡1` is full-dimensional (increasing angles in a semicircle). So the real closure and the arithmetic closure really
differ. The inert-prime generalization (`p|t_ij` means equal directions mod `p`) is also correct.

### 3.3 "Ambient identities also hold for angle models": CORRECT, a remark

This confirms the coordinator's first worry.

## 4. III.4 What survives when the plane or the parent moves

### 4.1 Uniform `O_ε((2T)^ε)` count on an arbitrary rational plane: CORRECT

This is the main claim of the block. I checked every case.

* **dim 1.** Pairwise coprimality forces the multiplier to be 1, so there is at most one sample.
* **dim 2, basis.** Points are `(p,q)` in a basis of `U∩Z^10`. `gcd(p,q)=1` holds, because any common factor
  divides all ten norms.
* **Proportional forms.** Two coprime positive values fixed by a rational ratio confine the samples to one affine
  line. The triangle sextic for a pair separating a nonconstant norm is a polynomial of degree ≤6 in the line
  parameter. It is not identically zero, because at the real zero of that norm `D_ij=0`, so `LHS=-A^2≤0<RHS`.
  Hence at most 6 samples. (Trivial subcase that Astra skips: if all norms are constant on the line there is one
  sample.)
* **Thue case.** `F(p,q)=h=4(t_ij t_ik t_jk)^2`, and each primitive irreducible factor takes a divisor value. For a
  factor of degree `d≥3`, Bombieri–Schmidt gives at most `c·d^{1+ω(b)}` primitive solutions, with `c` absolute,
  for an irreducible integral binary form of degree ≥3. I verified the statement through secondary sources: the
  MaRDI entry for Invent. Math. 88 (1987) 69–81 and later papers that quote it (arXiv 1508.03602, 2206.13653). I
  have not read the original text. Summing over divisors gives `O_ε(T^ε)`, independent of the coefficients. This is
  exactly what makes the bound uniform in the plane: Bombieri–Schmidt has no dependence on the form's height or
  discriminant. That answers the block's question "what theorem gives uniformity in the plane?".
* **Two nonproportional factors of degree ≤2.** The form `b_2^{d_1}G_1^{d_2}-b_1^{d_2}G_2^{d_1}` is nonzero by
  unique factorization and has degree ≤4. That gives at most 4 slopes, each with 2 primitive representatives.
* **`cL^6`.** At most two lines `L=±β`. On each, some `L_e` is nonconstant, since otherwise all forms are
  proportional to `L`. The real-zero argument then gives at most 6 per line.
* **`cQ^3`.** If some necessary equation fails on the conic, Bézout gives at most 12 samples. Otherwise the III.2
  converse realizes real stars on the positive arc through the sample, and the identities extend to the whole
  conic.
  * Ellipse: every line through the origin meets it, and at the zero of `L_e` the identity `D_ij=s^2+t^2` gives
    `0≥t^2>0`. This uses that `s_ij` is a real polynomial in `n`, which it is.
  * Hyperbola: I checked that the II.6 argument needs only one sample. The positivity `a_e,b_e>0` comes from the
    no-real-zero property of every `n_e` on the branch, which (II.10) gives once the identities extend, and not from
    samples at both ends. Irrational asymptotes give `a_eb_e≠0`, and pairwise nonproportionality gives distinct
    `b_e/a_e`. The cocycle and the five-star lemma are real arguments. ref2_plane found II.6 CORRECT.

  So containment is impossible and at most 12 samples remain.

The cases are exhaustive. The constant depends only on ε. **Verdict: CORRECT.** Two small corrections:
* The bound is `O_ε(T^{ε})` with an absolute dependence on the degree ≤6 hidden in the constant, as stated.
* The phrase "ordered ratio classes up to an absolute factor" should read "at most the number of norm vectors"
  (uniqueness in III.3).

The standing cautions are respected. Nothing uses a fixed-S Subspace theorem. The Thue bound is uniform in the
primes and the coefficients, and the plane may move.

### 4.2 Specified singleton norms: CORRECT

* With `S_i=a_i∏n_ij`, `c'_i=a_ic_i`, `D_ij=a_ia_j∏...` and `A_ijk` multiplied termwise by `a_i,a_j,a_k`, the
  right side stays `4(t_ij t_ik t_jk)^2`.
* Bombieri–Schmidt is coefficient-uniform.
* The real-zero, ellipse and hyperbola arguments only use `a_ia_j>0`.

So the count is `O_ε(T^ε)` per *specified* singleton vector and per plane. As Astra says, the number of choices is
not controlled.

### 4.3 `O_ε(K L^{1+ε})` over K planes, and the interpolation degree (III.9 line 905): CORRECT, trivial

The count of signed 10-tuples with `∏|t|≤L` is `2^{10}Σ_{N≤L}τ_{10}(N)=O(L^{1+ε})`, and each contributes
`O(L^ε)`. Interpolation needs `binom(d+4,4)>s`, so `d=O(1+s^{1/4})`, which gives
`O_ε(1+(AK)^{1/4}L^{(1+ε)/4})`. The count of rational 2-planes of covolume ≤P by `O(P^{10}log P)` is a crude but
valid upper bound; Schmidt's asymptotic is `~cP^{10}`. The condition `α+10κ+ℓ<4` for degree `o(m)` is correct
arithmetic.

This is a restatement of "few points imply a low-degree interpolating polynomial". All its content is in the
hypotheses on `A,K,L`, none of which the growth problem supplies. In the growth setting `L` is at least `m^{c}` for
an unspecified `c` (cofactor log-size `O(log m)`), so `ℓ` is not known to be <4.

### 4.4 Four-to-five count: CORRECT (for the pure private ten-pair model)

* `H_5=H_4`: every 2|3 cut stays nonconstant after deleting a vertex.
* `J=c∏(q_i-q_j)/∏q_i^2=τ/H^3` with `|c|=1/1024` (my derivation:
  `q_j-q_i=2it_ijn_ij/(\bar P_i\bar P_j)`), confirmed exactly on 40 genuine tuples.
* As a function of the fifth mark the equation is a nonzero quartic: the true fifth mark is a root in 40/40 cases.
* The six internal `|t|` are determined by the four marks, because a primitive Gaussian `A` with given `A/\bar A` is
  unique up to ±1.

The count `8⌊L_5⌋` and the min-vertex bound `L^{2/5}`, giving `O(m^4L^{2/5})`, are correct. The phrase "at most four
completions per fixed cofactor product" in summary item 11 is correct as stated: fixed τ, ≤4.

Significance check, which Astra makes himself: the Walsh argument of II.1 (`E f_J^2≤2L/m` for `|J|=2,4`, then
Chebyshev) shows that **pure ten-pair tuples are O_{C,ε}(m^4) on any growing grid**. I checked that
derivation: singleton patterns should carry weight ≈1/16 each, and pure profiles carry 0. So everything in III.2,
III.3, 4.1 and 4.4 about the *pure* ten-pair model concerns a set that the growth problem already discards.

### 4.5 Fifteen cuts (III.4 formula and the doubleton cone): CORRECT

* Each singleton enters `∏D_ij` four times and each doubleton six times, and
  `|∏(q_i-q_j)|=2^{10}|τ|/(σ^2∏n^3)=2^{10}|τ|σ/H^3` with `H=σ∏n`. Checked by hand.
* `H_5=H_4a_5`.
* The doubleton closure is the 7-dimensional cone over `Gr(2,5)` (rescaled). The free singletons absorb the
  kernel-vector condition, and the cubic identity is automatic. Dimension 9 with the `a_i`.

This is a correct and useful negative statement: in the realistic fifteen-cut setting even the fixed-`t` sextic
disappears.

### 4.6 Conditional ratio polynomial, area formula, bridge to 3-spaces

* **Conditional ratio polynomial: CORRECT, trivial.** Its hypothesis is a list of `o(W/log W)` values of
  `|τ|σ/H^3`, which is not available.
* **Area formula `|det| = abc/(2√N) = 4(N/H_klm)|t t t|`: CORRECT.** Circumradius formula, then
  `|z_k-z_l|=2√N|t_kl|/√D_kl` and `√(D_klD_lmD_mk)=H_klm`. In the pure model a split triangle has exactly two
  crossing pairs for each split cut, and singletons appear twice. Astra himself says it transfers nothing new.
* **Bridge to 3-spaces: CORRECT for fixed `U,t`.**
  * The Veronese parametrization has a bounded gcd.
  * The degree-12 binary form reduces to `L^{12}` or an indefinite `Q^6`.
  * The splitting over `Q(i)` (an irreducible `n_e` would divide both `s±it` and hence `2it`) is right.
  * The final step to a "rational Gaussian linear template" uses the strict linear system of III.5–III.8, outside
    this block. It is consistent with II.6's cocycle, but I did not check it against III.8's definitions.
    PLAUSIBLE.

## 5. Genus-four curve (III.9 lines 935-941): CORRECT, peripheral

`z^3=(c/λ)∏q_i∏(q_i-q_j)` is the cubic-lift condition: multiply `J/λ=H^{-3}` by `(∏q_i)^3`. Correct. With four marks
fixed it becomes `z^3=k·x∏_{i≤4}(x-q_i)`, a quintic with simple roots `0,q_1..q_4`. Since 5 is not divisible by 3,
`∞` is also totally ramified. Riemann–Hurwitz: `2g-2=3(-2)+6·2=6`, so `g=4`. The curve is irreducible because the
right side is not a cube. Faltings over `Q(i)` gives finiteness for fixed anchors and a fixed cube class. As Astra
says, this is not uniform in anchors or twists, which is exactly the uniformity the problem needs. Peripheral.

## 6. Novelty and significance

* **Novelty.** The fixed-`t` fourfold (dimension 4, degree 12, irreducible, dominant on ratios), the parity
  obstruction, the coefficient-uniform plane count through Bombieri–Schmidt, and the four-to-five quartic are not in
  the record. I grepped `docs/claude_review`: round 27 has "fixed cofactors" only in the split Cox-ring setting, and
  the cubic cover `v^3=J/b` first appears in Astra's own Part II. The coordinator had guessed the dimension 4.
* **Significance.** Nothing here changes the growth rate, the uniform bound, or the status of Question A, Question B
  or route 5.
  * The most consequential item is negative: dominance of a *single* `t` over ratio space kills the "fixed
    cofactors give a proper subvariety, then Schwartz–Zippel" route proposed in `reply3.txt`.
  * In the fifteen-cut setting the doubleton norms carry no fixed-`t` equation beyond the Gr(2,5) cone.
  * The plane-uniform count is a correct and neat use of Bombieri–Schmidt. However, it applies to the pure
    ten-pair model, which II.1's Walsh bound already makes `O(m^4)`. The singleton extension needs control of the
    number `A·K` of choices, which nothing supplies.
  * The interpolation bound is a triviality conditional on those same missing inputs.

  Peripheral overall. It is useful mainly as a clean closure of one proposed line.
