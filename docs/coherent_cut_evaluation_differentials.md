# Differential constraints on coherent cut evaluations

The exact Specht filtration gives additional linear dependencies among the
cut evaluations arising from one coherent configuration. The dependencies
hold at every active depth. They sharpen the count of independent coherent
basepoint conditions, but do not put short numerical relations into the
next kernel. A separate comparison below shows why the aggregate
successive-minima count cannot supply that missing implication.

## 1. The coherent evaluation polynomial

Let `m=2q`, `1<=k<=floor(m/6)`, and put

```text
s=q-k,       p=q-3k,       N_k=binom(m,s).
```

Use `K_k` from the [exact Specht filtration](kernel_specht_filtration.md).
Fix distinct rational numbers `b_1,...,b_m`. For a cut `S` of size `s`,
write `F_(S,Q)` for the restriction which sets rows in `S` to `(0,1)`
and rows outside it to `(1,z_i)`. Package the scalar evaluations as

```text
E_(Q,b)(x)=sum_(|S|=s) F_(S,Q)(b_outside) product_(i in S) x_i.
                                                               (1)
```

Thus `E_(Q,b)` is squarefree and homogeneous of degree `s`. Define

```text
D_1=sum_i partial_(x_i),       D_b=sum_i b_i partial_(x_i).
```

**Differential identity.** For every `Q in K_k` and all scalars `a,c`,

```text
(a D_1+c D_b)^(p+1) E_(Q,b)=0.                         (2)
```

Equivalently, the restriction of `E_(Q,b)` to every affine plane parallel
to `span{(1,...,1),(b_1,...,b_m)}` has degree at most `p`. At the terminal
depth `m=6k`, both `D_1 E` and `D_b E` vanish.

Here is a proof which does not require a graph spanning assertion for
`K_k`. In the spin-one model write

```text
v_i=(1,b_i,b_i^2),       u=(0,0,1),
E_(Q,b)=[lambda^s] Q_tilde(v_i+lambda x_i u).          (3)
```

The multilinear polynomial `Q_tilde` corresponds to the original
degree-two-each binary invariant. Formula (3) follows by multilinearity.
The constituent `(2a,2b,2c)` has maximum collapse degree `q-c`, so only
the layer `c=k` can contribute to (3).

For that layer, start with the highest-weight polynomial

```text
Delta_1^(a-b) Delta_2^(b-k) Delta_3^k,
a+b+k=q.
```

After any mixing of its vector slots, let the first two coordinate
columns be `A,B`, and let `X` be the correspondingly mixed column of
the `x_i`. The degree-two leading part of `Delta_3` in `lambda` is
the fixed form determinant times `det(A,B,X)^2`. Translating the
original `x` by `alpha*1+beta*b` replaces `X` by
`X+alpha*A+beta*B`, so this factor does not change. The other leading
factors have total `x` degree

```text
(a-b)+2(b-k)=q-3k=p.
```

Consequently the top coefficient has degree at most `p` along that
plane. This survives the multilinear weight projection: introduce row
scalars `z_i` before projecting. The same `z_i` multiply each original
first coordinate, second coordinate, and `x_i`, so the determinant
identity still holds before taking the coefficient of `product z_i`.
The highest-weight orbit spans the constituent; linearity proves (2)
for all of `K_k`.

## 2. An explicit matrix for the coherent rank bound

Let `rho_(b,k)` be the rank of `Q |-> E_(Q,b)` on `K_k`. Define the
rational matrix `M_(b,k)` with columns indexed by `s`-subsets `S` and
rows indexed by pairs `(U,j)`, where

```text
|U|=2k-1,       0<=j<=p+1.
```

Its entries use the elementary symmetric polynomial `e_j`:

```text
M_((U,j),S) = e_j(b_(S minus U))   if U is a subset of S,
             0                   otherwise.             (4)
```

Indeed the coefficient of `x_U` in
`D_1^(p+1-j) D_b^j E` is `(p+1-j)! j!` times the corresponding
matrix entry applied to the coefficient vector of `E`. Hence

```text
rho_(b,k) <= N_k-rank_Q M_(b,k)
          <= N_k-binom(m,2k-1).                          (5)
```

For the second inequality, the `j=0` rows are, up to one factorial,
the map `D_1^(p+1)` from squarefree degree `s` to degree `2k-1`.
Each one-step downward map is onto at these degrees: all are at most
`m/2`, by the elementary `SL_2` up/down calculation in the
[all-cut hierarchy](all_cut_invariant_relation_hierarchy.md).
Their composition is therefore onto.

There is also an exact description of the first upper bound in (5).
For `R=Q[x_1,...,x_m]/(x_1^2,...,x_m^2)`, it is the degree-`s`
dimension of

```text
R / (sum_i x_i, sum_i b_i x_i)^(p+1).                   (6)
```

This follows by dualizing multiplication and differentiation in the
squarefree monomial basis. No general Hilbert-function formula or
surjectivity of the coherent evaluation map onto this space is asserted.
For the rectangular terminal case, the later
[odd-block balancing theorem](coherent_nagata_balancing.md) proves
surjectivity and gives the exact rank at every distinct configuration.

## 3. An evaluation identity in the rectangular terminal case

If `m=6k`, the sole terminal constituent is `(2k,2k,2k)`. Its generic
spin polynomials transform under a common linear map `g` of the three
spin coordinates by `det(g)^(2k)`. This is immediate for `Delta_3^k`
and is preserved by mixing vector slots and by weight projection.
Taking `g=diag(1,1,1+lambda)` in (3) gives

```text
E_(Q,b)(b_1^2,...,b_m^2)=Q((1,b_1),...,(1,b_m)).         (7)
```

Thus, when `m=6k`, vanishing of all coherent cut evaluations
already implies numerical vanishing at the same configuration `b`.
The evaluation equation at `b` is redundant; (7) supplies its explicit
linear combination of cut equations.

## 4. Actual numerical zeros with coherent interior basepoints

Let `P` be any configuration of distinct rational binary directions,
and let `Lambda_P` be its numerical-zero space in `K_k`. For a fixed
distinct rational configuration `b`, put

```text
U_(P,b,k)=ker(E_(-,b)) intersect Lambda_P.
```

Every depth-`k` restriction of this space has a common interior image
point: evaluate its matching coordinates at the outside `b` values.
The point is defined because a product of `2k` disjoint differences is
nonzero there, and `4k<=q+k` at an active depth. All points come from
the same globally labelled `b`. With `r_k=dim K_k`,

```text
dim U_(P,b,k) >= r_k-rho_(b,k)-1
              >= r_k-N_k+binom(m,2k-1)-1.              (8)
```

The space contains `Lambda_P intersect K_(k+1)`. If `m=6k` and
`P=((1,b_i))`, (7) improves the first lower bound in (8) to
`r_k-rho_(b,k)`.

These are actual rational numerical-zero spaces in the invariant system.
Clearing denominators makes their bases integral. Neither this step nor
the dimension assertion gives small coefficient heights or a full core
profile. They therefore do not counterexample an arithmetic descent
theorem with the required height hypothesis.

## 5. The all-depth successive-minima deficit

The [higher-kernel height coefficient](higher_kernel_evaluation_heights.md)
has the following useful probability form. For
`X~Binomial(m,1/2)`,

```text
A_k/2^m = m/4-E max(k,|X-m/2|).                         (9)
```

This follows from
`g_k(s)=s-m/2+max(k,|s-m/2|)` and the symmetry of `X`.
Since `k<=m/6` and `max(c,x)<=c+x^2/(4c)` for `c>0,x>=0`,
the variance `E(X-m/2)^2=m/4` gives the uniform bound

```text
A_k >= 2^m (m/12-3/8).                                (10)
```

Consequently, at every active depth,

```text
A_k>N_k                  for even m>=6,
A_k>2N_k                 for even m>=8.                (11)
```

For completeness, `2^m/binom(m,m/2)` increases with even `m`.
At `m=10`, (10) divided by the central binomial coefficient already
exceeds one, and at `m=12` it exceeds two. The remaining exact values
are `A_1(6)=16>N_1(6)=15`, `A_1(8)=162>2*56`, and
`A_1(10)=1048>2*210`.

Suppose the aggregate successive-minima estimate is
`sum log lambda_j<=A_k w+o(w)`. At coefficient exponent `h>0`, this
guarantees at least `r_k-1-A_k/h-o(1)` independent relations below
`exp(hw)`, subject to integer rounding. For any `h<=1`, (11) makes
this guaranteed count strictly smaller than `r_k-N_k-1` whenever
that bound is informative. Thus even the weaker one-equation-per-cut
capacity can accommodate the count; the sharper coherent dependencies
in (5) increase the available dimension further.

This comparison is a limitation of an inference from aggregate height
and dimension alone. It does not place the actual short relations in
one of the spaces (8), or obstruct a theorem using their arithmetic
coefficient structure. Such a theorem remains necessary for descent.

## 6. Integral circuits when `m=6k` and their coefficient heights

Return to `m=6k`, and now use actual integral binary rows `(X_i,Y_i)`.
In the centrally truncated half-angle normalization put

```text
P_i=X_i+iY_i=K_i product_(T containing i) H_T.
```

These `P_i` are half-angle numerators. Their norms generally differ;
they are not the original equal-radius circle points. In particular,
equal-radius dot-product identities cannot be applied to the circuits
below without an additional justified transformation.

For `|S|=2k`, define the Gaussian integer

```text
e_S=Q((0,1) on S; (P_j,Y_j) outside S).
```

The polynomial before evaluation is divisible by every outside `P_j`:
setting one more outside `P_j` to zero gives the next, identically
vanishing collapse restriction. Divide these monomial factors to write

```text
e_S=(product_(j outside S) P_j) r_S.                    (12)
```

The quotient polynomial has degree one in each of the `4k` outside
binary rows. It has diagonal weight zero and is invariant under the
common shear `Y_j -> Y_j+tP_j`. In a finite-dimensional `SL_2` module,
a weight-zero vector killed by the raising operator lies in the trivial
summand; hence this is a full `SL_2` invariant. Integral matching
straightening expresses it as an integer linear combination of perfect
matchings on the outside rows. Each matching bracket is

```text
Delta_ij=P_iY_j-P_jY_i=X_iY_j-X_jY_i in Z.
```

Thus every `r_S` is an ordinary integer. Specialization and removal of
the monomial factor in (12) do not increase the polynomial coefficient
sum, so the quotient has coefficient sum at most `C_Q`. Expressing it
in a fixed integral matching basis costs at most a constant `L_m`
depending only on `m`.

Put `E(x)=sum e_S x_S` and `R_Q(x)=sum r_S x_S`. The homogeneous-row
version of (2) and (7), with the first two spin coordinates now
`P_i^2,P_iY_i`, gives

```text
sum_i P_i^2 partial_i E=0,   sum_i P_iY_i partial_i E=0,
E(Y_1^2,...,Y_m^2)=Q(actual rows).                       (13)
```

Since `E(x)=(product_i P_i) R_Q(x_i/P_i)`, (13) implies the
denominator-free integral circuit identities

```text
sum_i X_i partial_i R_Q=0,   sum_i Y_i partial_i R_Q=0.
                                                               (14)
```

Equivalently, for every `(2k-1)`-set `U`,

```text
sum_(i outside U) (X_i,Y_i) r_(U union {i})=(0,0).        (15)
```

The exact reconstruction in these coordinates is
`sum_S r_S product_(j outside S)P_j product_(i in S)Y_i^2=Q(actual)`.
It is not the substitution `R_Q(Y^2)`.

The reconstruction holds as a polynomial identity and gives the stronger
coefficient accounting

```text
sum_(|S|=2k) ||r_S||_(monomial,1)=C_Q.                  (15a)
```

The summands have disjoint monomial supports: the zero-`P` set of a
monomial in the `S` summand is exactly `S`. Equivalently, the terminal
spin polynomial has `GL_3` diagonal weight `(2k,2k,2k)`, so every
binary monomial has exactly `2k` rows of each `P` exponent `0,1,2`.
The convention here is `v_i=(P_i^2,P_iY_i,Y_i^2)`; no factor of two
is inserted in the middle spin coordinate. Identity (15a) concerns
raw polynomial coefficients, before conversion to a matching basis.
The separate [integral circuit checker](check_terminal_integral_circuit_audit.py)
verifies reconstruction, the raw coefficient identity, matching-basis
conversion constants in two finite cases, and prime-power cancellation.

Use the exact core factorization
`|Delta_ij|=b_ij product_(T containing i,j)n_T`, with
`n_T=Norm(H_T)` and `B=product_(i<j)b_ij`. A perfect matching of
the outside `4k` rows has at least
`max(0,|T minus S|-2k)` edges wholly in `T`. Hence

```text
product_T n_T^max(0,|T minus S|-2k) divides r_S,
D_U=product_T n_T^max(0,|T minus U|-2k-1)
    divides every r_(U union {i}), i outside U.          (16)
```

These are universal content floors, not exact contents of the array
belonging to a particular numerical-zero polynomial. They are sharp
in the following precise sense. Let `G_U` be the gcd of all the
values `r_(U union {i})(Q)` as `i` varies and `Q` varies over the
saturated integral lattice `K_k`, without imposing numerical vanishing.
Then

```text
D_U | G_U | D_U B.                                      (17)
```

For the upper bound, work prime by prime. Choose `i in T minus U`
when this set is nonempty, so the exponent in (16) is minimal, and
choose an outside perfect matching attaining that minimum. Attach
each of the `2k` inside vertices to one matching edge to form `2k`
disjoint triangles. The resulting terminal graph invariant has exactly
that matching, up to sign, as its `r_S`. Each matching edge occurs
once, so the residual valuation is at most `v_p(B)`. The core norm
supports are disjoint. This proves (17); it does not bound accidental
content for a fixed numerical zero.

In a fixed fair profile `log n_T=w+O(eta w)`, the floor in (16) has

```text
log D_U=J_k w+O_m(eta w),
J_k=2^(2k-2)[(4k+1)binom(4k,2k)-2^(4k)].               (18)
```

Indeed a subset `T` chooses its intersection with `U` in `2^(2k-1)`
ways, while the remaining `4k+1` labels contribute
`sum_a binom(4k+1,a)max(0,a-2k-1)`.
Each `r_S` has the upper bound

```text
log |r_S| <= log C_Q+2k*2^(6k-2)w
             +O_m(eta w)+log B+O_m(1),                 (19)
```

when it is nonzero. A zero coordinate may simply be omitted from a
maximum norm. If the circuit at `U` is nonzero, divide its entries
by their actual positive gcd `g_U`. Since `g_U>=D_U`,

```text
log max_i |r_(U+i)/g_U|
 <= log C_Q+H_k w+O_m(eta w)+log B+O_m(1),
H_k=2^(2k-2)[(2k+1)2^(4k)-(4k+1)binom(4k,2k)].         (20)
```

The first three `H_k` are `18,2600,266560`. This is an upper bound,
not an equality for the primitive height of a specified relation.

## 7. Actual coherent basepoints are excluded at every short depth

There is a positive arithmetic conclusion at every active depth, not
only the terminal one. Return to `m=2q`, `s=q-k`, and an outside set
`O` of size `n=q+k`. The same specialization argument as in (12)
gives a residual integral multilinear polynomial

```text
e_S=(product_(j in O)P_j) r_S(P,Y),
r_S=sum_(J subset O, |J|=2k) c_J Y_J P_(O minus J),
sum_J |c_J|<=C_Q.                                      (21)
```

Here `Y_J=product_(j in J)Y_j`. Unless `m=6k`, this
residual polynomial need not be an invariant: its `P` degree is `s`
and its `Y` degree is `2k`.

**Actual-value nonvanishing theorem.** Assume primitive rows and the
full core profile, with `log Norm(H_T)>=(1-eta)w`. If

```text
log C_Q < 2^s(1-eta)w-2 sum_(j in O)log|K_j|,           (22)
```

then `r_S(actual)=0` implies `r_S=0` as a polynomial, and hence
`F_(S,Q)=0`.

To prove this, fix a nonzero coefficient `c_J` and aggregate

```text
H_J^O=product_(T intersect O=J) H_T.
```

There are exactly `2^s` factors. The set `J` is nonempty and proper
in `O`, so none of these factors is affected by an empty/full-pattern
convention. Modulo `H_J^O`, every monomial in (21) except the `J`
monomial vanishes: if `J'!=J` and `|J'|=|J|`, some `j in J minus J'`
contributes the factor `P_j`. Primitivity and odd split support make
the `Y_j`, `j in J`, units modulo this aggregate block. The outside
core factors in `P_(O minus J)` are coprime to it. Therefore an exact
numerical zero gives

```text
H_J^O divides c_J product_(j in O minus J)K_j,
Norm(H_J^O)/Norm(gcd_G(H_J^O,product_(O minus J)K_j))
    divides c_J in Z.                                  (23)
```

Although the value of `r_S` may be Gaussian, it is exactly zero and
the coefficient `c_J` is an ordinary integer. This is precisely the
setting of the real norm-divisibility rule. A nonzero coefficient
thus forces

```text
log C_Q >= 2^s(1-eta)w
           -2 sum_(j in O minus J)log|K_j|,             (24)
```

contrary to (22). No determinant-residue loss occurs.

At the terminal depth `m=6k`, the residual polynomial is a
degree-one-each invariant on `4k` rows. Complementary balanced
coefficients agree up to sign. Applying (23) to `J` and `O minus J`
therefore doubles the disjoint norm divisors and gives the stronger
sufficient threshold

```text
log C_Q < 2^(2k+1)(1-eta)w
          -2 sum_(j in O)log|K_j|.                      (25)
```

This is also the aggregated form of the existing
[degree-one balanced threshold](invariant_degree_variation.md).
The uniform correction losses in (22) and (25) are respectively
`2(q+k)sigma` and `8k sigma`.

For a collection of relations satisfying (22), a common basepoint
equal to the actual outside matching point at a given cut forces
every polynomial restriction there to vanish. If this holds at every
cut, the collection lies in `K_(k+1)`. For `m=6k`, (25)
therefore excludes the actual coherent configuration as a common
basepoint for any nonzero short relation.

This is an arithmetic exclusion of the actual basepoint. Chow-form
vanishing supplies some common point, which need not be that point.
Identifying those points, or otherwise proving all polynomial
restrictions zero, remains the missing descent step. The independent
prime-power fixtures are in
[check_actual_cut_nonvanishing.py](check_actual_cut_nonvanishing.py).

## 8. The first covolume bound for the actual coherent kernel

The new actual-value exclusion can be compared with a direct lattice
bound for the kernel of the entire coherent evaluation matrix. When
`m=6k`, divide its scalar row `r_S` by the universal row content
from (16). In a fixed integral basis of `K_k`, every normalized row
has Euclidean height at most

```text
exp(Hcoord_k w+o(w)),
Hcoord_k=k*2^(2k)[2^(4k-1)-binom(4k,2k)].               (26)
```

Choose `rho` independent normalized rows. The covolume of their
saturated integer kernel is the Euclidean exterior norm divided by
the gcd of the maximal minors. Bounding that gcd below by one and
applying Hadamard gives the first-minimum upper exponent

```text
rho Hcoord_k/(r_k-rho),       if rho<r_k.                (27)
```

Using the universal bound `rho<=binom(6k,2k)-binom(6k,2k-1)` from
(5), the ratio of (27) to the exclusion exponent `2^(2k+1)` is
approximately `85.2941,217.4053,674.2263,1750.7616` for
`k=2,3,4,5`, respectively. Thus this direct product-of-row-heights
bound does not close those cases. The calculation discards the
maximal-minor content; a sharper arithmetic estimate for that
content could change the comparison. No impossibility claim for
such a refinement is made here.

The finite checks are in
[check_coherent_cut_evaluation_differentials.py](check_coherent_cut_evaluation_differentials.py).
