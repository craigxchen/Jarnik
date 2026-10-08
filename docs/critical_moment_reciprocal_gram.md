# Reciprocal magnitudes force an exact rank-three cross Gram matrix

This is a necessary condition for the **full** punctured orthogonal
eight-row moment problem. It does not assert an eight-row solution or
give a bound for arbitrary endpoint circle arcs. Unlike the auxiliary
coefficients in [the scalar-product relaxation](critical_moment_exact_crossproduct_relaxation.md),
the quantities here are computed from the actual magnitudes.

Use the notation of [the common-jet identities](critical_moment_common_jet_products.md):
`n=4s+2`, `s>=1`, positive distinct magnitudes `a_h`, groups `A,B` of
four rows, common odd moments through degree `2s-1`, and actual middle
coefficients `c_i`. Define

```
beta_i = sum_h S_ih/a_h,
gamma = sum_h a_h^-2,
C_ik = sum_h S_ih S_kh/a_h^2,     i in A, k in B.
```

Then all eight `beta_i` are distinct and

```
(beta_i-beta_k)^2 = 4 sum_(h in D_ik) a_h^-2,
C_ik = gamma - (beta_i-beta_k)^2/2.                    (1)
```

Consequently the `4 x 4` weighted cross Gram matrix `C` has rank
**exactly three**. For any ordered triples `I` of rows in `A` and `J`
of rows in `B`, write `V(beta_I)=product_(r<t)(beta_(I_t)-beta_(I_r))`.
Every three-by-three minor is explicitly

```
det C_(I,J) = -(1/4) V(beta_I) V(beta_J).              (2)
```

In increasing reciprocal-sum order all these minors are strictly
negative. If `H=I_4-(1/4)11^t`, then

```
H C H = (H beta_A)(H beta_B)^t                        (3)
```

has rank exactly one. These are reciprocal-magnitude constraints beyond
the scalar signed-product matrix of the common-jet note. No assertion
of their logical independence from all the full polynomial identities
is intended: they follow from those identities near zero.

## Proof from the residual polynomials

For a cross pair, the differing-root polynomial has the form

```
R_ik(X)=F_ik(X)+r_ik,   F_ik odd, deg F_ik=2s+1,
r_ik=(c_i-c_k)/2 != 0.
```

Its roots are `z_h=S_ih a_h`, `h in D_ik`. Put
`p_l=sum_(h in D_ik) z_h^-l`. Logarithmic differentiation at zero gives

```
R'(0)/R(0)=-p_1,
R''(0)/R(0)=p_1^2-p_2=0.
```

Since `beta_i-beta_k=2p_1`, the first identity in (1) follows. It also
shows these cross-pair reciprocal sums are different. The second follows
by writing `C_ik=gamma-2 sum_(D_ik) a_h^-2`.

For a within-group pair the residual has the form

```
R_ij(X)=Q_ij(X^2)+r_ij X,
r_ij=(c_i-c_j)/2 != 0,       Q_ij(0) != 0.
```

The same logarithmic derivative now gives
`beta_i-beta_j=-2r_ij/Q_ij(0) != 0`. Thus all eight reciprocal sums
are distinct.

Let `V_A` have rows `(1,beta_i,beta_i^2)`, and define `V_B` similarly.
Equation (1) factors the entire matrix as

```
C = V_A K V_B^t,
K = [[gamma,0,-1/2], [0,1,0], [-1/2,0,0]],
det K=-1/4.
```

Both Vandermonde matrices have column rank three, proving exact rank
three and (2). Centering removes every term of (1) except
`beta_A beta_B^t`, proving (3). Distinctness makes both centered vectors
nonzero.

## Reciprocal order and the signs of within-group distances

Let `u_i=(S_ih/a_h)_h`, so `beta_i=u_i dot 1` and `||u_i||^2=gamma`.
For any pair put

```
L_ij = ||u_i-u_j||^2 - (beta_i-beta_j)^2.
```

Cross pairs have `L_ik=0`. For a within-group residual write
`Q(Y)=q_0+q_2 Y+...`, where `q_2` is the coefficient of `X^2`
in `R(X)`. Logarithmic differentiation gives exactly

```
L_ij=-8 q_2/q_0.                                     (4)
```

The polynomial `R''` is even, of degree `2s`, with positive leading
coefficient and simple real roots by Rolle's theorem. None can be zero:
an even polynomial would then have a repeated zero. Hence
`sign(q_2)=(-1)^s`. In successor type A the differing roots have
`s+1` negative and `s+1` positive signs, so
`sign(q_0)=(-1)^(s+1)` and `L_ij>0`. In successor type B their counts
are `s,s+2` in either order, so `sign(q_0)=(-1)^s` and `L_ij<0`.

Similarly, `F'` for a cross pair is even with `2s` simple real roots
and positive leading coefficient. Thus `sign(F'(0))=(-1)^s`.
Together with the logarithmic derivatives this proves

```
pair type        sign(beta_i-beta_j)/sign(c_i-c_j)
cross                         (-1)^(s+1)
successor A                   (-1)^s
successor B                   (-1)^(s+1).             (5)
```

For the existing even-`s`, `q=2` orientation with
`c_6<c_7<c_0<c_1<c_2<c_3<c_4<c_5`, (5) therefore requires

```
beta_4<beta_5<beta_0<beta_1<beta_2<beta_3<beta_6<beta_7. (6)
```

## The fourth-order cross Gram has a rank-one quotient

There is a further simultaneous condition from the next even reciprocal
coefficient. Set

```
eta_i=sum_h S_ih/a_h^3,
gamma_4=sum_h a_h^-4,
C^(4)_ik=sum_h S_ih S_kh/a_h^4       (i in A,k in B).
```

On a cross support, `R(X)=F(X)+r` gives `e_2(1/z)=e_4(1/z)=0`.
For degree three, `e_4=0` is interpreted as the elementary symmetric
function of order greater than the number of variables. Newton's
identities give `p_1^2=p_2` and `4p_1p_3=p_2^2+3p_4`. Substituting
`2p_1=beta_i-beta_k`, `2p_3=eta_i-eta_k`, and
`2p_4=gamma_4-C^(4)_ik` yields

```
C^(4)_ik=gamma_4
 -(2/3)(beta_i-beta_k)(eta_i-eta_k)
 +(beta_i-beta_k)^4/24.                               (7)
```

Let `P_A` and `P_B` be the Euclidean orthogonal projectors onto the
complements of `span(1,beta_A)` and `span(1,beta_B)`, respectively.
Both have rank two. Applying them to (7) cancels every term involving
`eta`, as well as all terms of the fourth power except its mixed
quadratic term. Therefore

```
P_A C^(4) P_B=(1/4)(P_A beta_A^2)(P_B beta_B^2)^t     (8)
```

has rank **exactly one**. Here vector powers are coordinatewise.
Neither projected square vector vanishes, because a nonzero quadratic
cannot agree with an affine function at four distinct `beta` values.
Thus (8) eliminates the unknown third reciprocal sums entirely.

An explicit formulation needs no orthogonal projectors. For any triple
`I` in one group define normalized second divided-difference weights

```
lambda_(I,i)=1/product_(j in I,j!=i)(beta_i-beta_j).
```

These weights annihilate `1,beta` and have pairing one with `beta^2`.
For every triple `I` in `A` and triple `J` in `B`, (8) gives the same
strictly positive value

```
sum_(i in I,k in J) lambda_(I,i) C^(4)_ik lambda_(J,k)
 =1/4.                                               (9)
```

Equivalently, the `6 x 6` bordered matrix

```
B_4=[[C^(4), 1_A, beta_A],
     [1_B^t, 0, 0],
     [beta_B^t, 0, 0]]
```

has rank **exactly five**. Each `5 x 5` minor formed by retaining
three ordinary rows `I`, three ordinary columns `J`, and both border
rows and columns equals

```
det (B_4)_(I,J,borders)=(1/4)V(beta_I)V(beta_J).        (10)
```

To see the rank assertion directly, a kernel vector has ordinary
component `y` orthogonal to `1_B,beta_B`. Projection of its upper
equations by `P_A`, followed by (8), also makes `y` orthogonal to
`beta_B^2`. This leaves a one-dimensional space; for each such `y`
the two border coordinates are determined uniquely. Expanding the
two border rows and columns of a `5 x 5` minor gives (9) multiplied
by the two Vandermonde factors, proving (10). These minors are positive
in increasing `beta` order, whereas the unbordered second-order minors
in (2) are negative in that order.
There is no sign contradiction here: the matrices in these two minor
identities have different entries and different borders.

These are coupled fourth-order necessary conditions on the actual
magnitudes. They are not an assertion that the entire `C^(4)` has
rank one, nor a contradiction for eight rows. They follow from the
full cross residual identities; no logical independence from those
identities is claimed.

## The within-group third coefficient is a divided-difference formula

For a within-group pair, the residual `Q(X^2)+rX` has zero cubic
coefficient, so `e_3(1/z)=0`. Newton's identity
`p_1^3-3p_1p_2+2p_3=0` gives, with `u=beta_i-beta_j!=0`,

```
eta_i-eta_j=-u^3/8+(3u/4)(gamma-C_ij),
C_ij=gamma-(4/3)(eta_i-eta_j)/u-u^2/6.                (11)
```

Here `C_ij=sum_h S_ih S_jh/a_h^2` extends the earlier notation to
within-group entries. Thus the off-diagonal same-group Gram entries
are determined by the divided differences of `eta` as a function of
`beta`. The within-group distance quantity becomes

```
L_ij=(8/3)(eta_i-eta_j)/(beta_i-beta_j)
     -(2/3)(beta_i-beta_j)^2.
```

Consequently the known successor-A sign requires the divided difference
to be strictly greater than `(beta_i-beta_j)^2/4`; successor B requires
it to be strictly smaller. This is another joint necessary condition,
not a proof that these matrices cannot be positive semidefinite.

## Scope and exact checks

Cross-null distances, positive within-group distances on each sheet,
and negative distances across sheets alone are not contradictory in
a space with one negative square. In `R^(4,1)`, take four points with
spatial parts `e_1,-e_1,e_2,e_3` and time zero. Take four further points
with spatial part `y e_4` and time `t`, where
`(y,t)=(0,1),(3/4,5/4),(0,-1),(3/4,-5/4)`. They have exactly those
distance signs. This is only a distance-sign example: it lacks the
common coordinate magnitudes, common spatial norms, and moment
conditions of the required sign-matrix model.

The [exact checker](check_critical_moment_reciprocal_gram.py) checks
the matrix factorization, all sixteen three-by-three minors, centering,
the fourth-order formula on pair fixtures, its sixteen divided differences
and bordered minors on abstract coefficient data, three rational
pair-residual examples, the orientation table, and this
distance-sign example. Pair examples and abstract matrix data are not
presented as a full eight-row moment solution. The full system and its
connection to arbitrary endpoint arcs remain unresolved.

Further consequences are now proved separately. The
[prefix obstruction](critical_moment_reciprocal_prefix_obstruction.md)
uses one cross rectangle to exclude every magnitude tuple of the
prescribed large-shift scalar-cross-product construction. The
[cross-support Hankel condition](critical_moment_cross_support_hankel.md)
combines reciprocal and forward moments to give an exact projection
norm on each cross support, with strict lower-degree inequalities.
The [reciprocal concentration bound](critical_moment_reciprocal_concentration.md)
uses the fourth reciprocal coefficient to bound the effective support
of each reciprocal-square measure independently of the degree.
None excludes arbitrary full templates or gives a circle-arc bound.

The [reciprocal relaxation construction](critical_moment_reciprocal_relaxation_feasibility.md)
now satisfies the punctured orthogonal code, common first moment, all
cross Gram identities, distinct reciprocal sums and the required
within-group distance signs simultaneously. Its cross supports fail the
critical root-sign order and the fourth reciprocal coefficient identity.
Thus the stated weaker conditions do not by themselves imply the full
moment system.

The full residual identity also gives
[uniform reciprocal-square tail tightness](critical_moment_uniform_reciprocal_tail.md),
independent of the degree after the smallest magnitude is normalized.
This strengthens total-mass bounds but does not bound the number of
columns or establish a lattice-circle estimate.
