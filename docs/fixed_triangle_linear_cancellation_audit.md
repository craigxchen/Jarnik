# Fixed linear triangle cancellation: a universal leading-term obstruction

A nonzero fixed linear combination of triangle determinants cannot cancel
the cubic small-arc term on every ordered shape. If that leading term
vanishes universally, the combination is already an exact affine area
identity for arbitrary planar points. This is a restricted geometric audit,
not an obstruction to cancellation on a special arithmetic family.

The [higher-degree Gram-jet barrier](higher_degree_gram_jet_barrier.md)
concerns prescribed source products and explicitly leaves general signed
combinations open. The elementary injectivity statement here concerns
linear combinations of actual triangle determinants instead.

## 1. The precise universal statement

Fix `m>=3` and fixed real coefficients `lambda_ijk`, for `i<j<k`.
For planar points write

```text
D_ijk=det(P_j-P_i,P_k-P_i),
F(P)=sum_(i<j<k) lambda_ijk D_ijk.
```

Let `E_ij=D_0ij` for `1<=i<j<m`. The affine area relation

```text
D_ijk=E_ij-E_ik+E_jk        (0<i<j<k)
```

reduces `F` to a fixed sum `sum mu_ij E_ij`. This reduction is valid
for every planar configuration, without concyclicity or integrality.

Now take actual real circle points

```text
P_0=R(1,0),       P_i=R(cos(epsilon a_i),sin(epsilon a_i)),
N=R^2,           0<a_1<...<a_(m-1)<1.
```

For a fixed ordered shape `a` and `epsilon -> 0+`,

```text
E_ij=(N epsilon^3/2) a_i a_j(a_j-a_i)+O_a(N epsilon^5),
F=(N epsilon^3/2) P_mu(a)+O_(a,lambda)(N epsilon^5),
P_mu(a)=sum_(1<=i<j<m) mu_ij a_i a_j(a_j-a_i).         (1)
```

The cubic polynomials `a_i a_j(a_j-a_i)` are linearly independent:
the monomial `a_i a_j^2`, with `i<j`, occurs in exactly the polynomial
indexed by that pair. Therefore the following are equivalent:

* `P_mu` vanishes on the open ordered shape simplex;
* every `mu_ij` is zero;
* `F` is identically zero for all planar point configurations.

The implication from vanishing on the simplex uses the elementary fact
that a real polynomial vanishing on an open set is the zero polynomial.
Consequently a universal improvement from cubic order to
`o(N epsilon^3)`, or to `O(N epsilon^5)`, cannot be obtained from a
nonzero fixed-coefficient linear form using circle order alone. This
statement fixes the coefficients and quantifies over the full open set
of shapes; it is not an assertion along one selected sequence of shapes.

## 2. Rational circle realizations and the essential caveat

If `P_mu` is nonzero, there is an ordered positive rational shape `a`
with `P_mu(a)!=0`, by density and continuity. Instead use the rational
circle phases with half-angle parameters `t_i=epsilon a_i`, `t_0=0`.
For rational positive `epsilon`, they admit actual integral circle
realizations by Gaussian denominator clearing and full gcd removal.
In every such primitive realization the exact identity is

```text
E_ij/N = 4 epsilon^3 a_i a_j(a_j-a_i)
           /[(1+epsilon^2 a_i^2)(1+epsilon^2 a_j^2)].   (2)
```

Its angular span is `Delta=2 arctan(epsilon a_(m-1))`. Consequently

```text
F/(N Delta^3) -> P_mu(a)/(2 a_(m-1)^3) != 0.
```

The radius `N` here is whatever the exact primitive clearing produces;
common Gaussian division leaves the displayed ratio unchanged. This gives
actual primitive integer realizations at arbitrarily small angular width.
There is no claim that `N^(1/4) Delta` remains bounded, so these examples
do not refute an endpoint-specific arithmetic cancellation theorem.

Conversely, a nonzero fixed form can cancel its cubic term on a proper
shape subvariety. An explicit example is

```text
F=D_012-D_123,
(t_0,t_1,t_2,t_3)=(0,epsilon,2epsilon,3epsilon).
```

Using (2) and its three-index version gives the nonzero exact value

```text
F/N = 72 epsilon^5
       /[(1+epsilon^2)(1+4epsilon^2)(1+9epsilon^2)].    (3)
```

Extra ordered rows can be retained without changing this identity.
Thus even bounded integer coefficients, after passing to a subsequence
on which they are fixed, do not yield a contradiction: the shapes may
converge to or lie in the zero set of that particular cubic. The theorem
excludes universal cancellation on an open set, not cancellation on
special shape loci, nor unbounded or arithmetic-dependent coefficients.

## 3. Anchored Pluecker triples are the same primitive Ptolemy triples

For five points in their order on a minor circle arc, put
`l_ij=|z_i-z_j|` and `E_ij=D_0ij>0`, for `1<=i<j<=4`.
The anchored determinant relation is

```text
E_12 E_34 + E_14 E_23 = E_13 E_24.
```

Since `E_ij=l_0i l_0j l_ij/(2R)`, all three terms have the common
positive factor `product_(j=1)^4 l_0j/(4R^2)`. After its cancellation
the relation is exactly

```text
l_12 l_34+l_14 l_23=l_13 l_24.
```

In particular, divide the three integer determinant products by their
full ordinary gcd. The resulting positive coprime integer triple is
exactly the primitive Ptolemy triple: it is the unique such triple
proportional to the three chord matching products. No residual source
factor appears merely by changing from chords to triangle determinants.
This includes all accidental contents and shared prime powers.

The [positive-pentagon content note](positive_pentagon_primitive_content.md)
already retains the simultaneous primitive cancellation of these five-point
Ptolemy systems. The calculation here identifies the proposed anchored
Pluecker input with that existing input; it supplies no new conductor
estimate. Other nonlinear functions of the triangle data, or arithmetic
restrictions forcing special shape loci, are outside this audit.

## Verification

The [checker](check_fixed_triangle_linear_cancellation_audit.py) verifies
the affine reduction and independent cubic coefficients symbolically,
checks the nonzero exceptional-shape identity on fully primitive rational
circle realizations, and compares full gcd-normalized determinant triples
with independently reduced Gaussian chord-product triples. Finite checks
supplement the universal polynomial argument above. No new uniform
circle-point or endpoint height bound is claimed.
