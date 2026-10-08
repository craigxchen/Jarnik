# A product identity for the genuinely different source-anchor targets

Reanchor a rational circle configuration at each actual source point and
then apply the **same** diagonal map. The resulting least squared target
radii need not agree. They nevertheless obey an exact product relation
whose cancellation, away from the source conductor, is charged only to
source primitive pair residues. There is no diagonal-parameter or
discriminant factor in this charge.

This couples the family introduced in
[the adaptive diagonal note](adaptive_diagonal_radius_distortion.md).
It differs from the earlier
[all-anchor normalization](mobius_all_anchor_descent_audit.md), which
rewrites one fixed physical map in different frames and keeps its target
radius fixed. It does not produce an anchor with a bounded endpoint
constant, or prove an upper bound for the minimum full target radius.

## 1. Rational phases and the anchor-indexed maps

Let `q_0=1,q_1,...,q_(m-1)` be distinct rational unit-circle phases. Choose
ordinary coordinate-primitive Gaussian half-angle rows `H_i` with

```text
q_i=H_i/bar(H_i),       H_0=1.
```

Odd-odd rows are allowed. Define

```text
e_ik=Norm(gcd_G(H_i,H_k)),
b_ik=|det(H_i,H_k)|/e_ik,          B=product_(i<k) b_ik,
N=least primitive squared radius of all the source phases.
```

The residues `b_ik` are positive integers. At odd primes they are the
usual primitive source pair residues. All least primitive squared circle
radii here are odd; the half-angle parity reduction is included when
computing them.

More explicitly, the ordinary primitive relative half-angle row is
`H_i bar(H_k)/e_ik`, up to sign, and its imaginary coordinate has absolute
value `b_ik`. Its norm is the corresponding two-point denominator norm
or twice that norm. In particular `Norm(H_i)` divides `2N`. This includes
the ratio i represented by `1+i`; no Gaussian-conjugate-primitivity or
common-row-unit assumption is imposed on these half-angle rows.

Fix coprime positive integers `a,d`, put

```text
U=a+d,             V=a-d,
F(q)=(U q+V)/(V q+U),
Q_ij=F(q_i/q_j)=(U q_i+V q_j)/(V q_i+U q_j).
```

This is exactly the effect of `diag(a,d)` on half-angle rows. Indeed, for
`w=H_i bar(H_j)=x+iy`, its transformed row is `a x+i d y`, and

```text
U w+V bar(w)=2(a x+i d y).
```

No numerator or denominator vanishes: `a,d>0` gives `|U|>|V|`, while
the q's have modulus one. The key exact symmetry is

```text
Q_ji=Q_ij^(-1),              Q_jj=1.                     (1)
```

Let `n_ij=n_ji` be the least squared radius of the two phases `{1,Q_ij}`.
Let `N_j` be the least primitive squared radius of the entire target
`{Q_ij:0<=i<m}`. This is the source-reanchor-then-diagonal radius;
the different `N_j` are not assumed equal.

## 2. The exact integral quotient

For an odd split prime `p=pi bar(pi)`, use the additive valuation on
`Q(i)` and set

```text
t_ij=v_pi(Q_ij),
t_ji=-t_ij,                  t_jj=0.
```

The unit-circle identity makes `v_bar(pi)(Q_ij)=-t_ij`. Thus

```text
v_p(n_ij)=|t_ij|,
v_p(N_j)=max_i t_ij-min_i t_ij.                          (2)
```

The second identity is the full signed-valuation width, including both
Gaussian orientations. In particular a star lcm alone must not be
substituted for `N_j`.

Define

```text
A = product_(i<j) n_ij^2 / product_j N_j.
```

Then A is a positive integer, with the exact local formula

```text
v_p(A)=sum_j L_j(p),
L_j(p)=sum_i |t_ij|-(max_i t_ij-min_i t_ij) >=0.          (3)
```

Indeed, since `t_jj=0`, subtracting the range removes precisely the
largest positive magnitude and the largest negative magnitude, if present.
All other magnitudes remain nonnegative. Antisymmetry (1) makes
`sum_(i,j)|t_ij|=2sum_(i<j)|t_ij|`, proving (3). Odd inert primes and
two do not divide any of these least primitive squared radii, so they
cause no omitted valuation in the integrality assertion.

## 3. No parameter loss away from the source conductor

For a positive integer X, write `X^[N]` for its full part supported on
primes not dividing `2N`. The result is

```text
A^[N] divides (B^[N])^(m-2).                            (4)
```

Consequently the product comparison, involving actual different target
radii, is

```text
product_(i<j) (n_ij^[N])^2 / (B^[N])^(m-2)
  <= product_j N_j^[N]
  <= product_(i<j) (n_ij^[N])^2.                         (5)
```

In particular, outside the primes dividing `2NB` there is the exact
prime-by-prime identity

```text
product_j N_j = product_(i<j) n_ij^2.                    (6)
```

Here (6) asserts equality of the parts supported outside `2NB`, not
necessarily equality of the complete integers. Neither (4) nor its
exceptional prime set depends on `a,d`.

To prove (4), fix an odd split prime `p` not dividing N. Since `q_0=1`,
the source signed widths imply that all source q's are pi-adic units.
Every half-angle norm is also a p-adic unit. For distinct i,k, direct
subtraction therefore gives

```text
v_pi(q_i-q_k)
 =v_pi(2i det(H_i,H_k)/(bar(H_i)bar(H_k)))
 =v_p(b_ik)=:beta_ik.                                    (7)
```

The Gaussian source gcd has no p-part here, and two is a unit. At odd p,
U and V cannot both be divisible by p. If exactly one is divisible by p,
every numerator and denominator defining Q is a unit, so every `t_ij=0`.

Otherwise both U and V are units. At a fixed anchor j, write

```text
A_i=U q_i+V q_j,          D_i=V q_i+U q_j,
t_ij=v_pi(A_i)-v_pi(D_i).
```

All four types of valuations are nonnegative. If `t_ij,t_kj>0`, then

```text
min(t_ij,t_kj)
 <=min(v_pi(A_i),v_pi(A_k))
 <=v_pi(A_i-A_k)=v_pi(U(q_i-q_k))=beta_ik.                (8)
```

If both are negative, subtract the two D's to obtain

```text
min(|t_ij|,|t_kj|)<=beta_ik.                              (9)
```

This proof also covers primes dividing `ad`, where both the numerator
and denominator may vanish modulo pi. Such simultaneous cancellation
can only decrease the signed valuation relative to the numerator or
denominator bound used in (8)--(9). There is no determinant-unit
hypothesis. If `a=d`, then the map is the identity, and the displayed
bound is still valid.

For each sign group at anchor j, choose an index of largest magnitude.
Each remaining magnitude in that sign group is its pairwise minimum
with the chosen one, and so is bounded by the corresponding beta in
(8) or (9). These chosen source pairs are distinct and avoid j. Hence

```text
L_j(p)<=sum_(i<k, i,k!=j) beta_ik.
```

Summing over anchors counts each source pair exactly `m-2` times:

```text
v_p(A)<= (m-2) sum_(i<k) beta_ik=(m-2)v_p(B).
```

This proves (4). Inert primes contribute zero automatically. The prime
two has explicitly been excluded from the pair-content argument.

## 4. Exact fixtures, including a nontrivial correction

For the order-of-operations example

```text
H=(1, 2+i, 3+i),              (a,d)=(2,1),
```

the shared edge norms are

```text
(n_01,n_02,n_12)=(17,37,197),
(N_0,N_1,N_2)=(629,3349,7289),
A=1.
```

Thus the changed source anchor has genuinely changed the target radius,
but their complete product is the square of `17*37*197`.

The correction in (4) can be nontrivial and sharp. Take

```text
H=(1,4+i,9+i),              (a,d)=(1,2).
```

Here the source has `N=697`, `B=5`, and

```text
(n_01,n_02,n_12)=(5,85,1469),
(N_0,N_1,N_2)=(85,7345,124865),
A=5.
```

The prime five does not divide the source radius. It divides the
primitive source residue between the last two rows, and its exponent
attains `(m-2)v_5(B)=1`. Thus deleting the source-residue correction
would make the asserted product identity false even for three actual
rational circle phases.

## Verification and scope

The [exact checker](check_adaptive_diagonal_anchor_radius_product.py)
computes every reanchored target from its literal transformed half-angle
rows, obtains its radius from all target pairs with the parity correction,
and checks (3)--(6). It includes primes dividing U or V, primes dividing
ad, odd-odd source and target rows, and both displayed fixtures.
It passes 532 exact configurations, including 94 with a nontrivial
source-residue correction away from N. Sol independently checked the
signed-valuation and parity conventions; the root agent independently
read the proof and reran the checker, confirming both fixtures and the
parameter-prime cases.

The new relation is the parameter-independent content charge for the
product of **different** anchor-indexed target radii. It explains why new
diagonal primes outside the source radius and residues cannot disappear
from that product. Product information permits very uneven individual
radii, and (4) deliberately leaves the primes dividing N untreated.
No existence of a good source anchor, bounded-endpoint target, or general
growth improvement is inferred.
