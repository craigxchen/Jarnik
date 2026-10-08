# Gale radicals and the existing circle squareclasses

The squareclasses of the seven integral Gale metric values do not give an
independent parity invariant. Their pair products are exactly the norm
squareclasses of the full rational half-angle ratios of the original
circle points. Consequently the existing six-wise separation theorem
forces degree 64 in the seven-point extracted endpoint case. The earlier
degree-at-least-eight conclusion is valid but weaker than this existing
restriction.

## 1. Identify the full half-angle norm, including units

Let the positive rational binary metric be

```text
q(x,y)=a x^2+2bxy+c y^2,  ac-b^2=r^2>0,  r in Q,  a>0.
```

Put `ell(x,y)=a x+b y+i r y`. Then `Norm(ell(v))=a q(v)`.
The associated rational-circle phase is `ell(v)/conjugate(ell(v))`.
A common rotation cancels in phase ratios, and conjugating every phase
does not change the squareclasses below.

For two marked vectors set `Q_i=q(v_i)` and
`beta=ell(v_i) conjugate(ell(v_j))`. Their phase ratio is

```text
z_i/z_j=beta/conjugate(beta),
Norm(beta)=a^2 Q_i Q_j.                              (1)
```

Choose any nonzero integral Gaussian numerator `h_ij` representing the
full half-angle ratio, so `z_i/z_j=h_ij/conjugate(h_ij)`. Equality with
(1) implies `beta/h_ij` is real and belongs to `Q(i)`, hence is rational.
Therefore, in `Q*/Q*^2`,

```text
[Norm(h_ij)]=[Q_i Q_j].                              (2)
```

Rational contents in `h_ij` only change its norm by a rational square.
The *full* ratio is important: a ratio of Gaussian units equal to `i`
has numerator `1+i` and contributes the squareclass of two. Omitting
that unit bit would make (2) false for an allocation-only numerator.
In the common-unit extraction used by the six-wise theorem, the full
and allocation half-angle norm classes coincide.

## 2. Identify the primitive Gale metric values

With `D_i=product_(j!=i)det(v_i,v_j)`, the seven-point integral radical
normalization is

```text
tau=gcd_Q(Q_i^5/D_i^2 : i=1,...,7),
c_i=Q_i^5/(D_i^2 tau).
```

For every pair,

```text
c_i c_j/(Q_i Q_j)
 = [Q_i^2 Q_j^2/(D_i D_j tau)]^2.                  (3)
```

Thus the integral metric values satisfy exactly

```text
[c_i c_j]=[Q_i Q_j]=[Norm(h_ij)].                   (4)
```

This argument retains the rational gcd `tau`; it need not be a rational
square in an arbitrary initial metric normalization.

## 3. Repeated metric squareclasses are precisely square phase ratios

For any nonzero `beta in Q(i)`,

```text
beta/conjugate(beta) is a square in Q(i)*
 if and only if Norm(beta) is a square in Q*.       (5)
```

If its norm is `s^2` with positive rational `s`, then
`beta/conjugate(beta)=(beta/s)^2`. Conversely, if the ratio is `gamma^2`,
then `Norm(gamma)=1`, since this positive rational norm has square one.
It follows that `beta/gamma` equals its conjugate and is rational;
its square is `Norm(beta)`. Combining (4) and (5),

```text
c_i c_j is a rational square
 if and only if z_i/z_j is a square in Q(i)*.       (6)
```

Two integral Gaussian points in one squareclass can be written
`z_i=epsilon b w_i^2` with the same integral squarefree representative
`b` and unit class, by unique factorization. Their consistent roots
have a common modulus and lie on an arc of length at most
`C/(2 sqrt(|b|))`. The equal-norm Gaussian spacing `sqrt(2)` then gives
the exact per-class count

```text
1+C/(2 sqrt(2|b|)).                                (7)
```

Consequently an endpoint arc with `C<2 sqrt(2)` has pairwise distinct
`c_i` squareclasses, with no fair-profile assumption. See
[the root packing proof](luna_fresh_squareclass_root_descent.md).

## 4. The degree of the squared radical

Suppose `c_0,...,c_6` have pairwise distinct rational squareclasses and
let `T=sum_i a_i sqrt(c_i)` with all `a_i` nonzero rational numbers.
The seven radicals are linearly independent over `Q`: in their
multiquadratic field they belong to distinct character subspaces.
In particular `T` is nonzero.

Let `K=Q(sqrt(c_0),...,sqrt(c_6))` and
`K_0=Q(sqrt(c_i c_0):1<=i<=6)`. Clearly `T^2` belongs to `K_0`.
An automorphism of `K` fixes `T^2` exactly when it sends `T` to `T`
or `-T`. Linear independence makes this equivalent to giving the same
sign to all seven radicals. That is exactly the subgroup fixing every
pair radical. The fixed fields therefore coincide:

```text
Q(T^2)=K_0.                                       (8)
```

Seven distinct relative squareclasses already require the generating
squareclass group to have dimension at least three over `F_2`, so the
degree-at-least-eight conclusion follows from root packing alone.

More strongly, in the common-unit endpoint configurations supplied by
the central extraction, the
[six-wise norm-squareclass theorem](norm_squareclass_separation.md)
excludes every nonempty product of at most six primitive anchor norms
from being a rational square. For seven selected points there are only
six nonanchor labels. By (4) these are exactly `[c_i c_0]`, so all six
are independent. Equations (4) and (8) give

```text
[Q(T^2):Q]=[K_0:Q]=64.                            (9)
```

The hypotheses of the six-wise theorem, including its full cut support,
small correction budget, common unit class, and nonzero-word check, are
essential. Equation (9) is a statement about those extracted endpoint
configurations, not about arbitrary seven rational points or arbitrary
positive metric data.

## 5. Consequence for the research route

The positivity proof for the alternating Gale radical remains an exact
identity and its norm inequality remains valid. But its exclusion of
degree at most four does not exclude an additional extracted endpoint
case: the earlier six-wise theorem already forces the full degree 64.
Similarly, detecting a repeated pair squareclass recovers root descent.
A successful continuation must address the full independent-class case
or obtain information beyond these parity restrictions.

The [exact checker](check_gale_radical_squareclass_identification.py)
verifies 168 pair comparisons in rational metrics, all four unit ratios,
and a rational seven-node example whose squared radical has all 64
Galois conjugates. The example is not asserted to be a fair endpoint
configuration; it checks the field and normalization calculations.
