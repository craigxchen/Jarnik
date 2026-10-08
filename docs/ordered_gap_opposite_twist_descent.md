# Opposite-twist descent for the ordered pair--triple gap

The endpoint-versus-interior gap factor has an exact Gaussian
factorization which is more oriented than the radius-denominator formula.
Removing the gap does not rotate all five points together: it rotates the
outer pair and the interior triple in opposite directions.  After this
descent, geometric interlacing is only a one-parameter phase condition.
This explains precisely why primitive outer-pair reflection, by itself,
does not add a second height inequality to the primitive interior-triangle
bound.

Let `z_0,...,z_4` be distinct Gaussian integers of common norm `N`, primitive
as a full tuple.  Put `O={0,4}` and `I={1,2,3}`.  At an odd split rational
prime `p=pi conjugate(pi)`, write

```text
e=v_p(N),                 t_j=v_pi(z_j),
E=[min(t_0,t_4),max(t_0,t_4)],
I_p=[min(t_1,t_2,t_3),max(t_1,t_2,t_3)].
```

These are all prime factors of `N`.  If `2|N`, then `1+i` divides every
`z_j`; if an inert prime `q=3 mod 4` divides `N`, then `q` divides every
`z_j`.  Either conclusion contradicts full Gaussian primitivity.  Hence
`N` is odd and supported only on split primes.

Full primitivity says that the convex hull of `E union I_p` is `[0,e]`.
Let `h_p` be the distance between these two intervals.  The ordered gap
factor is

```text
K=product_p p^h_p.
```

This is the same `K` as in
[the pair--triple radius denominator](ordered_pair_triple_gap_denominator.md).

## Exact oriented descent

For every split prime with a nonzero gap choose its factor as follows:

```text
E lies above I_p: include pi^h_p in alpha;
I_p lies above E: include conjugate(pi)^h_p in alpha.
```

Then

```text
Norm(alpha)=K,
alpha divides z_0,z_4,
conjugate(alpha) divides z_1,z_2,z_3.                 (1)
```

Consequently the five Gaussian integers

```text
x_0=z_0/alpha,       x_4=z_4/alpha,
y_j=z_j/conjugate(alpha)       (j=1,2,3)             (2)
```

all have norm

```text
M=N/K.                                                   (3)
```

They are primitive as a full five-tuple, and their endpoint--interior gap
factor is one.

Division preserves every angular difference inside either group.  It
divides the physical arc lengths inside both groups by `sqrt(K)`.  Thus an
original group span at most `C N^(1/4)` becomes at most

```text
C N^(1/4)/sqrt(K) = C K^(-1/4) M^(1/4).              (3a)
```

The two groups need not lie in one short descended arc because their
relative phase has changed.

These assertions are primewise.  Suppose first that `E` lies above `I_p`.
Because their convex hull is `[0,e]`, write

```text
I_p=[0,v],       E=[v+h_p,e].
```

Division in (2) changes the allocation intervals, at the new norm exponent
`e-h_p`, to

```text
I'_p=[0,v],      E'=[v,e-h_p].                         (4)
```

Their union still has minimum zero and maximum `e-h_p`, and the two
intervals meet.  If `I_p` lies above `E`, then

```text
E=[0,b],         I_p=[b+h_p,e]
```

and division by the conjugate choices in (1) gives

```text
E'=[0,b],        I'_p=[b,e-h_p].                       (5)
```

Again the intervals meet and fill both extrema.  Primes with `h_p=0` are
unchanged.  Equations (4)--(5) prove the common norm, full primitivity, and
gap removal simultaneously.  They also show why ordinary division of all
five rows by one common Gaussian integer cannot perform this descent.

The subgroup gcds retain the same information.  If

```text
G_O=gcd_G(z_0,z_4),       G_I=gcd_G(z_1,z_2,z_3),
```

then (1) gives `alpha|G_O` and `conjugate(alpha)|G_I`.  Dividing the `x`
pair by `G_O/alpha` makes it a primitive equal-norm pair, so its two entries
are conjugate up to a Gaussian unit.  Its reflection axis is therefore a
coordinate axis or a diagonal after a unit rotation.  This reflection
belongs to the independently normalized outer pair; it is not a symmetry
of the three `y` rows.

## Interlacing is one phase interval

Let `phi=arg(alpha)`.  Choose compatible real lifts of the descended angles

```text
theta_0 < theta_4,                  arg(x_0),arg(x_4),
psi_1 < psi_2 < psi_3,              arg(y_1),arg(y_2),arg(y_3).
```

Here the lifts are inherited from the original minor arc and chosen on its
common semicircle branch; in particular `theta_4-theta_0<pi`.

Multiplication in (2) sends endpoint angles to `theta+phi` and interior
angles to `psi-phi`.  The original cyclic order

```text
z_0 < z_1 < z_2 < z_3 < z_4
```

on one lifted minor arc is therefore equivalent to

```text
psi_3-theta_4 < 2 phi < psi_1-theta_0.                (6)
```

The order within each group is already fixed and is unaffected by `phi`.
In particular the middle interior point supplies no further inequality on
the twist.  The existence of the phase interval in (6) is equivalent to
the outer descended span exceeding the interior descended span.

There is also an exact tangent-half-angle form.  Write

```text
alpha=a+ib,       alpha^2=A+iB,
A=a^2-b^2,        B=2ab.
```

The construction chooses at most one of `pi,conjugate(pi)` over every
rational prime and `K` is odd.  Hence `a,b` are coprime of opposite parity,
and

```text
A^2+B^2=K^2,       gcd(A,B)=1,
lambda=alpha/conjugate(alpha)=(A+iB)/K.               (7)
```

Thus the gap is the hypotenuse of the primitive Pythagorean phase which
oppositely rotates the two descended groups.  For `y*conjugate(x)=C+iD`,
the orientation comparison between the original endpoint `alpha*x` and
interior point `conjugate(alpha)*y` is the sign of

```text
Im(conjugate(alpha)^2 y conjugate(x))=A D-B C.         (8)
```

On this fixed lifted semicircle branch, all interlacing tests are homogeneous
linear inequalities in the same primitive vector `(A,B)`.  After the orders
inside the two groups are fixed, their intersection is exactly the angular
interval (6); the inequalities from the middle point are redundant.

There is a literal integer separation at either boundary of (6), but it is
only at the full-circle lattice scale.  For an endpoint `x` and an interior
point `y`, both of norm `M`, put `lambda=alpha/conjugate(alpha)`.  Then

```text
lambda-y/x
  =(alpha*x-conjugate(alpha)*y)/(conjugate(alpha)*x).  (9)
```

The numerator is a nonzero Gaussian integer whenever the corresponding
original points are distinct, while the denominator has absolute value
`sqrt(KM)=sqrt(N)`.  Hence

```text
abs(lambda-y/x) >= N^(-1/2).                          (10)
```

This is an exact chord separation between two unit phases.  At the target
arc scale `Delta` of order `N^(-1/4)`, (10) does not constrain `Norm(alpha)`:
it is much smaller than the available phase window.  Thus outer-pair
reflection plus order gives one interval for the same twist and the usual
`N^(-1/2)` rational phase resolution, rather than a second independent
congruence capable of improving the exponent in the known bound
`K <= N Delta_I^3/16`.

This is a limitation of this particular coupling.  It does not rule out an
additional arithmetic invariant involving all five descended points.

For overlapping windows there is an additional exact restriction:
[crossing endpoint pairs have coprime oriented factors](ordered_overlapping_gap_factors.md).
On adjacent windows their product divides the common interior-pair gcd.
The phases can still cancel, and its resulting size inequality does not
improve the individual triangle bounds at arc constant `1/2`.

## Actual ordered fixture

For

```text
(z_0,...,z_4)=
((4,-33),(9,-32),(12,-31),(23,-24),(24,-23)),
N=1105,        K=5,
```

take `alpha=-1+2i`.  Equations (2) give

```text
(x_0,y_1,y_2,y_3,x_4)=
((-14,5),(11,10),(10,11),(5,14),(-14,-5)).
```

Every descended point has norm `221`, their full Gaussian gcd is a unit,
and their new pair--triple gap factor is one.  The descended outer pair is
the conjugate pair `-14+/-5i`; the descended triple lies in a separate
sector.  The opposite rotations by `alpha` and `conjugate(alpha)` put them
back into the displayed strict interlacing order.  This concrete tuple
shows the obstruction without using an allocation-only model.

The [checker](check_ordered_gap_opposite_twist_descent.py) exhausts the
primewise allocation statement through exponent eight, tests random
multi-prime Gaussian tuples with nested powers and units, and verifies the
ordered fixture exactly.
