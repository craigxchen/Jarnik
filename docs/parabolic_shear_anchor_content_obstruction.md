# A critical parabolic shear requires a large-content anchor

This note sharpens the discriminant charge for the exact equal-diagonal
case of the [critical triangular reduction](mobius_finite_sample_critical_form.md).
It uses the actual source anchor, including its ordinary coordinate gcd.
The result excludes that case at small-content anchors. It neither
constructs a descent map nor excludes the general unequal-diagonal case,
and does not prove the uniform endpoint bound.

Let `z_0,...,z_(m-1)` be a primitive Gaussian integral circle tuple, with
common norm `N`. Thus `N` is odd and supported on split primes. Fix an
actual anchor `j`, use the phases `q_i=z_i/z_j`, and put

```text
g_j=gcd_Z(Re z_j, Im z_j).
```

Suppose that in these source coordinates, and after anchoring the image
at its corresponding point, the half-angle map has primitive representative

```text
M = ((a,b),(0,a)),       a>0, b!=0, gcd(a,b)=1.
```

These are exact hypotheses. In particular, equality of the two diagonal
entries is not a consequence of the near-shear reduction. Write `N'` for
the least squared radius of the full image and `Nclip` for the source
conductor clipped to the coefficient intervals.

## 1. An exact improvement over the discriminant

The new divisibility is

```text
Nclip | |b|(b^2+4a^2) g_j.                            (1)
```

The general discriminant bound would only give
`Nclip | b^2(b^2+4a^2)`. The additional factor of `b` can be charged to
the actual anchor content instead.

To prove (1), write the map on complex half-angle vectors as

```text
Mz=alpha z+beta conjugate(z),
U=2alpha=2a-ib,             V=2beta=ib.
```

At an odd split prime `p=pi conjugate(pi)`, let
`e=v_p(N)` and `e_j=v_pi(z_j)`. The source phase valuations range from
`-e_j` to `e-e_j`. The coefficient interval has endpoints

```text
v_pi(V)-v_pi(U),       v_pi(conjugate(U))-v_pi(conjugate(V)).
```

If `p|b`, coprimality of `a,b` makes `U` a unit at both orientations.
With `t=v_p(b)`, the coefficient interval is exactly `[-t,t]`.
Consequently its clipped width is

```text
min(t,e_j)+min(t,e-e_j)
 <= t+min(e_j,e-e_j)
  = v_p(b)+v_p(g_j).                                 (2)
```

If `p|(b^2+4a^2)`, then `p` does not divide `b`. The two orientations
of `U` cannot both have positive valuation: a common divisor would divide
`4a` and `2ib`, contrary to coprimality at an odd prime. The coefficient
interval therefore has length exactly `v_p(b^2+4a^2)`. At every other
odd prime its length is zero. These two prime sets are disjoint, since
their common odd prime would divide `a` and `b`. There is no source
conductor at two. Multiplying (2) and these remaining bounds proves (1).

The factor `g_j` is necessary. Take an odd split prime `p=pi bar(pi)`,
an integer `t>=1`, and the primitive physical tuple

```text
(p^t, pi^(2t), bar(pi)^(2t)),       N=p^(2t).
```

Anchor at `p^t` and take `a=1,b=p^t`. The source valuation interval
and the coefficient interval are both `[-t,t]`. Hence `v_p(Nclip)=2t`,
whereas `v_p(b)=t`. The missing depth is exactly `v_p(g_j)=t`.
This is an exact local sharpness example, not a large endpoint cluster.

## 2. Radius and condition number force the anchor cost

Assume now that both endpoint constants are at most two, `N'<=N`, and
`m>=8`. Retain `q_m`, `u_m`, `g_m`, `b_m`, and `D_m` from the critical
triangular note. To avoid confusing its coefficient `b_m` with this
matrix's shear entry `b`, put

```text
h_m=g_m+b_m,
eta_m=1-q_m-(3/2)u_m-3h_m,
K_m=40 D_m^3 2^(m q_m+(3/2)u_m).
```

If this actual anchor also satisfies the small-diagonal bound
`a<=D_m N^h_m`, then

```text
g_j >= K_m^(-1) N^eta_m,
eta_m=1/4-61/(4m)+O(m^(-2)).                         (3)
```

In particular `eta_m` is positive for all sufficiently large `m`.
The small-diagonal hypothesis holds at the anchor selected by the
critical triangular theorem; exact diagonal equality is the additional
condition being tested here.

Indeed, the condition number `kappa` obeys

```text
kappa+kappa^(-1)=2+b^2/a^2,
|b|<=a sqrt(kappa),
b^2+4a^2<=5a^2 kappa.
```

The finite-sample upper bound `kappa<=4(2N)^u_m` therefore gives

```text
|b|(b^2+4a^2)<=40a^3(2N)^((3/2)u_m).
```

Combine this with (1) and the erased-conductor estimate
`Nclip>=2^(-m q_m) N^(1-q_m)` to obtain (3). This proof does not assume
integrality of the image cotangents at the original source residue scale,
and it includes all new primes in the actual least image radius.

More generally, if `a<=A_0 N^h`, the same proof has exponent
`1-q_m-(3/2)u_m-3h` and constant
`40 A_0^3 2^(m q_m+(3/2)u_m)`. Thus a subpower diagonal and a subpower
anchor content cannot coexist with a large equal-diagonal endpoint map
as `m` tends to infinity. This is a conditional map obstruction.

There is also a uniform version with coarse absolute constants. For
`m>=256`, all the hypotheses of (3) imply

```text
g_j >= 2^(-57) N^(1/8).                               (4)
```

Here is an elementary verification, so that no limit hides dependence on
`m`. Write `m=4k+r`, `0<=r<=3`, and `s=m-2k`. The definitions give

```text
q_m = m(m-1)/[2k((m-1)(s-1)+1)] <= 8/m.
```

The last inequality follows from
`16k(s-1)-m^2=16k(k-1)+8kr-r^2>0` for `k>=2`.
Also `F_m>=m(m-2)/4` and `Q_m<=m^2/4`, so

```text
h_m <= 2/(m-2)+9/[4(m-1)] <= 7/m       (m>=8).
```

The exact formula for `u_m` gives `u_m<=1/2+2/m`. Thus
`eta_m>=1/4-32/m>=1/8`. Finally
`m(m-1)/F_m<=5` and `m b_m<=m h_m<=7`, giving `D_m<=2^14`;
`m q_m<=8` and `(3/2)u_m<=1` give `K_m<2^57`. These prove (4).

## 3. Only a small set of anchors can pay this cost

For source endpoint constant at most `sqrt(2)`, the established
[interior-allocation estimate](endpoint_allocation_rounding.md) says

```text
sum_j log g_j <= (sqrt(m)/2) log N.                   (5)
```

Let `S` be the set of source anchors at which there exists a map satisfying
all the hypotheses of (3), with the same `m,N` and small-diagonal bound.
The map and its image may vary with the anchor. If
`eta_m>0` and `N>=K_m^(2/eta_m)`, then (3) and (5) imply

```text
|S| <= sqrt(m)/eta_m.                                (6)
```

Indeed every such anchor contributes at least
`(eta_m/2)log N` to the nonnegative sum in (5). This makes the eligible
fraction tend to zero as `m` grows. It does not assert that every anchor
admits a small-diagonal representative with equal diagonals.
In particular `m>=256` and `N>=2^912` give `|S|<=8 sqrt(m)` directly
from (4)--(5).

For a literal binary allocation at every source prime, every physical
row has `g_j=1`. In that class (3) yields the finite radius bound
`N<=K_m^(1/eta_m)` whenever its exponent is positive. For example,

```text
eta_128=29828621/231197785 > 1/8.
```

The uniform version (4) gives `N<=2^456` for every such binary-source
map with `m>=256`. This radius bound is independent of both the source
radius and the number of source points.

For the integer shear `((1,b),(0,1))`, the small-diagonal hypothesis
holds automatically at its stated anchor; no selection theorem is
needed. Thus the anchor count (6) applies to the existence of any such
integer shear based at a source point. On a binary source with `m>=256`
and `N>2^456`, none can have both a radius-nonincreasing image and the
stated endpoint constants. Squarefree `N` is one natural binary case.

This family corollary retains the map-existence and diagonal hypotheses;
it is not a radius bound for arbitrary binary-profile endpoint tuples.

## 4. Scope and verification

The ordinary-anchor content bound (5) is already known. The new step is
its exact appearance in the shear clipping divisor (1), which replaces
one full factor of the shear height and gives (3)--(6).

The general critical representative can have `a!=d`. Then the two
discriminant factors are `b^2+(a+d)^2` and `b^2+(a-d)^2`, and neither
is forced to be a scalar square. Argument (2) does not apply to them.
No equality of diagonals, small content at the selected anchor, or
radius-decreasing map is presumed. The uniform count remains unproved.

The [exact checker](check_parabolic_shear_anchor_content.py) tests the
clipped widths, independently computed Gaussian coefficient intervals,
ordinary anchor contents, source and image least norms, and the
sharpness example. It also checks the rational exponent calculations
and the condition-number inequalities using rational certificates.
A separate agent audited the local divisibility.
