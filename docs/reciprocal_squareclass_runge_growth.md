# Reciprocal squareclasses improve the conditional growth rate

Suppose a primitive configuration of `m` lattice points on squared radius
`N` lies in an arc of length at most `C N^(1/4)`. Assume there exists a
nonconformal rational projective map of its half-angle coordinates whose
full primitive image has squared radius `N'` and lies in an arc of
length at most `C'(N')^(1/4)`. For any fixed positive `C,C'`, as `N`
tends to infinity,

```text
m = O_(C,C')(sqrt(log N / loglog N))
  = O_(C,C')(sqrt(log R / loglog R)), R=sqrt(N).       (1)
```

The constant can be an absolute constant times
`ceil(C/2)ceil(C'/2)`. No ordering of the two radii is required.
Sections 1--3 first prove the estimate when `C,C'<=2` and `N'<=N`;
Section 5 removes both restrictions. This improves the earlier conditional
bound `O(log R/(loglog R logloglog R))`. It uses integer square equations
in addition to the divisor count. It does not assert that an appropriate
map exists for an arbitrary endpoint configuration, so it does not improve
the current general bound or prove the requested uniform count.

## 1. Exact integer square equations at both actual anchors

Assume throughout Sections 1--3 that both endpoint constants are at most
two and `N'<=N`.

Use the actual-anchor triangular representative and the notation of
[the reciprocal stretch theorem](mobius_reciprocal_stretch_grid.md):

```text
M=((a,b),(0,d)),       a,d>0,       delta=ad,
T=a^2+b^2+d^2,       Wgram=a^2-b^2-d^2+2iab,
C=Wgram z_0,       Q=Nclip,
E=N/Q,       E'=N'/Q,
H=4delta^2 E E',       A=TE,       B=2delta E.
```

A negative determinant is handled by conjugating the entire target, as
in that theorem. All occurrences of `delta` refer to this triangular
integral representative, including the cost of actual target anchoring.

The normalized squared stretches give positive divisors

```text
n_i=2delta E x_i
   =TE+Re(conjugate(C) z_i)/Q,       n_i | H.          (2)
```

Each divisor value occurs at most twice. The conductor theorem gives
`Q|gcd(Re C,Im C)`, so the additional coordinate

```text
Y_i=Im(conjugate(C) z_i)/Q
```

is an integer as well. The physical norm identity is

```text
|conjugate(C) z_i/Q|^2
 = (T^2-4delta^2) E^2.
```

Substituting the real part from (2) yields

```text
Y_i^2=2A n_i-n_i^2-B^2.                            (3)
```

No pointwise primitivity, fixed Gaussian unit, or half-angle parity is
assumed here. The common projection divisibility supplies both integer
coordinates at once.

## 2. A squareclass relation gives a monic polynomial square

Define one integer parameter and one integer root for each divisor label:

```text
Z=2E'A=2EE'T,
c_n=E'n+E(H/n) in Z_(>0).                           (4)
```

Equation (3) becomes exactly

```text
E'Y_i^2=n_i(Z-c_(n_i)).                             (5)
```

Two distinct labels have the same root if and only if

```text
c_n=c_v  <=>  nv=(E/E')H=B^2.                       (6)
```

Every root consequently has at most two divisor labels, and every
label accounts for at most two physical points. Select one label for
each distinct root. The resulting set has size `s>=ceil(m/4)`.
This removes the complementary-pair degeneracy; no further extraction
or change of radius takes place.

Let `r=omega(H)`, the number of distinct prime divisors of `H`. Encode
a label `n` by its prime-exponent parities and one final coordinate equal
to one. These vectors lie in `F_2^(r+1)`. If `s>r+1`, at most `r+2`
selected vectors have a nonempty linear dependence. The final coordinate
forces the dependence to have even size, say `2k<=r+2`; its other
coordinates say that the product of its labels is an integer square.

Multiplying (5) for these `2k` labels gives

```text
P(Z)=product_(j=1)^(2k)(Z-c_(n_j))
    =(E'^k product_j Y_j / sqrt(product_j n_j))^2.   (7)
```

The displayed square root is rational, and `P(Z)` is an integer. A
rational number whose square is integral is an integer. Hence (7) is
an integer square value of a monic polynomial with distinct positive
integer roots. This step retains the entire denominator; no assertion
that each individual `n_i` is a square is needed.

Since `H=4delta^2 EE'` and all three factors are positive integers,

```text
H>=4,       E,E'<=H/4,
c_n<=H(E+E')<=H^2.                                (8)
```

If `Z>H^2`, the
[elementary even-product lemma](runge_even_product_square_bound.md)
with `D=H^2` gives `Z<=(8kH^2)^(2k+1)`. If `Z<=H^2`, the same
upper bound is immediate. Also `r<=log_2 H` and `k<=r+1<=H`.
Thus, with `L_H=log H`,

```text
log Z <= (r+3)(log 8+3L_H) <= 5(r+3)L_H.            (9)
```

The last inequality uses `H>=4`. We have proved the finite arithmetic
alternative

```text
either  m<=4(r+1),
or      log Z<=5(r+3)log H.                         (10)
```

## 3. Applying the critical conductor budget

Put `W=log N`. The already established actual-anchor bounds, collected
in [the previous conditional growth audit](reciprocal_divisor_conditional_growth_audit.md),
give, uniformly for `m>=256`,

```text
log H <= 74log 2 + 44W/m.                           (11)
```

The [condition-number lower bound](mobius_endpoint_condition_number.md)
has `kappa>=N^(z_m)/512`, with `z_m>=1/4` in this range. Since
`T=delta(kappa+kappa^(-1))>=kappa` and `E,E'>=1`, (4) gives

```text
log Z >= W/4-8log 2.                               (12)
```

Consider any sequence of these admissible configurations with `m`
unbounded. The original divisor count gives `m<=2tau(H)<=2H`, so
`log H` tends to infinity. Equation (11) then implies `W/m` tends
to infinity. Along this sequence, eventually

```text
log H <=88W/m,       log Z>=W/8.                    (13)
```

In the first branch of (10), `r>=m/4-1`. In the second branch,
(9) and (13) give

```text
W/8 <=440(r+3)W/m,
r+3 >=m/3520.                                      (14)
```

Consequently `r` is at least a fixed positive multiple of `m` for
all sufficiently large members of any such sequence.

Order the distinct primes of `H` increasingly. Its `j`-th prime is
at least `j+1`, so the elementary factorial bound gives

```text
H>=(r+1)!,
log H >= log((r+1)!) >> r log(r+1).                 (15)
```

Combining (13)--(15) proves

```text
m^2 log m = O(W).                                  (16)
```

This implies (1). For example, if `m>W^(1/4)` then
`log m>(log W)/4`, and (16) gives (1) directly; the complementary
case is smaller than the asserted bound for large `W`. Bounded values
of `m` are absorbed by the absolute asymptotic constant.

The reasoning uses neither a prime number theorem nor an ineffective
height estimate for a moving curve. In particular, it does not depend
on the earlier maximal-order divisor-function bound.

## 4. Scope, degeneracies, and verification

The [six-point Pell sharpness family](full_clipping_six_point_pell_sharpness.md)
has `E=E'=delta=1`, `H=4`, and labels `1,2,4`. Their roots in (4)
are `5,4,5`. The repeated root is removed in Section 2, leaving two
labels in a parity-augmented space of dimension two. There is no
dependence to which the distinct-root lemma could be applied. The
known family is therefore compatible with this theorem.

The [polynomial checker](check_reciprocal_squareclass_runge.py) verifies the
truncation and its repeated-root exception. The independent
[projection checker](check_reciprocal_squareclass_projection.py) verifies
the integer square equations on actual mapped tuples and the passage
from even squareclass relations to monic integer square values.
Root, Astra and Sol independently audited the projection normalization,
parity-augmented dependence, and asymptotic deduction; Luna independently
proved the elementary polynomial lemma.

The remaining general gap is still map existence. The direct cotangent
offset kernel has logarithmic height at least
`(1/2-O(1/m))log N-2log L-O(1)` for small common residue modulus `L`, as shown in
[the cut and residue dictionary](integer_cotangent_offset_height_cut_residue_dictionary.md).
It does not supply the small `H=N^(O(1/m))` used in (11). No passage
from an arbitrary endpoint tuple to the present conditional class has
been proved.

The [direct monic bridge](direct_projection_monic_height_bridge.md) supplies
an unconditional integer parameter `2N` and roots equal to the original
source's radial deficits. Those roots have height `O(sqrt(N))`, and their
exact common normalization removes only the actual anchor content, up to
two. Thus existence of a monic model itself is no longer the missing step;
the small root and conductor heights in (8), (11) are what still require
the nonconformal map.

## 5. Arbitrary radius ordering and fixed endpoint constants

The [arbitrary-unit retention theorem](mobius_arbitrary_unit_retention_addendum.md)
does not assume radius descent. For a nonconformal endpoint map with both
constants at most two, and `m>=8`, it and its actual-target-anchored inverse
give

```text
N' >= N^(1-4/m)/16,
N  >= (N')^(1-4/m)/16,
log N' <= (log N+4log 2)/(1-4/m).                   (17)
```

Apply Sections 1--3 with the larger realization as source. Equation (17)
bounds its logarithmic norm by `(1+O(1/m))log N+O(1)` in terms of
either original norm. This proves (1), with an absolute constant when
`C,C'<=2`, regardless of which radius is larger. Bounded point counts
are absorbed in the asymptotic statement.

Now allow arbitrary fixed `C,C'>0`. Partition the original source and
target arcs into respectively

```text
J=ceil(C/2),       J'=ceil(C'/2)
```

pieces of their own endpoint constant at most two. Assign boundary points
to just one piece. One intersection of the source and target label
classes has `q>=m/(JJ')` corresponding points. Divide the full Gaussian
gcd from each of these two selected tuples separately. Their new least
squared radii `S,S'` satisfy `S<=N` and `S'<=N'`. The normalized arc
constant decreases under such division, so both are still at most two.
Explicitly, division by a Gaussian gcd of modulus `g>=1` changes
`L` to `L/g` and `N` to `N/g^2`; hence `L/N^(1/4)` is divided
by `sqrt(g)`.

The induced map remains rational and nonconformal: division of an entire
tuple leaves each ratio to its actual anchor unchanged, and changes
between actual anchor frames are rational conformal transformations. If the
selected set is small, its cardinality is already harmless; otherwise
apply the radius-unordered result to these primitive tuples. This gives

```text
q <= K sqrt(log S/loglog S)
```

for an absolute constant and sufficiently large `S`. Bounded `S` permits
only boundedly many lattice points and is absorbed by enlarging `K`.
The function `sqrt(log S/loglog S)` is increasing for `S>exp(e)`, so
`S<=N` then implies

```text
m <= JJ' q = O(JJ' sqrt(log N/loglog N)).            (18)
```

No primitive subset radius is replaced by a larger hypothetical one, and
no comparison of the two inherited parent radii is needed in this step.
The estimate also applies in terms of `N'` by symmetry. Root and Luna
independently audited the radius comparison and subset normalization.

The existence of any such nonconformal map to another fixed-constant
endpoint configuration now suffices for the conditional bound. A proof
that every large endpoint tuple has one is still missing.
