# Runge's method on the six-point circle cover: integrality audit

The genus-five equation has integer homogeneous representatives, but these
are not integral points on its usual affine hyperelliptic chart. The exact
Gaussian lift does not supply that missing affine integrality: every
rational point on the allowed cover already has an integral Gaussian lift
after the specified homogeneous construction and primitive division.

Removing the ten triangle-degeneration divisors creates more poles, but
also introduces denominators at the corresponding parameter values. An
explicit good-prime calculation below describes this cost. It is a scoped
applicability obstruction, not a theorem excluding every Runge function,
covering construction, or exceptional pole-splitting case.

## 1. The theorem whose hypotheses must be checked

For a fixed smooth projective curve over a number field, a rational function
with `r` Galois orbits of poles gives the standard Runge finiteness criterion
`r>|S|` for points whose function value is `S`-integral; `S` includes the
archimedean places. The Bombieri--Sprindzuk version permits the field and
the set `S` to vary while this inequality holds, but keeps the underlying
curve and function fixed. These are Theorems 2 and 3 in
[Aaron Levin, *Variations on a theme of Runge*](https://arxiv.org/pdf/0805.1345).

Thus varying prime names are not by themselves fatal. The missing input
can instead be a bound on their number, an integrality statement for the
chosen functions, or uniform control of the auxiliary functions as the
curve coefficients vary. Effective finiteness for each fixed curve is not
automatically a bound polynomial in a moving coefficient height.

## 2. Exact denominators of the ordinary affine chart

Use the controlled integer frame of
`six_point_moving_frame_radius_height.md`. Its squarefree binary form
`Delta(s,t)` has degree twelve and coefficient height `U^O(1)`. A rational
cover point can be written

```text
gcd(s,t)=1,       d in Z,       d^2=Delta(s,t).
```

On the chart `t!=0`, the affine coordinates are

```text
x=s/t,       y=d/t^6,
y^2=f(x),    f(x)=Delta(x,1).                      (1)
```

In particular the reduced denominator of `x` is exactly `|t|`. Its
minimal finite denominator support contains every prime dividing `t`.
Being an integer triple `(s,t,d)` on the weighted homogeneous equation
does not make `x,y` integral or give a point on a fixed affine integral
model.
The equation in three affine coordinates is a surface (the weighted
cone); obtaining the curve requires quotienting the scaling
`(s,t,d)->(lambda s,lambda t,lambda^6 d)`.

Let `a=Delta(1,0)`. If `p|t` and `p` does not divide `a`, then

```text
d^2 = a s^12 mod p,
v_p(d)=0,
v_p(x)=-v_p(t),       v_p(y)=-6v_p(t).             (2)
```

For an odd such prime, (2) requires `a` to be a square modulo `p`; this
is a splitting restriction when `a` is nonsquare, not a finite support
restriction. In the split-leading case the chosen local branch at
infinity can depend on the prime.

At the fiber `t=0`:

* If `a` is a nonzero rational square, `x` has two rational pole points,
  hence two pole orbits over `Q`.
* If `a` is nonzero and nonsquare, those two points form one orbit with
  residue field `Q(sqrt(a))`.
* If `a=0`, squarefreeness makes infinity a simple branch point. There
  is one pole point, with multiplicity two, hence one orbit.

Pole multiplicity does not increase the orbit count. In the first case
ordinary Runge for (1) requires `1+omega(t)<2`, so it only applies when
`t=+1` or `-1`, before including any other required places. In the other
two cases the strict criterion already fails at the archimedean place.

Adjoining a square root of a positive nonsquare `a` makes the two points
rational, but introduces two real archimedean places. It does not create
a strict Runge inequality by itself. Negative leading coefficients require
a separate real-locus discussion; for real rational points the affine
coordinate is then bounded at large absolute value directly from `f(x)<0`.

The same analysis applies after a fixed rational Mobius change of
parameter. Its denominator is a fixed linear form in `(s,t)`, up to fixed
coefficient denominators. A point-dependent unimodular change can send
the primitive parameter to denominator one, but its coefficients and the
transformed degree-twelve equation then depend on that point's height.
Such a change cannot silently be included in the controlled `U^O(1)` frame.

Even ordinary affine integrality is not intrinsic to the permitted frame.
For a frame built from a rational seed at `(s:t)=(1:0)`, the fixed
unimodular parameter change `(s',t')=(s,2s+t)` moves that same seed to
`(1:2)`. It changes coefficient heights by only an absolute factor, keeps
the homogeneous integer lift, and gives affine parameter `s'/t'=1/2`.

## 3. Gaussian integrality does not remove these denominators

The raw Gaussian lift is

```text
Z_i=d(2A x_i+B y_i+c1)
       + i(d^2 y_i-Bc1+2A c2),                    (3)
```

where here `x_i(s,t),y_i(s,t)` are the homogeneous moving marked rows,
not the affine cover coordinates `x,y` in (1). All coefficients and
evaluated terms in (3) are integers in the cleared frame. Dividing their
common Gaussian gcd yields the least integral circle realization.

Consequently this integrality holds for every rational cover point in the
nondegenerate circle locus. It imposes no additional assertion that `t`
is a unit, or that `s/t` is integral away from a specified fixed set. In
particular, at a prime in (2), the homogeneous lift is still integral.
If that prime also avoids `K(1,0)` and the content exceptions, it does not
even divide the primitive circle radius: affine chart denominators and
circle conductor primes need not coincide.

The exact radius-content identity is

```text
N=|K(s,t)|/C(s,t,d),
C=h N_G(e) N_G(f),
1<=C<=E_V^2 |E_R|^(3/2),                          (4)
```

with `E_V,E_R=U^O(1)` in the controlled frame. A bounded common content
does not make the individual triangle minors units. Their product is
exactly `-K`, and large valuations in this product normally remain in
`N` rather than disappear from a reciprocal minor.

## 4. The ten triangle divisors and their exact pole orbits

Write the ten chosen balanced triangle quadratics as `M_I`; their roots
are simple, pairwise disjoint, and disjoint from `Delta`. Thus each has
four distinct geometric points above its two roots on the double cover.

The splitting-field theorem `six_point_triangle_splitting_fields.md` gives

```text
F_I=Q(sqrt(D_I)),
D_I=-(sum_(i in I)u_i) product_(i in I)u_i.         (5)
```

At either root the conic is two distinct nonparallel lines, each determined
by its three marked points over `F_I`. Its quadratic discriminant is a
negative square over that field. Consequently the cover points above
the roots have residue field `F_I(i)`.

There are therefore:

```text
one Q-orbit per partition if D_I has neither squareclass 1 nor -1;
two Q-orbits per partition if its squareclass is 1 or -1.          (6)
```

In the second case, if the base roots are rational each has two conjugate
points over `Q(i)`; if `F_I=Q(i)`, all four points are defined over that
quadratic field and form two orbits. If `b` partitions fall in this
exceptional squareclass case, the union of the ten divisors has

```text
r_Q=10+b,       r_(Q(i))=20+2b.                    (7)
```

Thus the generic count is ten over `Q`, not forty. No claim that `b=0`
for every allowed central vector is made.
For the two explicit central vectors in the circle-cover note,
`(1,-20,250,-1000,1445,-676)` and
`(-45,616,-3850,15895,-27500,14884)`, direct integer square tests give
`b=0` in both cases: none of the ten numbers `|D_I|` is a square.

## 5. Exact pole-function values and the good-prime qualification

Choose an integer linear form `ell(s,t)` nonzero at every projective root
of `K`. A choice among finitely many small integer linear forms suffices,
since there are only twenty roots. The function

```text
psi=ell(s,t)^20/K(s,t)                             (8)
```

is homogeneous of degree zero, so is a rational function on the cover.
It has exactly the forty geometric poles above `K=0`, with the orbit
counts (7). Globally its valuation is exactly

```text
v_p(psi(P))=20v_p(ell(s,t))-sum_I v_p(M_I(s,t)).     (9)
```

Formula (9) retains specialization cancellations.

There is an explicit nonzero exception integer `E(U)` containing the
frame denominators, the relevant contents and discriminants, and the
resultants of distinct `M_I`, of `K` with `Delta`, and of `K` with `ell`.
Because all degrees are fixed and coefficients have height `U^O(1)`,
it may be chosen of size `U^O(1)`.

If `p` does not divide `E(U)` and a primitive parameter has
`v_p(M_I(s,t))=e>0`, then every other minor, `ell(s,t)`, and `Delta(s,t)`
is a unit at `p`. Thus

```text
v_p(psi(P))=-e.                                   (10)
```

This prime must belong to any `S` witnessing integrality of the function
(8). It also satisfies `p=1 mod 4` for an actual rational cover point:
at that reduced two-line conic, `Delta` is a nonzero negative square,
while `d^2=Delta` is a square. Over `Q(i)` the rational prime splits,
and (10) holds at both primes above it.

If the integral model instead uses a fixed multiple `c psi`, its valuation
is `v_p(c)-e`. Such a rescaling can absorb bounded denominator valuations;
it must be retained in the exception integer or in the explicit depth
comparison `e>v_p(c)`. A constant allowed to depend on the chosen frame
still cannot be treated as independent of `U` without a height estimate.

For example, **conditional on the presence of ten distinct such good
primes, one assigned to each partition**, in the generic `b=0` case
ordinary Runge has at least eleven places against ten pole orbits over
`Q`, or at least twenty-one against twenty over `Q(i)`. It fails its
strict inequality. The circle's real locus stays away from `K=0`, so a
tubular refinement might avoid charging its archimedean place; even
then these finite-place counts are ten versus ten, or twenty versus
twenty, with no strict surplus.

The qualification is essential. The full fair profile alone has not been
shown here to provide these ten good primes. An integer of size `U^O(1)`
may contain a substantial part of the relevant prime mass when
`U=exp(7w+o(w))`. At primes dividing `E(U)`, the exact formula (9) must
be used: different minor roots can meet, the numerator can vanish, and
valuations can cancel. Their valuations and number cannot be discarded
by calling the exception set finite.

## 6. Other poles and polynomial dependence on U

A function with poles at the coalescence divisor `d=0` has twelve
geometric pole points, but their number of rational orbits is the number
of irreducible factors of `Delta`. The verified degree-twelve irreducible
witnesses have only one such orbit. For example `ell(s,t)^6/d` has
valuation `6v_p(ell(s,t))-v_p(d)`; making these points into poles also
introduces the corresponding denominator condition. The endpoint bound
`|d|<=C_end U^a H` is not a fixed-support or bounded-prime-count condition.

There is a genuine polynomial Runge bound in the restricted affine
integral, split-leading case: for a squarefree integer polynomial `f`
of degree twelve with square leading coefficient and coefficient size
at most `B`, its integer solutions have `|x|<=B^O(1)` (with harmless
absolute constants). One elementary proof takes the polynomial part
`S(x)` of the Laurent expansion of `sqrt(f(x))`. Its rational coefficients
and a common denominator are bounded by powers of `B`, and
`f-S^2` is a nonzero polynomial of degree at most five. For sufficiently
large integer `|x|`, an integer multiple of `y-S(x)` has absolute value
less than one and must vanish; then `f-S^2` vanishes, whose integer roots
are bounded by the same coefficient estimates. The same argument works
after clearing a common parameter denominator of polynomially bounded
size. That common denominator is the missing input here.

In fact, with the additional endpoint restriction, bounded parameter
denominators can be handled directly by growth of the degree-twelve
polynomial or by the branch-approximation estimates already proved in
`six_point_endpoint_parameter_gap.md`.

The present identities supply neither such bounded denominators nor the
required small `S` for the multiple-pole alternatives. They also do not
exclude exceptional splitting cases or more elaborate covers satisfying
an independently proved integral-point hypothesis. A valid Runge
continuation must specify those functions, their exact evaluation
denominators, the actual pole orbits over the chosen field, and the
coefficient dependence of its effective bound.
