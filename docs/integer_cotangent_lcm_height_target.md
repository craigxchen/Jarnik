# An exact ordinary-integer lcm formulation of the endpoint target

The Gaussian denominator bookkeeping can be eliminated completely if
all pairs are retained. This note proves the resulting radius formula
and states an equivalent integer height target. The height target is
**unproved**. This reformulation does not improve the established general
growth bound.

## 1. The integer data and their edge norms

Let `L>0` and let `X_1,...,X_k` be distinct integers such that

```text
Q_ij=(X_i X_j+L^2)/(X_i-X_j) is an integer (i<j).
```

Adjoin the anchor `q_0=1` to the rational circle points

```text
q_i=(X_i+iL)/(X_i-iL),                 1<=i<=k.
```

Set `Q_0i=X_i`. For every edge `e` of these `k+1` labels, define

```text
g_e=gcd(Q_e,L),
a_e=Q_e/g_e,       b_e=L/g_e,
epsilon_e=2 if a_e,b_e are both odd, and 1 otherwise,
n_e=(a_e^2+b_e^2)/epsilon_e.
```

All these quantities are ordinary integers and `n_e>=1`. In particular,
the definitions allow `Q_e=0` and retain common integer factors and the
prime two. Put

```text
N=lcm_(all edges e) n_e.                              (1)
```

**Exact radius formula.** The least radius of a Gaussian integral
realization of this rational circle tuple is

```text
R=sqrt(N).                                            (2)
```

Here a realization multiplies all `q_i` by the same complex number.
The primitive realization is unique up to a Gaussian unit. This is
the same notion of least radius as in
[the Gaussian normalization note](integer_cotangent_gcd_reduction.md).

## 2. Proof of the radius formula

First take any coprime integers `a,b`, with `b>0`. The Gaussian gcd of
`a+ib` and `a-ib` has norm one, except when `a,b` are both odd, when
it has norm two. An odd Gaussian prime dividing both would divide
`2a` and `2ib` and contradict integer coprimality. At `1+i`, the norm
of `a+ib` is exactly divisible by two when both are odd; otherwise
that Gaussian prime is absent. Thus the norm of the denominator in
the reduced Gaussian fraction

```text
(a+ib)/(a-ib)
```

is exactly `(a^2+b^2)/epsilon`. This proves that each `n_e` is the
norm of the reduced Gaussian denominator of `q_j/q_i`; reversing an
edge conjugates the ratio and leaves its denominator norm unchanged.

Construct a primitive integral tuple by taking a Gaussian lcm of the
reduced denominators of the `q_i`, as in the normalization note.
Write it as `z_0,...,z_k`, and put `N_0=|z_i|^2`.

No inert rational prime or ramified Gaussian prime divides this
primitive tuple's common norm. Such a prime would have the same
positive Gaussian valuation in every coordinate, contradicting
primitivity. For a split rational prime `p=pi bar(pi)`, put

```text
s=v_p(N_0),       e_i=v_pi(z_i).
```

Then `v_bar(pi)(z_i)=s-e_i`. Primitivity at both orientations gives

```text
min_i e_i=0,       max_i e_i=s.
```

The norm of the reduced denominator of `z_j/z_i` has `p`-valuation
`|e_j-e_i|`. Its maximum over all pairs is therefore exactly `s`.
Taking the ordinary integer lcm over all pair denominator norms
recovers every prime power of `N_0`. Hence `N=N_0`, proving (2).

The use of **all edges** is essential. For example, `L=12` and
`(X_1,X_2)=(4,36)` give anchor edge norms `5,5`, but the third edge
has `Q_12=-9` and norm `25`. A primitive realization is

```text
5,       -4+3i,       4+3i,
```

whose radius squared is `25`, not the anchor-only ordinary lcm `5`.
The missing factor comes from opposite Gaussian orientations.

## 3. Exact angular scale and its invariances

Suppose now that all `X_i>0`, and put `A=min_i X_i`. The arguments of
the displayed points lie in the arc from zero to
`2 arctan(L/A)<pi`, with both ends attained. Consequently its shortest
containing arc, on the primitive circle, has normalized length

```text
C_* = arc length/sqrt(R)
    = 2 N^(1/4) arctan(L/A).                           (3)
```

Multiplying every `X_i` and `L` by a common positive integer leaves
each `n_e`, `N`, `A/L`, and `C_*` unchanged. Comparing `R/A^2` alone
would fail this elementary invariance. The relevant asymptotic ratio
as `A/L` tends to infinity is

```text
C_*^4 ~ 16 N/(A/L)^4.                                 (4)
```

Swapping the anchor with the point `A`, followed by reflection, sends
the other finite coordinates to

```text
X -> A+(A^2+L^2)/(X-A),
infinity -> A.
```

Every new coordinate is integral by the original divisibility.
This transformation preserves the multiset of all edge cotangents up
to signs, and therefore `N` and `L`. Since the other coordinates are
strictly greater than `A`, their images are also strictly greater
than `A`; the new minimum is again `A`. Hence this particular
complement symmetry preserves the exact quantity in (3).

## 4. A height inequality equivalent to the uniform goal

The following is the remaining **unproved integer target**:

> There exist an integer `k>=1` and a constant `c>0` such that, for
> every `L>0` and every `k` distinct integers `X_i>=L` satisfying
> the pair divisibilities, the all-edge lcm in (1) satisfies
>
> `N >= c (min_i X_i/L)^4`.                           (H)

It is equivalent to the requested existence of a point-count bound
depending only on `C` for every arc of length `C sqrt(R)`.

**From (H) to uniformity.** If `A>=L`, the elementary inequality

```text
arctan u >= (pi/4)u       (0<=u<=1)
```

and (3) give `C_* >= (pi/2)c^(1/4)`. Put

```text
c_0=(pi/4) min(1,c^(1/4)).
```

An arc of length at most `c_0 sqrt(R)` cannot contain `k+1` points.
To verify every normalization in this assertion, divide any proposed
integral tuple by its Gaussian gcd. Its normalized arc length can
only decrease, and its primitive radius is at least one. Its angular
span is then at most `c_0<=pi/4`. Anchor at one endpoint. The half-angle
cotangents of all remaining ratios are positive rational numbers at
least one. Clear their denominators and the denominators of every
pair half-angle cotangent to obtain exactly the integer system above
with `X_i>=L`. Reanchoring preserves the primitive least radius.
Inequality (H) would force normalized length at least `2c_0`, a
contradiction.

For `C>0`, partitioning an arbitrary arc into at most `ceil(C/c_0)` pieces of
that length therefore gives the uniform bound

```text
M(C) <= k ceil(C/c_0).                                (5)
```

Radii below one contain no nonzero lattice point; the radius-zero
case, if included, has just the origin and can be handled separately.

**From uniformity to (H).** Suppose the radius-independent bound at
`C=1` is an integer `M_1`, and take `k=M_1`. A clique with `k` finite
coordinates and its anchor has `M_1+1` points on the primitive circle,
so (3) must exceed one. Since `arctan(L/A)<=L/A`, it follows that

```text
N > (1/16)(A/L)^4.
```

Thus (H) follows with `c=1/16`. Neither direction assumes that the
common residue parameter `L` is bounded.

## 5. An explicit obstruction at three finite coordinates

The established
[four-point Pell family](four_point_bonus_counterexample.md) supplies
the integer cotangent triple

```text
L=24,
X=(90UV+30V^2+3, 96UV, 160V^2+128UV+16),
U+V sqrt(5)=(9+4 sqrt(5))^n,       n=1 mod 10.
```

Here `A=96UV`. Its exact all-edge lcm is the primitive squared radius
`N=5abc` from that family, with

```text
a=1+10V^2+2UV,
b=1+10V^2-2UV,
c=13+130V^2+38UV.
```

As `V` tends to infinity, `N` has order `V^6` while `(A/L)^4` has
order `V^8`. Hence (H) is false for `k=3`, even with fixed `L`.
This is an exact existing circle family, not a formal divisor profile.

The later [five-point cyclotomic subendpoint family](cyclotomic_five_row_subendpoint_family.md)
also disproves (H) for `k=4`. Its normalized containing arc length tends
to zero. Anchor at an endpoint and clear all cotangent denominators
to obtain four finite coordinates as above. The primitive radius is
at least one, so the angular span tends to zero and eventually `A>=L`.
If (H) held for `k=4`, equations (3)--(4), or the arctangent inequality
in Section 4, would give a fixed positive lower bound on that normalized
length. This is a contradiction. Thus any affirmative universal (H)
must use at least five finite coordinates. This consequence does not
affect the equivalence with the full uniform theorem.

The [endpoint-swap calculation](endpoint_swap_arithmetic_sensitivity.md)
now gives nonconformal images of the four-point Pell family at the same exponent:
every fixed positive rational alignment has a bounded normalized image
width, and one alignment has exact width below fifteen. A nearby moving
alignment instead has unbounded width. The result is specific to this
four-point source and does not prove (H) for larger cliques.

For comparison, the audited
[six-factor five-point family](six_factor_primitive_content_bound.md)
does satisfy (H) for its own four-cotangent realizations, with
`c=1/1600`. Indeed its proved lower bound `C_*>=1/sqrt(10)` and
`arctan(L/A)<=L/A` imply exactly that inequality. The family theorem
does not cover arbitrary four-cotangent cliques.

## 6. Verification and scope

[check_integer_cotangent_normalization.py](check_integer_cotangent_normalization.py)
checks the ideal-index gcd formula against an independent Gaussian
Euclidean algorithm, including its signed quotient. It compares the
ordinary all-edge lcm with the directly constructed primitive Gaussian
tuple and checks every possible reanchoring in each sampled clique.
It includes an anchor-only lcm counterexample, larger cliques, and the
Pell triple above. There were 13,599 signed pair checks, 1,863 clique
checks, and three exact Pell comparisons. A separate agent independently
audited both directions of the height equivalence and checked four Pell
instances. These finite checks supplement the proofs.

The local cardinality bounds in
[integer_cotangent_local_bounds.md](integer_cotangent_local_bounds.md)
depend on `L`. They do not prove (H). No estimate in this note removes
the dependence on radius from the general point-count bound.

The [intrinsic clearing-scale identity](intrinsic_cotangent_scale_triangle_divisibility.md)
now removes arbitrary enlargement of `L`: its minimal value is the lcm
of the reduced edge cotangent denominators. Their triple products equal
the triangle half-determinants divided by the exact Gaussian triangle
content and a parity factor. The resulting bound on the minimal scale
still grows with the primitive radius and does not prove (H).
