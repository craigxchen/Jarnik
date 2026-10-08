# Pair-selected and fixed alignments fail on one five-point endpoint family

There is a fixed family of five primitive circle points on
`5 sqrt(R)` arcs for which all three interior-pair-selected projective
maps have unbounded target normalized width. The two extreme source
points are kept as the endpoints. Each map preserves the four selected
points up to a conformal reflection, but its fifth point forces a larger
Gaussian denominator lcm. On the same source, every fixed positive
rational alignment except the identity also fails; a three-point target
subset already proves this second assertion.

This excludes the rule “choose a pair of actual interior source points”
as a universal endpoint-scale repair with these fixed endpoints. It does
not exclude other source-dependent parameters, other endpoint choices, or a
method restricted to configurations with more than five points. It gives
no new general point-count bound. It also does not exclude a theorem
restricted to sufficiently small source endpoint constants: this source
constant tends to `2 sqrt(5)`.

## 1. The all-row map and its exact residual denominator

Let `B=P+iQ`, where `P,Q>0`, and write each interior half-angle row as

```text
H_h=u_h+i v_h,      d_h=Q u_h-P v_h>0,
Q H_h=d_h+v_h B.
```

For distinct interior indices `i,j`, the pair-selected map, after a
conformal reflection, has the half-angle representatives

```text
1, B, F_h,
F_h=d_i d_j v_h+v_i v_j d_h B.                         (1)
```

Indeed the endpoint-fixing alignment scales the coefficient of `B` by
`lambda=d_i d_j/(v_i v_j Norm(B))`. Reflection sends the transverse
coordinate `w=d/v` to `Norm(B)/w`; the composition therefore sends it
to `d_i d_j/(v_i v_j w)`, giving (1). It exchanges the endpoints and
preserves the full intervening interval. In particular

```text
F_i=d_i v_i Q H_j,       F_j=d_j v_j Q H_i.             (2)
```

Thus the selected four-point subset is unchanged as a set after this
conformal adjustment. Conformal adjustment does not change the least
primitive squared radius: equivalently, use the ordinary lcm of all
relative-phase denominator norms.

For any half-angle representative `H`, let `A(H)` be its reduced
Gaussian phase denominator. If `Lambda_ij` is the Gaussian denominator
lcm of `1,B,H_i,H_j`, then the full target denominator is exactly

```text
Lambda'_ij=Lambda_ij lcm_G,h!=i,j
                 [A(F_h)/gcd_G(A(F_h),Lambda_ij)].     (3)
```

The same identity with the original `H_h` gives the source denominator.
There is no implication from equality of the selected subset radii to
control of the residual lcm in (3).

## 2. A genuine bounded-endpoint source

Use the [five-point cubic family](five_point_affine_shape_cubic_family.md),
with its parameter here called `t`, on positive multiples of `600`.
Its extreme points are `z_0,z_4`. Anchor at `z_0` and reflect the relative
phases so that all half-angle ratios are positive. Exact primitive
polynomial half-angle representatives over `Q(i)[t]` are

```text
H_0=1,
H_1=t^3+13t+12i,
H_2=t^3+31t+30i,
H_3=t^5+45t^3+444t+i(40t^2+720),
B=H_4=t^4+25t^2-36+60it.                              (4)
```

Their phase fractions are exactly the reflected source relative phases;
polynomial coordinate contents have already been removed. Here

```text
d_1=48(t^4+10t^2+9),
d_2=30(t^4+37t^2+36),
d_3=20(t^6+49t^4+504t^2+1296).                         (5)
```

These are positive, and so are all imaginary parts. Thus every interior
point lies strictly between `1` and `B`. The endpoint angle is

```text
Delta=2 arctan[60t/(t^4+25t^2-36)] ~ 120/t^3.
```

The known exact least squared source radius is

```text
N=product_(a=1)^6(t^2+a^2)/720^2 ~ t^12/720^2.
```

Consequently `Delta N^(1/4) -> 2 sqrt(5)`, and these are actual
bounded-endpoint sources.
The three selected alignment parameters tend respectively to `4`, `2`
and `1/2`; all three maps are nonconformal for sufficiently large `t`.

## 3. Exact polynomial denominator calculation

All gcds and lcms in this section are monic in `Q(i)[t]`. Work with
conjugates of the phase denominators, which leaves every radius norm
unchanged. The full source denominator lcm is

```text
L=t^6+it^5+45t^4+85it^3+404t^2+1164it-720.
```

For each pair `i,j` from `{1,2,3}`, the lcm of the selected subset
`1,B,H_i,H_j` is already `L`. Applying (1) to the remaining row and
removing its full polynomial coordinate content gives the following
monic, conjugate-coprime polynomial `K_ij`. The table also gives its
exact gcd `G_ij` with `L` and its residual factor `E_ij=K_ij/G_ij`.

| Pair | `K_ij` | `G_ij` | `E_ij` |
|---|---|---|---|
| `1,2` | `t^5+21t^3+40t+(20i/3)(t^2+4)` | `t+i` | `t^4-it^3+20t^2-(40i/3)t+80/3` |
| `1,3` | `t^5+29t^3+240t+20i(t^2+18)` | `t^2+5it-6` | `t^3-5it^2+10t-60i` |
| `2,3` | `t^5+(140/3)t^3+544t+(160i/3)t^2+960i` | `t^2-4it+12` | `t^3+4it^2+(56/3)t+80i` |

Thus the complete target polynomial denominator lcms are exactly

```text
L'_12=L E_12,       deg L'_12=10,
L'_13=L E_13,       deg L'_13=9,
L'_23=L E_23,       deg L'_23=9.                       (6)
```

This calculation retains the fifth row and every denominator overlap.
In particular it does not infer a target height from the selected
four-row subset.

## 4. Specialization retains the degree, including all Gaussian contents

For completeness, the polynomial degrees in (6) determine the true
primitive specialized radii up to fixed factors. Given the reduced
half-angle polynomials `J_s` and their monic Gaussian lcm `A`, the
polynomials

```text
Z_s=bar(A) J_s/bar(J_s)
```

have equal norm `A bar(A)` and polynomial Gaussian gcd one. At a
prime polynomial dividing `bar(A)`, a denominator attaining its largest
valuation has numerator coprime to that prime; hence the corresponding
`Z_s` has valuation zero. This proves the gcd assertion with all
multiplicities. The anchor row is included.

Clear the finitely many rational coefficients by one fixed integer.
A Bezout identity for the `Z_s` then bounds their common Gaussian
divisor at every integer specialization by a fixed nonzero Gaussian
integer. Dividing that common content therefore changes the squared
radius by only a bounded factor. This includes ordinary contents,
ramified primes, shared split primes and exceptional specialization
gcds; no coprimality assumption at each value of `t` is needed.

It follows from (6) that the **least primitive** target squared radii
obey

```text
N'_12=Theta(t^20),       N'_13=Theta(t^18),
N'_23=Theta(t^18).
```

All maps preserve the exact endpoint angle `Delta`. Their normalized
widths consequently satisfy

```text
C'_12=Theta(t^2),       C'_13=Theta(t^(3/2)),
C'_23=Theta(t^(3/2)).                                 (7)
```

Every pair-selected candidate therefore fails, even if the pair is
chosen adaptively at each source. Equality of the selected four-point
radius does not provide the missing all-row height control.

## 5. Every fixed nonidentity positive rational alignment also fails

Fix a positive rational `lambda!=1`, independently of `t`. The
endpoint-fixing alignment has representatives

```text
1, B, C_h=d_h+lambda v_h B.
```

It is conformally equivalent to the corresponding endpoint-swapping
alignment. Consider only the target subset `1,B,C_3`. Write

```text
d_3=20(t^2+4)(t^2+9)(t^2+36),
v_3=40(t^2+18),
B=(t+i)(t+2i)(t+3i)(t-6i).
```

The polynomial `C_3` has degree six, with leading coefficient
`20+40lambda>0`, and is coprime to its conjugate. To see the latter,

```text
C_3-bar(C_3)=120i lambda t v_3.
```

A common root would therefore have `t=0` or `t^2=-18`. But
`C_3(0)=25920(1-lambda)!=0`, while at `t^2=-18` one has
`C_3=d_3=45360!=0`. Moreover the exact monic anchor gcd is

```text
gcd(B,C_3)=gcd(B,d_3)=(t+2i)(t+3i)(t-6i),
```

of degree three, for every `lambda`. The complete polynomial denominator
lcm for this three-point subset consequently has degree

```text
4+6-3=7.
```

The specialization argument of Section 4 gives its true primitive
squared radius `Theta_lambda(t^14)`. Its squared radius divides the full target
squared radius. Since the endpoints are retained, the full target
normalized width is at least `c_lambda sqrt(t)` for all sufficiently
large `t`, where `c_lambda>0` is fixed. Thus no fixed positive rational
alignment other than `lambda=1` preserves bounded endpoint scale on
this source family. The identity case has the original bounded width.

This is a complete success/failure classification for fixed positive
rational parameters on this family, without asserting exact full target
degrees for every parameter. The constants depend on the fixed parameter;
the assertion gives no bound for parameters varying with `t`, even when
they have a finite positive nonidentity limit.

For arbitrary varying parameters, the
[adaptive height bound](endpoint_adaptive_alignment_height_floor.md)
instead gives a necessary condition: bounded full target width requires
`p+q` at least a constant times `t^(2/5)` for reduced `p/q!=1`.
This does not exclude parameters above that height.

## Verification

The [exact checker](check_endpoint_pair_selected_five_point_obstruction.py)
derives (4) from the original Gaussian point polynomials, verifies all
three raw maps and polynomial content reductions, checks the full source
and target lcms and their degrees, and verifies polynomial gcd one for
the complete integral circle tuples. It also checks four integer
specializations using independent all-edge and Gaussian-LCM primitive
radius formulas. For the fixed-parameter classification it verifies the
factorizations, both exclusion identities, the exact anchor gcd and
four sample rational-parameter lcms. The symbolic exclusion argument,
rather than the samples, covers every positive rational nonidentity
parameter. The specialization argument above, rather than the finite
checks, proves the growth statements.
