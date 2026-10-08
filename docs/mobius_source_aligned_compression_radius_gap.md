# A source-aligned Möbius compression forces radius growth

A rational circle map whose expanding input singular axis is an **actual
source direction** has a precise angle–radius tradeoff. If the source
angular diameter is `Theta<pi/2`, its least transformed squared radius
`N'` satisfies

```text
N' >= kappa^2 cos^2(Theta/2)/Theta^2,                (1)
```

where `kappa>1` is the singular-value ratio of the same rational real
`2 x 2` half-angle map. The transformed angular diameter does not exceed
the source diameter. Thus a source endpoint arc with
`Theta<=C N^(-1/4)` gives

```text
N' >= (kappa^2/C^2) cos^2(Theta/2) N^(1/2).       (2)
```

In particular, choosing the critical compression `kappa` of order
`N^(1/2)` would require `N'` of order at least `N^(3/2)`. For `M>=14`
source points, `C<=2`, and unbounded primitive `N`, no nonconformal map
in this source-aligned class can have `N'<=N`: the established full-unit
endpoint condition theorem forces a larger `kappa` than (1) permits.
This excludes a natural source-derived descent choice, rather than an
arbitrary near-shear map whose singular axis misses the source arc.

## Exact source and transformed Gaussian conductors

Divide the **complete** common Gaussian gcd of the original equal-norm
rows first. Their relative phases and all projective maps remain the same;
the primitive source squared radius is `N`, and its normalized arc
constant only decreases. Choose the original source anchor and primitive
integer Gaussian half-angle columns `H_0=1,H_i=x_i+i y_i` with
`q_i=H_i/bar(H_i)`. Even-norm columns are allowed when row units differ.
For an integral representative `M` of the rational projective map put

```text
r_i = gcd_Z((M H_i)_1,(M H_i)_2),    B_i=M H_i/r_i.
```

The **least** transformed squared radius, retaining every row and the
transformed anchor, is the **ordinary integer lcm of the norms** of the
reduced Gaussian denominators of
`(B_i/bar B_i)/(B_j/bar B_j)` over all pairs. Taking a Gaussian lcm of
both pair orientations would double-charge conjugate factors. Equivalently,
the Gaussian lcm of the reduced denominators relative to **one fixed**
anchor has norm `N'`. At each odd
split prime `p=pi bar pi`, put

```text
s_i'(p)=v_pi(B_i)-v_barpi(B_i),
v_p(N')=max_i s_i'(p)-min_i s_i'(p).                 (3)
```

No inert or ramified prime contributes to a primitive squared circle
radius; parity only changes unit representatives of the phases. To see
that (3) includes a literal Gaussian anchor, choose row `j`, set
`t_i=s_i'-s_j'`, and assign the target anchor valuations
`v_pi(Z_j')=-min_i t_i`, `v_barpi(Z_j')=max_i t_i` at every split prime.
Then every `Z_j'(q_i'/q_j')` is Gaussian integral, and removing any factor
from this anchor would fail for a row attaining one extreme. This
construction proves both sufficiency and minimality of (3), including
ordinary contents `r_i` and any new primes in the images.

For a literal full-cut source with disjoint conjugate-primitive blocks,
the prime of cut `S` has the exact source signed profile

```text
s_i(p)=e_p(1_(i in S)-1_(0 in S)).                   (4)
```

If source correction factors are present, add their actual signed
valuations to (4). Writing the same real map as
`M z=alpha z+beta bar z`, `w=beta/alpha`, its phase action after a
common output rotation is

```text
q_i'=(q_i+w)/(bar w q_i+1),
s_i'(p)=v_pi(q_i+w)-v_pi(bar w q_i+1)
        up to a common shift independent of i.                 (5)
```

The numerator and denominator in (5) do not vanish on the unit circle
for a nonsingular nonconformal map (`Norm(w)!=1`). Formula (5) retains
endpoint cancellations and primes outside the source conductor. A
full-cut source need not remain full-cut: in an exact four-row fixture,
`M=diag(3,1)` erases **all seven old cut-prime widths**, while the least
target radius comes entirely from new primes and is much larger. That
fixture has no endpoint hypothesis and is only a conductor audit.

The particularly explicit source-derived map, for `H_j=x+i y`, is

```text
M_A = [ A x   A y ],       A integer >1,
      [ -y      x ],
det M_A=A Norm(H_j),       kappa(M_A)=A.                  (6)
```

It sends row `j` to the primitive output anchor `(1,0)`. Every other raw
image is

```text
(A dot(H_j,H_i), det(H_j,H_i)),
r_i=gcd_Z(A dot(H_j,H_i),det(H_j,H_i)).              (7)
```

Equations (3) and (7) give its exact least Gaussian denominator. The
integer gcd in (7) can erase old source cuts; it cannot erase the new
conductor without lowering the *actual* `N'` in (3).

## The geometric inequality and the endpoint corollary

Let `e` be the expanding right singular direction, equal projectively
to row `j`. In its orthonormal input coordinates, `M` is an output
orthogonal factor times `diag(kappa,1)` up to a common scalar. Since row
`j` belongs to an arc of diameter `Theta<pi/2`, every source
half-angle direction differs from `e` by `|theta_i|<=Theta/2<pi/4`
after choosing the short lift. The exact two-point determinant identity
gives

```text
|q_i'-q_j'|
  = |q_i-q_j|/sqrt(kappa^2 cos^2(theta_i)+sin^2(theta_i))
  <= Theta/(kappa cos(Theta/2)).                         (8)
```

Two distinct Gaussian integer points on a circle of squared radius
`N'` have normalized phase separation at least `1/sqrt(N')`, even
after the whole target tuple is divided by its common gcd. Comparing
with (8) proves (1). This lattice separation is why choosing a huge
Möbius matrix cannot bypass the radius cost through ordinary row gcds.

The lifted phase derivative is

```text
|d phi'/d phi|
 = kappa/(kappa^2 cos^2(phi/2)+sin^2(phi/2)) <=1           (9)
```

throughout that arc. The last inequality follows from
`cos^2(phi/2)>=1/2` and `(kappa^2+1)/2>=kappa` for `kappa>=1`.
Hence the image arc has diameter `Theta'<=Theta` and does not cross the
affine pole in the aligned half-angle chart. If `N'<=N`, its **actual**
normalized arc constant satisfies
`Theta' (N')^(1/4)<=Theta N^(1/4)<=C`. Thus the target endpoint
hypothesis required next is automatic for a proposed radius descent.

For `C<=2` and `N>=5`, `Theta<=2 N^(-1/4)<pi/2`.
The [full-unit endpoint condition theorem](mobius_endpoint_condition_number.md)
applies to the **same** nonconformal rational matrix `M`, with all `M`
source and target rows retained. Its singular-value ratio obeys

```text
kappa >= N^(z_M)/512,       z_M>1/4 for every M>=14,
z_14=139/552=1/4+1/552.                             (10)
```

The conformal identity and conjugate rotations have `kappa=1` and are
excluded from (10); they preserve least radius. From (1),
`cos(Theta/2)>=1/sqrt(2)` and `N'<=N` give
`kappa<=sqrt(2) C N^(1/4)`. Combining with (10),

```text
N^(z_M-1/4) <= 512 sqrt(2) C.                        (11)
```

The displayed `z_M` is the explicit cut exponent in the cited theorem.
Its global monotonicity can be checked directly: for `k=floor(M/4)`
and `D_M=(M-1)(M-2k-1)+1`, write
`q_M=M(M-1)/(2k D_M)` and
`z_M=(1-q_M)/2-M/(4F_M)`. At fixed `k`,
`q_M-q_(M+1)=M[(2k+1)(M-1)-2]/(2k D_M D_(M+1))>0`.
At the transition `4k+3 -> 4k+4`, the comparison after positive
denominators are cleared reduces to `2(k+1)^2>0`. Also `M/F_M`
decreases across both parities, since
`F_(2h)=h(h-1)` and `F_(2h+1)=h^2`.
Therefore `z_M` increases for `M>=8`, and in particular
`z_M>=z_14=139/552` for `M>=14`. Thus (11) bounds `N`
for this whole source-aligned class (for example
`N<=(512 sqrt(2) C)^552` uniformly for `M>=14`). The result makes no
claim that the singular axis must equal a source point. The following
extension locates the only surviving geometric regime.

## A necessary power approach of the expansion line

Let `d` be the projective angular distance in `[0,pi/2]` from the
**contracting** right singular line of a nonsingular real map to the
entire lifted source half-angle arc. This is the line where the circle
map has its largest local angular derivative, rather than an
affine-chart pole. For `d>0`, every source half-angle direction has
`|cos(theta)|>=sin(d)` in the expanding-axis coordinates. The exact
derivative therefore gives

```text
|d phi'/d phi| <= 1/(kappa sin^2(d)).               (12)
```

Order the `M` source rows along their short phase-arc lift. A
nonsingular projective circle map is monotone on this lift, even if
an affine cotangent chart passes through infinity. Its image lift has
`M-1` consecutive gaps. Each gap has circular arc length at least
the chord between two **distinct** primitive target lattice phases,
which is at least `1/sqrt(N')`. The image lift thus has length at least
`(M-1)/sqrt(N')`. Integrating (12) along the original arc gives the
exact finite-sample necessary condition

```text
kappa sin^2(d) <= Theta sqrt(N')/(M-1).             (13)
```

This argument uses the length of the traversed image lift, not its
shortest containing arc; it remains valid if the image crosses the
chosen angular seam. If `d=0`, (13) is vacuous, as expected.

If the source and full target both satisfy endpoint constants at most
`2` and `N'<=N`, combine (13) with the **same physical map** lower
condition bound (10):

```text
sin^2(d) <= [512 C/(M-1)] N^(1/4-z_M).            (14)
```

For `M>=14`, its right side tends to zero by a fixed power of `N`.
Thus any nonconformal endpoint radius-nonincreasing map in this range
must move its contracting singular line toward the source arc. A
source-aligned expanding axis has `d>=pi/2-Theta/2`, so (14) directly
recovers the exclusion above. Near-shear maps whose contracting line
lies in or very close to the arc remain possible under these
inequalities; their exact Gaussian conductor is the unresolved issue
in [the finite-sample critical-form audit](mobius_finite_sample_critical_form.md).

The [exact checker](check_mobius_source_aligned_compression.py) verifies
the Gaussian lcm, ordinary row contents, old-cut erasure and new-prime
radius on a literal full-cut source, and the determinant/angle and
finite-sample lift inequalities on an actual Pell endpoint source. The
Pell fixture has four rows and
does not assert the `M>=14` conclusion.
