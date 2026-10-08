# Eight-contact reconstruction and the frame resultant

The balanced K5 pencil of
[the star-pencil reduction](four_row_star_pencil_reduction.md) has a
particularly rigid missing ingredient: its degree-(4,3) complex frame.
For a *fixed* rational affine pencil with ten prescribed Gaussian-rational
contacts, a genuine frame is unique after monic normalization and is
automatically rational. Its existence is decided by one 8-by-8 linear
solve, two complex contact residuals, and seven real quadratic coefficient
tests. No further rational-root choice for the frame or the five private
roots is needed.

There is also an explicit 7-by-7 Sylvester resultant for constructing every
monic determinant-one frame in this degree chart. These are reconstruction
and denominator statements, not an exclusion of rational profiles or a
uniform lattice-arc bound. In particular, the certified real profile does
not provide the *rational input pencil* assumed in the descent statement.

## 1. The fixed input and the interpolation matrix

Let K be a subfield of R and normalize the affine pencil as

```text
A_a=t+alpha_a,       B_a=Y_a t+beta_a,       a=0,...,4.
```

Assume the Y_a are distinct and give ten distinct nonreal oriented points
z_ab in K(i), with their conjugates all distinct from these ten points.
For a<b impose the exact polynomial identities

```text
A_a B_b-B_a A_b = (Y_b-Y_a)(t-z_ab)(t-bar(z_ab)).       (1)
```

These are input conditions, not consequences of the interpolation below.
For an edge e=ab choose either endpoint a and define the row

```text
C_e = (A_a(z_e), A_a(z_e)z_e, ..., A_a(z_e)z_e^4,
       B_a(z_e), B_a(z_e)z_e, ..., B_a(z_e)z_e^3).     (2)
```

Equation (1) makes the contact equation independent of that choice:
A_a(z_e) is nonzero because alpha_a is real and z_e is nonreal, and
the two pencil columns are proportional at z_e. A frame

```text
h=sum_(j=0)^4 h_j t^j,        g=sum_(j=0)^3 g_j t^j
```

has every required contact exactly when its nine coefficients lie in
ker C. Throughout this note bar(h) conjugates coefficients and leaves t
fixed. Put

```text
W(t)=Im(bar(h)g)=(bar(h)g-h bar(g))/(2i).
```

## 2. Any eight contacts determine a genuine frame

**Uniqueness theorem.** Suppose h,g are coprime, deg h=4, deg g<=3,
and obey eight distinct contact rows (2). Then the kernel of those eight
rows is exactly the line K(i)(h,g), even if scalars are extended to C.

Indeed, a second solution h',g' makes hg'-h'g vanish at every one of
the eight contacts, since A_a(z_e) is nonzero. This polynomial has degree
at most seven, so it vanishes identically. Coprimality implies
h'=f h and g'=f g for a polynomial f. The degree-four bound on h'
forces f to be constant. Thus the eight rows have rank eight.

In particular, a frame with nonzero constant W is coprime: a common zero
of h and g would be a zero of W. For such a frame every choice of eight
of the ten contacts works. Since h_4 is nonzero, deleting the h_4 column
from these eight rows leaves an invertible 8-by-8 matrix. Consequently:

* The frame normalized by h_4=1 is reconstructed by Cramer's rule.
* If pencil and contacts are in Q and Q(i), the normalized frame is in
  Q(i)[t], even if its existence was first established over C.
* An identically singular selected 8-by-8 matrix rules out a genuine
  frame in this monic degree chart. It is not a request to try a different
  eight-contact subset.

The last assertion uses all the hypotheses: it does not apply to a
lower-degree h, a frame with a common factor, repeated contacts, or a
singular pencil direction.

## 3. A necessary and sufficient finite compatibility test

Select any eight rows, delete column h_4, and call the resulting matrix E.
Let v be the deleted column. If det E=0, reject this degree chart.
Otherwise solve

```text
E x=-v,       h_4=1.                                    (3)
```

Use the remaining two rows to test their two complex residuals exactly.
For 0<=m<=7 define the rational quadratic expressions

```text
w_m = sum_(j+k=m) (Re(h_j) Im(g_k)-Im(h_j) Re(g_k)).     (4)
```

The constant-frame tests are exactly

```text
w_1=...=w_7=0,          kappa=w_0 != 0.                  (5)
```

One can additionally require g_3!=0 for the exact (4,3) chart. These
are seven real quadratic equations in the reconstructed coefficients,
not fresh polynomial root-solving problems. Before elimination, the
contact equations are linear in those coefficients. Equivalently, write
the solution by the signed 8-by-8 cofactors of the selected eight rows;
the two remaining residuals are the corresponding 9-by-9 determinants.
Clearing the nonzero denominator det E and its conjugate in (4) gives
explicit polynomial compatibility equations in the input data.

When (3)--(5) hold, set

```text
Q_a=h A_a+g B_a.
```

Every Q_a is a monic quintic with all four incident contact roots.
Its remaining root is explicitly

```text
u_a=-h_3-alpha_a-Y_a g_3-sum_(b!=a) z_ab.               (6)
```

In particular u_a is Gaussian rational for rational input. One then
checks that all fifteen roots u_a,z_ab are nonreal, pairwise distinct,
and disjoint from their conjugates. These are separate nonvanishing
conditions. They are not entailed by (3)--(5). Likewise, rank four of the
constant column span is a separate minor condition if required.

The identity

```text
Im(bar(Q_a)Q_b)=kappa (Y_b-Y_a)(t-z_ab)(t-bar(z_ab))    (7)
```

is now automatic. Thus

```text
P_a=bar(Q_0) Q_a / ((t-z_0a)(t-bar(z_0a))),   a=1,...,4
```

are monic degree-eight polynomials with constant imaginary parts
kappa(Y_a-Y_0), and the separated fifteen factors give the full profile.
Conversely a monic degree-(4,3) frame for this input passes every test.
If determinant one is wanted, replace B_a by kappa B_a and g by
g/kappa. This is rational and preserves Q_a; no square root or norm
equation is introduced. If the original pencil is required to stay
fixed, one must instead test the additional equality kappa=1.

## 4. The explicit frame resultant

For a monic determinant-one frame write

```text
h=p+i q,       g=r+i s,       p monic of degree four.
```

The coefficient of t^7 in W=1 gives Im g_3=0. Hence, in the exact
(4,3) chart,

```text
deg r=3,       deg q<=3,       deg s<=2,
p s-r q=1.                                               (8)
```

Conversely, take *any* monic real quartic p and real cubic r with
Res(p,r)!=0. The seven unknown coefficients of q and s in (8) solve
the seven coefficient equations in degrees zero through six. Their
matrix has columns

```text
-r, -t r, -t^2 r, -t^3 r, p, t p, t^2 p,                (9)
```

so its determinant is the usual Sylvester resultant up to sign.
It is invertible and gives the unique q,s with the stated bounds.
This parametrizes this entire monic determinant-one frame chart by the
eight coefficients of p,r on the resultant-nonzero open set. In
particular it is a rational parametrization over Q, not a rational-point
problem for an additional conic.

If p,r are integral and p is monic, every denominator of q,s divides
Res(p,r), by the adjugate formula. More sharply,

```text
q,s are integral  <=>  Res(p,r)=+1 or -1.                (10)
```

The reverse implication is the integral adjugate formula. For the
forward implication, reduce (8) modulo p and multiply over its four
roots, counted with multiplicity. Monicity gives

```text
Res(p,r) Res(p,q)=1.
```

Both factors are integers, proving (10). Rational coefficients or
cleared nonmonic factors need not satisfy this unit-resultant condition.
It therefore cannot be used to exclude the rational Gaussian profiles
sought in the main problem; their fixed denominators can absorb a
nonunit resultant.

## 5. Scope and checks

The matrix descent eliminates a possible *second* arithmetic obstacle:
after a rational pencil and ten oriented rational contacts are found,
an existing genuine frame and the five private roots cannot require a
new number field. Finding such a compatible input remains unresolved.
The norm quadratics alone do not select oriented contacts, do not impose
the two last interpolation residuals, and do not impose (5).

The exact [checker](check_five_star_contact_frame_reconstruction.py)
verifies the reconstruction on rational data, rejects both a changed
ninth contact condition and a nonconstant-W interpolant, and verifies
the Sylvester construction and denominator phenomenon. Its interpolation
fixtures are ten independent scalar directions, not a rational K5
profile. It neither manufactures such a profile nor assigns rationality
to the stored real algebraic certificate.

The simpler recovery of a private root from the first two power sums is
already implicit in [the root-completion note](quartic_root_completion_and_quadratic_split_criterion.md).
It is not the new assertion here; the eight-contact uniqueness,
rational frame descent, and explicit frame resultant are.
