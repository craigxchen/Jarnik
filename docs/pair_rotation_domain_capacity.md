# A pair rotation has bounded integral-domain capacity on an endpoint arc

Let `Z` be a finite set of Gaussian integers on the origin-centered
circle of radius `R>0`, contained in an arc of length `L=C sqrt(R)`.
Choose distinct anchors `a,b in Z`, and consider the rotation

\[
\rho=b/a.
\]

Write `rho=A/B` in reduced Gaussian integers, so
`gcd_G(A,B)=1` and `Norm(A)=Norm(B)=H`. Its domain inside the selected
arc is

\[
S_{a\to b}=\{z\in Z:\rho z\in\mathbb Z[i]\}.
\]

Then, without any primitivity assumption on `Z`,

\[
\boxed{S_{a\to b}=Z\cap B\mathbb Z[i],\qquad H\mid R^2,}\tag{1}
\]

and

\[
\boxed{|S_{a\to b}|\le1+\frac{L}{\sqrt{2H}}
\le1+\frac{L|b-a|}{2R}
<1+\frac{C^2}{2}.}\tag{2}
\]

In particular, `|S_(a->b)|<=ceil(C^2/2)`. If `C<=sqrt(2)`, the
rotation is integral on exactly its starting anchor among the rows
of `Z`. This counts all rows with an integral image; the image need
not belong to `Z`.

The result complements `pair_reflection_domain_capacity.md`. It
rules out forcing a single pair rotation to be integral on a positive
fraction of a growing endpoint cluster. It does not bound the cluster
itself or a number of separately chosen rotations growing with it.

## 1. Reduced denominator and exact origin-centered division

Choose `g=gcd_G(a,b)` and take `A=b/g`, `B=a/g`. Since the anchors
have equal norms,

\[
H=\operatorname{Norm}(A)=\operatorname{Norm}(B)
=R^2/\operatorname{Norm}(g).
\]

Coprimality gives `Az/B in Z[i]` if and only if `B|z`, proving
the domain assertion. Also `B|a`, so `H|R^2`. The complete
integrality domain on the full source circle is exactly

\[
B\{u\in\mathbb Z[i]:\operatorname{Norm}(u)=R^2/H\}.
\tag{3}
\]

Dividing only `S_(a->b)` by `B` therefore produces integer points
on an **origin-centered** circle of radius `R/sqrt(H)`, in an arc
of length at most `L/sqrt(H)`. There is no fractional-center
normalization in this argument.

## 2. The anchor chord pays for the denominator

Two distinct equal-norm Gaussian integers have squared separation

\[
|u-v|^2=2\bigl(\operatorname{Norm}(u)-u\cdot v\bigr),
\]

a positive even integer. Thus their distance is at least `sqrt(2)`.
In particular `A!=B` gives `|A-B|>=sqrt(2)`. The exact anchor identity

\[
|b-a|=\frac{R}{\sqrt H}|A-B|
\]

therefore yields

\[
\boxed{\sqrt H\ge\frac{\sqrt2 R}{|b-a|}.}\tag{4}
\]

Order the divided domain points along their containing arc. Each
successive chord has length at least `sqrt(2)`, and the corresponding
arc segment is at least that long. Consequently

\[
(|S_{a\to b}|-1)\sqrt2\le L/\sqrt H.
\]

Combining this with (4) proves the first two inequalities in (2).
The last is strict because the positive anchor arc gap `ell_ab`
satisfies `|b-a|<ell_ab<=L`.

Retaining that gap gives the more precise statement

\[
\boxed{|S_{a\to b}|<1+
\frac{C^2}{2}\frac{\ell_{ab}}L,\qquad
|S_{a\to b}|\le
\left\lceil\frac{C^2\ell_{ab}}{2L}\right\rceil.}\tag{5}
\]

Since the starting anchor belongs to the domain, (2) proves the
singleton assertion even at `C=sqrt(2)`. More generally the domain
is a singleton whenever `C^2 ell_ab/(2L)<=1`.

## 3. Arbitrary nested source levels

At a split rational prime `p=pi bar(pi)` with `p^e || R^2`, put
`t_z=v_pi(z)`, so `v_bar(pi)(z)=e-t_z`. The reduced denominator has
exponents

\[
v_\pi(B)=(t_a-t_b)_+,\qquad
v_{\bar\pi}(B)=(t_b-t_a)_+.
\]

Hence the exact local domain condition is

\[
\boxed{(t_a-t_b)_+\le t_z\le e-(t_b-t_a)_+.}\tag{6}
\]

Inert and ramified source factors cancel from `b/a`, since their
valuations are fixed by the common norm. Formula (6) needs no
squarefree or two-level assumption.

For binary source levels, (6) says that a domain row must agree
with the starting anchor `a` at every coordinate separating `a`
from `b`; coordinates where the anchors agree impose no condition.
For comparison, the pair-reflection domain requires agreement
with the common anchor value on every coordinate where the anchors
agree, with separating coordinates unrestricted.

## 4. Both domain-capacity bounds admit arbitrarily large binary profiles

These two kinds of small domains do not imply a finite code size.
For every even `m=2n>=4`, use one binary coordinate for every
`n`-element subset of the labels, and let row `i` have value one
precisely on subsets containing `i`. Every coordinate is balanced.

Fix different labels `a,b,z`. There is a balanced subset containing
`a` and excluding both `b,z`: add `n-1` other labels to `a`, which
is possible because `n-1<=m-3`. At this coordinate the anchors
`a,b` differ and `z` disagrees with the starting anchor. Thus no
third row is in the rotation domain; it is exactly `{a}`.

There is also a balanced subset containing `a,b` and excluding `z`:
add `n-2` other labels. Here the anchors agree and `z` disagrees.
Thus the reflection domain is exactly `{a,b}`. These statements
hold simultaneously for every pair, for arbitrarily large `m`.

This is a valuation-profile obstruction to deriving cardinality
from the two domain bounds alone, not an assertion that such a
profile has an endpoint realization. Actual endpoint phase accuracy
is still required. A bounded collection of `q` pair rotations
touches at most `q ceil(C^2/2)` source rows, and intersecting their
domains cannot increase the retained cardinality. An adaptive
operation using a growing number of rotations is outside this claim.

The independent checker `check_pair_rotation_domain_capacity.py`
verifies reduced Gaussian denominators, exact domains, full-circle
domain counts, the equal-norm separation and anchor identities,
the arc capacity inequalities, and the balanced binary examples.
It passes on 21,664 ordered anchor-pair domains in 1,000 integer-circle
arcs and on the balanced codes through ten labels. Arithmetic and domain
checks are exact; the arc-length inequalities use floating-point angles
with a stated numerical tolerance in the checker.
