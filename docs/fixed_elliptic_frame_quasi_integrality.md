# Fixed elliptic functions have full-height denominators

## Status

Let `E/Q` be fixed and let `f in Q(E)` be nonconstant.  For every
sequence `Q_n in E(Q)` with canonical height tending to infinity, the
denominator of the reduced rational number `f(Q_n)` has logarithm

```text
deg(f) hat h(Q_n)+o(hat h(Q_n)).                       (1)
```

This is a direct application of Lang's classical arithmetic-approximation
theorem for curves, together with standard height comparison on an elliptic
curve.  It is stronger than qualitative finiteness of integral values: it
excludes values whose denominators omit any fixed positive proportion of the
full height.

As a consequence, a fixed rational frame with a nonconstant entry cannot be
evaluated along a height-diverging rational sequence and cleared by common
integer denominators of logarithm `o(hat h)`.  This simplifies the fixed-map
case of the existing elliptic frame obstruction.  It gives no uniform result
for rational functions that move with the evaluation point and proves no
global lattice-arc bound.

## 1. The arithmetic approximation theorem

We use Theorem 2 of Serge Lang,
[Integral points on curves](https://www.numdam.org/articles/10.1007/BF02698777/),
*Publications Mathematiques de l'IHES* **6** (1960), 27--43, at printed
page 40.  In Lang's multiplicative-height notation, its number-field case
says the following.

> Let `C/K` be a curve of genus at least one, let `x in K(C)` be
> nonconstant, and fix a place `v` of `K`.  For every `rho,c>0`, the
> `K`-rational points satisfying
>
> ```text
> |x(P)|_v <= c/H_C(P)^rho                              (2)
> ```
>
> have bounded height.

Lang explicitly includes finite and archimedean places in his convention.
The theorem is the power-rate Diophantine approximation input in his proof
of Siegel's theorem.  The quantifier `rho>0` is what is needed below;
qualitative finiteness of integral points by itself would not give the same
conclusion.

Apply (2) with `K=Q`, `C=E`, the real place, and `x=1/f`.  Fix any ample
projective multiplicative height `H_E` on `E`.  Outside a bounded-height set,

```text
|f(Q)| < H_E(Q)^rho                                    (3)
```

for each fixed `rho>0`, whenever `f(Q)` is finite and nonzero.  Zeros of `f`
are harmless for the upper bound, and the finitely many rational poles have
bounded height.  An ample logarithmic height on the fixed elliptic curve is
`O(hat h(Q)+1)`.  Since `rho` in (3) can be arbitrarily small, every sequence
with `hat h(Q_n)->infinity` therefore satisfies

```text
log^+ |f(Q_n)|=o(hat h(Q_n)).                          (4)
```

This holds for arbitrary height-diverging sequences in `E(Q)`, not only for
multiples of one point.  In particular, if `P` is nontorsion, then (4)
applies to `Q_n=nP`, for which `hat h(nP)=n^2 hat h(P)`.

## 2. Exact denominator identity

For a finite value write

```text
f(Q)=a_Q/b_Q,       a_Q,b_Q in Z,       gcd(a_Q,b_Q)=1,
b_Q>0.                                                   (5)
```

With the ordinary absolute logarithmic height

```text
h(a/b)=log max(|a|,b),                                  (6)
```

there is an exact identity

```text
log b_Q=h(f(Q))-log^+|f(Q)|.                            (7)
```

Indeed, if `|a_Q|<=b_Q`, the subtracted term is zero and
`h(f(Q))=log b_Q`.  If `|a_Q|>b_Q`, subtracting
`log(|a_Q|/b_Q)` from `log|a_Q|` again gives `log b_Q`.

Let

```text
D=f^*(infinity),       d=deg D=deg(f).                  (8)
```

Functoriality of Weil heights identifies `h(f(Q))`, up to `O(1)`, with a
height for `D`.  The usual elliptic height comparison gives

```text
h(f(Q))=d hat h(Q)+O_f(sqrt(hat h(Q))+1).               (9)
```

For completeness, the square-root error follows directly from the divisor
class.  If `T` is the group-law sum of the points of `D`, counted with
multiplicity, then

```text
D ~ (d-1)[O]+[T].                                      (10)
```

A height for `[T]` is `hat h(Q-T)+O(1)`, while a height for `[O]` is
`hat h(Q)+O(1)`.  Hence

```text
h_D(Q)
 =(d-1)hat h(Q)+hat h(Q-T)+O_f(1)
 =d hat h(Q)-2 <Q,T>+O_f(1).                           (11)
```

Cauchy--Schwarz for the Neron--Tate pairing bounds the last pairing by
`O_f(sqrt(hat h(Q)))`, proving (9).  The normalization here is the one in
which a degree-`d` map has leading height `d hat h`.

Combining (4), (7), and (9) proves

```text
log b_(Q_n)=d hat h(Q_n)+o(hat h(Q_n)).                 (12)
```

Thus asymptotically all of the value height is nonarchimedean denominator
height.  Equivalently, for every `epsilon>0`, all sufficiently large points
in any height-diverging sequence satisfy

```text
log b_Q >= (d-epsilon)hat h(Q).                         (13)
```

## 3. Fixed-frame consequence

Let

```text
U=(f_1 f_2; f_3 f_4),       f_j in Q(E),                (14)
```

be a fixed rational matrix map, and let `Q_n in E(Q)` have
`hat h(Q_n)->infinity`.  Suppose a positive integer `q_n` clears all four
finite values:

```text
q_n f_j(Q_n) in Z       for every j.                    (15)
```

If some entry `f_j` is nonconstant of degree `d_j`, its reduced denominator
divides `q_n`.  Formula (12) gives

```text
log q_n >= d_j hat h(Q_n)+o(hat h(Q_n)).                (16)
```

In particular

```text
log q_n=o(hat h(Q_n))                                  (17)
```

is impossible.  This conclusion does not use contact divisibilities, a
determinant-one identity, or an upper bound for the cleared numerators.

If every entry is constant, there is no denominator obstruction of this
kind.  Constant determinant-one frames certainly exist.  In the contact
application they are excluded only after adding the unbounded contact
moduli: fixed norm values, even after multiplication by a common denominator
of logarithm `o(hat h)` and subpower error factors, cannot be divisible by a
modulus whose logarithm is comparable with `hat h`.  More explicitly, the
logarithm of such a cleared fixed norm is
`2log q_n+o(hat h)+O(1)=o(hat h)`.  Determinant one alone does not exclude
the constant case.

## 4. Relation to the existing fixed-map obstruction

The theorem in
[Fixed rational frames cannot carry the elliptic contact profile](fixed_rational_frame_contact_height.md)
assumes a fixed rational frame and common clearing denominators satisfying
`log q_n=o(hat h(Q_n))`.  Under that explicit hypothesis, (16) already
excludes every frame with a nonconstant entry.  The all-constant case is then
excluded by its unbounded contact divisibilities.  Thus the classical
quasi-integrality theorem gives a shorter proof of the stated fixed-map
exclusion.

The contact argument remains a stronger geometric diagnosis: it shows that
the twenty-five contacts force a common pole divisor of degree at least
twenty-five and records the resulting projective frame-height gap.  Formula
(12) instead detects the earlier incompatibility between fixed
nonconstancy and subpower denominator clearing, independently of where the
poles occur.

This observation also strengthens the qualitative use of Siegel in
[Four fixed norm directions and the finite-cover barrier](disjoint_elliptic_first_frame_barrier.md):
for a fixed function, not only are integral values finite, but denominators
smaller than the full degree-height scale by a fixed proportion occur only
at bounded height.

## 5. Scope

Fixedness is essential.  Lang's bounded-height exceptional set depends on
the function `f`, as do the divisor `D`, its degree, and every constant in
the height comparison.  For a sequence of functions `f_n`, one may move a
zero or pole with `Q_n`, or change coefficients so that the values have
small denominators.  No uniform conclusion follows from (12) without
separate bounds on the degrees and coefficient heights.

Consequently this note does not replace the bounded-degree moving-frame
argument, does not constrain arbitrary pointwise choices of frames, and
does not advance the global endpoint uniformity problem beyond the fixed-map
case.

For an explicit moving-map result, including arbitrary rational-point
sequences and its required logarithmic coefficient-height margin, see
[the moving elliptic denominator estimate](moving_map_archimedean_denominator_exclusion.md).
