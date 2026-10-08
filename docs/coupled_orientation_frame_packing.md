# Coupled determinants and the number of short oriented frames

Fix the rational source vectors and rational contact norms in
[the oriented contact lattice](oriented_contact_lattice_rigidity.md).
Allow the Gaussian orientations of those same prime powers to vary.
Two determinant-one frames then satisfy a coupled determinant identity
that uses both agreeing and flipped orientations. It gives an exact
CRT condition, a small norm cofactor, and an exact theorem: if
`D>H^6`, there are at most two Gaussian-unit classes of frames with
Euclidean norm at most `H`, across all orientations.

The count is a count of frames for fixed source vectors and fixed
contact data, up to common Gaussian units. It is not a bound on the
number of points on the original circle. Arbitrary rotations or a
new median normalization generally change the primitive source vectors
and do not supply further frames for these fixed data.

## 1. Fixed unoriented data

For each contact `e`, fix a primitive integer coefficient vector
`lambda_e=(r_e,s_e)` and an odd positive integer `d_e` supported on
split rational primes. The integers `d_e` are pairwise coprime. Put

```text
D=product_e d_e.
```

An admissible frame is `x=(z,w) in Z[i]^2` with
`q(x)=Im(bar(z)w)=1`, such that for each `p^a || d_e`,
`r_e z+s_e w` is divisible by one of `pi^a,bar(pi)^a` above `p`.
Its orientation is allowed to depend on the prime power. Each prime
has a unique orientation for a given frame: the row is an integral
primitive vector because the frame lies in `SL_2(Z)`, so it cannot
vanish at both orientations of an odd split prime.

The coefficient vectors can be the repeated representatives of the
ten actual five-row clusters, or include the five private singleton
representatives as additional data. In the latter case their rational
norms must also be fixed when applying the results below.

## 2. The two determinants

Let `x=(z,w)` and `x'=(z',w')` be admissible frames, and define

```text
u=det_C(x,x')=z w'-w z',
v=det_C(x,bar(x'))=z bar(w')-w bar(z').
```

Direct expansion gives

```text
Norm(v)-Norm(u)=4q(x)q(x')=4.                          (1)
```

Partition the rational prime powers in `D` according to whether the
two frames agree or disagree on their orientation. Let `Delta_A` be
the product of the agreeing Gaussian prime powers, choosing the
orientation used by `x`, and let `Delta_F` be the product of the
flipped ones in that same reference orientation. Then

```text
Delta_A | u,       Delta_F | v,
D_A=Norm(Delta_A),  D_F=Norm(Delta_F),
D_A D_F=D,         gcd(D_A,D_F)=1.                    (2)
```

For agreement, the same real coefficient vector annihilates both
frames modulo the selected Gaussian prime power. For disagreement,
it annihilates `x` and `bar(x')` there. Their respective determinants
therefore vanish to the full prescribed exponent. This proof needs
the same source contact vector for the two frames; common rational
norms alone would not justify (2).

If the frames are not Gaussian-unit multiples, then `u!=0`, since
both are primitive and have `q=1`. Put `a=Norm(u)>0`. Equations
(1)--(2) imply

```text
a=0 mod D_A,       a=-4 mod D_F,
D | a(a+4).                                          (3)
```

The odd prime-power partition in (3) is exact: `gcd(a,a+4)` divides
four. Thus every odd prime in `D` must occur on precisely its assigned
side. No radical approximation is involved.

For a specified agreement/flipping pattern, (3) is one integer
residue class modulo `D`. If both frame norms are at most `H` and
`H^4<D`, it gives at most one possible value of `a`, since

```text
0<a<=H^4.
```

This is uniqueness of the determinant norm for that pattern, not
uniqueness of the frames from their determinant norm.

## 3. The norm cofactor and cross-orientation separation

There is an exact Gaussian-integer cofactor

```text
c=(u/Delta_A)(v/Delta_F) in Z[i],
a(a+4)=D Norm(c),
(a+2)^2-D Norm(c)=4.                                  (4)
```

This is a near-square equation with a sum-of-two-squares cofactor.
The cofactor is not asserted to be a square, fixed, or a unit; (4)
alone is not a contradiction or a fixed Pell-family classification.

Since the positive integer `a(a+4)` is at least `D`, every pair of
nonunit-equivalent admissible frames satisfies

```text
|det_C(x,x')| >= s(D),
s(D)=(sqrt(D+4)-2)^(1/2) ~ D^(1/4).                  (5)
```

Unlike the same-orientation bound `|det|>=sqrt(D)`, this bound holds
across all Gaussian orientation choices of the fixed rational data.

For the full `m`-row profile put `A_m=2^(m-3)` and suppose the common
frame-height bound is `H=exp(A_m w+o(w))`. If every paired nonconstant
cut is prescribed, including the private singleton/complement pairs,
then

```text
log D=(8A_m-2)w+o(w),
log Norm(c)<=2w+o(w).                                 (6)
```

If only the nonprivate paired contacts are prescribed, the corresponding
bounds are

```text
log D=(8A_m-2m-2)w+o(w),
log Norm(c)<=(2m+2)w+o(w).                            (7)
```

These follow from `a<=H^4` and (4). For five rows the two versions
have `log D=30w+o(w)` and `20w+o(w)`, respectively.

## 4. At most two frames when `D>H^6`

**Theorem.** For the fixed data of Section 1, if `H>=sqrt(2)` and
`D>H^6`, at most two Gaussian-unit classes of admissible frames have
`||(z,w)||_2<=H`. No lower bound on `|z|` is required.

Write the integer frame matrix and its Gram form as

```text
U=[[Re z,Re w],[Im z,Im w]],
F=U^t U=[[A,B],[B,C]],    AC-B^2=1.
```

The matrix is positive definite, and `A+C=||U||_F^2<=H^2`.
For every `p^e || d_j`, independently of the Gaussian orientation,

```text
A r_j^2+2B r_j s_j+C s_j^2=0 mod p^e.                 (15)
```

The evaluation vector `(r_j^2,2r_j s_j,s_j^2)` is primitive. It can
therefore be the last column of an integral unimodular matrix.
Multiplying any matrix of three coefficient rows `(A,B,C)` by that
matrix makes its last column divisible by `p^e`. Thus the determinant
of any three Gram coefficient rows is divisible by `D`, retaining
every prescribed exponent. But

```text
||(A,B,C)||_2 <= A+C <= H^2,
|det of any three coefficient rows| <= H^6 < D.       (16)
```

All such determinants vanish, so the Gram forms span at most a
2-dimensional real linear space. In dimension one there is only one
positive form with determinant one. In dimension two, choose one
positive form and use an `SL_2(R)` congruence to send it to the identity.
An orthogonal change diagonalizes a second independent symmetric form.
The whole span then consists of diagonal matrices, and its positive
determinant-one locus is

```text
F(s)=g^t diag(exp(s),exp(-s)) g,   det(g)=1.           (17)
```

Here `g` is fixed for this span. If `alpha,beta` are the squared norms
of its two rows, then `alpha beta>=1` by the determinant constraint.
Consequently

```text
tr F(s)=alpha exp(s)+beta exp(-s)
       =2 sqrt(alpha beta) cosh(s+s_0).
```

The trace bound confines all parameters `s` to an interval of length
at most `2R`, where `R=arcosh(H^2/2)`.

For two frames, apply the same right change of coordinates `U^(-1)`
to both Gaussian row vectors: they become `(1,i)` and
`(1,i)V`, where `V=U' U^(-1) in SL_2(R)`, while their complex
determinant is unchanged. Direct expansion gives
`Norm(det_C((1,i),(1,i)V))=||V||_F^2-2`, and cyclicity of trace gives
`||V||_F^2=tr(F^(-1)F')`; this identity holds for every pair of frames,
without a coplanarity assumption. In the parametrization (17), it gives

```text
tr(F(s)^(-1) F(t)) = Norm(det_C(x,x'))+2
                   = 2 cosh(s-t).                   (18)
```

The first equality holds for every pair of determinant-one frames,
without the plane hypothesis. Indeed `x=(1,i)U`; common right
multiplication by `U^(-1)` preserves the complex determinant and sends
the pair to `(1,i),(1,i)V`, where `V=U'U^(-1)`. If the entries of `V`
are `a,b,c,d`, its determinant norm is
`(b+c)^2+(d-a)^2=||V||_F^2-2`, while
`tr((U^tU)^(-1)U'^tU')=||V||_F^2`.

Equation (3) therefore separates distinct parameters by at least

```text
delta=arcosh(sqrt(D+4)/2)>arcosh(H^2/2)=R.             (19)
```

The strict inequality follows from `D>H^6>=H^4-4`.
Three distinct parameters would occupy length at least `2delta>2R`,
which is impossible. Finally equal Gram matrices mean
`U' U^(-1) in SO_2(Z)`, precisely multiplication by a common Gaussian
unit. This proves the stated count of frame classes.

With all paired nonconstant cut norms fixed, including the private
ones, the full profile has `D=exp((8A_m-2)w+o(w))` and
`H=exp(A_m w+o(w))`. The theorem applies for `A_m>1`, hence every
`m>=4`. For nonprivate contact data alone it applies when
`A_m>m+1`, hence every `m>=6`. For five rows it gives the constant
bound two with the private norms fixed. The error terms must be
uniform over all candidates counted at each `w`.

The constant two cannot be replaced by one for arbitrary fixed
contact systems. One example is

```text
x=(1+8i,1+9i),       x'=(-5+9i,1-2i),
H^2=147,
lambda_1=(1,11729),  pi_1=42-103i,  d_1=12373,
lambda_2=(1,5690),   pi_2=64+91i,   d_2=12377.
```

Both rational norms are prime. The frames agree on the first
orientation and disagree on the second; their determinant norms
are `12373` and `12377`, and `D=153140621>147^3`.
This example does not claim five-row profile weights or all ten
shared five-vector contacts. It only establishes sharpness for the
general fixed-contact theorem.

## 5. A finite packing bound below that threshold

Here is a bound with explicit finite hypotheses. Fix numbers `H>=h>0`
and consider admissible frames satisfying

```text
|z|>=h,       ||(z,w)||_2<=H.
```

Associate the complex slope `zeta=w/z`. The determinant-one equation
gives

```text
Im(zeta)=1/|z|^2,
|Re(zeta)|<=H/h.                                     (8)
```

The slope identifies a frame up to a Gaussian unit: equal slopes make
two primitive Gaussian frames proportional, and `q=1` makes the
proportionality factor a unit. For different slopes, (5) gives

```text
|zeta-zeta'|=|det_C(x,x')|/|z z'|>=s(D)/H^2.          (9)
```

If

```text
s(D)>=2(H/h)^2,                                      (10)
```

the difference of imaginary parts is at most `1/h^2`, hence at most
half the separation in (9). The real parts are consequently separated
by at least `s(D)/(2H^2)`. They lie in an interval of length `2H/h`.
The number of Gaussian-unit classes of admissible frames is therefore

```text
at most 1+4H^3/(h s(D)).                              (11)
```

This is one-dimensional packing because the determinant-one condition
places the slopes in a much thinner horizontal strip than their
pairwise separation. It is not a lattice-volume argument.

For an asymptotic counting statement, use one common height window
for all candidates at each `w`, with

```text
log h=A_m w+o(w),       log H=A_m w+o(w).              (12)
```

The error terms in this window are uniform over the candidates being
counted. Under all paired contact data, including private norms,
(10) holds for all sufficiently large `w`, and (11) becomes

```text
number of short frame classes <= exp(w/2+o(w)).        (13)
```

With the nonprivate contacts alone, it becomes

```text
number of short frame classes <= exp((m+1)w/2+o(w)),   (14)
```

provided `m>=4`, when (10) again holds eventually. For five rows this
is `exp(3w+o(w))`; adding fixed private norm data improves it to
`exp(w/2+o(w))`.

The rational vectors and rational norms may vary with `w`, but they
must be held fixed inside each set counted by (11)--(14). A different
primitive reanchoring or an arbitrary common rotation is not presumed
to produce another member of that same set. These bounds can coexist
with one exceptional short lift, and give no uniform lattice-arc
obstruction by themselves.

## 6. Exact arithmetic certificate

Run `python3 docs/check_coupled_orientation_frame_packing.py`. The
checker verifies the determinant/Gram identities, full-exponent CRT
constraints, Gaussian norm cofactor, and three-form divisibility.
It includes one fixed five-vector system with all ten median contact
assignments under several orientations, and the sharp two-frame
example above. These finite checks do not certify asymptotic profile
heights or substitute for the proof of the theorem.
