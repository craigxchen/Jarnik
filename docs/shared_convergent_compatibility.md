# Shared modular roots and the cost of Euclidean descent

## Status

The last-convergent parametrization in
[sparse_convergent_parametrization.md](sparse_convergent_parametrization.md)
is exact. This note examines its cross-row arithmetic and a proposed
Euclidean height descent.

The shared-root conditions are exactly the normalized determinant
divisibilities; general CRT adds no further compatibility equation.
Removing the small continued-fraction prefix decreases height by only
the logarithm of the tiny residue. Crossing the following large partial
quotient can change the height scale, but its integer transformation
already has large coefficients. The distinction is quantitative.

No uniform endpoint bound or new positive residue exponent is proved.

## 1. Shared roots are equivalent to normalized determinants

Let `P_i` be conjugate-primitive Gaussian divisors, with
`gcd(P_i,bar(P_j))=1` for every pair of indices. This is the
orientation condition supplied by independent fixed core blocks.
Write

```text
Q_i=N(P_i),
r_i^2=-1 mod Q_i,
P_i Z[i]={x+it:x=r_i t mod Q_i}.
```

Let `h_i=x_i+i t_i` be a conjugate-primitive multiple of `P_i`,
with `t_i!=0`. Then `gcd(t_i,Q_i)=1`, and hence

```text
r_i=x_i t_i^(-1) mod Q_i.                       (1)
```

Put `G_ij=gcd_G(P_i,P_j)`. The orientation assumption gives
`D_ij=N(G_ij)=gcd(Q_i,Q_j)`. The following conditions are
equivalent:

```text
r_i=r_j mod D_ij,
D_ij divides x_i t_j-x_j t_i.                   (2)
```

This follows immediately by multiplying (1) by the two invertible
denominators. In the actual fixed-block model the root equality is
already guaranteed by the common oriented divisor. Conversely,
the displayed determinant condition recovers that equality exactly.

General CRT says these pairwise compatibilities are necessary and
sufficient for a unique residue `r` modulo

```text
Q_*=lcm_i Q_i
```

having every prescribed `r_i` as a reduction. It automatically
satisfies `r^2=-1 mod Q_*`. The corresponding Gaussian ideal is
the least common multiple of the ideals `P_i Z[i]`, equivalently
the intersection of those ideals as sets. A generator is a
Gaussian lcm of the `P_i`, with norm `Q_*`.

Thus gluing all the oriented roots does not yield an additional
higher-order consistency equation. It encodes the Gaussian lcm
whose denominator cost is already retained by the exact
reconstruction theorem. This statement concerns the literal CRT
conditions; it does not rule out a theorem controlling the joint
distribution of their selected convergents.

## 2. Exact small-prefix transformation

Consider one row and write `P,Q,r` for its fixed divisor and root.
Suppose the sparse window has bounds `U,T>=1`, with

```text
2UT<|P|,      Q=|P|^2.
```

Choose the sign of its candidate so `t>0`. Let `k/t` be the
selected last convergent of `r/Q`, and `k_-/t_-` the preceding
convergent. For the initial convergent use `k_-=1,t_-=0`.
Set

```text
epsilon=k t_--k_- t in {1,-1},
x=r t-Q k,
A=t_- r-k_- Q.
```

The unimodular matrix whose columns are the convergent and its
predecessor has entries of absolute value at most `T`. Applying
its inverse to `(r,Q)` gives, up to the common sign `epsilon`,
the primitive integer pair `(A,-x)`. The exact identity is

```text
t A=epsilon Q+t_- x.                            (3)
```

Since `t_-<=t<=T`, `|x|<=|P|U`, and `T|P|U<Q/2`, equation
(3) gives

```text
Q/(2t) < |A| < 3Q/(2t),
|A|>|x|.                                       (4)
```

Non-strict inequalities in (4) would also suffice; the sparse
window is strict. The height of the transformed rational tail is
therefore

```text
log max(|A|,|x|)=log Q-log t+O(1).               (5)
```

In the extracted setting `log t=o(w)`, where `w` is an
individual block log norm. Thus stripping this entire small
continued-fraction prefix decreases logarithmic height by only
`o(w)`. Its small coefficients do not produce a descent by a
positive fraction of a block.

## 3. The first step that changes scale has large coefficients

Let `t_next` be the denominator of the next convergent. The
already proved exact gap gives

```text
t_next > |P|/U-T.
```

Consequently the unimodular matrix through that next convergent
has an entry at least `t_next`. Its logarithmic coefficient
height is at least

```text
log(|P|/U-T).                                  (6)
```

The next partial quotient also satisfies

```text
a_next > |P|/(UT)-2.                            (7)
```

For example, once `|P|/(UT)>=4`, its logarithmic size is at
least `log|P|-log U-log T-log 2`. Traversing this quotient
replaces the large tail numerator by a Euclidean remainder of
size at most `|x|`, so it can lower the scale. But the required
transformation is no longer a negligible-height coefficient
change: in the sparse block windows, `|P|/U` is exponentially
large on the block scale.

Equations (5)--(7) isolate the limitation of this descent route.
The available small-prefix transformation preserves a small
coefficient budget but gains only `log t+O(1)` in height. The
next transformation pays a large coefficient height. This is a
statement about these exact continued-fraction transformations,
not an impossibility claim for every possible descent method.

## 4. Projectivizing the tails discards the shared oriented root

Return to rows `i,j` and a common modulus `D=D_ij`. Modulo `D`,
the small-prefix transformed pair satisfies

```text
(A_i,-x_i)=r_i (t_(i,-),-t_i) mod D.             (8)
```

Here `r_i` and `t_i` are units modulo `D`. Thus its projective
reduction is simply

```text
[A_i:-x_i]=[t_(i,-):-t_i] mod D.                (9)
```

The shared root has disappeared as a unit scalar. Indeed,
already `[r_i:Q_i]=[1:0] mod D`; this projective reduction is
the same even for different unit roots. Accordingly one cannot
recover an extra oriented-root constraint by comparing only the
projective continued-fraction tails. Retaining their affine
scales restores the original equations (1)--(2).

This is the precise reduction issue in applying projective or
Euclidean reciprocity to the roots: the root information sits
in the unit scale at a modulus where the denominator vanishes.
The sparse archimedean gap does not make that unit scale small.

## 5. A common integer residue is available at negligible cost

There is a further exact arithmetic normalization that may be
useful for another route. Set

```text
ell=lcm_i |t_i|,
H_i=(ell/t_i)h_i=X_i+i ell.
```

All `X_i` are integers. The common divisor `G_ij` divides both
`H_i` and `H_j`, and their difference is real. Conjugate
coprimality of `G_ij` therefore gives

```text
N(G_ij) divides X_i-X_j.                         (10)
```

Moreover

```text
log ell <= sum_i log|t_i|.
```

In central truncation with `k` rows, this is
`O(k^2 W/sqrt(M))`, hence is negligible relative to every block
throughout `k=floor(c log_2 M)`, `c<1/2`.

This produces a common **integer** imaginary residue without
losing the per-block height budget. The new `H_i` need not be
primitive, so they must not be substituted into arguments that
require primitive coordinates without retaining their contents.
Equation (10) is another form of the existing normalized
determinant conditions; direct subtraction or a Euclidean step
on the `X_i` has not yielded a new global height lower bound.

The remaining issue is genuinely simultaneous: the unique tiny
convergents of many coherently shared Gaussian moduli must obey
the full cut-factor profile. Root gluing, small-prefix stripping,
and common-residue normalization preserve that problem, but do
not settle it.
