# Reanchoring, prime flips, and frozen-divisor rigidity

## Status and scope

The sparse reconstruction theorem supplies exact coordinate charts on
primitive phase configurations. This note checks whether changing charts
or applying a symmetry to their unfrozen factors produces an actual
smaller-radius cluster. Three conclusions are proved:

1. Reanchoring only changes the coordinates of the same configuration.
2. Flipping a prime orientation outside the frozen blocks preserves the
   parent radius but forces an exponentially larger primitive residue.
3. A rational projective map preserving the oriented frozen divisors
   either is conformal or has coefficient height comparable to their
   total norm height.

The stronger theorem in
[projective_radius_inflation.md](projective_radius_inflation.md) removes
the frozen-divisor preservation hypothesis: every sufficiently small
nonisometric rational projective map inflates the least radius. None of
these results proves endpoint uniformity. Higher-height projective maps
and other arithmetic constructions remain outside their scope.

## 1. The common arithmetic window

Use the central truncation and simultaneous pair extraction of
[endpoint_central_truncation.md](endpoint_central_truncation.md). There
are eight selected rows, numbered `0,...,7`, inherited norm height
`W=log(R^2)`, central correction height `s=log D`, and 128 oriented
core blocks indexed by `S subset {1,...,7}`. Write these blocks as
`H_S`, suppressing the prime used for core blocks in that note.

They are pairwise coprime, coprime to each other's conjugates, and
individually conjugate-primitive. In the extracted asymptotic regime,

```text
log N(H_S)=w+o(w) uniformly in S,
w=(W-2s)/128,
s=o(w),                 W=128w+o(w).             (1)
```

The actual primitive anchor numerators satisfy

```text
h_i=x_i+i t_i,    h_i/bar(h_i)=z_i/z_0,
gcd_G(h_i,bar(h_i))=1,
0<|t_i|=exp(o(w)),
log|h_i|=W/4+o(w)=32w+o(w),
h_i=K_i product_(S containing i) H_S.            (2)
```

Here `log|K_i|<=s`; the Gaussian units of the original selected points
are equal. The `o(w)` precision in (2) uses the separate pair-moment
selection. Merely knowing an error `o(W)` without comparison to the
individual block scale is insufficient for the arguments below.

Freeze the 49 nonempty blocks of the optimal whole-cut anchor-majority
design from
[sparse_phase_lattice_reconstruction.md](sparse_phase_lattice_reconstruction.md),
and put

```text
F = this collection of 49 cuts,
Q_i=product_(S in F, i in S) H_S.
```

Each row occurs in 33 or 34 of these cuts. Consequently

```text
Q_i divides h_i,
log N(Q_i)>=33w+o(w).                            (3)
```

The exact finite design and its restricted optimality statement are in
[fano_freezing_audit.md](fano_freezing_audit.md).

## 2. Reanchoring is an exact coordinate change

Changing the anchor from row zero to row `a` gives the canonical primitive
transition

```text
h_ai = h_i bar(h_a)/N(gcd_G(h_i,h_a)),
h_ai/bar(h_ai)=z_i/z_a.                          (4)
```

This is the all-anchor identity proved in
[all_anchor_primitive_transition.md](all_anchor_primitive_transition.md).
It includes the full Gaussian cancellation rather than an uncharged
rational denominator.

On the core cuts, reanchoring has the following concrete effect. A cut
`S` is represented by the side not containing the current anchor. If
`a` is outside `S`, its representative and Gaussian block stay the same.
If `a` is inside `S`, replace the side by its complement in the full
eight-row set and replace `H_S` by `bar(H_S)`. All block norms stay the
same. The central exponents

```text
r_p=max_i min(a_ip,e_p-a_ip)
```

do not depend on the anchor, so neither the correction modulus `D` nor
the core radius changes. Reanchoring preserves every actual angular
distance and every actual point.

The least radius of the relative configuration is also anchor-invariant:
any integral realization with one anchor is already a realization with
any other selected point designated as its anchor. Taking the minimum
radius over these identical sets of realizations gives the assertion.

Thus overlapping reconstruction charts compose to the identity on their
common configurations. They do not by themselves construct a second
configuration on a smaller circle.

## 3. Unfrozen prime flips force a large residue

Choose a nonempty set of rational split primes belonging only to
unfrozen nonconstant core blocks. At each chosen prime interchange the
two conjugate Gaussian prime factors throughout every original row;
keep all Gaussian units fixed. This is an actual integral operation:

```text
a_ip maps to e_p-a_ip at every chosen prime.      (5)
```

All points keep the same norm. Every absolute allocation difference
`|a_ip-a_jp|` is unchanged. The new canonical primitive numerators
`h_i'=x_i'+i t_i'` therefore satisfy

```text
|h_i'|=|h_i|,          Q_i divides h_i'.           (6)
```

At least one chosen prime has a nonconstant selected core cut, hence
has a row `i` with nonzero signed allocation relative to row zero.
That prime's signed valuation in `h_i'` differs from the one in `h_i`.
For this row `h_i'` cannot equal `h_i` or `-h_i`.

Since both numbers have equal modulus, a zero real determinant would
force exactly one of these two equalities. Therefore

```text
0 != Im(h_i' bar(h_i)) is in N(Q_i) Z.
```

Taking absolute values and using (6) proves the exact inequality

```text
|t_i'| >= N(Q_i)/|h_i| - |t_i|.                  (7)
```

Combining (1)--(3) with (7) gives

```text
|t_i'| >= exp((1-o(1))w).                        (8)
```

This is an obstruction to the actual proposed flip, not only to an
attempted recovery of its factors. Let `Delta'` be the shortest angular
interval containing the flipped tuple, including its anchor. For the
row above, the angular distance of `h_i'/bar(h_i')` from one has half-angle
sine `|t_i'|/|h_i'|`. Hence

```text
Delta' >= 2|t_i'|/|h_i'|,
Delta' sqrt(R) >= exp((1-o(1))w).                (9)
```

So the fixed-radius flip loses the endpoint arc scale by an exponential
factor in one block's height.

Even removing the selected tuple's common Gaussian divisor afterward
does not rescue it. If that divisor is `g`, its norm height obeys

```text
log N(H_empty) <= log N(g)
              <= log N(H_empty)+2s=w+o(w).       (10)
```

Indeed, the core contributes exactly its all-common empty block; at
each prime the remaining contribution is at most `2r_p log p`.
The flip preserves the norm of this common divisor, because it only
interchanges its two conjugate valuations. Primitive normalization
reduces `sqrt(R)` by the factor `N(g)^(1/4)`. Consequently the normalized
arc length on the primitive circle still satisfies

```text
Delta' sqrt(R/|g|) >= exp((3/4-o(1))w).          (11)
```

Flipping an all-common cut has no relative effect, and is excluded by
the nonconstant-cut hypothesis. Flips touching frozen blocks are not
covered by this section.

## 4. Exact projective rigidity from frozen divisors

This statement needs no asymptotic assumptions. Let `h_i=x_i+i t_i` be
conjugate-primitive Gaussian integers, and let pairwise coprime,
nonunit, conjugate-primitive Gaussian blocks `H_S` divide all incident
`h_i`. Every frozen block must be incident to at least one nonanchor
row; an all-common empty block is therefore excluded.

Represent a rational projective map by a primitive integral matrix

```text
A=[[a,b],[c,d]],       det A != 0,
H=max(|a|,|b|,|c|,|d|)>=1.
```

After applying it to `(x_i,t_i)`, reanchor at the image `(a,c)` of the
old anchor `(1,0)`. Define

```text
q=a^2+c^2>0,        r=ab+cd,        s_A=ad-bc,
g_i=q x_i+(r+i s_A)t_i,
2g_i=(q+s_A-ir)h_i+(q-s_A+ir)bar(h_i).           (12)
```

Suppose a primitive output representative still contains every incident
oriented frozen divisor. Then every such divisor also divides the raw
`g_i`: removing coordinate content or a ramified factor can only remove
factors. As it divides `h_i` and is coprime to `bar(h_i)`, (12) forces

```text
H_S divides L(A)=q-s_A+i r.
```

Pairwise coprimality now gives the global exact restriction

```text
P_F=product_(S in F) H_S divides L(A).           (13)
```

There is a useful factorization

```text
L(A)=(a-i c)((a-d)+i(b+c)),      |L(A)|<=4H^2.   (14)
```

Since the first column is nonzero, `L(A)=0` is equivalent to
`a=d` and `b=-c`. This is precisely a conformal matrix; after reanchoring
it leaves every relative phase unchanged. Otherwise, (13)--(14) imply

```text
log H >= (1/4) sum_(S in F) log N(H_S)-log 2.    (15)
```

For the 49-block window of Section 1 this is

```text
log H >= (49/4-o(1))w.                          (16)
```

Thus a coefficient height `exp(o(w))` is incompatible with any
nonconformal projective change preserving these frozen oriented factors.
All rational matrix denominators were charged by clearing them and
choosing the primitive integral representative before defining `H`.

A reflection generally conjugates the frozen divisors instead of
preserving their orientation, and so need not satisfy the preservation
hypothesis. It preserves the least radius, just as a rotation does.
Equation (15) alone is a matrix-height obstruction; it does not claim
that every large matrix increases the circle radius.

## 5. What the projective radius theorem adds

The subsequent calculation in
[projective_radius_inflation.md](projective_radius_inflation.md) applies
to the same raw map (12) even when it changes all frozen data. In its
notation, `u=q+s_A-ir`, `v=q-s_A+ir`, and

```text
gcd_G(g_i,g_j) divides q u v t_ij.              (17)
```

The proof keeps prime-power valuations and all coordinate contents.
When `uv!=0`, old large shared factors have disappeared from (17).
For `m` nonanchor rows, primitive lcm reconstruction then gives

```text
log R_new=mW/4+O(m^2(E+log H+1)),
log(Delta_new sqrt(R_new))
   =(m-2)W/8+O(m^2(E+log H+1)),                 (18)
```

where `E` bounds the errors in `log|h_i|=W/4+O(E)` and the logarithms
of the old primitive residues. In the extracted regime the error in
(18) is `o(W)`. For `m>=3` the transformed least radius exceeds even
the original parent-circle radius, and its normalized arc length
diverges. The exceptions `uv=0` are exactly rotations and reflections.

Equations (13) and (17) describe two versions of the same obstruction:
a small nonconformal projective map cannot retain the large old common
divisors, and losing them makes the reconstructed radius larger. The
remaining question requires a construction beyond these small
projective changes or prime-flip chart symmetries.

## Verification

The root and uniformity-audit agents independently checked the flip
inequality, including the higher-prime-power common-divisor formula and
the factor `N(g)^(1/4)` in primitive radius normalization. The root also
checked the exact projective divisibility and factorization in Section 4.
The fresh-algebraic agent separately checked Sections 2--4, including
the primitive normalization exponent and the frozen-block height bound.

For the complementary projective radius theorem, an independent exact
Gaussian-integer enumeration checked (17) and the raw coordinate-content
bound in 13,920 nonexceptional matrix/pair cases: all invertible integral
matrices with entries in `{-1,0,1}`, and conjugate-primitive inputs
`x+i t` with `1<=x<=8` and `0<|t|<=4`. All checks passed. This finite
verification supplements the prime-valuation proof; it is not used as a
substitute for that proof.
