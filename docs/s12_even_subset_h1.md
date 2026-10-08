# The even-subset module has one nontrivial affine `S_12` torsor

Let

```text
V={even subsets of {0,...,11}}/<complement>.
```

This is a ten-dimensional `F_2` representation of `S_12`.  The exact
calculation below gives

```text
dim H^1(S_12,V)=1.                                  (1)
```

The nonzero class is the affine action on odd subsets modulo complement.
Consequently an affine group with full linear image `S_12` has only the
linear, odd-subset, and full-translation fiber patterns described below.

## 1. Coordinates and the Coxeter certificate

Represent a class in `V` by the unique even 12-bit vector whose bit 11 is
zero.  Bits 0 through 9 are independent, bit 10 is their parity, and bit 11
is zero.  Let `s_i=(i,i+1)`, for `0<=i<=10`, and let `A_i` be its induced
ten-dimensional matrix: permute the 12 bits and, if bit 11 becomes one,
xor the all-ones vector to return to the canonical representative.

Use the cocycle convention

```text
f(gh)=f(g)+g f(h),       X_i=f(s_i).
```

The standard type-`A_11` Coxeter presentation of `S_12` is complete:

```text
s_i^2=1,
s_i s_j=s_j s_i                    if |i-j|>1,
s_i s_(i+1) s_i=s_(i+1) s_i s_(i+1).
```

Therefore an assignment of the eleven `X_i` extends to a cocycle on
`S_12` if and only if it kills these relators.  Expanding them gives exactly

```text
(I+A_i)X_i=0,
(I+A_j)X_i+(I+A_i)X_j=0                         (|i-j|>1),
(I+A_i A_j+A_j)X_i+(A_i+I+A_j A_i)X_j=0         (j=i+1).  (2)
```

There are respectively `11*10`, `45*10`, and `10*10` scalar equations,
so (2) is a `660 x 110` matrix over `F_2`.  Order its columns by `(i,r)`,
meaning coordinate `r` of `X_i`, with column number `10i+r`.  Exact reduced
row elimination has rank 99.  Its eleven free columns are

```text
1,12,23,34,45,56,67,78,89,99,109.                  (3)
```

Equivalently they are

```text
(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),
(6,7),(7,8),(8,9),(9,9),(10,9).
```

Thus

```text
dim Z^1(S_12,V)=110-99=11.                          (4)
```

The checker constructs every row directly from (2), verifies the Coxeter
relations for the `A_i`, asserts (3), and hashes the 99 reduced rows.  With
each row stored as a 14-byte little-endian integer, their SHA-256 digest is

```text
f49378980e1deaebda82fd425ad983aff39f3aec87b253208e0f8e274f8ca68e.
```

This is a finite reproducible row-reduction certificate rather than a rank
reported by an external algebra package.

## 2. Coboundaries and the odd-subset class

For `b in V`, the coboundary is

```text
X_i=(A_i+I)b.                                        (5)
```

The ten coboundaries obtained from the coordinate basis of `V` are linearly
independent.  This also follows without computation: an element fixed by
every transposition has a representative `S` for which a transposition
changes `S` either by zero or by the all-ones vector.  A transposition
changes at most two bits, so the latter is impossible.  Thus every
transposition fixes `S`, making its 12-bit vector constant and its class
zero modulo complements.  Hence `V^(S_12)=0` and

```text
dim B^1(S_12,V)=10.                                  (6)
```

For the remaining class, regard odd subsets modulo complement as a torsor
under `V` and choose the base point `{11}`.  Every `s_i` with `i<10` fixes
the base point, while `s_10` moves it by the class of `{10,11}`.  Thus its
cocycle has

```text
X_i=0  (i<10),       X_10=[{10,11}].                (7)
```

In the canonical ten coordinates, the last value is the all-ones vector,
because the representative with bit 11 zero is the complement
`{0,...,9}`.  Substitution verifies all equations (2).  The eleven vectors
consisting of the ten basis coboundaries and (7) have rank 11.  Combining
this with (4) proves that every cocycle has the unique form

```text
X_i=(A_i+I)b+epsilon*1_(i=10)[{10,11}],
b in V,       epsilon in F_2,                        (8)
```

and proves (1).

## 3. Translation kernels and fiber degrees

The `S_12`-module `V` is simple.  Indeed, take a nonzero class and represent
it by a nonconstant even subset `S`.  Choose `a in S` and `b` outside `S`.
Then

```text
[S]+(a b)[S]=[{a,b}].
```

All pair classes are conjugate and generate the even-subset space, hence
their images generate `V`.  Every nonzero invariant subspace is therefore
all of `V`.

Let `G` be a subgroup of the affine group `V semidirect S_12` whose linear
projection is all of `S_12`.  Its translation kernel `G intersect V` is an
invariant subspace, so it is either zero or `V`.

If the kernel is zero, `G` is the graph of a cocycle.  Translation-conjugate
graphs differ by a coboundary, so (1) leaves two actions.  The zero class is
the linear action on even subsets modulo complement, with orbit sizes

```text
1, 66, 495, 462
```

corresponding to subset weights `0,2,4,6` (at weight six, complements are
identified).  The nonzero class is the action on odd subsets modulo
complement, with orbit sizes

```text
12, 220, 792
```

corresponding to weights `1,3,5`.  If the translation kernel is `V`, its
translations act transitively and the fiber degree is

```text
|V|=1024.                                            (9)
```

Thus the possible orbit-degree patterns are exactly

```text
linear:  {1,66,495,462},
odd:     {12,220,792},
full:    {1024}.
```

This classifies the affine fiber action under the hypothesis that the
linear image is the full `S_12`.  It does not assert that a particular
arithmetic cover has full linear image or determine which translation
kernel it realizes.

## Verification

Run `python3 docs/check_s12_even_subset_h1.py`.  Besides the rank certificate,
it enumerates all 1024 points for both complement classes and verifies the
displayed orbit sizes exactly.
