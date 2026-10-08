# Below-half-height Walsh characters must be coset characters

For the repeated Walsh model with at least four copies of each label and
one distinct physical-column flip per row, every integer row certificate
below the unweighted half-height threshold is a signed character on an
affine subspace. This classifies all integer certificates, including
ones with large row coefficients; it is not a support-size cutoff.

Combined with the archived nonlinear 32-row assignment, the classification
gives a concrete saturated character lattice whose every nonzero character
has actual prime-weighted height strictly above half the full primitive
height. It does not assert that those rows lie on a short arc.

## An elementary Fourier gap

Let `G=F_2^t`, `M=|G|`, and `f:G->{0,1,-1}` be nonzero. Normalize its
Fourier transform by

```text
fhat(a) = (1/M) sum_x f(x)(-1)^(a dot x),
L(f) = sum_a |fhat(a)|.
```

Exactly one of the following holds:

* `f` is a signed character on an affine subspace, and `L(f)=1`.
* `L(f)>=3/2`.

Here a signed character on `x_0+V` means
`f(x_0+v)=epsilon (-1)^ell(v)` on that coset and zero elsewhere.
The subspace and its character may initially be trivial.

To prove the gap, translate a support point to zero and multiply by a
global sign so that `f(0)=1`. Translation, sign change, and multiplication
by a character preserve `L(f)`. Fourier inversion gives
`sum_a fhat(a)=1`; if `N` is the total absolute negative Fourier mass,
then `L(f)=1+2N`.

If the support is not a subspace, there are support points `x,y` with
`x+y` outside the support. They are nonzero and distinct, hence linearly
independent. Multiply by a character taking the prescribed signs at
`x,y` to obtain `f(x)=f(y)=1` and `f(x+y)=0`. Partition the Fourier
coefficients by the two character values at `x,y`. The four signed
cell sums, by Fourier inversion, are

```text
(1+1+1+0)/4 = 3/4,
(1+1-1-0)/4 = 1/4,
(1-1+1-0)/4 = 1/4,
(1-1-1+0)/4 = -1/4.
```

Thus `N>=1/4` and `L(f)>=3/2`. If the support is a subspace but the
signs are not a character, choose `x,y` violating multiplicativity.
The same normalization gives `f(x+y)=-1`, and the last cell has sum
`-1/2`, an even stronger gap. Otherwise the support is a subspace and
the signs form a character. Direct Fourier inversion then gives norm
one. Undo the initial translation and sign change.

## Classification in the physical-column model

There are `b>=4` physical copies of every nonzero Walsh label, so
`r=b(M-1)`. Each row flips one distinct physical column; no physical
column is flipped twice. For an integer row certificate `c!=0` with
`sum_x c_x=0`, let

```text
v_j = (1/2) sum_x c_x s_xj,
A = max_x |c_x|,
H(c) = sum_(physical j) |v_j|.
```

All `v_j` are integers. The unflipped coefficient sum is
`bM L(c)/2`, with the same normalized Fourier norm as above. Flipping
the column assigned at row `x` changes its coefficient by a number of
absolute value `|c_x|`. Distinctness of the flipped columns therefore
gives

```text
H(c) >= bM L(c)/2 - sum_x |c_x|.                         (1)
```

Fourier inversion implies `L(c)>=A`. If `A>=2`, (1) gives

```text
H(c) >= M A (b/2-1) >= (b-2)M >= bM/2.                 (2)
```

If `A=1` and `c` is not a signed coset character, the Fourier gap gives

```text
H(c) >= (3b/4-1)M >= bM/2.                             (3)
```

Consequently `H(c)<bM/2` forces `c` to be a signed character on an
affine subspace. Its zero row sum forces its restriction character
`ell` to be nonzero, so the support size `h` is at least two.

For such a certificate on `P=x_0+V`, let `k` be the number of selected
rows whose assigned flip label restricts to `ell` on `V`. Before the
flips, exactly `M/h` labels have restriction `ell`, each with absolute
coefficient `h/2`. A matching selected flip lowers this coefficient
by one, while a nonmatching selected flip raises a zero coefficient
to one. Thus the exact physical coefficient sum is

```text
H(c) = bM/2+h-2k.                                      (4)
```

In particular, the strict primitive half-height condition is equivalent
to the majority-gain condition

```text
H(c)<r/2  if and only if  2k-h>b/2.                     (5)
```

Equations (2)--(5) concern all integer row certificates. They do not
assume their support is small, their coefficients are bounded in advance,
or the flip assignment is affine.

## A full saturated-lattice obstruction at M=32

Use five copies and the nonlinear assignment

```text
[28,21,11,30,24,23,13,3,12,7,25,4,31,1,22,3,
 19,18,10,9,2,16,17,27,26,6,15,5,14,8,20,29].
```

Every label occurs once, except label 3 which occurs twice, so every
label has at least three unflipped physical copies. The
[exact nonlinear coset audit](walsh_nonlinear_coset_height_obstruction.md)
enumerates all affine subspaces and nonzero restriction characters.
For support sizes `h=2,4,8,16,32`, the largest matching counts are
`k=2,3,4,3,2`. Formula (4) gives respective minima

```text
H(c)=78,78,80,90,108.
```

All other certificates have `H(c)>=bM/2=80` by (2)--(3).
Therefore the exact physical minimum over every nonzero zero-sum
integer certificate is `78`, strictly above `r/2=155/2`.

This also covers the entire saturated rational-row character lattice.
Indeed an unflipped copy of each Walsh label is present. If
`lambda` is rational with zero row sum and `lambda^T S` integral,
comparison of the unflipped label and the column flipped at row `x`
forces `2lambda_x` integral. Conversely every zero-sum integer `c`
gives integral `c^T S/2`. Thus the saturation is exactly

```text
{lambda^T S: lambda rational, sum lambda=0} intersect Z^r
 = {c^T S/2: c integer, sum c=0}.
```

For actual distinct split primes with
`ell<=log p_j<ell+1/r`, `ell>1`, put `W_0=sum_j log p_j`.
Every such nonzero character has exact conjugate-primitive Gaussian
height `V=sum_j |v_j|log p_j` and therefore

```text
V-W_0/2 > 78ell-(155ell+1)/2 = (ell-1)/2 > 0.           (6)
```

The [nearby-prime construction](strict_obtuse_prime_box_countermodels.md)
provides actual prime choices of this form; this is not a substitution
of arbitrary real weights for prime logarithms. The statement concerns
the primitive height `W_0`. Common content would increase the source
height and is not included in this obstruction. No endpoint arc
realization is claimed.

The [checker](check_walsh_all_character_half_height_classification.py)
exhausts every `{0,1,-1}` function on `F_2^3`, verifies the Fourier
gap and its equality class, and invokes the existing exhaustive
32-row coset audit. The proof above, rather than a bounded coefficient
search, covers all integer coefficients.
