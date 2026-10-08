# A degree-22 odd-column code beyond the short Newton certificate

This note records one finite boundary example.  It is a binary sign code,
not an odd-moment solution and not a Gaussian endpoint construction.

Take the normalized tensor Hadamard matrix \(H_{12}\otimes H_2\), where
\(H_{12}\) is the normalized Paley matrix of order twelve.  Number tensor
rows and columns in row-major order.  Select rows

```text
0, 1, 2, 4, 6, 9, 19, 23
```

and delete the constant column 0 and balanced column 1.  The retained
length-22 code has minority-size profile

```text
(n1,n2,n3,n4) = (2,3,2,15).
```

Thus its defect is \(D=2\cdot9+3\cdot4+2=32\), while its coefficient-one
excess is

```text
D + n4 - 2m = 32 + 15 - 44 = 3.
```

It has odd columns.  Its parity column is the deleted balanced
column, and restoring that column and the constant column recovers the
orthogonal order-24 row restriction.  In the padded-incidence notation this
has \(L=24\), \(K=0\), and therefore survives that check exactly.

The short weighted-Newton certificate does not exclude this code.  Define


\[
 \mathcal L=\{S^{\mathsf T}\mu:\mu\in\mathbf Q^8,
 \ \mathbf1^{\mathsf T}\mu=0,\ S^{\mathsf T}\mu\in\mathbf Z^{22}\}.
\]

There is no nonzero \(c\in\mathcal L\) with \(\lVert c\rVert_1\le10\), and
the bound is sharp: taking

```text
lambda = (-1,0,0,0,0,0,0,1), q = 2
```

gives \(c=\lambda^{\mathsf T}S/q\) with eleven entries equal to \(-1\) and
the rest zero.

Here is the complete finite proof of the lower bound.  Write
\(\mu=(x_0,\ldots,x_6,-\sum x_i)\), so \(c=B^{\mathsf T}x\), where row \(i\)
of \(B\) is row \(i\) of \(S\) minus row 7.  Seven retained coordinates,
at indices

```text
0, 1, 2, 3, 4, 6, 10
```

form an invertible coordinate map.  If their values are
\(y\in\mathbf Z^7\), rational inversion reconstructs every coordinate of
\(c\), with common denominator six.  Any \(c\) of mass at most ten must
have \(\lVert y\rVert_1\le10\).  There are exactly 433,905 such signed
integer vectors, including zero.  Exact rational reconstruction finds no
nonzero integral \(c\) of total mass at most ten.  Conversely, every
rational zero-sum \(\mu\) is represented: clearing its denominators gives
integer \(\lambda\) and positive integer scale \(q\).  Hence this enumeration is
exactly the certificate search, rather than a bounded search in
\(\lambda\).

The accompanying checker also scans all
\(\binom{24}{8}=735{,}471\) normalized row subsets.  There are 1,815
positive row subsets after the prescribed deletions.  Exactly 1,320 have
an odd column, and every one has full 24-column profile

```text
(n0,n1,n2,n3,n4) = (1,2,3,2,16).
```

This aggregate count is only for normalized row subsets of this tensor
matrix; it is not a classification of degree-22 codes.

The full moment system was subsequently excluded for all 21,120 balanced
deletions in the odd-positive family by the stronger
[critical root-order argument](tensor24_odd_column_family_critical_order.md).
The lattice minimum above remains a valid limitation of the shorter
linear-certificate method; it does not imply moment feasibility.

Run

```text
python3 docs/check_tensor24_odd_column_certificate_gap.py
```

to reconstruct the matrix, repeat the complete row-subset scan, verify the
padded identities, and audit the exact certificate lattice.
