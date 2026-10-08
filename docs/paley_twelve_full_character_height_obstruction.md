# A full-character obstruction outside the Walsh family

The Walsh Fourier-gap argument does not extend to arbitrary Hadamard
matrices. In the normalized Paley matrix of order twelve, a balanced
four-row certificate already has normalized transform norm `4/3`, below
the Walsh gap `3/2`. Nevertheless, a different exact calculation shows
that five copies with at most two row-assigned flips per label have
**every** nonzero saturated integer character above primitive half-height.
This is a finite obstruction to the character method, not an endpoint
configuration or a theorem for all Paley orders.

## The exact matrix and finite transform classification

Let `q=11`, `M=12`, and let `chi` be the quadratic character modulo
eleven, with `chi(0)=0`. Index the first row and column separately;
the remaining rows and columns are indexed by `i,j in F_11`. Put

```text
H_(0,a)=H_(x,0)=1,
H_(i+1,j+1)=-1                  if i=j,
H_(i+1,j+1)=-chi(i-j)           if i!=j.
```

Direct integer multiplication verifies `H H^T=12I`. Every nonconstant
column has zero row sum. For a zero-sum integer vector c, put

```text
T_a=sum_x c_x H_(x,a),       L(c)=(1/M) sum_a |T_a|.
```

The all-ones column contributes zero. Fourier inversion is replaced here
by ordinary Hadamard inversion, giving `L(c)>=max_x|c_x|`.

There are 73,788 nonzero zero-sum vectors in `{0,1,-1}^12`, or
36,894 after identifying a vector with its negative. The
[exact checker](check_paley_twelve_full_character_height_obstruction.py)
enumerates all of them. It gives the following integer results:

| Support h | Minimum sum abs(T) | Number with sum abs(T)=12, up to sign | Minimum above 12 |
| --- | --- | --- | --- |
| 2 | 12 | 66 | none |
| 4 | 16 | 0 | 16 |
| 6 | 12 | 110 | 20 |
| 8 | 20 | 0 | 20 |
| 10 | 20 | 0 | 20 |
| 12 | 12 | 11 | 28 |

The norm-one vectors are exactly these explicitly generated vectors:

* differences of two row-coordinate unit vectors;
* `(H_a+H_b)/2` and `(H_a-H_b)/2` for distinct nonconstant columns;
* single nonconstant columns.

The enumeration compares the full set of norm-one vectors with this
generated set, not only their counts. In particular, it supplies both
the classification and the support-specific gaps used next.

## Five copies and arbitrary bounded-multiplicity flip assignments

Take `b=5` physical copies of each of the eleven nonconstant columns.
At every row flip one physical entry, using distinct physical columns.
Assume each original column label is assigned at most twice. There are
`r=55` physical primes, and each label retains at least three unflipped
copies. For a zero-sum integer c let

```text
v_j=(1/2) sum_x c_x s_xj,       B(c)=sum_j |v_j|.
```

These coefficients are integers. Each selected flip changes one
coefficient by a number of absolute value `|c_x|`, so

```text
B(c) >= (b/2) sum_a |T_a| - sum_x |c_x|.                  (1)
```

For amplitude `A=max|c_x|>=2`, Hadamard inversion gives

```text
B(c) >= M A (b/2-1)=18A>=36.                             (2)
```

For amplitude one and transform norm greater than one, the last column
of the table and (1) give bounds `36,44,42,40,58` at respective supports
`4,6,8,10,12`. These all exceed `55/2`.

For the classified norm-one vectors, let k count selected rows whose
assigned flip label is in the support of their nonzero transform
coefficients. At those rows the signs agree with the corresponding
baseline coefficient, so a matching flip lowers its magnitude by one.
A nonmatching selected flip raises a zero coefficient to one. Hence

```text
B(c)=bM/2+h-2k=30+h-2k.                                 (3)
```

A row-pair certificate has `h=2` and `k<=2`, so `B(c)>=28`.
A half-sum or half-difference of columns has `h=6`, only two transform
labels, and thus `k<=4` by the assignment multiplicity, again giving
`B(c)>=28`. A single-column certificate has `h=12` and `k<=2`, giving
`B(c)>=38`. Combining these cases with (2) proves

```text
B(c)>=28>55/2       for every nonzero zero-sum integer c.  (4)
```

This holds for **every** assignment with the stated multiplicity,
without searching through assignments. The canonical assignment
`a_x=1+(x mod 11)` attains 28, so the universal lower bound is sharp.

## Saturation, actual prime heights, and scope

The bound covers all rational-row characters with integral physical
coefficients. Comparing an unflipped copy of each label with the column
flipped at row x shows that integrality of `lambda^T S` forces
`2lambda_x` integral. Conversely every zero-sum integer c has
`c^T S/2` integral. Thus

```text
{lambda^T S: lambda rational, sum lambda=0} intersect Z^r
 = {c^T S/2: c integer, sum c=0}.
```

Choose actual distinct split rational primes with
`ell<=log p_j<ell+1/r`, `ell>1`, and put `W0=sum_j log p_j`.
Such clusters follow from the fixed-modulus prime count and pigeonhole
argument in the [nearby-prime construction](strict_obtuse_prime_box_countermodels.md).
Before flipping, every off-diagonal physical Gram equals `-5` at unit
weights; the two row-assigned flips can increase it by at most four.
The actual weighted Gram is therefore at most `-ell+1<0`.

For every nonzero saturated character, the associated conjugate-primitive
Gaussian composite has exact height `V=sum_j |v_j|log p_j`, so (4) gives

```text
V-W0/2 > 28ell-(55ell+1)/2=(ell-1)/2>0.                  (5)
```

All physical primes are distinct and every composite is nonunit by (4).
No claim is made about the actual rows' angular span. Common content is
not included in the primitive-height obstruction (5).

The example exposes a new limitation in extending the Walsh work:
general Hadamard cores need not have the Walsh transform gap or the
affine-subspace classification. The finite result here does not establish
an analogous obstruction at unbounded Paley orders, or exclude other
simultaneous arithmetic arguments. The uniform circle-arc goal remains open.
