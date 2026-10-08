# Small full-support stresses satisfy the gradient core constraints

The exact gradient divisibility and the three infinitesimal `SL_2`
identities do not, by themselves, force a five-row gradient to vanish.
For every actual source tuple in the full-cut setting, they admit an
integer vector with all five coordinates nonzero and residual exponent
six, below the exponent eight permitted by the elementary gradient bound.

This is a conditional construction on any such actual tuple. It does
not construct an unbounded endpoint family, nor does it realize these
vectors as gradients of one small invariant. That last simultaneous
realization is an additional arithmetic condition that remains open.

## 1. The exact source and gradient constraints

Use the actual conjugate-primitive rows

```text
P_i=X_i+iY_i=K_i product_(U containing i) H_U,  i=1,...,5,
n_U=Norm(H_U),       log|K_i|<=sigma,
(1-eta)w<=log n_U<=(1+eta)w.
```

The nonempty cuts have disjoint odd split-prime supports, including
their conjugates. Put `Delta_ij=X_i Y_j-Y_i X_j`. The actual primitive
pair residues imply

```text
|Delta_ij|=b_ij product_(U containing i,j) n_U,
1<=b_ij<=exp(beta),       beta=4sigma+log T.             (1)
```

The `b_ij` are positive integers. This keeps the ordinary norm, the
Gaussian modulus, and the actual residue separate; no coprimality of
different `b_ij` is required.

The companion [gradient and discriminant calculation](five_row_gradient_discriminant_core_count.md)
proves that an integer invariant `Q`, of degree two in each row and
vanishing at this tuple, has integer gradient multipliers

```text
grad_i Q=lambda_i(-Y_i,X_i),
G_i | lambda_i,
G_i=product_U n_U^g_i(U),
g_i(U)=max(0,2|U\{i}|-4),                              (2)
```

and

```text
sum_i lambda_i X_i^2=sum_i lambda_i X_iY_i
                      =sum_i lambda_i Y_i^2=0.         (3)
```

The divisor in (2) has no correction-factor loss. One direct proof is
that the discriminant of `Q` in row `i` is a four-row invariant of
degree four in each remaining row and evaluates to `lambda_i^2`.
In the common determinant-one chart `(P,Y)`, every monomial has total
first-coordinate degree eight. At a cut with `a=|U\{i}|`, it therefore
contains at least `max(0,4a-8)` first coordinates from that cut.
Thus `H_U^max(0,4a-8)` divides the discriminant value. Since
`lambda_i` is an ordinary integer and the exponent is even, this
already gives the ordinary divisor `n_U^max(0,2a-4)` in (2).
Explicitly, if `pi^m` is a prime-power factor of `H_U`, the displayed
Gaussian divisibility says `2v_p(lambda_i)>=m max(0,4a-8)`;
an ordinary integer has the same valuation at `pi` and `bar(pi)`.
This gives the claimed power of `p^m=Norm(pi^m)` in `lambda_i`.
Disjoint supports allow multiplication over all cuts.

Writing `lambda_i=G_i ell_i`, the companion calculation gives the
necessary upper bound

```text
|ell_i|<=4 C_Q exp((8+72eta)w+4beta+sigma).             (4)
```

Here `C_Q` is the coefficient sum in the original polynomial
coordinates. The exact coefficient conversion and the two incident
derivative terms are included in the constant four. The following
construction only concerns (2)--(4), not the common discriminant
polynomials of a particular `Q`.

## 2. The normalized three-by-three minors

Form the integer matrix with five columns

```text
A_i=G_i (X_i^2, X_iY_i, Y_i^2)^t.
```

Every three columns are independent, since the five directions are
distinct. For `I={i,j,k}`, in increasing order,

```text
det A_I=G_i G_j G_k Delta_ij Delta_ik Delta_jk.          (5)
```

For a cut `U` of size `s`, set `a=|I intersect U|`. Its forced exponent
in (5) is

```text
d_U(I)=sum_(i in I) g_i(U)+binomial(a,2).
```

The minimum over the ten triples, and the remaining exponent, are

| `s` | Minimum `l_s` | When `d_U(I)-l_s=1` |
| --- | ---: | --- |
| 1 | 0 | never |
| 2 | 0 | `U` is contained in `I` |
| 3 | 3 | `|U intersect I|=1` |
| 4 | 9 | never |
| 5 | 15 | never |

All remaining exponents are zero. Consequently the positive integer

```text
L=product_(|U|=3) n_U^3
  product_(|U|=4) n_U^9 n_{12345}^15                  (6)
```

divides every three-by-three minor. More precisely,

```text
|det A_I|/L
 = b_ij b_ik b_jk
   product_(|U|=2, U subset I) n_U
   product_(|U|=3, |U intersect I|=1) n_U.             (7)
```

There are exactly three factors in each of the last two products.
Thus, with `B=max_I |det A_I|/L`,

```text
1<=B<=exp(6(1+eta)w+3beta).                            (8)
```

Canceling the exact common factor before estimating is essential:
the error in (8) is `6eta w`, rather than the sum of errors for the
raw minors and their common divisor. The fair logarithms of the raw
core minors and `L` are respectively `96w` and `90w`.

## 3. An integer kernel vector with every coordinate nonzero

Take the cofactor kernel vector on columns `1,2,3,4`, divide its
coordinates by `L`, and insert zero in position five. Call it `u`.
Do the same on columns `1,2,3,5`, with zero in position four, to
obtain `v`. Equations (5)--(6) show that both are integer vectors,
both lie in `ker A`, and every coordinate in their four-column
supports is nonzero. Also `||u||_infinity,||v||_infinity<=B`.

For each of positions one through three, at most one value of the
integer `t` makes `u_i+t v_i=0`. Among `t=1,2,3,4`, at least one
therefore avoids all three vanishings. Position four is always
`u_4!=0`, and position five is `t v_5!=0`. The resulting vector

```text
ell=u+t v in Z^5,
A ell=0,       all ell_i!=0,
max_i |ell_i|<=5 exp(6(1+eta)w+3beta)                 (9)
```

is the promised full-support certificate. Defining
`lambda_i=G_i ell_i` satisfies (2)--(3) exactly.

For fixed small errors, its residual exponent is six, below the
exponent eight in (4). This remains true with the displayed
correction and primitive-residue charges. In particular, proving that
all the actual gradient multipliers are nonzero, as in the
[ramification height gap](five_row_singular_point_height.md), does
not remove this obstruction.

The vectors in (9) have the ordinary affine-stress interpretation as
well. Set `c_i=lambda_i Norm(P_i)` and `z_i=P_i/bar(P_i)`. Then

```text
sum_i c_i=0,       sum_i c_i z_i=0.                    (10)
```

These are integer stresses of the five distinct points on the unit
circle. Common rotation and scaling of the circle preserve them.
Thus any sign restriction implied solely by circular convexity is
already satisfied by these exact vectors; it cannot exclude (9).

The missing input is a restriction on which of these small vectors
can occur simultaneously as the five row derivatives of a **single**
small invariant at its actual zero. The cofactor construction neither
supplies nor disproves such a restriction.

## Verification

[check_five_row_gradient_cofactor_stress.py](check_five_row_gradient_cofactor_stress.py)
checks all 310 cut/triple exponent pairs, the six-block formula,
the cofactor construction, and literal Gaussian source tuples both
without corrections and with shared correction primes. It checks
the actual Gaussian gcd residues, the ordinary stress identities,
and nonvanishing of all five constructed multipliers. These are
finite checks of the displayed identities, not a uniform arc proof.
