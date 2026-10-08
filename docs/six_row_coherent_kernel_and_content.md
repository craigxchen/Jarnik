# The exact six-row coherent kernel and its selected-minor content

The entire six-row coherent cut matrix has rank four at every configuration
of distinct binary directions. In the established gradient basis, its
one-dimensional kernel is the matching-coordinate vector of those same
six directions. The primitive kernel therefore has height `18w+o(w)` in
the full core profile. Exact content in four selected rows improves the
crude `32w` lattice estimate to `18w`, but the actual-cut exclusion
threshold is only `8w+o(w)`.

The calculation below is symbolic. The polynomial certificate is
[check_six_row_coherent_kernel_and_content.py](check_six_row_coherent_kernel_and_content.py).
It verifies all fifteen kernel equations, all five selected cofactors,
and all sixty-four core cut orders, not merely ranks at sampled points.

## 1. Bases and the complete matrix

For six binary rows write `Delta_ij=X_iY_j-X_jY_i`, and use the
matching coordinates from [the Segre gradient calculation](segre_gradient_arithmetic.md):

```text
A=Delta_12 Delta_34 Delta_56,
B=Delta_12 Delta_36 Delta_45,
C=Delta_14 Delta_23 Delta_56,
D=Delta_16 Delta_23 Delta_45,
E=Delta_16 Delta_25 Delta_34,
Phi=BCE-AD(A+B+C+D+E).
```

The five polynomials `G=(G_A,G_B,G_C,G_D,G_E)=grad Phi` form the
fixed basis of `K_1(6)` used here. In particular,

```text
G_A=-D(2A+B+C+D+E),   G_B=CE-AD,   G_C=BE-AD,
G_D=-A(A+B+C+2D+E),   G_E=BC-AD.
```

For a two-element cut `S`, specialize its rows to `(0,1)` and factor
out one `P_j` from each of the four outside rows, as in
[the integral coherent circuit construction](coherent_cut_evaluation_differentials.md).
Let `R_S` be the resulting row of five numerical degree-one matching
values for these five basis polynomials. Thus `R` is a `15 by 5`
matrix, and a coefficient vector `c` satisfies `R c=0` exactly when
all its normalized coherent cut evaluations vanish.

The following table specifies every entry. For each cut let its outside
labels be `a<b<c<d`, and put

```text
U_S=Delta_ab Delta_cd,       V_S=Delta_ac Delta_bd.
R_S=U_S u_S+V_S v_S.
```

The vectors in the table are in the column order `A,B,C,D,E`.

| `S` | `u_S` | `v_S` |
|---|---|---|
| 12 | `(1,1,0,0,0)` | `(-1,0,0,0,0)` |
| 13 | `(-1,-1,-1,-1,0)` | `(1,0,1,0,0)` |
| 14 | `(1,1,1,1,1)` | `(0,0,-1,0,0)` |
| 15 | `(-1,0,-1,-1,-1)` | `(0,0,1,1,0)` |
| 16 | `(0,0,0,1,1)` | `(0,0,0,-1,0)` |
| 23 | `(0,0,1,1,0)` | `(0,0,-1,0,0)` |
| 24 | `(-1,0,-1,-1,-1)` | `(1,0,1,0,0)` |
| 25 | `(1,1,1,1,1)` | `(-1,-1,-1,-1,0)` |
| 26 | `(-1,-1,0,-1,-1)` | `(0,1,0,1,0)` |
| 34 | `(1,0,0,0,1)` | `(-1,0,0,0,0)` |
| 35 | `(-1,-1,0,-1,-1)` | `(1,1,0,0,0)` |
| 36 | `(1,1,1,1,1)` | `(0,-1,0,0,0)` |
| 45 | `(0,1,0,1,0)` | `(0,-1,0,0,0)` |
| 46 | `(-1,-1,-1,-1,0)` | `(1,1,0,0,0)` |
| 56 | `(1,0,1,0,0)` | `(-1,0,0,0,0)` |

Let `M=(A,B,C,D,E)^t` denote the matching values of the actual rows.
Direct Pluecker simplification gives the polynomial identity

```text
R M=0.                                                    (1)
```

The checker obtains `R` independently from the ten disjoint-triangle
products and their integral conversion to the gradient basis, and
verifies (1) in six independent affine variables. Each equation is
separately homogeneous in the binary rows, so the identities extend
to all binary directions.

## 2. Four rows prove rank four everywhere on the distinct locus

Use, in this order, the cuts `12,45,23,34`. Their rows simplify to

```text
R_12=(-B,A,0,0,0)/Delta_12,
R_45=(0,-D,0,B,0)/Delta_45,
R_23=(0,0,-D,C,0)/Delta_23,
R_34=(-E,0,0,0,A)/Delta_34.                              (2)
```

These displayed divisions are exact polynomial cancellations. Define

```text
F=Delta_12 Delta_56 Delta_36 Delta_16 Delta_45.
```

The signed vector of the five `4 by 4` minors of (2), with sign
`(-1)^j` when zero-based column `j` is deleted, is exactly

```text
c_tree=F M.                                               (3)
```

For example, deleting column `A` gives
`A^2 B D/(Delta_12 Delta_45 Delta_23 Delta_34)=F A`.
All factors are nonzero when the directions are distinct. Equations
(1)--(3) therefore prove, without a generic-position qualification,

```text
rank R=4,       ker_Q R=Q M.                              (4)
```

With `g=gcd(A,B,C,D,E)>0`, the primitive integer coefficient vector
in these coordinates is `M/g`, up to sign. The corresponding formal
polynomial is nonzero and vanishes at the actual configuration, since
Euler's identity gives `M dot grad Phi(M)=3Phi(M)=0`.

In the triangle basis `T_123,T_124,T_125,T_134,T_135`, its unnormalized
coefficient vector is

```text
(2A+B+C+D+E, -A-B-C-D, A+C, A+B, -A).                    (5)
```

The change of basis is integral with integral inverse, so it introduces
no varying denominator or projective-height loss. Other fixed coefficient
norms change logarithmic heights only by `O(1)`.

## 3. Exact content of the four selected normalized rows

Use the actual half-angle core factorization

```text
|Delta_ij|=b_ij product_(T containing i,j)n_T,
B_res=product_(i<j)b_ij,
log n_T=w+O(eta w).
```

The norm blocks have pairwise disjoint rational prime support. The
known matching-coordinate gcd satisfies

```text
G_match=product_T n_T^max(0,|T|-3),
G_match | g | G_match B_res.                             (6)
```

For a two-element cut put

```text
D_S=product_T n_T^max(0,|T minus S|-2).
```

Every entry of `R_S` is divisible by `D_S`. Divide each of the four
rows in (2) by its respective `D_S`. The signed cofactor vector is

```text
c_tree_normalized=q_tree (M/g),
q_tree=F g/(D_12 D_45 D_23 D_34) in Z minus {0}.          (7)
```

Integrality follows either from the integral minors and primitivity
of `M/g`, or from the primewise description below. In particular the
positive gcd of these five selected minors is exactly `|q_tree|`.

For a core cut `T`, the exponent in the compulsory part of (7) is

```text
e(T)=#{ij in {12,56,36,16,45}:i,j in T}
      +max(0,|T|-3)
      -sum_(S in {12,45,23,34})max(0,|T minus S|-2).      (8)
```

It is one on precisely the following fourteen subsets and zero on
the other fifty:

```text
12,16,36,45,56,
124,136,245,356,
1234,1236,1245,2345,3456.
```

Let `G_tree` be the product of their fourteen norm blocks. Because
the five factors of `F` use distinct edges, (6)--(8) give the exact
divisibility bounds

```text
G_tree | |q_tree| | G_tree B_res^2.                       (9)
```

Thus the selected normalized cofactor content has logarithm
`14w+O(eta w)+O(log B_res)`. This is a statement about the fixed
four independent rows in (2). It is not a claim that the gcd of
all `4 by 4` minors of the overdetermined `15 by 5` matrix has
the same content.

## 4. The exact kernel height and the remaining gap

Each matching coordinate has three brackets and raw logarithmic
size `48w+O(eta w)+O(log B_res)`. The sum of the core gcd exponents
in (6) is

```text
sum_(a=0)^6 binom(6,a)max(0,a-3)=30.
```

Consequently the primitive generator of the coherent kernel has

```text
log ||M/g||=18w+O(eta w)+O(log B_res)+O(1).              (10)
```

This agrees with the existing exact matching-height computation.
For comparison, each row of `R` has raw height `32w`, while its
mandatory row content has height `24w`. The four normalized row
heights therefore gave the crude exterior bound `32w`. Equation
(9) removes the remaining `14w` and gives exactly `18w`.

The short actual-cut nonvanishing theorem at `m=6,k=1` excludes a
nonzero coherent-kernel vector only below `8w+o(w)`. Thus the exact
minor content materially sharpens the lattice calculation but does
not reach that exclusion threshold. No uniform endpoint bound or
short-relation descent follows from this six-row calculation.
