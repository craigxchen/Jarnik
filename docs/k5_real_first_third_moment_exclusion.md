# Real K5 first-and-third-moment exclusion

The later [zero-edge boundary argument](five_row_nine_factor_moment_exclusion.md)
extends this exclusion to one zero edge with the other nine absolute
values distinct. Together with parity extension, it excludes every
five-row first-and-third-moment code of length at most nine.

**Theorem.** Assign a nonzero real number `x_ij=x_ji` to every edge of the
complete graph on five vertices. Suppose the five vertex sums

```text
sum_(j!=i) x_ij
```

are equal, and the five vertex cube sums

```text
sum_(j!=i) x_ij^3
```

are equal. Then two edge numbers have the same absolute value.

This holds for arbitrary real ratios. It removes the prescribed-ratio
hypothesis of the earlier rational ideal computation. The proof below
uses a complete finite order classification and four integer polynomial
certificates; its [standalone checker](check_k5_real_moment_exclusion.py)
uses only the Python standard library.

This is a theorem about the `K5` moment model, not a uniform lattice-circle
bound. It does not classify other ten-column sign matrices, allow repeated
absolute values without further analysis, or justify replacing arbitrary
Gaussian factors by a fixed-coefficient asymptotic model. In particular,
the [ten-factor real witness](tenfactor_real_outercut_witness.md) has a
different column pattern, including an identical-column pair, and is not
contradicted.

## Exhaustive reduction to four orders

Assume, for contradiction, that all ten absolute values are distinct.
For any two vertices, subtract their sum and cube-sum equations. Their
common edge cancels, leaving six distinct-magnitude, nonzero signed numbers
whose first and third power sums vanish. The
[six-root successor theorem](critical_odd_moment_root_order.md) says that
their signs, in increasing absolute-value order, must be one of

```text
+--++-,   -++--+,   +++--+,   ---++-.
```

The sign in this test is the edge's actual sign multiplied by its
incidence difference for the vertex pair. All ten vertex pairs must use
the same global edge order and actual signs.

Permute the vertices so that the edge of smallest absolute value is `01`,
and negate every edge if necessary so that its sign is positive. These
operations preserve both moment hypotheses. A depth-first enumeration
then considers each remaining edge and both possible signs. A prefix is
removed exactly when, for some vertex pair, its nonzero incidence
subsequence is not a prefix of any of the four displayed patterns.
Every possible solution survives this necessary test.

The numbers of surviving prefixes at depths zero through ten are

```text
1, 1, 18, 108, 192, 216, 216, 192, 240, 48, 48.
```

The twelve vertex permutations preserving the unordered edge `01` divide
the 48 full orders into four orbits of size twelve. Representatives are:

| Orbit | Edges in increasing absolute-value order, with signs |
| --- | --- |
| 0 | `+01 -02 -03 +12 -24 -14 +34 -13 -23 -04` |
| 1 | `+01 -02 +12 -03 -24 -14 +34 -13 -23 -04` |
| 2 | `+01 +02 +03 -12 +24 +14 -34 +13 +23 +04` |
| 3 | `+01 +02 -12 +03 +24 +14 -34 +13 +23 +04` |

The checker regenerates the full search and its orbit quotient. Neither
the counts nor the exclusion relies on a bounded search over coefficient
magnitudes.

## Complete positive parametrizations of the linear equations

In one representative, let `r_0<...<r_9` be the positive absolute values.
Set `g_0=r_0` and `g_i=r_i-r_(i-1)` for `i>=1`, so every gap is strictly
positive. The four vertex-sum differences give an integer matrix equation
`B g=0`, where `B` has rank four.

The checker lists six explicit integer columns `R_0,...,R_5` for each
representative. It verifies

```text
B R_j=0,        rank(R)=6.
```

They therefore span the whole six-dimensional kernel, so every solution
has a unique expression `g=sum_j lambda_j R_j`. Each parameter is a literal
gap coordinate, as follows (all indices start at zero):

| Orbit | `(lambda_0,...,lambda_5)` |
| --- | --- |
| 0 | `(g_0,g_3,g_8,g_2,g_1,g_9)` |
| 1 | `(g_3,g_0,g_8,g_2,g_1,g_9)` |
| 2 | `(g_3,g_8,g_0,g_2,g_1,g_9)` |
| 3 | `(g_3,g_8,g_0,g_2,g_1,g_9)` |

In particular every `lambda_j` is strictly positive. This is an exact
parametrization of all first-moment solutions in the order cone, not an
assumption that they are positive combinations of numerically sampled
rays. The ranks, kernel identities, and coordinate inverses are checked
with exact arithmetic.

## Four positive polynomial contradictions

After this substitution let

```text
F_i(lambda)=sum_(e incident to i) x_e^3
            -sum_(e incident to 0) x_e^3,       i=1,2,3,4.
```

The cube-sum hypothesis would make all four cubics zero. In orbits zero
and two respectively, the following cubic combinations have only
nonnegative coefficients and are not zero:

```text
orbit 0: -F_3+3F_4,
orbit 2: -F_1-F_2+F_3-F_4.
```

For orbit one, let the four rows of the following matrix be the
coefficients of linear forms `L_i(lambda)`:

```text
[   0,   0,   0,   36,    36,   36 ]
[   0,   0,  36,    0,   346,   36 ]
[ -18, -36,  27, -364, -1266, -216 ]
[  18,  36, -27,   18,   920,  216 ]
```

Then `sum_i L_i F_i` is a nonzero quartic with nonnegative coefficients.
For orbit three, the same conclusion holds with matrix

```text
[    0, -2340, -11172, -3756, -5586,  -312 ]
[    0,  2028,   4308,  3444,  5598,  -312 ]
[  156, -1494,    208,  3276,  2871,  1872 ]
[ -156,  3834,   6032,  -156,  2079, -1872 ]
```

These are integer identities obtained by expanding the actual substituted
edge cubes. The checker reconstructs them from incidence, the gap columns,
and the displayed multiplier coefficients. The four positive polynomials
have respectively 49, 78, 54, and 95 nonzero monomials; their smallest
nonzero coefficients are respectively 3, 108, 3, and 624.

Since all six parameters are positive, each polynomial is strictly
positive. But each is a polynomial combination of the four equations
`F_i=0`, so it would have to vanish. This excludes every representative,
and therefore every possible real solution with distinct absolute values.

The linear-programming searches used to discover two certificates are not
part of the verification. In particular, no numerical feasibility verdict,
exceptional-ratio assumption, or rational-function specialization is used
in this proof. The independent
[generic rational-function computation](k5_variable_ratio_qq_generic_note.md)
is preserved with its own narrower specialization scope.
