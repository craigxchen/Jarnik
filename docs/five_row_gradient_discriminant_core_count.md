# The five-row gradient discriminant core count

At an actual zero of a five-row degree-two-each integer invariant, its
row-gradient multipliers have a large forced ordinary divisor.  The
discriminant of the invariant as a quadratic in row `i` avoids dividing
individual graph derivatives by `P_i`: for every row,

```text
G_i=product_U n_U^max(0,2|U\{i}|-4) divides lambda_i,   (1)
sum_U max(0,2|U\{i}|-4)=24.
```

The available graph upper bound is about `32w`, while (1) is about
`24w`.  This exact count alone therefore gives no coefficient-height
gap and does not exclude a smooth irreducible five-row section.  The
finite polynomial-basis and literal derivative checks are in
[check_five_row_gradient_discriminant_core_count.py](check_five_row_gradient_discriminant_core_count.py).

## Actual rows and the integral gradient

Use the inherited five-row core from
[odd_support_relation_rigidity.md](odd_support_relation_rigidity.md) and
the pair-bracket scale of
[five_row_del_pezzo_arithmetic.md](five_row_del_pezzo_arithmetic.md):

```text
P_i=X_i+iY_i=K_i product_(U containing i) H_U,
n_U=Norm(H_U),
(1-eta)w<=log n_U<=(1+eta)w,
1<=|K_i|<=exp(sigma),
|Delta_ij|=b_ij product_(U containing i,j) n_U,
1<=b_ij<=exp(beta).                                    (2)
```

The core blocks and their conjugates have disjoint odd split-prime
support.  The inherited full cut is retained; the empty cut has no
incidence and may be omitted.  With primitive pair residues bounded by
`T`, one may take `beta=4sigma+log T`.  Every row lies in sixteen of the
thirty-two labelled cuts, so (2) gives

```text
|P_j|/|P_i|<=exp(16eta w+sigma).                         (3)
```

Let `Q` be a nonzero integer simultaneous `SL_2` invariant of degree
two in each of the five rows, with `Q(P)=0`.  Conjugate primitivity of
`P_i` implies `gcd(X_i,Y_i)=1`.  Euler's identity in row `i` yields a
unique integer `lambda_i` satisfying

```text
(Q_(X_i)(P),Q_(Y_i)(P))=lambda_i(-Y_i,X_i).             (4)
```

For `Q=a_i X_i^2+b_i X_iY_i+c_i Y_i^2`, its coefficients are integral
polynomials of degree two in each of the other four rows.  From
`Q(P)=0` and (4), the ordinary row discriminant obeys

```text
Disc_i(Q)(P)=b_i(P)^2-4a_i(P)c_i(P)=lambda_i^2.         (5)
```

This holds even when a multiplier is zero.  Simultaneous `SL_2`
invariance also gives the three exact stress equations

```text
sum_i lambda_i X_i^2=0,
sum_i lambda_i X_iY_i=0,
sum_i lambda_i Y_i^2=0.                                 (6)
```

They follow by applying the three infinitesimal `SL_2` motions to
`Q` and then using (4).  The normalized-chart ramification meaning of
`lambda_i=0` is established in
[five_row_singular_point_height.md](five_row_singular_point_height.md).

## Integral graph expansion and the forced divisor

The twenty-two degree-two-each bracket graph monomials span the
six-dimensional five-row invariant space.  Six of them have a selected
`6 x 6` polynomial coefficient minor of determinant `-1`; its inverse
has maximum column `l1` sum `2`.  Thus every **integer** `Q` has an
integral expansion in these six graphs, with graph coefficient sum

```text
H_G<=2C_Q,                                             (7)
```

where `C_Q` is the ordinary polynomial coefficient sum.  In row zero,
each of the twenty-one diagonal or cross terms in its discriminant
expands with **integer** coefficients in the four-row bracket graphs of
degree four in each row.  The checker constructs the literal
polynomials, finds a five-element four-row graph basis, and verifies
all twenty-one complete polynomial identities; the common denominator
is one.  Row permutation proves the same statement for every `i`.
This finite certificate is stronger than merely citing rational
`SL_2` invariance: it prevents any hidden fixed denominator in (1).

For a cut `U`, put `a=|U\{i}|`.  A four-row degree-four bracket graph
has at least

```text
g(a)=max(0,4a-8)                                       (8)
```

edges internal to the `a` rows in `U`.  This follows by counting their
`4a` incident half-edges: at most `4(4-a)` can cross to the other
rows.  The fifteen literal four-row multigraphs attain the minima

```text
a:       0  1  2  3  4
g(a):    0  0  0  4  8.
```

Each internal bracket is divisible by `n_U` exactly as in (2).
Distinct cut blocks have disjoint norm-prime support, so every graph
term, and hence the integral discriminant, is divisible by
`product_U n_U^g(a)`.  Since (5) is a square and every `g(a)` is even,
ordinary prime valuations give (1), with **no `K_i` denominator loss**.
Among thirty-two inherited cuts, eight have `a=3` and two have `a=4`;
their contributions are `8*2+2*4=24`.  In particular, if
`lambda_i!=0`,

```text
log|lambda_i|>=24(1-eta)w.                              (9)
```

The literal derivative route is weaker before cancellation.  Let
`L_i=Q_(X_i)+iQ_(Y_i)=iP_i lambda_i`.  Differentiating a graph term at
an edge `ij` replaces that bracket by `iP_j`, up to sign.  At a cut
`U`, after division by `P_i`, the **formal core-block exponents** in
its two Gaussian orientations have minima shown here, with
`a=|U\{i}|`.  The table suppresses correction and residue factors,
which can add valuations, including a denominator from `K_i`:

```text
                  a=0   a=1   a=2   a=3   a=4
i outside U, pi:    0     0     0     2     4
i outside U, bar:   0     0     0     1     3
i inside U, pi:    -1    -1     0     2     4
i inside U, bar:    0     0     0     2     4.
```

The negative entries are legitimate valuations of **individual**
derivative terms divided by `P_i`; only their sum is the integer
`lambda_i`.  The discriminant argument supplies the missing conjugate
powers without assuming termwise integrality or cancellation.

## Archimedean budget and limit of this route

Each graph has five pair brackets and two edges incident to row `i`.
Each term of its differentiated graph has four remaining brackets and
one factor `P_j/P_i`.  Equations (2), (3), and (7) therefore give the
explicit bound

```text
|lambda_i|<=4C_Q exp((32+48eta)w+4beta+sigma).         (10)
```

The `4` charges two derivative terms and the fixed graph coefficient
conversion; `beta` retains the actual pair correction and residue
cost.  Comparing (9) and (10) leaves an `8w` exponent margin before
errors.  It supplies no positive lower bound on `log C_Q`; a new
relation among the five multipliers or a shorter gradient estimate is
needed.  The Veronese stress equations (6) are exact simultaneous
compatibilities. The [cofactor construction](five_row_gradient_cofactor_stress_obstruction.md)
shows that they admit vectors satisfying (1) with every multiplier
nonzero and reduced exponent six, within the exponent-eight budget.
That construction need not be the gradient of any one small invariant.

The checker uses a separate integer primitive-row fixture to verify
(4)--(6) with all five multipliers nonzero.  That fixture is an algebraic
test of the identities, not an endpoint arc example.
