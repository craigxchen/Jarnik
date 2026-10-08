# Generic twelve-row coherent kernel and a six-row core audit

The terminal coherent evaluation kernel at twelve rows is generated, at
one explicit distinct rational configuration and on a nonempty Zariski-open
set, by six-row polar relations localized on six labels.  This is an exact
generic structural statement.  It supplies neither a coefficient-height
bound for sums of those generators nor a statement at **every** distinct
configuration, so it does not prove the radius-uniform arc bound.  A
separate six-row valuation audit shows that the gcd of *all* normalized
maximal minors has no compulsory factor from an isolated pure core; it
must not be confused with the large cofactor content of a selected row
set.

The later [odd-block balancing theorem](coherent_nagata_balancing.md)
does prove that the full coherent map has rank 121 at every distinct
twelve-row configuration. The subsequent
[determinantal generation theorem](distinct_configuration_determinantal_kernel_generation.md)
now proves that this particular localized family spans at every distinct
configuration. The finite calculations below remain generic certificates;
the new scheme argument supplies the all-point result. Neither argument
yet controls integral presentation heights.

The deterministic [checker](check_generic_terminal_coherent_syzygies.py)
uses exact arithmetic modulo `101` and completes in a few seconds.  A
nonzero minor modulo `101` is a nonzero integral minor and therefore gives
a characteristic-zero rank lower bound.  The checker uses NumPy only for
array operations; every reduction is integer arithmetic modulo `101`.

## 1. Twelve-row local relations

Let `K_2(12)` be the terminal invariant constituent.  By the
[Specht filtration](kernel_specht_filtration.md),

```text
dim K_2(12)=dim Specht(4,4,4)=462.
```

Fix distinct binary directions `P_i=(1,b_i)`, and let `E_P` send a
terminal invariant to its coherent four-label cut array.  For a sorted
six-set `I`, use the five integral six-row gradient invariants
`G_A,...,G_E` from [the Segre calculation](segre_gradient_arithmetic.md).
Let `M(P_I)=(A,B,C,D,E)` be its five matching values, and set

```text
Q_I(P)=A(P_I)G_A(I)+B(P_I)G_B(I)+C(P_I)G_C(I)
       +D(P_I)G_D(I)+E(P_I)G_E(I).
```

The exact [six-row coherent-kernel identity](six_row_coherent_kernel_and_content.md)
says that every two-label cut of `Q_I(P)` evaluates to zero at `P_I`.
For each of the five gradient invariants `G_j(I^c)` on the complementary
six-set, the product `Q_I(P)G_j(I^c)` belongs to `K_2(12)` and to
`ker E_P`.  Indeed, a nonzero four-label cut of a product of two
two-triangle invariants must select exactly two labels from each
six-set.  The first factor then has zero coherent evaluation.  There
are `binom(12,6)*5=4620` displayed generators; their span is denoted
`L_P`.

Take the explicit configuration `b_i=i`, `0<=i<12`.  The checker gives

```text
rank_F101 E_P = 121,
dim_F101 span{Q_I(P)G_j(I^c)} = 341.
```

For the first rank, it forms 462 products of four disjoint ternary
triangle determinants indexed by increasing `3×4` tableaux, then
computes their `495×462` coherent cut matrix.  Each column is a
member of `K_2(12)`.  The same 462 deterministic ternary test
configurations used below give these tableau products an evaluation
matrix of rank 462 modulo `101`.  Thus they are an integral source
basis over `Q`, using the known source dimension 462, without invoking
a separate standard-monomial basis theorem.  A coherent row for a
four-set `S` is zero unless each
triangle contains one label of `S`; in the nonzero case it is the
product of the four signed differences of the two outside directions.
For the second rank, it evaluates all 4620 localized generators at 462
deterministically generated ternary-vector configurations.  The
resulting polynomial-evaluation matrix has rank 341 modulo `101`.
Those ternary configurations are rank-separating test functionals, not
putative circle realizations.

The modular ranks imply characteristic-zero lower bounds 121 and 341.
Since `L_P⊂ker E_P` and `dim K_2(12)=462`, their sum already fills the
source.  Thus both ranks are exact over `Q`, and

```text
ker E_P = L_P,       rank E_P=121                 (b_i=i).
```

The two nonzero minors remain nonzero on a nonempty Zariski-open subset
of distinct configurations.  The same dimension argument proves the
equality there.  It does not prove constancy on the whole distinct-row
locus for the **localized-generator span**, where special arithmetic
configurations may lie.  A separate
[Cox–Nagata balancing theorem](coherent_nagata_balancing.md) proves that
the coherent *image* has rank 121 and equals the full joint
differential kernel at every distinct twelve-row configuration.  That
all-distinct image theorem does not establish that these 4620
localized generators span its source kernel at every configuration.
That additional assertion is now proved by the determinantal theorem
linked above.

As a separate finite check, the joint squarefree lowering operators
`D_1=Σ∂/∂x_i` and `D_b=Σb_i∂/∂x_i` on degree four in twelve variables
have common-kernel dimension 121 modulo `101` at this configuration.
The modular result makes the characteristic-zero common-kernel
dimension **at most** 121.  Since the 121-dimensional coherent image
lies in it, equality follows at this configuration and generically.
The cited balancing theorem supplies the stronger all-distinct
statement independently.

## 2. One-row pencil

Keep `b_1,...,b_11=1,...,11` fixed and let `b_0=t`.  Every coherent
cut coefficient is at most linear in `t`, so `E(t)=E_0+tE_1`.  The
checker computes, modulo `101`,

```text
rank E_0=121,       rank E_1=121,
rank [E_0;E_1]=176.
```

Here the semicolon means vertical stacking of the two `495×462`
matrices.  There is also an exact upper bound at this fixed set of
eleven directions.  Write the coherent polynomial as

```text
R_Q(x,t)=A_0(x_1,...,x_11)+t A_1(x_1,...,x_11)
         +x_0 B(x_1,...,x_11).
```

The terminal differential identities `D_1R_Q=D_bR_Q=0` give

```text
D_1' A_1=D_b' A_0=0,
C:=D_b' A_1=D_1' A_0=-B.
```

The degree-three polynomial `C` is killed by both primed operators.
The kernel of the map `(A_0,A_1,B)↦C` consists of two independent
degree-four joint-harmonic polynomials.  Direct modular reductions of
the two lowering matrices on the eleven fixed labels give common-kernel
dimensions `55` in degree four and `66` in degree three.  Those are
characteristic-zero *upper* bounds, so the stacked image has dimension
at most `2·55+66=176`.  The modular lower bound proves its exact
characteristic-zero rank 176.  The same argument applies on a nonempty
Zariski-open set.

There is an exact determinant-line consequence for this generic
one-row pencil.  Its constant polynomial kernel is

```text
C={v:E_0v=E_1v=0},       dim C=462−176=286.
```

Every localized generator is at most linear in `t`: when label `0`
belongs to its six-set, exactly one bracket in each six-row matching
value involves that label; otherwise its coefficient is constant.
The generators span the kernel over `Q(t)` by Section 1, so choose 55
linear ones independent modulo the field-span of `C`.  Their leading
coefficients are independent modulo `C`.  Otherwise a constant linear
combination of those generators, after subtracting `t` times a vector
of `C`, would be a constant kernel vector, contrary to independence.
Their wedge with a basis of `C` therefore has degree exactly 55.

That wedge has no nonconstant common polynomial factor.  Extend the
constant field to split a putative factor.  If the wedge then
vanished at a finite `t=alpha`, a constant
combination of the selected degree-at-most-one vectors and the basis
of `C` would equal `(t−alpha)v` for a constant vector `v`.  Since it is
in the polynomial kernel, `v` would also lie in `C`, again contrary to
the selected independence, which persists under field extension.
The primitive Plücker coordinates of the
kernel map to `Gr(341,462)` consequently have degree **55** in this
row parameter.  Equivalently, its polynomial-kernel minimal indices
are 286 zeros and 55 ones.  This is a generic determinant-line
degree, not a bound for the arithmetic height of an integral kernel
vector or the kernel covolume.  It does not assert the same rank
behavior at every distinct arithmetic configuration.  At a parameter
where the actual evaluation rank drops, the nonvanishing wedge gives
the extension of the generic Grassmann kernel, which may be smaller
than the full specialized kernel.  The argument with the other eleven
directions fixed as `1,...,11` persists on a Zariski-open set of their
configurations; label symmetry gives the same generic degree in each
row parameter.

## 3. Six-row pure-core prime valuation

For comparison, use the exact four-row tree and cofactor from
[the six-row kernel note](six_row_coherent_kernel_and_content.md):

```text
rows S = 12,45,23,34,
cofactor vector = ±F M,
F=Delta_12 Delta_56 Delta_36 Delta_16 Delta_45.
```

At a split Gaussian prime `pi` belonging to one core indexed by
`T⊂{1,...,6}`, assume the **pure-core** valuations

```text
v_pi(Delta_ij)=1 if {i,j}⊂T, and 0 otherwise.
```

This is the leading core contribution before possible pair-correction
factors.  Put `r=|T|`.  The exact integral matching-basis identities
show that the five components of `M` have common `pi`-order
`g=max(0,r−3)`: among all fifteen perfect matchings the minimum
number of edges internal to `T` is `g`, and the five matching
coordinates generate the fifteen integrally.

Divide each scalar cut row by its mandatory content
`D_S=pi^max(0,|T\S|−2)`.  If `f` is the number of factors of `F`
internal to `T` and `d` is the sum of the four selected row-content
exponents, the normalized cofactor vector has minimum order
`f−d+g`.  By permuting the tree labels, this order is zero for every
core size:

| `r` | tree-label permutation for `T={1,...,r}` | `f` | `d` | `g` |
|---:|:---|---:|---:|---:|
| 0, 1 | identity | 0 | 0 | 0 |
| 2 | `(1,3,2,4,5,6)` | 0 | 0 | 0 |
| 3 | identity | 1 | 1 | 0 |
| 4 | `(1,2,3,5,4,6)` | 1 | 2 | 1 |
| 5 | identity | 2 | 4 | 2 |
| 6 | identity | 5 | 8 | 3 |

Permuting labels acts by an integral unimodular change of the
six-row gradient basis.  Therefore a normalized `4×4` minor is a
`pi`-unit for each `T`.  The gcd of *all* normalized rank-four minors
has no mandatory pure-core factor at six rows.  If edge corrections
`b_ij` are `pi`-divisible, the selected permuted tree instead bounds
the extra valuation by at most twice the sum of their fifteen
valuations; no all-prime exact gcd formula is claimed.

This full-matrix gcd differs from the cofactor content of one chosen
four-row matrix.  The latter can be large and is the content relevant
when estimating the covolume of the kernel from that chosen matrix.
The six-row result therefore does not refute a refined selected-row
method at twelve rows.  The twelve-row generic spanning statement
likewise contains no coefficient-height control over cancellations
between localized generators.

## 4. Twelve-row pure-core initial ranks

There is a finite analogue of the six-row no-forced-factor result.
Work over a discrete valuation parameter `t` and fix a core
`T⊂{0,...,11}`.  In local projective binary coordinates, choose the
rows in `T` with slopes `t a_i` and the other rows with slopes `b_i`.
The leading bracket order is one on an edge wholly inside `T` and
zero otherwise, with generic nonzero leading coefficients.  In a
four-label cut `S`, a nonzero four-triangle graph restricts to a
perfect matching of the eight outside labels.  Its mandatory row
content has order

```text
d_S=max(0,|T\S|−4).
```

The checker constructs the matrix of leading coefficients after
dividing each row by `t^d_S`.  For `T={0,...,r−1}`, it uses
`a_i=i+1`, `b_i=i+17` modulo `101`; all required differences and
noncore slopes are units.  The normalized `495×462` initial matrix
has rank **121 modulo 101 for every `r=0,...,12`**.  Label symmetry
covers every core subset of a given size.  Thus at each core size a
maximal minor has a nonzero leading-coefficient polynomial in the
residue parameters: no extra core power divides *every* normalized
maximal minor as a polynomial identity.  The rank 121 holds on a
nonempty Zariski-open set of residues in each stratum.  The upper
bound 121 follows because all 122-minors of the generic coherent
matrix vanish identically by Section 1, and row normalization by
powers of `t` preserves that identity before specializing `t=0`.

An optional extended checker run,
`python3 docs/check_generic_terminal_coherent_syzygies.py --extended-core`,
also tests the localized syzygies at these boundary fixtures.  If a
six-set `I` meets `T` in `r_I` labels, every matching component of
`Q_I(t)` is divisible by `t^max(0,r_I−3)`.  Divide by that factor and
take its leading coefficient; the resulting six-row relation times
each complementary gradient invariant lies in the kernel of the
normalized initial cut map by specialization of the exact identities
for `t≠0`.  The 4620 resulting polynomial vectors have rank **341
modulo 101 for every core size `0,...,12`** at the displayed residues.
Together with initial-map rank 121 and source dimension 462, they
span the initial kernel there.  The nonzero-minor argument gives the
same equality on a nonempty Zariski-open set of residues in each
stratum.  This extended check takes about half a minute on the stated
runtime; the default checker does not run it.

This finite certificate does **not** say that every actual extracted
Gaussian residue tuple lies in those open sets.  It also concerns
the gcd over all maximal minors, not the selected-row cofactor
content needed for a kernel covolume estimate.  No radius exponent
gain follows from this test alone.
