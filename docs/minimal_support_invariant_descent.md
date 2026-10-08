# Factoring exact relations and projecting their row support

A small exact invariant relation can be replaced by an irreducible
invariant factor that still vanishes, with a bounded coefficient cost
for a fixed number of rows. Dropping unused rows then amplifies the
inherited block scale by an exact power of two. This is a valid route
to induction if a theorem for the resulting smaller support and
multidegree is available.

It does not, by itself, close the proper-restriction-image gap. An
explicit irreducible five-row invariant has actual rational zeros,
and its multiples give a six-dimensional numerical-zero module on
ten rows with proper first-smaller images everywhere. This example
is not an asymptotic full-profile configuration.

The exact algebra is checked by
[check_minimal_support_invariant_descent.py](check_minimal_support_invariant_descent.py).

## 1. Irreducible factors remain multihomogeneous invariants

Let `Q` be a nonzero integer simultaneous `SL_2` invariant, with
degree at most two in each binary row. Suppose that it vanishes at
the actual rows. Factoring over `Q` and clearing contents gives an
irreducible primitive integer factor `G` with

```text
G(P_1,...,P_m)=0.                                       (1)
```

The factor `G` is multihomogeneous. For each row separately, the
largest and smallest homogeneous degrees of a product are the sums
of the corresponding degrees of its factors. The product has only
one degree, so every factor has only one degree. Write these row
degrees as `d_i`; they belong to `{0,1,2}`.

The factor is also an `SL_2` invariant. The connected group cannot
nontrivially permute the finitely many irreducible factors, so it
preserves each factor up to a scalar character. There is no nontrivial
character of `SL_2`. Equivalently, apply its upper and lower unipotent
one-parameter groups. Each fixes the factor up to a polynomial scalar
`c(t)` satisfying `c(t+u)=c(t)c(u)` and `c(0)=1`; this forces `c=1`.
Those subgroups generate `SL_2`. The argument works over an algebraic
closure as well, and hence for the rational irreducible factor.

Here is an explicit, deliberately generous coefficient bound:

```text
C_G <= 2^(3^m-1) C_Q.                                  (2)
```

To see it, dehomogenize each row by `y_i=1`. Multihomogeneity means
the coefficient sum is unchanged. Substituting
`x_i=t^(3^(i-1))` is injective on the monomials, since each exponent
is at most two. It therefore preserves coefficient sums for `Q`,
`G`, and the remaining integer factor, and gives a univariate
factorization with degree at most `3^m-1`.

For an integer univariate factor `g` of `f`, write
`g=a product(t-alpha_j)`. Its Mahler measure is
`M(g)=|a| product max(1,|alpha_j|)`, and

```text
||g||_1 <= 2^deg(g) M(g) <= 2^deg(g) M(f)
         <= 2^deg(g) ||f||_1.
```

The middle inequality follows from multiplicativity and `M(h)>=1`
for the remaining nonzero integer factor. The last follows by
averaging `log|f|` on the unit circle: the average contribution of
`t-alpha` is `log max(1,|alpha|)`, and `|f|<=||f||_1` there.
This proves (2), without assuming a short generating set for the
relation ideal. In logarithmic height the cost is `O_m(1)`.

At distinct directions a bracket is nonzero; hence any product of
brackets, including a triangle product, is nonzero. Such a factor
can never be the numerical zero selected in (1). More generally,
the argument only selects a factor whose actual value is zero; it
does not presume that every other possible factor is nonzero.

## 2. Exact support projection and the four-row boundary

Let `J` be the support of `G`, with `s=|J|`. Keep the same actual
numerators `P_i` for `i in J`. Define new blocks by

```text
B_U=product_(T: T intersect J=U) H_T,       U subset J.  (3)
```

For every nonempty `U`, there are exactly `2^(m-s)` factors in
(3), and

```text
P_i=K_i product_(U containing i) B_U,
(1-eta)w_s <= log N(B_U) <= (1+eta)w_s       (U nonempty),
w_s=2^(m-s)w.                                         (4)
```

The new blocks and their conjugates remain pairwise coprime. The
individual corrections are still `K_i`, so their height bound stays
`sigma`, and all primitive pair numerators and residues are unchanged.
The block `B_J` divides every retained numerator; it is deliberately
retained. Dividing out that common Gaussian factor is unnecessary
and would change the displayed profile. The empty-intersection
product may be recorded as well, but never divides a retained row.
If the original block convention omits `H_empty`, it has only
`2^(m-s)-1` factors; no weight assertion for `B_empty` is needed or
made. Neither common row content nor a radius renormalization is
silently being introduced.

Thus a factor of coefficient height `h+O_m(1)` is even smaller
relative to its inherited block scale when `s<m`. Its row degrees
can mix ones and twos. It is not automatically a relation of degree
two in each retained row, and no theorem for that different
multidegree may be applied without justification.

An irreducible numerical zero must involve at least four rows.
For at most three rows, a fixed multidegree invariant is a scalar
multiple of a bracket monomial. For three rows its exponents are
`(d_1+d_2-d_3)/2`, `(d_1+d_3-d_2)/2`, and
`(d_2+d_3-d_1)/2`; if these are not nonnegative integers there is
no invariant. Every nonzero such monomial has nonzero actual value.

For exactly four supported rows, the irreducible numerical zero
factor has multidegree `(1,1,1,1)`. The only positive degree patterns
bounded by two and of even total degree are, up to permutation,

```text
(1,1,1,1),       (2,2,1,1),       (2,2,2,2).
```

Every invariant in the middle pattern has the bracket between the
two heavy rows as a factor, leaving the first pattern. In the last
pattern, the invariant is a homogeneous quadratic in two independent
matching invariants, for example
`X=Delta_12 Delta_34` and `Y=Delta_13 Delta_24`. At actual distinct
integer directions `X/Y` is a nonzero rational number. A numerical
zero of that quadratic gives a rational root, and therefore a
rational linear factor in `X,Y`. Irreducibility again reduces to
the first pattern.

The balanced evaluation map for four-row degree-one invariants is
injective. The complementary real divisors, with one correction
power per outside row, then give

```text
log C_G >= 2(1-eta)w_4-8sigma
        = 2^(m-3)(1-eta)w-8sigma.                       (5)
```

These degree and balanced-map facts are also recorded in
[invariant_degree_variation.md](invariant_degree_variation.md).
Combining (2) and (5), a numerical relation with

```text
log C_Q+(3^m-1)log 2 < 2^(m-3)(1-eta)w-8sigma
```

cannot have a numerical-zero irreducible factor supported on four
or fewer rows. This is a legitimate minimal-support conclusion,
not an induction covering supports of five or more.

## 3. An irreducible five-row relation with an actual rational zero

Consider the two complementary five-cycles in the complete graph
on five rows, and set

```text
G=12 Delta_12 Delta_23 Delta_34 Delta_45 Delta_51
   +Delta_13 Delta_35 Delta_52 Delta_24 Delta_41.         (6)
```

It has degree two in each row and coefficient sum at most `416`.
On the chart of directions `(infinity,0,1,a,b)`, it is

```text
g(a,b)=-12a^2+(12+13b-b^2)a-12b.                        (7)
```

This is irreducible even over `C[a,b]`. Indeed, its discriminant
as a quadratic in `a` is

```text
(12+13b-b^2)^2-576b.
```

If this were a square `S(b)^2`, then
`(12+13b-b^2-S)(12+13b-b^2+S)=576b`. At least one factor has
degree two, and the other is nonzero, so the product cannot have
degree one. A polynomial that is a square in `C(b)` is already a
polynomial square up to a constant; thus the discriminant is not
a square in the fraction field either.

Chart irreducibility also gives irreducibility of the original
invariant. Any factor invisible on this chart would have its zero
set on the normalization boundary and thus be a product of bracket
factors. No bracket divides (6): the two cycles partition all ten
edges, so on a generic collision at any one edge exactly one
monomial vanishes and the other does not. Equivalently, use the
open chart with its invertible normalizing brackets and the
multihomogeneous invariance of factors from Section 1.

Nevertheless `g(2,3)=0`. Applying the common projective map `t->2t`
gives the actual integer binary rows

```text
P_1=1,       P_2=i,       P_3=2+i,       P_4=4+i,
P_5=6+i,                                               (8)
```

at which `G=0`. Every row in (8) is conjugate-primitive with odd
norm, and their directions are distinct. They can be realized on
a common Gaussian lattice circle by clearing the denominators of
`P_i/bar(P_i)`. This is a fixed configuration, not a short-arc
family with uniform full block weights.

## 4. Proper images can persist after irreducible descent

Adjoin any five additional distinct rational directions, for example
the primitive rows `8+i,10+i,12+i,14+i,16+i`. Let `L` range over
all degree-two-each invariants on these five new rows. That space
has dimension six. The products

```text
Q_L=G L                                                (9)
```

are six independent actual numerical zero relations on the ten
rows. Multiplication by the nonzero polynomial `G` is injective;
independence here is independence of invariant polynomials, not
their numerical values. A basis of `L` consists of six suitable
triangle products times the squared bracket on the remaining pair,
so these products have coefficient sum at most `416*32=13312`.

At a first-smaller cut on ten rows, there are four inside rows.
If `a` of them lie in the first five rows, the two restrictions in
(9) have homogeneous degrees `5-2a` and `2a-3`. A nonzero
restriction requires both degrees to be nonnegative, so `a=2`.
Both factors then restrict linearly. The three outside rows of the
last factor have a two-dimensional translation-invariant linear
space. Consequently

```text
rank{F_(S,Q_L):L} <= 2 < 9=dim W_6                     (10)
```

at every one of the `210` first-smaller cuts. The checker verifies
this rank statement by exact substitution and also verifies the
six-dimensional multiplier space and all numerical zeros.

Thus a proper-image condition alone does not identify a removable
bracket or triangle factor: a genuinely irreducible numerical-zero
factor can underlie it. Applying a smaller-support arithmetic theorem
to that factor remains a viable route, because (4) amplifies its
profile scale, but such a theorem for the five-row pattern in (6)
is not proved here. Nor does (8)--(10) satisfy or refute the full
asymptotic profile hypotheses. The additional global arithmetic
needed for the uniform bound remains open.

Independent audit: the root agent checked the factor-height argument,
the support projection including the empty-block convention, the
four-row reduction, the irreducible five-cycle construction, and the
proper-image module, and reran the exact checker successfully.
